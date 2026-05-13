#' Validate Model for Markov Simulation
#'
#' Internal helper that checks whether a model is compatible with Markov
#' simulation. For vglm models, verifies that the family is cumulative with
#' reverse = TRUE.
#'
#' @param object A fitted model object
#' @return NULL (invisibly) if valid, otherwise throws an error
#' @keywords internal
validate_markov_model <- function(object) {
  # Unwrap robcov_vglm if needed
  model_chk <- if (inherits(object, "robcov_vglm")) {
    object$vglm_fit
  } else {
    object
  }

  if (model_uses_offset(model_chk)) {
    stop_unsupported_offset()
  }

  # Check if it's a vglm/vgam model
  if (inherits(model_chk, c("vglm", "vgam"))) {
    # Extract family object
    fam <- model_chk@family

    # Check if it's cumulative family
    if (!grepl("cumulative", fam@vfamily[1], ignore.case = TRUE)) {
      stop(
        "For vglm/vgam models, only the cumulative family is supported.\n",
        "Current family: ",
        fam@vfamily[1],
        "\n",
        "Please refit your model with: family = cumulative(reverse = TRUE, ...)"
      )
    }

    # Check reverse coding
    # The reverse parameter is stored in the model's call, not in the family object
    model_call <- model_chk@call
    reverse_param <- NULL

    # Try to extract reverse from the family call
    if ("family" %in% names(model_call) && is.call(model_call$family)) {
      fam_args <- as.list(model_call$family)[-1] # Remove function name
      if ("reverse" %in% names(fam_args)) {
        reverse_param <- tryCatch(
          eval(fam_args$reverse),
          error = function(e) NULL
        )
      }
    }

    # If reverse was not explicitly set in the call, we cannot verify it
    # In this case, we should warn and reject (safer to be conservative)
    if (is.null(reverse_param)) {
      stop(
        "Cannot determine 'reverse' parameter from vglm model call.\n",
        "For safety, please explicitly specify: family = cumulative(reverse = TRUE, ...)\n",
      )
    }

    # Check that reverse = TRUE
    if (!isTRUE(reverse_param)) {
      stop(
        "For vglm/vgam models, the cumulative family must use reverse = TRUE.\n",
        "Current setting: reverse = ",
        reverse_param,
        "\n",
        "Please refit your model with: family = cumulative(reverse = TRUE, ...)"
      )
    }
  }

  # orm models are always compatible (they use reverse coding by default)
  # blrm models inherit from orm, so also OK
  # Other model types will be caught by the existing class check in soprob_markov

  invisible(NULL)
}

predict_vglm_response_markov <- function(object, newdata) {
  if (
    !isTRUE(attr(object, "markov_vglm")) &&
      !isTRUE(attr(object, "markov_split_assign"))
  ) {
    return(VGAM::predict(object, newdata, type = "response"))
  }

  tt <- stats::delete.response(stats::terms(object))
  X <- stats::model.matrix(
    tt,
    newdata,
    contrasts.arg = if (length(object@contrasts)) object@contrasts else NULL,
    xlev = object@xlevels
  )
  attr(X, "assign") <- object@assign

  lm2vlm <- utils::getFromNamespace("lm2vlm.model.matrix", "VGAM")
  X_vlm <- lm2vlm(
    X,
    Hlist = object@constraints,
    M = object@misc$M,
    xij = object@control$xij,
    Xm2 = NULL
  )

  eta <- matrix(
    X_vlm %*% VGAM::coefvlm(object),
    nrow = nrow(X),
    ncol = object@misc$M,
    byrow = TRUE
  )

  object@family@linkinv(eta = eta, extra = object@extra)
}

predict_orm_response_markov <- function(object, newdata) {
  Gamma <- get_effective_coefs(object)
  X <- orm_model_matrix(object, newdata, include_intercept = TRUE)
  X <- X[, colnames(Gamma), drop = FALSE]
  probs <- lp_to_probs(X %*% t(Gamma), nrow(Gamma))
  colnames(probs) <- as_state_labels(object$yunique)
  probs
}

select_posterior_draws <- function(model, n_draws = 100L, seed = NULL) {
  n_available <- nrow(model$draws)
  if (is.null(n_available) || n_available < 1) {
    stop("`blrm` model does not contain posterior draws.")
  }

  if (is.null(n_draws)) {
    n <- n_available
  } else {
    if (!is.numeric(n_draws) || length(n_draws) != 1 || is.na(n_draws) || n_draws < 1) {
      stop("`n_draws` must be a positive integer or NULL.")
    }
    n <- min(as.integer(n_draws), n_available)
  }

  if (n == n_available) {
    return(seq_len(n_available))
  }

  if (!is.null(seed)) {
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (old_seed_exists) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  sample.int(n_available, n)
}

blrm_design_matrix <- function(model, newdata, second = FALSE) {
  if (!requireNamespace("rms", quietly = TRUE)) {
    stop("Package 'rms' is required for `blrm` prediction.")
  }
  X <- rms::predictrms(model, newdata = newdata, type = "x", second = second)
  normalize_blrm_design_names(model, X, second = second)
}

normalize_blrm_design_names <- function(model, X, second = FALSE) {
  if (second) {
    return(X)
  }

  design <- model$Design
  rms_names <- design$colnames
  mm_names <- design$mmcolnames
  if (
    !is.null(rms_names) &&
      !is.null(mm_names) &&
      length(rms_names) == ncol(X) &&
      length(mm_names) == ncol(X) &&
      identical(colnames(X), mm_names)
  ) {
    colnames(X) <- rms_names
  }

  X
}

resolve_blrm_id_var <- function(model, newdata, id_var = NULL) {
  if (!is.null(id_var)) {
    return(id_var)
  }

  cluster_name <- model$clusterInfo$name
  if (!is.null(cluster_name) && length(cluster_name) == 1 && nzchar(cluster_name)) {
    return(cluster_name)
  }

  "id"
}

get_blrm_gamma_draws <- function(model, draw_indices) {
  if (!is.null(model$gamma_draws)) {
    gamma <- model$gamma_draws
    return(gamma[draw_indices, , drop = FALSE])
  }

  if (!requireNamespace("rmsb", quietly = TRUE)) {
    stop("Package 'rmsb' is required for `include_re = TRUE`.")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required for `include_re = TRUE`.")
  }

  stan_fit <- rmsb::stanGet(model)
  gamma <- tryCatch(
    {
      posterior::as_draws_matrix(stan_fit$draws(variables = "gamma"))
    },
    error = function(e) {
      mat <- as.matrix(stan_fit)
      mat[, grep("^gamma\\[", colnames(mat)), drop = FALSE]
    }
  )

  if (ncol(gamma) == 0) {
    stop("No random-effect draws named `gamma` were found in the `blrm` fit.")
  }

  gamma[draw_indices, , drop = FALSE]
}

cache_blrm_gamma_draws <- function(model, draw_indices, include_re) {
  if (!isTRUE(include_re)) {
    return(NULL)
  }

  gamma <- get_blrm_gamma_draws(model, draw_indices)
  rownames(gamma) <- as.character(draw_indices)
  gamma
}

subset_cached_draw_matrix <- function(x, all_draw_indices, draw_indices) {
  if (is.null(x)) {
    return(NULL)
  }

  idx <- match(draw_indices, all_draw_indices)
  if (anyNA(idx)) {
    stop("Cached posterior draw matrix is missing requested draw indices.")
  }

  x[idx, , drop = FALSE]
}

blrm_random_effect_matrix <- function(
  model,
  newdata,
  id_var,
  draw_indices,
  gamma_draws = NULL
) {
  cluster <- model$clusterInfo$cluster
  if (is.null(cluster) || length(cluster) == 0) {
    stop(
      "`include_re = TRUE` requires a `blrm` model fitted with `cluster()` ",
      "and stored cluster information."
    )
  }
  if (!id_var %in% names(newdata)) {
    stop("`newdata` must contain `", id_var, "` when `include_re = TRUE`.")
  }

  patient_ids <- as.character(newdata[[id_var]])

  gamma <- gamma_draws %||% get_blrm_gamma_draws(model, draw_indices)
  if (nrow(gamma) != length(draw_indices)) {
    stop("Random-effect draw matrix does not match requested posterior draws.")
  }

  gamma_names <- colnames(gamma)
  gamma_index <- NULL
  if (!is.null(gamma_names)) {
    cleaned_gamma_names <- sub("^gamma\\[([^]]+)\\]$", "\\1", gamma_names)
    if (all(patient_ids %in% gamma_names)) {
      gamma_index <- match(patient_ids, gamma_names)
    } else if (all(patient_ids %in% cleaned_gamma_names)) {
      gamma_index <- match(patient_ids, cleaned_gamma_names)
    }
  }
  if (is.null(gamma_index)) {
    cluster_levels <- levels(as.factor(cluster))
    gamma_index <- match(patient_ids, cluster_levels)
  }

  if (anyNA(gamma_index)) {
    missing_ids <- unique(patient_ids[is.na(gamma_index)])
    stop(
      "`newdata` contains IDs not present in the fitted `blrm` clusters: ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) " ..." else ""
    )
  }

  if (max(gamma_index) > ncol(gamma)) {
    stop("Random-effect draw matrix has fewer columns than fitted clusters.")
  }

  gamma[, gamma_index, drop = FALSE]
}

predict_blrm_response_markov <- function(
  object,
  newdata,
  include_re = FALSE,
  id_var = NULL,
  draw_indices = NULL,
  gamma_draws = NULL
) {
  if (is.null(draw_indices)) {
    draw_indices <- seq_len(nrow(object$draws))
  }
  manual_error <- NULL

  manual <- tryCatch(
    {
      X <- blrm_design_matrix(object, newdata = newdata, second = FALSE)
      draws <- object$draws[draw_indices, , drop = FALSE]
      ndraws <- nrow(draws)
      ns <- object$non.slopes
      K <- ns + 1L
      cn <- colnames(draws)
      tauinfo <- object$tauInfo
      tau_names <- if (!is.null(tauinfo) && length(tauinfo$name)) tauinfo$name else character()

      intercept_draws <- draws[, seq_len(ns), drop = FALSE]
      beta_names <- setdiff(cn, c(cn[seq_len(ns)], tau_names))
      if (length(beta_names) == 0) {
        xb <- matrix(0, nrow = ndraws, ncol = nrow(newdata))
      } else {
        beta_draws <- draws[, beta_names, drop = FALSE]
        missing_beta <- setdiff(beta_names, colnames(X))
        if (length(missing_beta) > 0) {
          stop(
            "Could not match `blrm` posterior coefficients to prediction ",
            "design columns. Missing columns: ",
            paste(utils::head(missing_beta, 10), collapse = ", "),
            if (length(missing_beta) > 10) " ..." else "",
            "."
          )
        }
        X <- X[, beta_names, drop = FALSE]
        xb <- beta_draws %*% t(X)
      }

      pppo <- object$pppo %||% 0L
      has_npo <- pppo > 0
      if (has_npo) {
        cppo <- object$cppo
        if (is.character(cppo)) {
          cppo <- eval(parse(text = cppo))
        }
        if (!is.function(cppo)) {
          stop("Only constrained partial proportional odds `blrm` models are supported.")
        }
        Z <- blrm_design_matrix(object, newdata = newdata, second = TRUE)
        tau_draws <- draws[, tau_names, drop = FALSE]
        zt <- tau_draws %*% t(Z)
        cppos <- vapply(object$ylevels[-1], cppo, numeric(1))
      }

      u_draws <- matrix(0, nrow = ndraws, ncol = nrow(newdata))
      if (include_re) {
        id_var <- resolve_blrm_id_var(object, newdata, id_var)
        u_draws <- blrm_random_effect_matrix(
          object,
          newdata,
          id_var,
          draw_indices,
          gamma_draws = gamma_draws
        )
      }

      ynam <- if (!is.null(object$yname)) {
        paste(object$yname, object$ylevels, sep = "=")
      } else {
        as_state_labels(object$ylevels)
      }
      out <- array(
        NA_real_,
        dim = c(ndraws, nrow(newdata), K),
        dimnames = list(draw_indices, rownames(newdata), ynam)
      )

      base_eta <- xb + u_draws
      cum_probs <- array(
        NA_real_,
        dim = c(ndraws, nrow(newdata), ns)
      )
      for (k in seq_len(ns)) {
        ep <- sweep(base_eta, 1L, intercept_draws[, k], "+")
        if (has_npo) {
          ep <- ep + cppos[k] * zt
        }
        cum_probs[, , k] <- stats::plogis(ep)
      }

      out[, , 1] <- 1 - cum_probs[, , 1]
      if (K > 2) {
        out[, , 2:(K - 1)] <- cum_probs[, , 1:(K - 2), drop = FALSE] -
          cum_probs[, , 2:(K - 1), drop = FALSE]
      }
      out[, , K] <- cum_probs[, , K - 1]
      out[out < 0] <- 0

      out
    },
    error = function(e) {
      manual_error <<- e
      NULL
    }
  )

  if (!is.null(manual)) {
    return(manual)
  }
  if (include_re) {
    stop(manual_error$message, call. = FALSE)
  }

  stop(
    "Manual `blrm` prediction failed: ",
    manual_error$message,
    call. = FALSE
  )
}

orm_model_matrix <- function(model, newdata, include_intercept = FALSE) {
  if (!inherits(model, "orm")) {
    stop("orm_model_matrix() requires an orm model.")
  }
  if (is.null(model$Design) || is.null(model$sformula)) {
    stop(
      "Fast orm prediction requires a model fitted by rms::orm() with ",
      "stored Design metadata."
    )
  }

  design <- model$Design
  newdata <- normalize_orm_prediction_data(model, newdata)

  oldopts <- options(
    contrasts = c(factor = "contr.treatment", ordered = "contr.treatment"),
    Design.attr = design
  )
  on.exit({
    options(contrasts = oldopts$contrasts)
    options(Design.attr = oldopts$Design.attr)
  }, add = TRUE)

  formulano <- add_rms_formula_helpers(model$sformula)
  Terms <- stats::terms(formulano, specials = "strat")
  attr(Terms, "response") <- 0L
  attr(Terms, "intercept") <- 1L

  strata_terms <- attr(Terms, "specials")$strat
  Terms_ns <- if (length(strata_terms)) {
    Terms[-strata_terms]
  } else {
    Terms
  }

  mf <- stats::model.frame(Terms, newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(Terms_ns, mf)[, -1L, drop = FALSE]

  expected_cols <- design$colnames
  if (length(expected_cols) != ncol(X)) {
    stop(
      "Could not reconstruct the orm design matrix. Expected ",
      length(expected_cols),
      " columns, got ",
      ncol(X),
      "."
    )
  }
  colnames(X) <- expected_cols

  if (include_intercept) {
    X <- cbind("(Intercept)" = 1, X)
  }

  X
}

normalize_orm_prediction_data <- function(model, newdata) {
  design <- model$Design
  categorical <- design$name[design$assume.code %in% c(5L, 8L)]

  for (nm in categorical) {
    if (!nm %in% names(newdata)) {
      next
    }

    levels_nm <- as.character(design$parms[[nm]])
    values <- newdata[[nm]]
    coerced <- factor(as.character(values), levels = levels_nm)
    bad <- is.na(coerced) & !is.na(values)
    if (any(bad)) {
      stop(
        "Values in `",
        nm,
        "` are not among fitted orm levels: ",
        paste(utils::head(unique(as.character(values[bad])), 5), collapse = ", "),
        if (length(unique(values[bad])) > 5) " ..." else ""
      )
    }
    newdata[[nm]] <- coerced
  }

  newdata
}

as_state_labels <- function(x) {
  as.character(x)
}

coerce_previous_state_values <- function(values, prototype, pvarname) {
  labels <- as_state_labels(values)

  if (is.factor(prototype)) {
    return(factor(
      labels,
      levels = levels(prototype),
      ordered = is.ordered(prototype)
    ))
  }

  if (is.integer(prototype)) {
    out <- suppressWarnings(as.integer(labels))
    if (anyNA(out) && any(!is.na(labels))) {
      stop(
        "`ylevels` must be integer-compatible when `", pvarname,
        "` is an integer previous-state variable."
      )
    }
    return(out)
  }

  if (is.numeric(prototype)) {
    out <- suppressWarnings(as.numeric(labels))
    if (anyNA(out) && any(!is.na(labels))) {
      stop(
        "`ylevels` must be numeric-compatible when `", pvarname,
        "` is a numeric previous-state variable."
      )
    }
    return(out)
  }

  if (is.character(prototype)) {
    return(labels)
  }

  labels
}

make_previous_state_column <- function(states, prototype, n, pvarname) {
  values <- coerce_previous_state_values(states, prototype, pvarname)
  rep(values, each = n)
}

normalize_previous_state_column <- function(values, prototype, pvarname) {
  coerce_previous_state_values(values, prototype, pvarname)
}

get_model_factor_levels <- function(model, varname) {
  model_chk <- if (inherits(model, "robcov_vglm")) {
    model$vglm_fit
  } else {
    model
  }

  if (inherits(model_chk, c("vglm", "vgam"))) {
    xlevels <- tryCatch(
      model_chk@xlevels,
      error = function(e) NULL
    )
    if (!is.null(xlevels[[varname]])) {
      return(as.character(xlevels[[varname]]))
    }
  }

  if (inherits(model_chk, c("orm", "blrm"))) {
    parms <- tryCatch(
      model_chk$Design$parms,
      error = function(e) NULL
    )
    if (!is.null(parms[[varname]])) {
      return(as.character(parms[[varname]]))
    }
  }

  NULL
}

get_time_info <- function(model, data, tvarname) {
  values <- if (!is.null(tvarname) && tvarname %in% names(data)) {
    data[[tvarname]]
  } else {
    NULL
  }
  model_levels <- if (!is.null(model)) {
    get_model_factor_levels(model, tvarname)
  } else {
    NULL
  }

  levels <- NULL
  ordered <- FALSE
  if (is.factor(values)) {
    levels <- model_levels %||% levels(values)
    ordered <- is.ordered(values)
  } else if (is.null(values)) {
    levels <- model_levels
  }

  list(
    is_factor = !is.null(levels),
    levels = levels,
    ordered = ordered
  )
}

default_sop_times <- function(data, tvarname, t_covs, default) {
  if (default == "max_sequence") {
    return(seq_len(max(data[[tvarname]], na.rm = TRUE)))
  }

  if (default == "fast") {
    if (!is.null(t_covs)) {
      return(seq_len(nrow(t_covs)))
    }
    if (!is.null(tvarname) && tvarname %in% names(data)) {
      return(sort(unique(data[[tvarname]])))
    }
    return(1)
  }

  if (is.null(tvarname) || !tvarname %in% names(data)) {
    stop("`times` must be specified if `tvarname` is not in data.")
  }

  sort(unique(data[[tvarname]]))
}

coerce_visit_times <- function(times, time_info, tvarname) {
  labels <- as.character(times)
  bad <- is.na(match(labels, time_info$levels))
  if (any(bad)) {
    bad_values <- unique(labels[bad])
    stop(
      "Requested visit values are not among fitted factor levels for `",
      tvarname,
      "`: ",
      paste(utils::head(bad_values, 5), collapse = ", "),
      if (length(bad_values) > 5) " ..." else ""
    )
  }

  factor(
    labels,
    levels = time_info$levels,
    ordered = isTRUE(time_info$ordered)
  )
}

resolve_sop_times <- function(
  model,
  data,
  times,
  tvarname,
  t_covs = NULL,
  default = c("unique", "max_sequence", "fast")
) {
  default <- match.arg(default)
  time_info <- get_time_info(model, data, tvarname)

  if (is.null(times)) {
    times <- if (time_info$is_factor) {
      time_info$levels
    } else {
      default_sop_times(data, tvarname, t_covs, default)
    }
  }

  if (time_info$is_factor) {
    times <- coerce_visit_times(times, time_info, tvarname)
  }

  if (!is.null(t_covs) && nrow(t_covs) != length(times)) {
    stop("`t_covs` must have one row per requested time point.")
  }

  list(times = times, time_info = time_info)
}

assign_sop_time <- function(data, tvarname, value, time_info) {
  if (isTRUE(time_info$is_factor)) {
    data[[tvarname]] <- factor(
      as.character(value),
      levels = time_info$levels,
      ordered = isTRUE(time_info$ordered)
    )
  } else {
    data[[tvarname]] <- value
  }
  data
}

validate_factor_gap <- function(gap, t_covs, time_info) {
  if (is.null(gap) || !isTRUE(time_info$is_factor)) {
    return(invisible(NULL))
  }

  if (
    is.null(t_covs) ||
      !gap %in% names(t_covs) ||
      !is.numeric(t_covs[[gap]])
  ) {
    stop(
      "Factor visit time with `gap` requires numeric gap values in ",
      "`t_covs[[\"", gap, "\"]]`; elapsed gaps are not inferred from ",
      "factor visit labels."
    )
  }

  invisible(NULL)
}

assign_sop_gap <- function(data, gap, times, index, t_covs, time_info) {
  if (is.null(gap)) {
    return(data)
  }

  validate_factor_gap(gap, t_covs, time_info)
  if (!isTRUE(time_info$is_factor) && (is.null(t_covs) || !gap %in% names(t_covs))) {
    data[[gap]] <- if (index == 1L) {
      times[index]
    } else {
      times[index] - times[index - 1L]
    }
  }

  data
}

assign_sop_t_covs <- function(data, t_covs, index) {
  if (is.null(t_covs)) {
    return(data)
  }

  for (nm in names(t_covs)) {
    data[[nm]] <- t_covs[index, nm]
  }
  data
}

assign_sop_visit <- function(
  data,
  tvarname,
  times,
  index,
  t_covs,
  gap,
  time_info
) {
  data <- assign_sop_time(data, tvarname, times[index], time_info)
  data <- assign_sop_gap(data, gap, times, index, t_covs, time_info)
  assign_sop_t_covs(data, t_covs, index)
}

#' Calculate State Occupation Probabilities for First-Order Markov Models
#'
#' Estimates state occupation probabilities over time by iterating a transition
#' matrix derived from a fitted model object (e.g., `vglm`, `rms`, or `rmsb`).
#' This function supports linear time iteration as well as complex, non-linear
#' time specifications (e.g., splines) via a covariate lookup table.
#'
#' @param object A fitted model object. Supported classes include:
#'   \code{"lrm"}, \code{"orm"} (from package `rms`),
#'   \code{"blrm"} (from package `rmsb`), and
#'   \code{"vglm"}, \code{"vgam"} (from package `VGAM`).
#' @param data A data frame containing the baseline covariates for the prediction.
#'   Rows represent unique patients. Columns must contain baseline covariates and
#'   initial values for time-varying variables.
#' @param times Visit-scale time points to iterate over. For numeric time
#'   variables this is a numeric vector. For factor-valued visit indices,
#'   values are matched to the fitted factor levels; if `NULL`, all fitted
#'   visit levels are used.
#' @param ylevels A character vector defining the names of the outcome levels (states).
#'   These must match the levels used in the fitted `object`.
#' @param absorb (Optional) A character vector of absorbing states (states from which
#'   transitions out are impossible). Defaults to \code{NULL}.
#' @param tvarname A character string specifying the name of the time variable in the
#'   model formula. Defaults to \code{"time"}.
#' @param pvarname A character string specifying the name of the previous state variable
#'   used in the model formula. Defaults to \code{"yprev"}.
#' @param p2varname Optional character string specifying the second previous
#'   state variable. `NULL` uses first-order recursion; any non-`NULL` value
#'   uses second-order recursion with histories `(p2varname, pvarname)`.
#' @param gap (Optional) A character string specifying the name of the variable representing
#'   the time gap (delta time) between observations, if used in the model. Defaults to \code{NULL}.
#' @param t_covs (Optional) A data frame used for explicit time-basis columns
#'   or other time-varying covariates.
#'   \itemize{
#'     \item **Structure:** The number of rows in \code{t_covs} must exactly match the length of \code{times}.
#'     \item **Columns:** Column names must match the specific basis variables used in the model formula (e.g., \code{t1}, \code{t2}).
#'     \item **Usage:** At step \code{i}, the values from the \code{i}-th row of \code{t_covs} are injected into the prediction data.
#'   }
#'   Inline model terms such as \code{rms::rcs(time, 4)} do not require
#'   \code{t_covs}; prediction reuses the fitted model's stored transform
#'   metadata.
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions. This does not
#'   integrate over the random-effects distribution; it uses the fitted effects
#'   for known IDs.
#' @param id_var Character ID column used when `include_re = TRUE`. If `NULL`
#'   for `blrm`, inferred from `model$clusterInfo$name` when available,
#'   otherwise `"id"`.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible `blrm` draw sampling.
#' @param ... Reserved for internal use.
#'
#' @details
#'
#' \strong{1. Data Expansion:}
#' We construct a "long" expansion dataset at every time point. For \eqn{N} patients
#' and \eqn{K} non-absorbing states, the expansion dataset contains \eqn{N \times K} rows.
#' \itemize{
#'   \item Rows \eqn{1 \dots N}: All patients assuming \eqn{y_{prev} = \text{State } 1}
#'   \item Rows \eqn{(N+1) \dots 2N}: All patients assuming \eqn{y_{prev} = \text{State } 2}
#'   \item ... and so on.
#' }
#' This allows a single \code{predict()} call to generate transition probabilities for the entire
#' cohort for all possible previous states in one step.
#'
#' \strong{2. Element-wise weighted sum}
#' We use an element-wise weighted sum approach based on the
#' Law of Total Probability:
#' \deqn{P(S_t = k) = \sum_{j} P(S_{t-1} = j) \times P(S_t = k | S_{t-1} = j)}
#'
#' where:
#' \itemize{
#'   \item \eqn{S_t}: State occupied at time \eqn{t}.
#'   \item \eqn{S_{t-1}}: State occupied at time \eqn{t-1}.
#'   \item \eqn{k}: The target state at time \eqn{t}.
#'   \item \eqn{j}: The origin state at time \eqn{t-1}.
#' }
#'
#' Let \eqn{\mathbf{v}_{prev, j}} be a vector of length \eqn{N} containing the probability that each
#' patient was in state \eqn{j} at time \eqn{t-1}.
#' Let \eqn{\mathbf{M}_{trans, j \to k}} be a vector of length \eqn{N} containing the transition
#' probability from \eqn{j} to \eqn{k} for each patient (derived from the batched prediction).
#'
#' The probability of being in state \eqn{k} at time \eqn{t} is updated as:
#' \deqn{\mathbf{v}_{curr, k} = \sum_{j} (\mathbf{v}_{prev, j} \odot \mathbf{M}_{trans, j \to k})}
#' where \eqn{\odot} denotes element-wise multiplication.
#'
#' \strong{3. Absorbing States:}
#' If an absorbing state \eqn{a} is present, the update logic handles the accumulation of probability mass:
#' \deqn{P(S_t = a) = P(S_{t-1} = a) + \sum_{j \neq a} P(S_{t-1} = j) \times P(S_t = a | S_{t-1} = j)}
#'
#' @return
#' An array of state probabilities.
#' \itemize{
#'   \item **Frequentist fit:** An array of dimension \code{[n_patients x n_times x n_states]}.
#'   \item **Bayesian fit:** An array of dimension \code{[n_draws x n_patients x n_times x n_states]}.
#' }
#'
#' @export
soprob_markov <- function(
  object,
  data,
  times = NULL,
  ylevels,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  p2varname = NULL,
  gap = NULL,
  t_covs = NULL,
  include_re = FALSE,
  id_var = NULL,
  n_draws = 100L,
  seed = NULL,
  ...
) {
  # --- 1. Initial Checks & Setup ---
  dots <- list(...)
  unknown_dots <- setdiff(names(dots), c(".draw_indices", ".gamma_draws"))
  if (length(unknown_dots) > 0) {
    stop("Unused arguments: ", paste(unknown_dots, collapse = ", "))
  }
  draw_indices_arg <- dots$.draw_indices
  gamma_draws_arg <- dots$.gamma_draws

  cl <- if (inherits(object, "blrm")) {
    "blrm"
  } else if (inherits(object, "robcov_vglm")) {
    "robcov_vglm"
  } else if (inherits(object, "orm")) {
    "orm"
  } else if (inherits(object, "vglm")) {
    "vglm"
  } else if (inherits(object, "vgam")) {
    "vgam"
  } else {
    class(object)[1]
  }
  ftypes <- c(
    orm = "rms",
    blrm = "rmsb",
    vglm = "vgam",
    vgam = "vgam",
    robcov_vglm = "robcov"
  )
  ftype <- ftypes[cl]

  if (is.na(ftype)) {
    stop("Object class not supported")
  }

  # Validate model is compatible with Markov simulation
  validate_markov_model(object)

  if (!pvarname %in% names(data)) {
    stop("Previous-state variable `", pvarname, "` not found in `data`.")
  }
  if (!is.null(p2varname) && !p2varname %in% names(data)) {
    stop("Second previous-state variable `", p2varname, "` not found in `data`.")
  }
  time_res <- resolve_sop_times(
    object,
    data,
    times,
    tvarname,
    t_covs = t_covs,
    default = "unique"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap, t_covs, time_info)

  draw_indices <- NULL
  if (ftype == "rmsb") {
    draw_indices <- draw_indices_arg %||% select_posterior_draws(object, n_draws, seed)
  }

  # Define prediction function
  prd <- switch(
    ftype,
    rms = function(obj, d) predict_orm_response_markov(obj, d),
    vgam = function(obj, d) predict_vglm_response_markov(obj, d),
    rmsb = function(obj, d) {
      predict_blrm_response_markov(
        obj,
        d,
        include_re = include_re,
        id_var = id_var,
        draw_indices = draw_indices,
        gamma_draws = gamma_draws_arg
      )
    },
    robcov = function(obj, d) {
      if (is.null(obj$vglm_fit)) {
        stop(
          "robcov_vglm object does not contain the original vglm fit. ",
          "Please re-run robcov_vglm() with the latest version of the package."
        )
      }
      predict_vglm_response_markov(obj$vglm_fit, d)
    }
  )

  # Prepare dimensions
  n_pat <- nrow(data)
  n_times <- length(times)
  ylevel_names <- as_state_labels(ylevels)
  absorb_names <- as_state_labels(absorb)
  n_states <- length(ylevel_names)
  yna <- ylevel_names[!ylevel_names %in% absorb_names] # Non-absorbing states
  n_yna <- length(yna)

  # Check Bayesian draws
  nd <- if (ftype == "rmsb" && length(object$draws)) length(draw_indices) else 0

  # Initialize Output Array
  # Structure: [Patients, Time, States]
  if (nd == 0) {
    P <- array(
      0,
      dim = c(n_pat, n_times, n_states),
      dimnames = list(rownames(data), times, ylevel_names)
    )
  } else {
    P <- array(
      0,
      dim = c(nd, n_pat, n_times, n_states),
      dimnames = list(draw_indices, rownames(data), times, ylevel_names)
    )
  }

  # --- 2. Time 1 Initialization (Start) ---
  # Update time variables for T1
  data <- assign_sop_visit(
    data,
    tvarname = tvarname,
    times = times,
    index = 1L,
    t_covs = t_covs,
    gap = gap,
    time_info = time_info
  )

  # Predict probabilities at T1
  p_t1 <- prd(object, data) # Returns [n_pat x n_states] or [nd x n_pat x n_states]

  if (nd == 0) {
    P[, 1, ] <- p_t1
  } else {
    P[,, 1, ] <- p_t1
  }

  if (!is.null(p2varname)) {
    return(soprob_markov_second_order_run(
      P = P,
      object = object,
      data = data,
      prd = prd,
      nd = nd,
      times = times,
      ylevel_names = ylevel_names,
      absorb_names = absorb_names,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      time_info = time_info
    ))
  }

  # --- 3. Prepare Expansion Data for Transitions ---
  # We create a long dataframe where every patient is repeated for every possible PREVIOUS state.
  # Order: Patient 1 (State A), Patient 2 (State A)... Patient 1 (State B)...
  # This ordering is crucial for the vectorized update logic later.

  edata_base <- data[rep(1:n_pat, times = n_yna), , drop = FALSE]

  # Assign the previous state variable
  # We repeat each state n_pat times

  edata_base[[pvarname]] <- make_previous_state_column(
    states = yna,
    prototype = data[[pvarname]],
    n = n_pat,
    pvarname = pvarname
  )

  # --- 4. Iterate Through Time (Vectorized over Patients) ---
  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      tvarname = tvarname,
      times = times,
      index = it,
      t_covs = t_covs,
      gap = gap,
      time_info = time_info
    )

    # Get Transition Probabilities
    # This returns matrix: [ (n_pat * n_yna) x n_states ]
    # Rows 1:N are transitions from State 1, Rows N+1:2N from State 2, etc.
    trans_probs <- prd(object, edata_base)

    # --- 5. The Update Step (Markov) ---
    # Formula: P(S_t = k) = Sum_over_j [ P(S_t-1 = j) * P(S_t=k | S_t-1=j) ]

    if (nd == 0) {
      # Initialize current time probabilities with 0
      p_current <- matrix(0, nrow = n_pat, ncol = n_states)
      colnames(p_current) <- ylevel_names

      # Add contribution from Non-Absorbing Previous States
      for (i in seq_along(yna)) {
        prev_state_name <- yna[i]

        # Extract prob of being in this previous state for all patients
        # Vector of length n_pat
        prob_prev <- P[, it - 1, prev_state_name]

        # Extract transition probs GIVEN this previous state
        # Rows correspond to the block for this state in edata_base
        row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
        probs_transition <- trans_probs[row_indices, ]

        # Weighted sum: multiply column-wise by the probability of being in prev state
        p_current <- p_current + (probs_transition * prob_prev)
      }

      # Add contribution from Absorbing States (if any)
      # Absorbing states transition to themselves with Prob 1
      if (!is.null(absorb)) {
        for (a_state in absorb) {
          # If you were in absorb state at t-1, you are in absorb state at t
          p_current[, a_state] <- p_current[, a_state] + P[, it - 1, a_state]
        }
      }

      P[, it, ] <- p_current
    } else {
      # trans_probs is [draw x (patient * previous state) x state].
      p_current <- array(0, dim = c(nd, n_pat, n_states))
      dimnames(p_current) <- list(dimnames(P)[[1]], rownames(data), ylevel_names)

      for (i in seq_along(yna)) {
        prev_state_name <- yna[i]
        prob_prev <- P[, , it - 1, prev_state_name]
        row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
        probs_transition <- trans_probs[, row_indices, , drop = FALSE]

        for (k in seq_len(n_states)) {
          p_current[, , k] <- p_current[, , k] +
            probs_transition[, , k] * prob_prev
        }
      }

      if (!is.null(absorb)) {
        for (a_state in absorb) {
          # total dead at t = new deaths + already dead
          p_current[, , a_state] <- p_current[, , a_state] +
            P[, , it - 1, a_state]
        }
      }
      P[, , it, ] <- p_current
    }
  }

  return(P)
}

match_state_indices <- function(values, ylevel_names, varname) {
  idx <- match(as.character(values), ylevel_names)
  if (anyNA(idx)) {
    bad <- unique(as.character(values[is.na(idx)]))
    stop(
      "Values in `", varname, "` are not among `ylevels`: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..." else ""
    )
  }
  idx
}

soprob_markov_second_order_run <- function(
  P,
  object,
  data,
  prd,
  nd,
  times,
  ylevel_names,
  absorb_names,
  tvarname,
  pvarname,
  p2varname,
  gap,
  t_covs,
  time_info
) {
  n_pat <- nrow(data)
  n_times <- length(times)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)
  prev_idx <- match_state_indices(data[[pvarname]], ylevel_names, pvarname)

  if (nd == 0) {
    joint_prev <- array(
      0,
      dim = c(n_pat, n_states, n_states),
      dimnames = list(rownames(data), ylevel_names, ylevel_names)
    )
    for (i in seq_len(n_pat)) {
      if (prev_idx[i] %in% absorb_idx) {
        joint_prev[i, prev_idx[i], prev_idx[i]] <- 1
        P[i, 1, ] <- 0
        P[i, 1, prev_idx[i]] <- 1
      } else {
        joint_prev[i, prev_idx[i], ] <- P[i, 1, ]
      }
    }
  } else {
    joint_prev <- array(
      0,
      dim = c(nd, n_pat, n_states, n_states),
      dimnames = list(dimnames(P)[[1]], rownames(data), ylevel_names, ylevel_names)
    )
    for (i in seq_len(n_pat)) {
      if (prev_idx[i] %in% absorb_idx) {
        joint_prev[, i, prev_idx[i], prev_idx[i]] <- 1
        P[, i, 1, ] <- 0
        P[, i, 1, prev_idx[i]] <- 1
      } else {
        joint_prev[, i, prev_idx[i], ] <- P[, i, 1, ]
      }
    }
  }

  if (n_times < 2) {
    return(P)
  }

  pair_grid <- expand.grid(
    h = seq_len(n_states),
    j = seq_len(n_states),
    KEEP.OUT.ATTRS = FALSE
  )
  n_pairs <- nrow(pair_grid)
  predictable_pair <- pair_grid$h %in% non_absorb_idx &
    pair_grid$j %in% non_absorb_idx
  block_rows <- lapply(seq_len(n_pairs), function(i) {
    ((i - 1) * n_pat + 1):(i * n_pat)
  })

  edata_base <- data[rep(seq_len(n_pat), times = n_pairs), , drop = FALSE]
  edata_base[[p2varname]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$h],
    prototype = data[[p2varname]],
    n = n_pat,
    pvarname = p2varname
  )
  edata_base[[pvarname]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$j],
    prototype = data[[pvarname]],
    n = n_pat,
    pvarname = pvarname
  )

  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      tvarname = tvarname,
      times = times,
      index = it,
      t_covs = t_covs,
      gap = gap,
      time_info = time_info
    )

    predict_rows <- unlist(block_rows[predictable_pair], use.names = FALSE)
    trans_probs <- prd(object, edata_base[predict_rows, , drop = FALSE])
    cursor <- 1L

    if (nd == 0) {
      joint_current <- array(0, dim = c(n_pat, n_states, n_states))
      for (pair_i in seq_len(n_pairs)) {
        h <- pair_grid$h[pair_i]
        j <- pair_grid$j[pair_i]
        prob_prev <- joint_prev[, h, j]
        if (!any(prob_prev > 0)) {
          if (predictable_pair[pair_i]) {
            cursor <- cursor + n_pat
          }
          next
        }
        if (j %in% absorb_idx) {
          joint_current[, j, j] <- joint_current[, j, j] + prob_prev
        } else if (predictable_pair[pair_i]) {
          rows <- cursor:(cursor + n_pat - 1L)
          transition <- trans_probs[rows, , drop = FALSE]
          for (l in seq_len(n_states)) {
            joint_current[, j, l] <- joint_current[, j, l] + transition[, l] * prob_prev
          }
          cursor <- cursor + n_pat
        }
      }
      P[, it, ] <- apply(joint_current, c(1, 3), sum)
      joint_prev <- joint_current
    } else {
      joint_current <- array(0, dim = c(nd, n_pat, n_states, n_states))
      for (pair_i in seq_len(n_pairs)) {
        h <- pair_grid$h[pair_i]
        j <- pair_grid$j[pair_i]
        if (j %in% absorb_idx) {
          joint_current[, , j, j] <- joint_current[, , j, j] +
            joint_prev[, , h, j]
        } else if (predictable_pair[pair_i]) {
          rows <- cursor:(cursor + n_pat - 1L)
          prob_prev <- joint_prev[, , h, j]
          if (any(prob_prev > 0)) {
            transition <- trans_probs[, rows, , drop = FALSE]
            for (l in seq_len(n_states)) {
              joint_current[, , j, l] <- joint_current[, , j, l] +
                transition[, , l] * prob_prev
            }
          }
          cursor <- cursor + n_pat
        }
      }
      P[, , it, ] <- apply(joint_current, c(1, 2, 4), sum)
      joint_prev <- joint_current
    }
  }

  P
}

#' Compute standardized state occupancy probabilities (Vectorized)
#'
#' Updates the standardize_sops wrapper to use the optimized
#' soprob_markov_vectorized function.
#'
#' @param model A fitted model object (orm, vglm, rmsb). For `vglm` models,
#'   the family must be `cumulative(reverse = TRUE, ...)`.
#' @param data A data frame containing patient trajectory data.
#' @param times Visit-scale time points to predict. Factor-valued visit
#'   indices use fitted factor levels when `times = NULL`.
#' @param ylevels States in the data.
#' @param absorb Absorbing state name.
#' @param varnames List of variable names: tvarname, pvarname, p2varname
#'   (optional second previous-state variable), id, tx, gap (optional).
#' @param t_covs Time-dependent covariate lookup table.
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions.
#' @param id_var Optional patient ID variable. Defaults to `varnames$id`.
#' @param n_draws Integer number of `blrm` posterior draws to sample, or `NULL`
#'   to use all stored draws.
#' @param seed Optional random seed for reproducible `blrm` draw sampling.
#'
#' @return A list with two components:
#'   \item{sop_tx}{Matrix (Time x State) or Array (Draws x Time x State) for Treatment}
#'   \item{sop_ctrl}{Matrix (Time x State) or Array (Draws x Time x State) for Control}
#'
#' @export
standardize_sops <- function(
  model,
  data = NULL,
  times = NULL,
  ylevels = factor(1:6),
  absorb = 6,
  varnames = list(
    tvarname = "time",
    pvarname = "yprev",
    p2varname = NULL,
    id = "id",
    tx = "tx",
    gap = NULL
  ),
  t_covs = NULL,
  include_re = FALSE,
  id_var = NULL,
  n_draws = 100L,
  seed = NULL
) {
  # --- 1. Setup & Data Extraction ---
  # Validate model compatibility
  validate_markov_model(model)

  if (
    !inherits(model, "orm") &&
      !inherits(model, "vglm") &&
      !inherits(model, "blrm") &&
      !inherits(model, "robcov_vglm")
  ) {
    stop("model must be an orm, rmsb, vglm, or robcov_vglm object.")
  }

  if (is.null(data)) {
    data <- model$x
    if (is.null(data)) {
      stop("No data provided and model was not fitted with x=TRUE")
    }
    data$y <- model$y
  }

  # Variable mapping
  id_var <- id_var %||% varnames$id
  tx_var <- varnames$tx
  tvar <- varnames$tvarname
  pvar <- varnames$pvarname
  p2var <- varnames$p2varname %||% NULL
  gap_var <- varnames$gap %||% NULL # Helper if varnames$gap is missing

  time_res <- resolve_sop_times(
    model,
    data,
    times,
    tvar,
    t_covs = t_covs,
    default = "max_sequence"
  )
  times <- time_res$times
  validate_factor_gap(gap_var, t_covs, time_res$time_info)

  if (!pvar %in% names(data)) {
    stop("Previous-state variable `", pvar, "` not found in `data`.")
  }

  draw_indices <- NULL
  gamma_draws <- NULL
  if (inherits(model, "blrm")) {
    draw_indices <- select_posterior_draws(model, n_draws, seed)
    gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  }

  # --- 2. Create Counterfactual Cohorts ---
  # Extract one row per patient (baseline)
  X_base <- data[!duplicated(data[[id_var]]), ]

  # Create Treated Cohort (Everyone tx=1)
  X_tx <- X_base
  X_tx[[tx_var]] <- 1

  # Create Control Cohort (Everyone tx=0)
  X_ctrl <- X_base
  X_ctrl[[tx_var]] <- 0

  # --- 3. Vectorized Prediction ---
  # These calls replace the loop.
  # Returns: [Patients x Time x States] (Freq) OR [Draws x Patients x Time x States] (Bayes)

  res_tx <- soprob_markov(
    object = model,
    data = X_tx,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvar,
    pvarname = pvar,
    p2varname = p2var,
    gap = gap_var,
    t_covs = t_covs,
    include_re = include_re,
    id_var = id_var,
    n_draws = n_draws,
    seed = seed,
    .draw_indices = draw_indices,
    .gamma_draws = gamma_draws
  )

  res_ctrl <- soprob_markov(
    object = model,
    data = X_ctrl,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvar,
    pvarname = pvar,
    p2varname = p2var,
    gap = gap_var,
    t_covs = t_covs,
    include_re = include_re,
    id_var = id_var,
    n_draws = n_draws,
    seed = seed,
    .draw_indices = draw_indices,
    .gamma_draws = gamma_draws
  )

  # --- 4. Marginalize (Average over Patients) ---

  # Helper to average over the patient dimension
  # Freq dims: [Pat, Time, State] -> Target: [Time, State] (Ave over dim 1)
  # Bayes dims: [Draw, Pat, Time, State] -> Target: [Draw, Time, State] (Ave over dim 2)

  calc_marginal <- function(arr) {
    ndims <- length(dim(arr))
    if (ndims == 3) {
      # Frequentist: Average over rows (Patients)
      # Result: [Time x States]
      return(apply(arr, c(2, 3), mean, na.rm = TRUE))
    } else if (ndims == 4) {
      # Bayesian: Average over dim 2 (Patients)
      # Result: [Draws x Time x States]
      return(apply(arr, c(1, 3, 4), mean, na.rm = TRUE))
    } else {
      stop("Unexpected array dimensions returned from vectorized function.")
    }
  }

  sop_tx_marg <- calc_marginal(res_tx)
  sop_ctrl_marg <- calc_marginal(res_ctrl)

  # Ensure column names
  if (length(dim(sop_tx_marg)) == 2) {
    colnames(sop_tx_marg) <- as.character(ylevels)
    colnames(sop_ctrl_marg) <- as.character(ylevels)
  } else {
    dimnames(sop_tx_marg)[[3]] <- as.character(ylevels)
    dimnames(sop_ctrl_marg)[[3]] <- as.character(ylevels)
  }

  return(list(
    sop_tx = sop_tx_marg,
    sop_ctrl = sop_ctrl_marg
  ))
}

standardize_time_map <- function(time_map) {
  if (is.null(time_map)) {
    stop("`time_map` must be supplied.")
  }

  if (is.data.frame(time_map)) {
    if (all(c("visit", "real_time") %in% names(time_map))) {
      visit <- time_map$visit
      real_time <- time_map$real_time
    } else if (all(c("visit", "time") %in% names(time_map))) {
      visit <- time_map$visit
      real_time <- time_map$time
    } else if (ncol(time_map) >= 2) {
      visit <- time_map[[1]]
      real_time <- time_map[[2]]
    } else {
      stop("`time_map` data frames must have at least two columns.")
    }
  } else if (is.atomic(time_map) && !is.null(names(time_map))) {
    visit <- names(time_map)
    real_time <- time_map
  } else {
    stop(
      "`time_map` must be a named numeric vector or a data frame with ",
      "visit and real-time columns."
    )
  }

  real_time <- suppressWarnings(as.numeric(real_time))
  visit <- as.character(visit)
  bad <- is.na(visit) | !nzchar(visit) | is.na(real_time) | !is.finite(real_time)
  if (any(bad)) {
    stop("`time_map` contains missing or non-finite visit/time values.")
  }
  if (anyDuplicated(visit)) {
    dup <- unique(visit[duplicated(visit)])
    stop("`time_map` contains duplicate visit entries: ", paste(dup, collapse = ", "))
  }

  data.frame(visit = visit, real_time = real_time, stringsAsFactors = FALSE)
}

map_sop_time_values <- function(times, time_map) {
  labels <- as.character(times)
  idx <- match(labels, time_map$visit)
  if (anyNA(idx)) {
    missing <- unique(labels[is.na(idx)])
    stop(
      "`time_map` is missing entries for visit value(s): ",
      paste(missing, collapse = ", ")
    )
  }
  time_map$real_time[idx]
}

validate_sop_xout <- function(xout, lower, upper) {
  if (!is.numeric(xout) || anyNA(xout) || any(!is.finite(xout))) {
    stop("`xout` must be a finite numeric vector.")
  }
  if (any(xout < lower | xout > upper)) {
    stop(
      "`xout` must stay within the supported time range [",
      lower,
      ", ",
      upper,
      "]."
    )
  }
  sort(unique(xout))
}

sop_measure_cols <- function(x) {
  intersect(c("estimate", "conf.low", "conf.high", "std.error", "draw"), names(x))
}

split_key <- function(data, cols) {
  if (length(cols) == 0) {
    return(rep("all", nrow(data)))
  }
  do.call(interaction, c(data[, cols, drop = FALSE], drop = TRUE, sep = "\r"))
}

coerce_state_like <- function(labels, prototype) {
  if (is.factor(prototype)) {
    return(factor(
      labels,
      levels = levels(prototype),
      ordered = is.ordered(prototype)
    ))
  }
  if (is.integer(prototype)) {
    return(as.integer(labels))
  }
  if (is.numeric(prototype)) {
    return(as.numeric(labels))
  }
  labels
}

rows_match_values <- function(data, values) {
  keep <- rep(TRUE, nrow(data))
  for (nm in names(values)) {
    keep <- keep & as.character(data[[nm]]) == as.character(values[[nm]])
  }
  keep
}

baseline_rows_for_anchor <- function(x, id_var) {
  newdata_orig <- attr(x, "newdata_orig")
  if (is.null(newdata_orig)) {
    stop(
      "`origin_time` requires SOP output with stored `newdata_orig` ",
      "attributes."
    )
  }
  if (!is.null(id_var) && id_var %in% names(newdata_orig)) {
    return(newdata_orig[!duplicated(newdata_orig[[id_var]]), , drop = FALSE])
  }
  newdata_orig
}

state_distribution_anchor <- function(
  x,
  combos,
  group_cols,
  filter_cols,
  states,
  origin_time,
  baseline,
  pvarname
) {
  n_out <- max(1L, nrow(combos)) * length(states)
  anchor <- x[rep(NA_integer_, n_out), , drop = FALSE]
  row <- 1L

  for (combo_i in seq_len(max(1L, nrow(combos)))) {
    combo <- if (length(group_cols) > 0) {
      combos[combo_i, group_cols, drop = FALSE]
    } else {
      data.frame()
    }

    baseline_i <- baseline
    use_filters <- intersect(filter_cols, names(combo))
    if (length(use_filters) > 0) {
      baseline_i <- baseline_i[
        rows_match_values(baseline_i, combo[use_filters]),
        ,
        drop = FALSE
      ]
    }
    if (nrow(baseline_i) == 0) {
      stop("No baseline rows are available for an empirical origin anchor.")
    }

    state_idx <- match(as.character(baseline_i[[pvarname]]), states)
    probs <- tabulate(state_idx, nbins = length(states)) / nrow(baseline_i)
    rows <- row:(row + length(states) - 1L)

    for (nm in group_cols) {
      anchor[[nm]][rows] <- combo[[nm]][1]
    }
    anchor$state[rows] <- coerce_state_like(states, x$state)
    anchor$estimate[rows] <- probs
    if ("conf.low" %in% names(anchor)) {
      anchor$conf.low[rows] <- probs
    }
    if ("conf.high" %in% names(anchor)) {
      anchor$conf.high[rows] <- probs
    }
    if ("std.error" %in% names(anchor)) {
      anchor$std.error[rows] <- 0
    }
    if ("draw" %in% names(anchor)) {
      anchor$draw[rows] <- probs
    }
    anchor$.sop_real_time[rows] <- origin_time
    row <- row + length(states)
  }

  anchor
}

empirical_baseline_anchor <- function(x, origin_time) {
  pvarname <- attr(x, "pvarname") %||% "yprev"
  if (!"state" %in% names(x)) {
    stop("SOP output must contain a `state` column.")
  }
  if (!"estimate" %in% names(x)) {
    stop("SOP output must contain an `estimate` column.")
  }

  states <- as_state_labels(attr(x, "ylevels") %||% unique(x$state))

  if (inherits(x, "markov_avg_sops")) {
    avg_args <- attr(x, "avg_args")
    variables <- names(avg_args$variables %||% list())
    by <- avg_args$by %||% character()
    group_cols <- unique(c(variables, by))
    combos <- if (length(group_cols) > 0) {
      unique(x[, group_cols, drop = FALSE])
    } else {
      data.frame(.dummy = 1L)
    }
    if (!pvarname %in% names(attr(x, "newdata_orig"))) {
      stop("Previous-state variable `", pvarname, "` not found in `newdata_orig`.")
    }
    baseline <- baseline_rows_for_anchor(x, attr(x, "avg_args")$id_var)
    return(state_distribution_anchor(
      x = x,
      combos = combos,
      group_cols = group_cols,
      filter_cols = by,
      states = states,
      origin_time = origin_time,
      baseline = baseline,
      pvarname = pvarname
    ))
  }

  if (pvarname %in% names(x)) {
    value_cols <- sop_measure_cols(x)
    meta_cols <- setdiff(names(x), c("time", value_cols, ".sop_real_time"))
    anchor_meta <- unique(x[, meta_cols, drop = FALSE])
    anchor <- x[rep(NA_integer_, nrow(anchor_meta)), , drop = FALSE]
    for (nm in meta_cols) {
      anchor[[nm]] <- anchor_meta[[nm]]
    }
    anchor$estimate <- as.numeric(
      as.character(anchor$state) == as.character(anchor[[pvarname]])
    )
    if ("conf.low" %in% names(anchor)) {
      anchor$conf.low <- anchor$estimate
    }
    if ("conf.high" %in% names(anchor)) {
      anchor$conf.high <- anchor$estimate
    }
    if ("std.error" %in% names(anchor)) {
      anchor$std.error <- 0
    }
    if ("draw" %in% names(anchor)) {
      anchor$draw <- anchor$estimate
    }
    anchor$.sop_real_time <- origin_time
    return(anchor)
  }

  by <- attr(x, "by") %||% character()
  group_cols <- by
  combos <- if (length(group_cols) > 0) {
    unique(x[, group_cols, drop = FALSE])
  } else {
    data.frame(.dummy = 1L)
  }
  baseline <- baseline_rows_for_anchor(x, NULL)
  if (!pvarname %in% names(baseline)) {
    stop("Previous-state variable `", pvarname, "` not found in `newdata_orig`.")
  }
  state_distribution_anchor(
    x = x,
    combos = combos,
    group_cols = group_cols,
    filter_cols = by,
    states = states,
    origin_time = origin_time,
    baseline = baseline,
    pvarname = pvarname
  )
}

interpolate_numeric_column <- function(time, value, xout) {
  ok <- !is.na(time) & !is.na(value)
  time <- time[ok]
  value <- value[ok]

  if (length(time) == 0) {
    return(rep(NA_real_, length(xout)))
  }

  agg <- stats::aggregate(value ~ time, FUN = mean)
  agg <- agg[order(agg$time), , drop = FALSE]
  if (nrow(agg) == 1L) {
    out <- rep(NA_real_, length(xout))
    out[xout == agg$time] <- agg$value
    return(out)
  }

  stats::approx(
    x = agg$time,
    y = agg$value,
    xout = xout,
    rule = 1,
    ties = "ordered"
  )$y
}

normalize_interpolated_sops <- function(x) {
  value_col <- if ("estimate" %in% names(x)) {
    "estimate"
  } else if ("draw" %in% names(x)) {
    "draw"
  } else {
    return(x)
  }

  key_cols <- setdiff(names(x), c("state", sop_measure_cols(x)))
  key <- split_key(x, key_cols)
  sums <- ave(x[[value_col]], key, FUN = function(z) sum(z, na.rm = TRUE))
  scale <- !is.na(sums) & sums > 0
  x[[value_col]][scale] <- x[[value_col]][scale] / sums[scale]
  x
}

sop_draw_attr_name <- function(x) {
  for (nm in c("simulation_draws", "bootstrap_draws", "draws")) {
    if (!is.null(attr(x, nm))) {
      return(nm)
    }
  }
  NULL
}

sop_has_uncertainty_cols <- function(x) {
  any(c("conf.low", "conf.high", "std.error") %in% names(x))
}

warn_missing_interpolation_draws <- function(x) {
  if (!sop_has_uncertainty_cols(x)) {
    return(invisible(NULL))
  }

  if (identical(attr(x, "method"), "posterior")) {
    warning(
      "`interpolate_sops()` found uncertainty columns but no stored posterior ",
      "draws. It will interpolate interval endpoints. For draw-level ",
      "interpolation, rerun `sops()` or `avg_sops()` with `return_draws = TRUE`.",
      call. = FALSE
    )
  } else {
    warning(
      "`interpolate_sops()` found uncertainty columns but no stored draws. ",
      "It will interpolate interval endpoints. For draw-level interpolation, ",
      "rerun `inferences(..., return_draws = TRUE)`.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

sop_draw_value_col <- function(draws) {
  if ("estimate" %in% names(draws)) {
    return("estimate")
  }
  if ("draw" %in% names(draws)) {
    return("draw")
  }
  NULL
}

repeat_anchor_for_draws <- function(anchor, draws, value_col) {
  if (is.null(anchor) || !"draw_id" %in% names(draws)) {
    return(NULL)
  }

  draw_ids <- sort(unique(draws$draw_id))
  anchor_value <- anchor$estimate
  anchor <- anchor[, intersect(names(anchor), c(names(draws), ".sop_real_time")), drop = FALSE]
  anchor[[value_col]] <- anchor_value

  bind_rows_fill(lapply(draw_ids, function(draw_id) {
    anchor_i <- anchor
    anchor_i$draw_id <- draw_id
    anchor_i
  }))
}

interpolate_sop_draws <- function(
  draws,
  time_map,
  xout,
  anchor,
  normalize
) {
  if (!is.data.frame(draws) || !all(c("draw_id", "time", "state") %in% names(draws))) {
    return(NULL)
  }
  value_col <- sop_draw_value_col(draws)
  if (is.null(value_col)) {
    return(NULL)
  }

  work <- as.data.frame(draws)
  work$.sop_real_time <- map_sop_time_values(work$time, time_map)
  draw_anchor <- repeat_anchor_for_draws(anchor, work, value_col)
  work <- bind_rows_fill(list(draw_anchor, work))

  group_cols <- setdiff(names(work), c("time", value_col, ".sop_real_time"))
  groups <- split(seq_len(nrow(work)), split_key(work, group_cols), drop = TRUE)
  pieces <- lapply(groups, function(idx) {
    group <- work[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    out <- meta[rep(1L, length(xout)), , drop = FALSE]
    out$time <- xout
    out[[value_col]] <- interpolate_numeric_column(
      time = group$.sop_real_time,
      value = group[[value_col]],
      xout = xout
    )
    out
  })

  result <- bind_rows_fill(pieces)
  result <- result[, intersect(names(draws), names(result)), drop = FALSE]
  if (isTRUE(normalize)) {
    result <- normalize_interpolated_sops(result)
  }
  result
}

interpolated_ci_from_draws <- function(draws, result, conf_level, conf_type) {
  value_col <- sop_draw_value_col(draws)
  if (is.null(value_col) || !"draw_id" %in% names(draws)) {
    return(result)
  }

  draws_for_ci <- draws
  if (value_col != "estimate") {
    draws_for_ci$estimate <- draws_for_ci[[value_col]]
  }
  group_cols <- setdiff(names(draws), c(value_col, "draw_id"))
  ci <- compute_ci_from_draws(
    draws_df = draws_for_ci,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )

  join_cols <- intersect(group_cols, names(result))
  ci <- ci[, c(join_cols, "conf.low", "conf.high", "std.error"), drop = FALSE]

  result$.sop_order <- seq_len(nrow(result))
  result_base <- result[
    ,
    setdiff(names(result), c("conf.low", "conf.high", "std.error")),
    drop = FALSE
  ]
  out <- left_join_preserve_order(result_base, ci, by = join_cols)
  out <- out[order(out$.sop_order), , drop = FALSE]
  out$.sop_order <- NULL
  rownames(out) <- NULL
  out
}

#' Interpolate SOPs from Visit Time to Real Time
#'
#' Maps visit-scale SOP output from [sops()] or [avg_sops()] to a real elapsed
#' time scale and linearly interpolates probabilities within patient, state,
#' draw, treatment, and strata groups.
#'
#' @param x A `markov_sops` or `markov_avg_sops` object.
#' @param time_map A named numeric vector mapping visit labels to real times, or
#'   a data frame with visit and real-time columns. Data frames may use columns
#'   `visit` and `real_time`, `visit` and `time`, or their first two columns.
#' @param xout Optional numeric real-time grid. If `NULL`, uses mapped visit
#'   times, plus `origin_time` when an empirical origin anchor is requested.
#' @param origin_time Optional real time for an empirical baseline anchor.
#' @param origin Origin handling. `"empirical_baseline"` adds an anchor from the
#'   stored `newdata_orig` and previous-state variable; `"none"` does not.
#' @param normalize Logical. If `TRUE`, normalize interpolated estimates so
#'   state probabilities sum to one within each time and group.
#'
#' @return An interpolated data frame with the same core columns as `x`, where
#'   `time` is on the real-time scale.
#'
#' @details
#' The recommended workflow for irregular assessment schedules is:
#' 1. Recode assessment days to visit indices, for example day 3, 7, 14, and 28
#'    to visit `1:4`.
#' 2. Fit the Markov model with the visit index as a factor.
#' 3. Compute visit-scale SOPs with [sops()] or [avg_sops()].
#' 4. Use `interpolate_sops()` or `time_in_state(..., time_map = ...)` for
#'    real-day summaries.
#'
#' With RCT standardization from [avg_sops()], the empirical baseline anchor is
#' shared across treatment counterfactual groups.
#'
#' If `x` contains stored simulation, bootstrap, or posterior draws, the draws
#' are interpolated too. Interval columns are then recomputed from the
#' interpolated draws, with the empirical origin anchor treated as fixed.
#' If interval columns are present but stored draws are unavailable,
#' `interpolate_sops()` warns and falls back to interpolating interval endpoints.
#'
#' @examples
#' \dontrun{
#' avg <- avg_sops(
#'   fit,
#'   newdata = baseline,
#'   variables = list(tx = c(0, 1)),
#'   times = NULL
#' )
#' interpolate_sops(
#'   avg,
#'   time_map = c("1" = 3, "2" = 7, "3" = 14, "4" = 28),
#'   xout = 0:28,
#'   origin_time = 0
#' )
#' }
#'
#' @export
interpolate_sops <- function(
  x,
  time_map,
  xout = NULL,
  origin_time = NULL,
  origin = c("empirical_baseline", "none"),
  normalize = TRUE
) {
  if (!inherits(x, c("markov_sops", "markov_avg_sops"))) {
    stop("`x` must be a `markov_sops` or `markov_avg_sops` object.")
  }
  if (!all(c("time", "state") %in% names(x))) {
    stop("`x` must contain `time` and `state` columns.")
  }

  origin <- match.arg(origin)
  time_map <- standardize_time_map(time_map)
  real_time <- map_sop_time_values(x$time, time_map)
  mapped_range <- range(time_map$real_time)

  use_origin <- !is.null(origin_time) && origin == "empirical_baseline"
  if (!is.null(origin_time)) {
    if (!is.numeric(origin_time) || length(origin_time) != 1 || is.na(origin_time)) {
      stop("`origin_time` must be a single finite numeric value.")
    }
    if (!is.finite(origin_time)) {
      stop("`origin_time` must be a single finite numeric value.")
    }
  }

  lower <- if (use_origin) origin_time else mapped_range[1]
  upper <- mapped_range[2]
  if (lower > upper) {
    stop("`origin_time` must not be greater than the largest mapped time.")
  }

  if (is.null(xout)) {
    xout <- sort(unique(c(if (use_origin) origin_time else NULL, real_time)))
  }
  xout <- validate_sop_xout(xout, lower, upper)

  work <- x
  work$.sop_real_time <- real_time

  anchor <- NULL
  if (use_origin) {
    anchor <- empirical_baseline_anchor(work, origin_time)
    work <- bind_rows_fill(list(anchor, work))
  } else {
    work <- as.data.frame(work)
  }

  value_cols <- sop_measure_cols(work)
  if (length(value_cols) == 0) {
    stop("`x` must contain an `estimate` or `draw` column to interpolate.")
  }
  group_cols <- setdiff(names(work), c("time", value_cols, ".sop_real_time"))
  groups <- split(seq_len(nrow(work)), split_key(work, group_cols), drop = TRUE)

  pieces <- lapply(groups, function(idx) {
    group <- work[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    out <- meta[rep(1L, length(xout)), , drop = FALSE]
    out$time <- xout
    for (nm in value_cols) {
      out[[nm]] <- interpolate_numeric_column(
        time = group$.sop_real_time,
        value = group[[nm]],
        xout = xout
      )
    }
    out
  })

  result <- bind_rows_fill(pieces)
  result <- result[, intersect(names(x), names(result)), drop = FALSE]
  if (isTRUE(normalize)) {
    result <- normalize_interpolated_sops(result)
  }

  draw_attr <- sop_draw_attr_name(x)
  if (is.null(draw_attr)) {
    warn_missing_interpolation_draws(x)
  }
  interpolated_draws <- if (!is.null(draw_attr)) {
    interpolate_sop_draws(
      draws = attr(x, draw_attr),
      time_map = time_map,
      xout = xout,
      anchor = anchor,
      normalize = normalize
    )
  } else {
    NULL
  }
  if (!is.null(interpolated_draws)) {
    result <- interpolated_ci_from_draws(
      draws = interpolated_draws,
      result = result,
      conf_level = attr(x, "conf_level") %||% 0.95,
      conf_type = attr(x, "conf_type") %||% "perc"
    )
  }

  for (a in names(attributes(x))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(x, a)
    }
  }
  if (!is.null(interpolated_draws)) {
    attr(result, draw_attr) <- interpolated_draws
  }
  attr(result, "time_map") <- time_map
  attr(result, "origin_time") <- if (use_origin) origin_time else NULL
  class(result) <- c("markov_interpolated_sops", class(x))
  result
}

trapezoid_auc <- function(time, value) {
  ok <- !is.na(time) & !is.na(value)
  time <- time[ok]
  value <- value[ok]
  if (length(time) < 2) {
    return(0)
  }
  agg <- stats::aggregate(value ~ time, FUN = sum)
  agg <- agg[order(agg$time), , drop = FALSE]
  if (nrow(agg) < 2) {
    return(0)
  }
  dt <- diff(agg$time)
  sum(dt * (head(agg$value, -1L) + tail(agg$value, -1L)) / 2)
}

time_in_state_tidy <- function(x, target_states, real_time = FALSE) {
  value_col <- if ("estimate" %in% names(x)) {
    "estimate"
  } else if ("draw" %in% names(x)) {
    "draw"
  } else {
    stop("Tidy SOP input must contain an `estimate` or `draw` column.")
  }

  keep <- as.character(x$state) %in% as.character(target_states)
  if (!any(keep)) {
    stop("Target states not found in SOP data frame.")
  }
  x <- x[keep, , drop = FALSE]
  measure_cols <- sop_measure_cols(x)
  group_cols <- setdiff(names(x), c("time", "state", measure_cols))

  agg_formula <- stats::as.formula(
    paste(value_col, "~", paste(c(group_cols, "time"), collapse = " + "))
  )
  by_time <- stats::aggregate(agg_formula, data = x, FUN = sum, na.rm = TRUE)

  if (!real_time) {
    if (length(group_cols) == 0) {
      out <- data.frame(total_time = sum(by_time[[value_col]], na.rm = TRUE))
      return(out)
    }
    total_formula <- stats::as.formula(
      paste(value_col, "~", paste(group_cols, collapse = " + "))
    )
    out <- stats::aggregate(total_formula, data = by_time, FUN = sum, na.rm = TRUE)
    names(out)[names(out) == value_col] <- "total_time"
    return(out)
  }

  groups <- split(seq_len(nrow(by_time)), split_key(by_time, group_cols), drop = TRUE)
  pieces <- lapply(groups, function(idx) {
    group <- by_time[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    meta$total_time <- trapezoid_auc(group$time, group[[value_col]])
    meta
  })
  bind_rows_fill(pieces)
}

time_in_state_bootstrap_df <- function(sops, target_states, real_time = FALSE) {
  if (!all(c("boot_id", "time", "tx") %in% colnames(sops))) {
    stop("Input data frame must contain 'boot_id', 'time', and 'tx' columns.")
  }

  state_cols_available <- grep("^state_", colnames(sops), value = TRUE)
  if (length(state_cols_available) == 0) {
    stop("Input data frame does not contain any columns starting with 'state_'.")
  }

  target_suffixes <- as.character(target_states)
  target_cols <- paste0("state_", target_suffixes)
  missing_cols <- setdiff(target_cols, colnames(sops))
  if (length(missing_cols) > 0) {
    stop(paste(
      "The following target state columns were not found in the bootstrap output:",
      paste(missing_cols, collapse = ", ")
    ))
  }

  prob_target <- if (length(target_cols) == 1) {
    sops[[target_cols]]
  } else {
    rowSums(sops[, target_cols, drop = FALSE])
  }

  agg_df <- sops[, c("boot_id", "tx", "time")]
  agg_df$prob <- prob_target

  if (real_time) {
    groups <- split(seq_len(nrow(agg_df)), split_key(agg_df, c("boot_id", "tx")))
    res_grouped <- bind_rows_fill(lapply(groups, function(idx) {
      group <- agg_df[idx, , drop = FALSE]
      data.frame(
        boot_id = group$boot_id[1],
        tx = group$tx[1],
        prob = trapezoid_auc(group$time, group$prob)
      )
    }))
  } else {
    res_grouped <- stats::aggregate(
      prob ~ boot_id + tx,
      data = agg_df,
      FUN = sum
    )
  }

  res_tx <- res_grouped[res_grouped$tx == 1, c("boot_id", "prob")]
  res_ctrl <- res_grouped[res_grouped$tx == 0, c("boot_id", "prob")]
  res_wide <- merge(
    res_tx,
    res_ctrl,
    by = "boot_id",
    suffixes = c("_tx", "_ctrl")
  )

  colnames(res_wide)[colnames(res_wide) == "prob_tx"] <- "SOP_tx"
  colnames(res_wide)[colnames(res_wide) == "prob_ctrl"] <- "SOP_ctrl"
  res_wide$delta <- res_wide$SOP_tx - res_wide$SOP_ctrl
  res_wide
}

#' Compute Total Time in Target State(s)
#'
#' Calculates the expected total time spent in specified target state(s).
#' By default this sums state occupancy probabilities over visit-scale unit
#' steps. When `time_map` or `origin_time` is supplied, it maps/interpolates to
#' real time and uses trapezoidal AUC. Already interpolated SOP output from
#' [interpolate_sops()] is integrated on its current real-time grid.
#' It automatically adapts to the input format:
#' \itemize{
#'   \item **Patient-Level:** If input is an array from \code{soprob_markov}, it returns time-in-state for each patient.
#'   \item **Bootstrap-Level:** If input is a data frame from \code{\link{bootstrap_standardized_sops}}, it returns mean time-in-state per treatment group and the difference for each bootstrap sample.
#' }
#'
#' @param sops Input object.
#'   \itemize{
#'     \item **Array:** `[Patients x Time x States]` (Frequentist) or `[Draws x Patients x Time x States]` (Bayes).
#'     \item **Data Frame:** Output from `bootstrap_standardized_sops` containing columns `boot_id`, `time`, `tx`, and `state_*`.
#'   }
#' @param target_states Vector of target state(s) to include in the time calculation.
#'   Can be integer indices or character names (e.g., \code{1} or \code{c("Home", "Rehab")}).
#'   For bootstrap outputs, these must match the suffix of the `state_` columns (e.g., if column is `state_1`, use `1`).
#' @param time_map Optional named numeric vector or data frame mapping visit
#'   labels to real elapsed times. Supplying this switches tidy SOP and array
#'   inputs to trapezoidal real-time AUC.
#' @param origin_time Optional real time for an empirical baseline anchor. For
#'   tidy outputs from [sops()] or [avg_sops()], this is passed to
#'   [interpolate_sops()].
#' @param xout Optional numeric real-time grid for tidy SOP outputs when
#'   `time_map` is supplied. This controls the interpolation grid used for AUC,
#'   for example `xout = 1:28` with `origin_time = 0` uses day 0 as an anchor
#'   but starts the AUC at day 1.
#' @param origin Origin handling for tidy SOP outputs when `origin_time` is
#'   supplied. See [interpolate_sops()].
#'
#' @return
#' \itemize{
#'   \item **Input = Array (Frequentist):** A named numeric vector of length \code{n_patients}.
#'   \item **Input = Array (Bayesian):** A matrix of dimension \code{[n_draws x n_patients]}.
#'   \item **Input = Bootstrap DF:** A tibble with columns:
#'     \itemize{
#'       \item `boot_id`: Bootstrap iteration
#'       \item `SOP_tx`: Mean time in target state (Treatment)
#'       \item `SOP_ctrl`: Mean time in target state (Control)
#'       \item `delta`: Difference (Treatment - Control)
#'     }
#' }
#'
#' @examples
#' \dontrun{
#' # --- Scenario 1: Patient-Level Estimates ---
#' sops_arr <- soprob_markov(model, data, times = 1:30, ylevels = 1:6)
#' days_home <- time_in_state(sops_arr, target_states = 1)
#'
#' # Real-time AUC after fitting on factor visit indices
#' real_days_home <- time_in_state(
#'   avg,
#'   target_states = 1,
#'   time_map = c("1" = 3, "2" = 7, "3" = 14, "4" = 28),
#'   origin_time = 0,
#'   xout = 1:28
#' )
#'
#' # --- Scenario 2: Bootstrap Inference ---
#' bs_res <- bootstrap_standardized_sops(model, data, n_boot=100)
#' # Get distribution of treatment effects
#' bs_effects <- time_in_state(bs_res, target_states = 1)
#' }
#'
#' @keywords time-in-state auc bootstrap
#' @export
time_in_state <- function(
  sops,
  target_states = 1,
  time_map = NULL,
  origin_time = NULL,
  xout = NULL,
  origin = c("empirical_baseline", "none")
) {
  origin <- match.arg(origin)
  use_real_time <- !is.null(time_map) || !is.null(origin_time)

  if (inherits(sops, "markov_interpolated_sops")) {
    if (!is.null(time_map) || !is.null(origin_time)) {
      stop(
        "`markov_interpolated_sops` is already on the real-time scale; ",
        "omit `time_map` and `origin_time`."
      )
    }
    if (!is.null(xout)) {
      xout <- validate_sop_xout(xout, min(sops$time), max(sops$time))
      missing_xout <- setdiff(xout, unique(sops$time))
      if (length(missing_xout) > 0) {
        stop(
          "`xout` for already interpolated SOPs must be contained in ",
          "the existing `time` values."
        )
      }
      sops <- sops[sops$time %in% xout, , drop = FALSE]
    }
    return(time_in_state_tidy(sops, target_states, real_time = TRUE))
  }

  if (inherits(sops, c("markov_sops", "markov_avg_sops"))) {
    if (!is.null(xout) && !use_real_time) {
      stop("`xout` requires `time_map` or an already interpolated SOP object.")
    }
    if (use_real_time) {
      if (is.null(time_map)) {
        stop("`time_map` must be supplied for real-time AUC.")
      }
      sops <- interpolate_sops(
        sops,
        time_map = time_map,
        xout = xout,
        origin_time = origin_time,
        origin = origin
      )
      return(time_in_state_tidy(sops, target_states, real_time = TRUE))
    }
    return(time_in_state_tidy(sops, target_states, real_time = FALSE))
  }

  # =========================================================================
  # BRANCH 1: Bootstrap Data Frame Input
  # =========================================================================
  if (inherits(sops, "data.frame")) {
    if (!is.null(xout)) {
      stop("`xout` is only supported for tidy SOP outputs.")
    }
    if (use_real_time) {
      if (is.null(time_map)) {
        stop("`time_map` must be supplied for real-time AUC.")
      }
      time_map <- standardize_time_map(time_map)
      sops$time <- map_sop_time_values(sops$time, time_map)
    }
    return(time_in_state_bootstrap_df(
      sops,
      target_states = target_states,
      real_time = use_real_time
    ))
  }

  # =========================================================================
  # BRANCH 2: Array Input (Original Logic)
  # =========================================================================

  if (!is.null(xout)) {
    stop("`xout` is only supported for tidy SOP outputs.")
  }

  # --- 1. Detect Dimensions ---
  dims <- dim(sops)
  if (is.null(dims)) {
    stop("Input 'sops' must be an array or data frame.")
  }

  ndim <- length(dims)
  dnames <- dimnames(sops)

  if (ndim == 3) {
    # Frequentist: [Patients, Time, States]
    state_dim <- 3
  } else if (ndim == 4) {
    # Bayesian: [Draws, Patients, Time, States]
    state_dim <- 4
  } else {
    stop("Input array must be 3D or 4D.")
  }

  # --- 2. Resolve Target States ---
  state_names <- dnames[[state_dim]]
  if (is.null(state_names)) {
    state_names <- as.character(1:dims[state_dim])
  }

  target_chars <- as.character(target_states)
  missing_states <- setdiff(target_chars, state_names)
  if (length(missing_states) > 0) {
    stop(paste(
      "Target states not found in input array:",
      paste(missing_states, collapse = ", ")
    ))
  }

  # --- 3. Filter & Aggregate States ---
  subset_args <- rep(list(TRUE), ndim)
  subset_args[[state_dim]] <- match(target_chars, state_names)

  sops_slice <- do.call("[", c(list(sops), subset_args, drop = FALSE))
  prob_in_target <- apply(sops_slice, (1:ndim)[-state_dim], sum)

  # --- 4. Integrate Over Time ---
  if (use_real_time) {
    if (is.null(time_map)) {
      stop("`time_map` must be supplied for real-time AUC.")
    }
    if (!is.null(origin_time) && origin == "empirical_baseline") {
      stop("Empirical `origin_time` anchoring is only supported for tidy SOP outputs.")
    }
    time_names <- dnames[[state_dim - 1L]]
    if (is.null(time_names)) {
      time_names <- seq_len(dims[state_dim - 1L])
    }
    time_map <- standardize_time_map(time_map)
    real_time <- map_sop_time_values(time_names, time_map)

    if (ndim == 3) {
      return(apply(prob_in_target, 1, function(z) trapezoid_auc(real_time, z)))
    }
    return(apply(prob_in_target, c(1, 2), function(z) trapezoid_auc(real_time, z)))
  }

  if (ndim == 3) {
    # Freq Input -> [Pat, Time] -> Apply over rows
    res <- rowSums(prob_in_target)
  } else {
    # Bayes Input -> [Draws, Pat, Time] -> Apply over Draws & Pat
    res <- apply(prob_in_target, c(1, 2), sum)
  }

  return(res)
}
#' Calculate Individual State Occupation Probabilities
#'
#' Computes individual-level state occupation probabilities (SOPs) for each row
#' in a dataset. Optionally aggregates results within strata defined by grouping
#' variables.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, or `blrm`). For
#'   `vglm` models, the family must be `cumulative(reverse = TRUE, ...)`.
#' @param newdata Optional. A data frame of new data for prediction. If NULL,
#'   uses the data used to fit the model.
#' @param times Visit-scale time points to estimate. For numeric time
#'   variables this is usually a numeric vector. For factor-valued visit
#'   indices, values are matched to fitted visit levels; if `NULL`, all
#'   fitted visit levels are used.
#' @param ylevels A vector of state levels. If NULL, attempts to infer from model.
#' @param absorb The absorbing state.
#' @param tvarname Name of the time variable in the model.
#' @param pvarname Name of the previous state variable in the model.
#' @param p2varname Optional second previous-state variable. `NULL` uses a
#'   first-order Markov recursion; a non-`NULL` column name uses a second-order
#'   recursion.
#' @param gap Name of the time gap variable (if used).
#' @param t_covs Optional time-varying covariate lookup table for explicit
#' basis columns. Inline terms such as `rms::rcs(time, 4)` can be used without
#' supplying `t_covs`.
#' @param by Optional character vector of variable names to stratify by. When
#'   provided, the function aggregates (averages) SOPs within each stratum
#'   defined by combinations of these variables. E.g., `by = "ecog"` aggregates
#'   within ECOG levels; `by = c("ecog", "age_group")` aggregates within each
#'   combination of ECOG and age group. NOTE: This is simple aggregation within
#'   observed strata, NOT G-computation standardization (use `avg_sops()` for that).
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs.
#' @param id_var Character ID column used when `include_re = TRUE`. For `blrm`,
#'   `NULL` is inferred from `model$clusterInfo$name` when available, otherwise
#'   `"id"`.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible random draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior SOP draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior SOP draws for
#'   extraction with `get_draws()`. Large draw objects are guarded to protect
#'   memory.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class `markov_sops` containing:
#'   \item{rowid}{Row identifier from newdata (omitted if `by` is specified)}
#'   \item{time}{Time point}
#'   \item{state}{State name}
#'   \item{estimate}{Probability of being in the state (individual-level or stratum average)}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     uncertainty summaries computed directly from SOP draws}
#'   \item{(by variables)}{Stratification variables if `by` is specified}
#'   Plus all columns from `newdata` (for individual-level results).
#'
#' @details
#' This function wraps `soprob_markov()` and converts its array output to a tidy
#' data frame. The output contains one row per patient-time-state combination
#' (or stratum-time-state if `by` is used).
#'
#' For `rmsb::blrm()` models, SOPs are computed on sampled posterior draws and
#' then summarized. The point estimate is the requested posterior summary of the
#' SOP draws, not a plug-in calculation from summarized model parameters.
#' State-wise medians and interval bounds are not constrained to sum to one
#' across states; use `posterior_summary = "mean"` when the displayed estimates
#' themselves need to preserve total probability. Draw-level probabilities
#' stored with `return_draws = TRUE` remain normalized within each draw.
#'
#' **Model Requirements:**
#'
#' For `vglm`/`vgam` models, only the `cumulative` family with `reverse = TRUE` is
#' supported. This is because the package's Markov simulation logic expects
#' higher-numbered states to represent worse outcomes, and uses reverse cumulative
#' probabilities to model the probability of being in state k or worse.
#'
#' **Stratification:**
#'
#' When `by` is specified, the function:
#' 1. Computes individual-level SOPs for all patients
#' 2. Groups results by time, state, and the specified variables
#' 3. Averages the SOPs within each group
#'
#' This aggregation preserves heterogeneity across observed strata without
#' creating counterfactual scenarios (unlike G-computation in `avg_sops()`).
#'
#' For computing marginal/standardized SOPs (G-computation), use `avg_sops()`
#' instead, which creates counterfactual datasets and averages over individuals.
#'
#' @seealso [avg_sops()] for marginal SOPs, [soprob_markov()] for the underlying
#'   computation.
#'
#' @examples
#' \dontrun{
#' # Individual-level SOPs
#' sops_ind <- sops(fit, newdata = baseline_data, times = 1:30, ylevels = 1:6)
#'
#' # Stratified by ECOG performance status
#' sops_ecog <- sops(fit, newdata = baseline_data, times = 1:30, ylevels = 1:6,
#'                   by = "ecog")
#'
#' # Stratified by multiple variables
#' sops_strat <- sops(fit, newdata = baseline_data, times = 1:30, ylevels = 1:6,
#'                    by = c("ecog", "age_group"))
#'
#' # Compatible with inferences()
#' sops_strat |> inferences(method = "simulation", n_sim = 1000)
#' }
#'
#' @importFrom stats model.frame
#' @export
sops <- function(
  model,
  newdata = NULL,
  times = NULL,
  ylevels = NULL,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  p2varname = NULL,
  gap = NULL,
  t_covs = NULL,
  by = NULL,
  include_re = FALSE,
  id_var = NULL,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  # --- 1. Setup & Defaults ---
  # Validate model compatibility
  validate_markov_model(model)

  if (is.null(newdata)) {
    # newdata <- model$x
    # if (is.null(newdata)) {
    #   newdata <- tryCatch(model.frame(model), error = function(e) NULL)
    # }
    # if (is.null(newdata)) {
    stop("Please provide `newdata` (should be baseline data).")
    # }
  }

  # Ensure rowid exists for tracking
  if (!"rowid" %in% names(newdata)) {
    newdata$rowid <- seq_len(nrow(newdata))
  }

  time_res <- resolve_sop_times(
    model,
    newdata,
    times,
    tvarname,
    t_covs = t_covs,
    default = "unique"
  )
  times <- time_res$times
  validate_factor_gap(gap, t_covs, time_res$time_info)

  if (is.null(ylevels)) {
    # Try to infer from model
    if (inherits(model, "vglm")) {
      # VGAM stores response levels
      ylevels <- model@extra$colnames.y
      if (is.null(ylevels)) stop("`ylevels` cannot be NULL")
    } else if (inherits(model, "robcov_vglm")) {
      # robcov_vglm stores extra slot in list
      ylevels <- model$extra$colnames.y
      if (is.null(ylevels)) stop("`ylevels` cannot be NULL")
    } else if (inherits(model, "blrm")) {
      ylevels <- model$ylevels
      if (is.null(ylevels)) stop("`ylevels` cannot be NULL")
    } else if (inherits(model, "orm")) {
      # rms orm stores levels
      ylevels <- model$yunique
      if (is.null(ylevels)) stop("`ylevels` cannot be NULL")
    } else {
      stop("`ylevels` cannot be NULL")
    }
  }

  if (inherits(model, "blrm")) {
    posterior_summary <- match.arg(posterior_summary)
    return(sops_blrm(
      model = model,
      newdata = newdata,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      by = by,
      include_re = include_re,
      id_var = id_var,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws
    ))
  }

  # --- 2. Compute SOPs (Vectorized) ---
  sops_array <- soprob_markov(
    object = model,
    data = newdata,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs
  )

  # --- 3. Validate 'by' Parameter ---
  if (!is.null(by)) {
    if (!is.character(by)) {
      stop("'by' must be a character vector of variable names.")
    }
    missing_vars <- setdiff(by, names(newdata))
    if (length(missing_vars) > 0) {
      stop(
        "Variables specified in 'by' not found in newdata: ",
        paste(missing_vars, collapse = ", ")
      )
    }
  }

  # --- 4. Tidy the Output ---
  # Array dims: [n_pat, n_times, n_states]
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]

  # Flatten array (Reminder for myself: R is column-major, iterates dim1 fastest)
  # as.vector order: [1,1,1], [2,1,1]... [N,1,1], [1,2,1]... [N,T,S]
  probs_flat <- as.vector(sops_array)

  # Construct indices matching as.vector order
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)

  # Build result by repeating newdata rows
  result <- newdata[idx_pat, , drop = FALSE]
  result$time <- idx_time
  result$state <- idx_state
  result$estimate <- probs_flat
  rownames(result) <- NULL

  # --- 5. Apply Stratified Aggregation if 'by' is Specified ---
  if (!is.null(by)) {
    # Group by time, state, and stratification variables, then average
    group_cols <- unique(c("time", "state", by))

    # Aggregate using base R for minimal dependencies
    agg_formula <- stats::as.formula(
      paste("estimate ~", paste(group_cols, collapse = " + "))
    )
    result <- stats::aggregate(
      agg_formula,
      data = result,
      FUN = mean,
      na.rm = TRUE
    )

    # Note: After aggregation, we lose individual patient information (rowid)
    # This is expected for stratified results
  }

  # --- 6. Store Attributes for Downstream Use ---
  attr(result, "model") <- model
  attr(result, "tvarname") <- tvarname
  attr(result, "pvarname") <- pvarname
  attr(result, "p2varname") <- p2varname
  attr(result, "ylevels") <- ylevels
  attr(result, "absorb") <- absorb
  attr(result, "gap") <- gap
  attr(result, "t_covs") <- t_covs
  attr(result, "by") <- by
  attr(result, "call_args") <- list(
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs,
    by = by
  )
  # Store original newdata for inferences() (needed for simulation-based inference)
  attr(result, "newdata_orig") <- newdata

  class(result) <- c("markov_sops", class(result))
  return(result)
}

sops_blrm <- function(
  model,
  newdata,
  times,
  ylevels,
  absorb,
  tvarname,
  pvarname,
  p2varname,
  gap,
  t_covs,
  by,
  include_re,
  id_var,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws
) {
  if (!is.null(by)) {
    if (!is.character(by)) {
      stop("'by' must be a character vector of variable names.")
    }
    missing_vars <- setdiff(by, names(newdata))
    if (length(missing_vars) > 0) {
      stop(
        "Variables specified in 'by' not found in newdata: ",
        paste(missing_vars, collapse = ", ")
      )
    }
  }

  draw_indices <- select_posterior_draws(model, n_draws, seed)
  n_pat <- nrow(newdata)
  n_times <- length(times)
  n_states <- length(ylevels)
  n_cells <- n_pat * n_times * n_states
  n_requested <- length(draw_indices)
  # `sops()` keeps individual-level draw summaries in memory. The default
  # guard is 50 million numeric cells, about 400 MB before data-frame overhead.
  max_draw_cells <- getOption("markov.misc.max_sops_draw_cells", 5e7)

  if (n_cells * n_requested > max_draw_cells) {
    stop(
      "`sops()` for `blrm` would require ",
      format(n_cells * n_requested, scientific = FALSE),
      " posterior draw cells, above the current limit of ",
      format(max_draw_cells, scientific = FALSE),
      ". This limit protects memory because individual-level `sops()` stores ",
      "draw x patient x time x state values. Lower `n_draws`, use ",
      "`avg_sops()` for marginal summaries, or increase option ",
      "`markov.misc.max_sops_draw_cells` if this allocation is intentional."
    )
  }

  draw_values <- matrix(NA_real_, nrow = n_requested, ncol = n_cells)
  rownames(draw_values) <- draw_indices
  chunks <- split_draw_indices(draw_indices)
  gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  row_cursor <- 1L
  for (chunk in chunks) {
    arr <- soprob_markov(
      object = model,
      data = newdata,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      id_var = id_var,
      n_draws = NULL,
      .draw_indices = chunk,
      .gamma_draws = subset_cached_draw_matrix(gamma_draws, draw_indices, chunk)
    )
    for (i in seq_along(chunk)) {
      draw_values[row_cursor, ] <- as.vector(arr[i, , , ])
      row_cursor <- row_cursor + 1L
    }
  }

  if (is.null(by)) {
    stats <- summarize_posterior_draw_matrix(
      draw_values,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    )
    idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
    result <- newdata[idx_pat, , drop = FALSE]
    result$time <- rep(rep(times, each = n_pat), times = n_states)
    result$state <- rep(ylevels, each = n_pat * n_times)
    result$estimate <- stats$estimate
    result$conf.low <- stats$conf.low
    result$conf.high <- stats$conf.high
    result$std.error <- stats$std.error
    rownames(result) <- NULL

    draws_df <- if (return_draws) {
      sops_draw_matrix_to_df(draw_values, result, draw_indices)
    } else {
      NULL
    }
  } else {
    draw_list <- vector("list", n_requested)
    for (i in seq_len(n_requested)) {
      arr_i <- array(draw_values[i, ], dim = c(n_pat, n_times, n_states))
      draw_i <- array_to_df_individual(arr_i, times, ylevels, newdata, by = by)
      draw_i$draw_id <- draw_indices[i]
      draw_list[[i]] <- draw_i
    }
    draws_df <- bind_rows_fill(draw_list)
    group_cols <- unique(c("time", "state", by))
    result <- summarize_posterior_draws_df(
      draws_df,
      group_cols = group_cols,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    )
    if (!return_draws) {
      draws_df <- NULL
    }
  }

  attr(result, "model") <- model
  attr(result, "tvarname") <- tvarname
  attr(result, "pvarname") <- pvarname
  attr(result, "p2varname") <- p2varname
  attr(result, "ylevels") <- ylevels
  attr(result, "absorb") <- absorb
  attr(result, "gap") <- gap
  attr(result, "t_covs") <- t_covs
  attr(result, "by") <- by
  attr(result, "call_args") <- list(
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs,
    by = by,
    include_re = include_re,
    id_var = id_var,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level
  )
  attr(result, "newdata_orig") <- newdata
  attr(result, "method") <- "posterior"
  attr(result, "engine") <- "posterior"
  attr(result, "n_draws") <- n_requested
  attr(result, "draw_ids") <- draw_indices
  attr(result, "conf_level") <- conf_level
  attr(result, "posterior_summary") <- posterior_summary
  if (!is.null(draws_df)) {
    attr(result, "draws") <- draws_df
  }

  class(result) <- c("markov_sops", class(result))
  result
}

split_draw_indices <- function(draw_indices, chunk_size = NULL) {
  if (is.null(chunk_size)) {
    chunk_size <- getOption("markov.misc.blrm_chunk_size", 10L)
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  split(draw_indices, ceiling(seq_along(draw_indices) / chunk_size))
}

summarize_posterior_draw_matrix <- function(draw_values, posterior_summary, conf_level) {
  alpha <- 1 - conf_level
  estimate <- switch(
    posterior_summary,
    mean = colMeans(draw_values, na.rm = TRUE),
    median = apply(draw_values, 2, stats::median, na.rm = TRUE)
  )
  conf.low <- apply(draw_values, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  conf.high <- apply(draw_values, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  std.error <- apply(draw_values, 2, stats::sd, na.rm = TRUE)
  data.frame(
    estimate = as.numeric(estimate),
    conf.low = as.numeric(conf.low),
    conf.high = as.numeric(conf.high),
    std.error = as.numeric(std.error)
  )
}

summarize_posterior_draws_df <- function(draws_df, group_cols, posterior_summary, conf_level) {
  alpha <- 1 - conf_level
  split_key <- interaction(draws_df[, group_cols, drop = FALSE], drop = TRUE)
  groups <- split(seq_len(nrow(draws_df)), split_key)
  out <- draws_df[vapply(groups, `[`, integer(1), 1L), group_cols, drop = FALSE]
  values <- lapply(groups, function(idx) draws_df$estimate[idx])
  out$estimate <- vapply(values, function(x) {
    switch(
      posterior_summary,
      mean = mean(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE)
    )
  }, numeric(1))
  out$conf.low <- vapply(values, stats::quantile, numeric(1), probs = alpha / 2, na.rm = TRUE)
  out$conf.high <- vapply(values, stats::quantile, numeric(1), probs = 1 - alpha / 2, na.rm = TRUE)
  out$std.error <- vapply(values, stats::sd, numeric(1), na.rm = TRUE)
  rownames(out) <- NULL
  out
}

sops_draw_matrix_to_df <- function(draw_values, result, draw_indices) {
  meta_cols <- setdiff(names(result), c("estimate", "conf.low", "conf.high", "std.error"))
  n_draws <- nrow(draw_values)
  n_cells <- ncol(draw_values)
  out <- result[rep(seq_len(n_cells), times = n_draws), meta_cols, drop = FALSE]
  out$draw_id <- rep(draw_indices, each = n_cells)
  out$estimate <- as.vector(t(draw_values))
  rownames(out) <- NULL
  out
}


#' Calculate Averaged State Occupation Probabilities (Marginal Effects)
#'
#' Computes standardized (marginal) state occupation probabilities using
#' G-computation. Creates counterfactual cohorts by setting all individuals
#' to each level of the treatment variable and averaging over the covariate
#' distribution.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, or `blrm`). For
#'   `vglm` models, the family **must** be `cumulative(reverse = TRUE, ...)`.
#' @param newdata Data frame for prediction. For **simulation inference**, pass
#'   baseline data only (one row per patient). For **bootstrap inference**, you
#'   must pass the full longitudinal dataset (all time points) since the model
#'   needs to be refit on bootstrap samples. If NULL, extracts from model.
#' @param variables A named list specifying the variable(s) to standardize over.
#'   E.g., `list(tx = c(0, 1))` creates counterfactual datasets for treatment
#'   and control.
#' @param by Optional character vector of additional variables to group by,
#' after standardization.
#' @param times Visit-scale time points. Numeric time variables use numeric
#'   values; factor-valued visit indices use fitted visit levels. If `NULL`,
#'   factor time uses all fitted visit levels and numeric time is inferred from
#'   `newdata`.
#' @param id_var Name of the patient ID variable. Required for bootstrap
#'   inference and for `blrm` random-effect prediction. If `NULL`, defaults to
#'   `"id"`; for `blrm` models with `include_re = TRUE`, it is first inferred
#'   from `model$clusterInfo$name` when available.
#' @param p2varname Optional second previous-state variable. `NULL` uses a
#'   first-order Markov recursion; a non-`NULL` column name uses a second-order
#'   recursion.
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible random draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior SOP draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior SOP draws for
#'   extraction with `get_draws()`.
#' @param ... Additional arguments passed to `sops()` (e.g., `ylevels`, `absorb`,
#'   `tvarname`, `pvarname`, `t_covs`). `t_covs` is only needed for explicit
#'   precomputed time-basis columns; inline terms such as `rms::rcs(time, 4)`
#'   can be used without it.
#'
#' @return A data frame of class `markov_avg_sops` with columns:
#'   \item{time}{Time point}
#'   \item{state}{State level}
#'   \item{(variables)}{Value of standardization variable (e.g., tx)}
#'   \item{estimate}{Average probability across individuals}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     uncertainty summaries computed directly from marginalized SOP draws}
#'
#' @details
#' This function implements G-computation (standardization) for Markov SOPs:
#'
#' 1. **Counterfactual Creation**: For each value in `variables`, creates a
#'    copy of `newdata` with that variable set to the specified value.
#'
#' 2. **Individual Prediction**: Calls `sops()` to compute individual-level
#'    SOPs for each patient under each counterfactual scenario.
#'
#' 3. **Marginalization**: Averages individual SOPs across patients within
#'    each time-state-treatment combination, yielding population-average
#'    (marginal) probabilities.
#'
#' The result represents the expected SOP if the entire population received
#' treatment vs. control, averaged over the observed covariate distribution.
#' For `rmsb::blrm()` models, patient-level posterior arrays are chunked and
#' marginalized by draw before summarizing to reduce memory use.
#' State-wise medians and interval bounds are not constrained to sum to one
#' across states; use `posterior_summary = "mean"` when the displayed estimates
#' themselves need to preserve total probability. Draw-level probabilities
#' stored with `return_draws = TRUE` remain normalized within each draw.
#'
#' @seealso [sops()] for individual-level SOPs, [inferences()] for bootstrap
#'   uncertainty, [standardize_sops()] for the underlying implementation.
#'
#' @examples
#' \dontrun{
#' # For simulation inference: use baseline data (one row per patient)
#' baseline_data <- data |> filter(time == 1)
#' result <- avg_sops(
#'   model = fit,
#'   newdata = baseline_data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6
#' ) |> inferences(method = "simulation", n_sim = 500)
#'
#' # For bootstrap inference: must use full data (all time points)
#' result_boot <- avg_sops(
#'   model = fit,
#'   newdata = data,  # Full longitudinal data
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6
#' ) |> inferences(method = "bootstrap", n_sim = 500)
#' }
#'
#' @export
avg_sops <- function(
  model,
  newdata = NULL,
  variables = NULL,
  by = NULL,
  times = NULL,
  id_var = NULL,
  p2varname = NULL,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  # --- 1. Input Validation ---
  # Validate model compatibility
  validate_markov_model(model)

  if (is.null(variables)) {
    stop(
      "`variables` is required for G-computation. ",
      "Specify the treatment variable, e.g., `variables = list(tx = c(0, 1))`."
    )
  }

  if (is.null(newdata)) {
    # this creates too many issues
    # newdata <- model$x
    # if (is.null(newdata)) {
    #   newdata <- tryCatch(model.frame(model), error = function(e) NULL)
    # }
    # if (is.null(newdata)) {
    stop("Provide newdata or ensure model stores data (x = TRUE).")
    # }
  }

  if (inherits(model, "blrm") && isTRUE(include_re)) {
    id_var <- resolve_blrm_id_var(model, newdata, id_var)
  } else {
    id_var <- id_var %||% "id"
  }

  # Validate id_var exists
  if (!id_var %in% names(newdata)) {
    stop("ID variable '", id_var, "' not found in data.")
  }

  # Validate variables exist in data
  missing_vars <- setdiff(names(variables), names(newdata))
  if (length(missing_vars) > 0) {
    stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }

  # --- 2. Extract Baseline Data (One Row Per Patient) ---
  # For standardization, we need unique patient baseline covariates
  # This matches standardize_sops() behavior
  baseline_data <- newdata[!duplicated(newdata[[id_var]]), ]

  # Could consider not checking ID, but taking tvarname as an argument and then
  # filtering for min.
  # newdata[newdata[, tvarname] == 1, ]

  # Ensure rowid exists
  if (!"rowid" %in% names(baseline_data)) {
    baseline_data$rowid <- seq_len(nrow(baseline_data))
  }

  # --- 3. Create Counterfactual Datasets ---
  # For each combination in variables, create a copy of baseline_data with
  # the variable(s) set to that value

  if (!is.list(variables)) {
    var_list <- list()
    for (i in seq_along(variables)) {
      var_list[[variables[i]]] <- unique(baseline_data[[variables[i]]])
    }
  } else {
    var_list <- variables
  }

  grid <- do.call(expand.grid, var_list)

  # expanded_data_list <- vector("list", nrow(grid))
  # for (i in seq_len(nrow(grid))) {
  #   dt_copy <- baseline_data
  #   for (v in names(grid)) {s
  #     dt_copy[[v]] <- grid[i, v]
  #   }
  #   expanded_data_list[[i]] <- dt_copy
  # }

  # newdata_expanded <- do.call(rbind, expanded_data_list)

  newdata_expanded <- create_counterfactual_data(baseline_data, grid, var_list)

  if (inherits(model, "blrm")) {
    posterior_summary <- match.arg(posterior_summary)
    result <- avg_sops_blrm(
      model = model,
      newdata_expanded = newdata_expanded,
      grid = grid,
      variables = var_list,
      by = by,
      times = times,
      id_var = id_var,
      p2varname = p2varname,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      ...
    )
    attr(result, "newdata_orig") <- newdata
    return(result)
  }

  # --- 4. Compute Individual SOPs ---
  sops_ind <- sops(
    model,
    newdata = newdata_expanded,
    times = times,
    p2varname = p2varname,
    ...
  )
  resolved_times <- attr(sops_ind, "call_args")$times

  # --- 5. Aggregate (Marginalize) ---
  # Group by time, state, and the variables used for standardization
  group_cols <- c("time", "state", names(var_list))

  # Add optional 'by' variables
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }

  # Validate grouping columns exist
  missing_groups <- setdiff(group_cols, names(sops_ind))
  if (length(missing_groups) > 0) {
    stop("Grouping variables missing: ", paste(missing_groups, collapse = ", "))
  }

  # Aggregate
  agg_formula <- stats::as.formula(
    paste("estimate ~", paste(group_cols, collapse = " + "))
  )
  result <- stats::aggregate(
    agg_formula,
    data = sops_ind,
    FUN = mean,
    na.rm = TRUE
  )

  # --- 5. Store Attributes ---
  # Copy attributes from sops_ind
  attr(result, "model") <- attr(sops_ind, "model")
  attr(result, "call_args") <- attr(sops_ind, "call_args")
  attr(result, "tvarname") <- attr(sops_ind, "tvarname")
  attr(result, "pvarname") <- attr(sops_ind, "pvarname")
  attr(result, "p2varname") <- attr(sops_ind, "p2varname")
  attr(result, "ylevels") <- attr(sops_ind, "ylevels")
  attr(result, "absorb") <- attr(sops_ind, "absorb")
  attr(result, "gap") <- attr(sops_ind, "gap")
  attr(result, "t_covs") <- attr(sops_ind, "t_covs")

  # Specific attributes for avg_sops/inferences
  attr(result, "avg_args") <- list(
    variables = var_list,
    by = by,
    times = resolved_times,
    id_var = id_var
  )
  # Store ORIGINAL newdata for bootstrap (not the expanded counterfactual)
  attr(result, "newdata_orig") <- newdata

  class(result) <- c("markov_avg_sops", class(result))
  return(result)
}

avg_sops_blrm <- function(
  model,
  newdata_expanded,
  grid,
  variables,
  by,
  times,
  id_var,
  p2varname,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  ...
) {
  args <- list(...)
  ylevels <- args$ylevels %||% model$ylevels
  absorb <- args$absorb %||% NULL
  tvarname <- args$tvarname %||% "time"
  pvarname <- args$pvarname %||% "yprev"
  gap <- args$gap %||% NULL
  t_covs <- args$t_covs %||% NULL

  time_res <- resolve_sop_times(
    model,
    newdata_expanded,
    times,
    tvarname,
    t_covs = t_covs,
    default = "unique"
  )
  times <- time_res$times
  validate_factor_gap(gap, t_covs, time_res$time_info)

  missing_by <- setdiff(by %||% character(), names(newdata_expanded))
  if (length(missing_by) > 0) {
    stop("Grouping variables missing: ", paste(missing_by, collapse = ", "))
  }

  draw_indices <- select_posterior_draws(model, n_draws, seed)
  chunks <- split_draw_indices(
    draw_indices,
    chunk_size = getOption("markov.misc.blrm_avg_chunk_size", 50L)
  )
  gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  n_cf <- nrow(grid)
  n_each <- nrow(newdata_expanded) / n_cf
  draw_results <- vector("list", length(draw_indices))
  out_i <- 1L

  for (chunk in chunks) {
    arr <- soprob_markov(
      object = model,
      data = newdata_expanded,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      id_var = id_var,
      n_draws = NULL,
      .draw_indices = chunk,
      .gamma_draws = subset_cached_draw_matrix(gamma_draws, draw_indices, chunk)
    )

    for (i in seq_along(chunk)) {
      arr_i <- arr[i, , , , drop = FALSE]
      arr_i <- array(arr_i, dim = dim(arr_i)[-1])
      if (is.null(by)) {
        draw_i <- marginalize_sops_array(
          sops_array = arr_i,
          grid = grid,
          times = times,
          ylevels = ylevels,
          variables = variables,
          n_cf = n_cf,
          n_each = n_each
        )
      } else {
        draw_i <- array_to_df_individual(
          sops_array = arr_i,
          times = times,
          ylevels = ylevels,
          newdata = newdata_expanded,
          by = NULL
        )
        group_cols <- unique(c("time", "state", names(variables), by))
        agg_formula <- stats::as.formula(
          paste("estimate ~", paste(group_cols, collapse = " + "))
        )
        draw_i <- stats::aggregate(
          agg_formula,
          data = draw_i,
          FUN = mean,
          na.rm = TRUE
        )
      }
      draw_i$draw_id <- chunk[i]
      draw_results[[out_i]] <- draw_i
      out_i <- out_i + 1L
    }
  }

  draws_df <- bind_rows_fill(draw_results)
  group_cols <- c("time", "state", names(variables))
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }
  result <- summarize_posterior_draws_df(
    draws_df,
    group_cols = group_cols,
    posterior_summary = posterior_summary,
    conf_level = conf_level
  )

  attr(result, "model") <- model
  attr(result, "call_args") <- list(
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs,
    include_re = include_re,
    id_var = id_var,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level
  )
  attr(result, "tvarname") <- tvarname
  attr(result, "pvarname") <- pvarname
  attr(result, "p2varname") <- p2varname
  attr(result, "ylevels") <- ylevels
  attr(result, "absorb") <- absorb
  attr(result, "gap") <- gap
  attr(result, "t_covs") <- t_covs
  attr(result, "avg_args") <- list(
    variables = variables,
    by = by,
    times = times,
    id_var = id_var
  )
  attr(result, "method") <- "posterior"
  attr(result, "engine") <- "posterior"
  attr(result, "n_draws") <- length(draw_indices)
  attr(result, "draw_ids") <- draw_indices
  attr(result, "conf_level") <- conf_level
  attr(result, "posterior_summary") <- posterior_summary
  if (return_draws) {
    attr(result, "draws") <- draws_df
  }

  class(result) <- c("markov_avg_sops", class(result))
  result
}


#' Inference for State Occupation Probabilities
#'
#' Adds confidence intervals to SOP objects using simulation-based or bootstrap
#' methods. The default method is simulation. For objects produced from
#' `rmsb::blrm()` models, posterior uncertainty is already computed by
#' `sops()`/`avg_sops()`, so `inferences()` ignores `method` and returns the
#' object unchanged.
#'
#' @param object A `markov_avg_sops` object from `avg_sops()` or a
#'   `markov_sops` object from `sops()`.
#' @param method Character. Inference method:
#'   \itemize{
#'     \item `"simulation"` (default): Uses simulation engines that do not
#'       refit models. Works for both individual and averaged SOPs.
#'     \item `"bootstrap"`: Resamples patients with replacement, refits model.
#'  Only works for `markov_avg_sops` objects.
#'   }
#' @param engine Character. Simulation engine used when `method = "simulation"`:
#'   \itemize{
#'     \item `"mvn"` (default): Draws coefficients from MVN(beta_hat, Sigma).
#'     \item `"score_bootstrap"`: Uses one-step score perturbation with
#'       cluster-level exponential multipliers. Requires a `robcov_vglm` model
#'       or an `orm` model with `cluster`, and `avg_sops()` objects.
#'   }
#' @param score_weight_dist Character. Cluster weight distribution for
#'   `engine = "score_bootstrap"`. Currently only `"exponential"` is supported.
#' @param n_sim Number of simulation draws (for simulation) or bootstrap
#'   iterations (for bootstrap). Default is 1000. For `blrm` SOP objects this
#'   argument is ignored; rerun `sops()`/`avg_sops()` with `n_draws` and `seed`
#'   to change posterior draws.
#' @param vcov Optional custom variance-covariance matrix. If provided,
#'   overrides the vcov extracted from the model.
#' @param cluster Optional row-level cluster vector for
#'   `engine = "score_bootstrap"` with `orm` models. Values should match
#'   `id_var` in the baseline data used by `avg_sops()`.
#' @param workers Number of parallel workers. If NULL (default) or 1, uses
#'   sequential processing. If > 1, uses parallel processing with that many
#'   workers.
#' @param conf_level Confidence level for intervals (default 0.95). For `blrm`
#'   SOP objects this argument is ignored; rerun `sops()`/`avg_sops()` with
#'   `conf_level` to change posterior intervals.
#' @param conf_type Type of confidence interval (simulation method only):
#'   \itemize{
#'     \item `"perc"` (default): Percentile-based intervals from the simulation
#'       distribution.
#'     \item `"wald"`: Uses simulation standard errors with normal quantiles.
#'   }
#' @param return_draws Logical. If TRUE, stores all individual simulation/bootstrap
#'   draws as an attribute. Extract with `get_draws()`. Default is TRUE
#' @param update_datadist Logical. Whether to update datadist for rms models
#'   during bootstrap. Default is TRUE.
#' @param use_coefstart Logical. Use original coefficients as starting values
#'   for bootstrap refitting. Default is FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object with added columns:
#'   \item{conf.low}{Lower confidence bound}
#'   \item{conf.high}{Upper confidence bound}
#'   \item{std.error}{Standard error from simulation/bootstrap}
#'
#'   If `return_draws = TRUE`, the object also has a `"simulation_draws"` or
#'   `"bootstrap_draws"` attribute containing all individual draws. Extract
#'   with `get_draws()`.
#'
#' @details
#' ## Simulation Method
#'
#' The simulation method works as follows:
#' 1. Extract coefficient vector beta_hat and (robust) variance-covariance Sigma
#' 2. Generate n_sim coefficient vectors via the selected engine:
#'    \itemize{
#'      \item `engine = "mvn"`: MVN draws from `N(beta_hat, Sigma)`.
#'      \item `engine = "score_bootstrap"`: one-step score perturbation with
#'        cluster-level exponential multipliers.
#'    }
#' 3. For each draw, replace model coefficients and compute SOPs
#' 4. Compute confidence intervals from the empirical distribution
#'
#' - Works for both individual-level (`sops()`) and averaged (`avg_sops()`) SOPs
#'   for `engine = "mvn"`
#' - `engine = "score_bootstrap"` supports `avg_sops()` with `robcov_vglm`
#'   models and with `orm` models when `cluster` is supplied.
#'
#' ## Bootstrap Method
#'
#' The bootstrap method resamples patients and refits the model:
#' 1. Sample patient IDs with replacement
#' 2. Refit model on bootstrap sample (handles missing states through releveling)
#' 3. Compute SOPs using G-computation
#' 4. Compute percentile-based confidence intervals
#' Bootstrap requires the full longitudinal dataset (all time
#' points) in the original `newdata` passed to `avg_sops()`, not just baseline
#' data. This is because the model must be refit on each bootstrap sample.
#'
#'
#' This design ensures consistency: the same vcov is used for both point
#' estimates and inference, regardless of how `inferences()` is called.
#'
#' @seealso [avg_sops()], [sops()], [get_draws()], [robcov_vglm()],
#'   [set_coef()]
#'
#' @examples
#' \dontrun{
#' # Step 1: Fit model and wrap with cluster-robust vcov
#' fit <- vglm(y ~ time + tx + yprev, family = cumulative(reverse = TRUE,
#' parallel = TRUE), data = data)
#' fit_robust <- robcov_vglm(fit, cluster = data$id)
#'
#' # Step 2a: Simulation inference (use baseline data only)
#' baseline_data <- data |> filter(time == 1)
#' result_sim <- avg_sops(
#'   model = fit_robust,
#'   newdata = baseline_data,  # Baseline only for simulation
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6",
#'   id_var = "id"
#' ) |>
#'   inferences(method = "simulation", n_sim = 1000)
#'
#' # Step 2a-bis: Score bootstrap simulation (requires robcov_vglm)
#' result_score <- avg_sops(
#'   model = fit_robust,
#'   newdata = baseline_data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6",
#'   id_var = "id"
#' ) |>
#'   inferences(method = "simulation", engine = "score_bootstrap", n_sim = 1000)
#'
#' # orm fits use rms::robcov() for full robust covariance matrices.
#' dd <- rms::datadist(data)
#' options(datadist = "dd")
#' fit_orm <- rms::orm(y ~ time + tx + yprev, data = data, x = TRUE, y = TRUE)
#' fit_orm_robust <- rms::robcov(fit_orm, cluster = data$id)
#'
#' avg_orm <- avg_sops(
#'   model = fit_orm_robust,
#'   newdata = baseline_data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6",
#'   id_var = "id"
#' )
#' result_orm <- avg_orm |>
#'   inferences(method = "simulation", engine = "mvn", n_sim = 1000)
#'
#' result_orm_score <- avg_orm |>
#'   inferences(
#'     method = "simulation",
#'     engine = "score_bootstrap",
#'     cluster = data$id,
#'     n_sim = 1000
#'   )
#'
#' # Step 2b: Bootstrap inference (must use full data)
#' result_boot <- avg_sops(
#'   model = fit_robust,
#'   newdata = data,  # Full longitudinal data for bootstrap
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6",
#'   id_var = "id"
#' ) |>
#'   inferences(method = "bootstrap", n_sim = 500)
#'
#' # Individual-level SOPs with simulation inference
#' ind_result <- sops(
#'   model = fit_robust,
#'   newdata = baseline_data,
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6") |>
#'   inferences(n_sim = 500)
#'
#' # Extract draws for custom analyses
#' get_draws(ind_result)
#'
#' # Compute treatment effect on time in state
#' library(dplyr)
#' library(tidyr)
#' treatment_effect <- draws |>
#'   filter(state == 1) |>
#'   group_by(draw_id, tx) |>
#'   summarise(total_time = sum(draw), .groups = "drop") |>
#'   pivot_wider(names_from = tx, values_from = total_time) |>
#'   mutate(effect = `1` - `0`)
#' quantile(treatment_effect$effect, c(0.025, 0.5, 0.975))
#' }
#'
#' @seealso [get_draws()] to extract individual draws from any inference method
#'
#' @export
inferences <- function(
  object,
  method = "simulation",
  engine = "mvn",
  score_weight_dist = "exponential",
  n_sim = 1000,
  vcov = NULL,
  cluster = NULL,
  workers = NULL,
  conf_level = 0.95,
  conf_type = "perc",
  return_draws = TRUE,
  update_datadist = TRUE,
  use_coefstart = FALSE,
  ...
) {
  # --- Input Validation ---
  if (!inherits(object, c("markov_avg_sops", "markov_sops"))) {
    stop(
      "inferences() requires a 'markov_avg_sops' or 'markov_sops' object. ",
      "Got: ",
      paste(class(object), collapse = ", ")
    )
  }

  if (inherits(attr(object, "model"), "blrm")) {
    return(object)
  }

  method <- match.arg(method, choices = c("simulation", "bootstrap"))
  engine <- match.arg(engine, choices = c("mvn", "score_bootstrap"))

  if (method == "bootstrap" && engine != "mvn") {
    stop("`engine` is only used when `method = \"simulation\"`.")
  }

  # --- Dispatch to Method-Specific Implementation ---
  if (method == "simulation") {
    inferences_simulation(
      object = object,
      engine = engine,
      score_weight_dist = score_weight_dist,
      n_sim = n_sim,
      vcov = vcov,
      cluster = cluster,
      conf_level = conf_level,
      conf_type = conf_type,
      workers = workers,
      return_draws = return_draws
    )
  } else if (method == "bootstrap") {
    inferences_bootstrap(
      object = object,
      n_boot = n_sim,
      workers = workers,
      conf_level = conf_level,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )
  }
}


# =============================================================================
# SIMULATION-BASED INFERENCE (MVN)
# =============================================================================

#' Simulation-Based Inference Using Multivariate Normal
#'
#' Internal function that implements MVN simulation-based inference for SOPs.
#' Draws coefficient vectors from a multivariate normal distribution centered
#' at the point estimates with the (robust) variance-covariance matrix.
#'
#' @param object A `markov_avg_sops` or `markov_sops` object.
#' @param engine Simulation engine (`"mvn"` or `"score_bootstrap"`).
#' @param score_weight_dist Cluster weight distribution for score bootstrap.
#' @param n_sim Number of simulation draws.
#' @param vcov Custom variance-covariance matrix (optional, overrides model vcov).
#' @param cluster Row-level cluster vector for orm score bootstrap.
#' @param conf_level Confidence level.
#' @param conf_type Type of confidence interval ("perc" or "wald").
#' @param workers Number of parallel workers. If NULL or 1, uses sequential
#'   processing. If > 1, uses parallel processing.
#' @param return_draws Store individual draws?
#'
#' @return The input object with added confidence intervals.
#'
#' @keywords internal
inferences_simulation <- function(
  object,
  engine,
  score_weight_dist,
  n_sim,
  vcov,
  cluster,
  conf_level,
  conf_type,
  workers,
  return_draws
) {
  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  call_args <- attr(object, "call_args")
  avg_args <- attr(object, "avg_args")

  tvarname <- attr(object, "tvarname")
  pvarname <- attr(object, "pvarname")
  p2varname <- attr(object, "p2varname")
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
  gap <- attr(object, "gap")
  t_covs <- attr(object, "t_covs")

  # For avg_sops objects
  is_avg <- inherits(object, "markov_avg_sops")

  if (is_avg) {
    variables <- avg_args$variables
    by <- avg_args$by
    times <- avg_args$times
    id_var <- avg_args$id_var
  } else {
    # For individual sops objects
    times <- call_args$times
    variables <- NULL
    by <- call_args$by %||% attr(object, "by")
    id_var <- NULL
  }

  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform simulation inference.")
  }

  # --- 2. Prepare Prediction Data (COMPUTED ONCE) ---
  if (is_avg) {
    # For avg_sops: create counterfactual datasets
    baseline_data <- newdata_orig[!duplicated(newdata_orig[[id_var]]), ]
    grid <- do.call(expand.grid, variables)
    newdata_pred <- create_counterfactual_data(baseline_data, grid, variables)
    n_cf <- nrow(grid)
    n_each <- nrow(baseline_data)
  } else {
    # For individual sops: use data directly
    newdata_pred <- newdata_orig
    grid <- NULL
    n_cf <- 1
    n_each <- nrow(newdata_pred)
  }

  # --- 3. Generate Coefficient Draws ---
  baseline_weights_draws <- NULL
  if (engine == "mvn") {
    # Check for mvtnorm package
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      stop(
        "Package 'mvtnorm' is required for MVN simulation inference.\n",
        "Install with: install.packages('mvtnorm')"
      )
    }

    # The vcov is extracted directly from the model object:
    # - For robcov_vglm: uses the stored cluster-robust vcov
    # - For orm with rms::robcov(): uses the stored robust vcov
    # - For plain vglm/orm: uses model-based vcov
    beta_hat <- get_coef(model)

    if (!is.null(vcov) && is.matrix(vcov)) {
      Sigma <- vcov
      Sigma <- validate_coef_vcov(beta_hat, Sigma, arg = "vcov")
    } else {
      Sigma <- get_vcov_robust(model)
      Sigma <- validate_coef_vcov(beta_hat, Sigma, arg = "model vcov")
    }

    beta_draws <- mvtnorm::rmvnorm(n_sim, mean = beta_hat, sigma = Sigma)
  } else if (engine == "score_bootstrap") {
    if (!is_avg) {
      stop(
        "`engine = \"score_bootstrap\"` currently supports only ",
        "'markov_avg_sops' objects."
      )
    }
    if (!is.null(vcov)) {
      stop(
        "`vcov` cannot be supplied when `engine = \"score_bootstrap\"` because ",
        "draws are generated from score perturbations."
      )
    }

    score_draws <- generate_score_bootstrap_draws(
      model = model,
      baseline_data = baseline_data,
      id_var = id_var,
      n_sim = n_sim,
      score_weight_dist = score_weight_dist,
      cluster = cluster
    )
    beta_draws <- score_draws$beta_draws
    baseline_weights_draws <- score_draws$baseline_weights
  } else {
    stop("Unknown simulation engine: ", engine)
  }

  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 4. Detect Fast Path Eligibility ---
  # The fast path uses pre-computed design matrix decompositions for supported
  # ordinal model backends.
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  use_fast_path <- inherits(model_chk, c("vglm", "orm")) && is.null(p2varname)

  if (use_fast_path) {
    # --- FAST PATH: Pre-build components once, then run efficient Markov loop ---
    components <- tryCatch(
      markov_msm_build(
        model = model_chk,
        data = newdata_pred,
        t_covs = t_covs,
        times = times,
        ylevels = ylevels,
        absorb = absorb,
        tvarname = tvarname,
        pvarname = pvarname,
        gap = gap
      ),
      error = function(e) {
        warning(
          "Fast path build failed, falling back to slow path: ",
          e$message
        )
        NULL
      }
    )

    if (!is.null(components)) {
      # Define fast analysis function
      analysis_fn <- function(i) {
        baseline_weights <- NULL
        if (is_avg && !is.null(baseline_weights_draws)) {
          baseline_weights <- baseline_weights_draws[i, ]
        }

        # Compute effective coefficients from beta draw
        Gamma_i <- get_effective_coefs(model_chk, beta = beta_draws[i, ])

        # Run fast Markov simulation
        sops_array <- tryCatch(
          markov_msm_run(components, Gamma_i, times, absorb),
          error = function(e) {
            warning("markov_msm_run failed in draw ", i, ": ", e$message)
            return(NULL)
          }
        )

        if (is.null(sops_array)) {
          return(NULL)
        }

        # Marginalize if needed (for avg_sops)
        if (is_avg) {
          result <- marginalize_sops_array(
            sops_array = sops_array,
            grid = grid,
            times = times,
            ylevels = ylevels,
            variables = variables,
            n_cf = n_cf,
            n_each = n_each,
            weights = baseline_weights
          )
        } else {
          # Individual-level: convert array to data frame
          result <- array_to_df_individual(
            sops_array,
            times,
            ylevels,
            newdata_pred,
            by = by
          )
        }

        result
      }
    } else {
      # Fast path build failed, fall back to slow path
      use_fast_path <- FALSE
    }
  }

  if (!use_fast_path) {
    # --- SLOW PATH: Use soprob_markov with coefficient replacement ---
    analysis_fn <- function(i) {
      baseline_weights <- NULL
      if (is_avg && !is.null(baseline_weights_draws)) {
        baseline_weights <- baseline_weights_draws[i, ]
      }

      # Replace coefficients
      model_i <- set_coef(model, beta_draws[i, ])

      # Compute SOPs
      sops_array <- tryCatch(
        soprob_markov(
          object = model_i,
          data = newdata_pred,
          times = times,
          ylevels = ylevels,
          absorb = absorb,
          tvarname = tvarname,
          pvarname = pvarname,
          p2varname = p2varname,
          gap = gap,
          t_covs = t_covs
        ),
        error = function(e) {
          warning("soprob_markov failed in draw ", i, ": ", e$message)
          return(NULL)
        }
      )

      if (is.null(sops_array)) {
        return(NULL)
      }

      # Marginalize if needed (for avg_sops)
      if (is_avg) {
        result <- marginalize_sops_array(
          sops_array = sops_array,
          grid = grid,
          times = times,
          ylevels = ylevels,
          variables = variables,
          n_cf = n_cf,
          n_each = n_each,
          weights = baseline_weights
        )
      } else {
        # Individual-level: convert array to data frame
        result <- array_to_df_individual(
          sops_array,
          times,
          ylevels,
          newdata_pred,
          by = by
        )
      }

      result
    }
  }

  # --- 6. Apply Across All Draws ---
  use_parallel <- !is.null(workers) && workers > 1

  if (use_parallel) {
    # Setup parallel
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("Package 'furrr' is required for parallel processing")
    }

    future::plan(future.callr::callr, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    # Build globals list dynamically based on fast/slow path
    globals_list <- c(
      "model",
      "model_chk",
      "beta_draws",
      "newdata_pred",
      "times",
      "ylevels",
      "absorb",
      "tvarname",
      "pvarname",
      "p2varname",
      "gap",
      "t_covs",
      "is_avg",
      "grid",
      "variables",
      "n_cf",
      "n_each",
      "baseline_weights_draws",
      "use_fast_path",
      "by"
    )
    if (use_fast_path) {
      globals_list <- c(globals_list, "components")
    }

    sim_results <- furrr::future_map(
      seq_len(n_sim),
      analysis_fn,
      .options = furrr::furrr_options(
        seed = TRUE,
        globals = globals_list,
        packages = c("rms", "VGAM", "stats", "markov.misc")
      ),
      .progress = FALSE
    )
  } else {
    sim_results <- lapply(seq_len(n_sim), analysis_fn)
  }

  # --- 6. Combine Results ---
  sim_results <- Filter(Negate(is.null), sim_results)

  if (length(sim_results) == 0) {
    stop("All simulation draws failed.")
  }

  # Add draw_id
  for (i in seq_along(sim_results)) {
    sim_results[[i]]$draw_id <- i
  }
  draws_df <- bind_rows_fill(sim_results)

  # --- 7. Compute Confidence Intervals ---
  if (is_avg) {
    group_cols <- c("time", "state", names(variables))
    if (!is.null(by)) {
      group_cols <- unique(c(group_cols, by))
    }
  } else {
    # For individual SOPs: include rowid unless 'by' is specified
    if (is.null(by)) {
      group_cols <- c("rowid", "time", "state")
    } else {
      group_cols <- unique(c("time", "state", by))
    }
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = draws_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )

  # --- 8. Merge with Original Object ---
  final_result <- merge(object, summary_stats, by = group_cols, all.x = TRUE)

  # Restore attributes
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(final_result, a) <- attr(object, a)
    }
  }
  class(final_result) <- class(object)

  # Add metadata
  attr(final_result, "n_sim") <- n_sim
  attr(final_result, "n_successful") <- length(sim_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "conf_type") <- conf_type
  attr(final_result, "method") <- "simulation"
  attr(final_result, "engine") <- engine
  if (engine == "score_bootstrap") {
    attr(final_result, "score_weight_dist") <- score_weight_dist
  }

  if (return_draws) {
    attr(final_result, "simulation_draws") <- draws_df
  }

  final_result
}

#' Generate One-Step Score-Bootstrap Draws
#'
#' Computes simulation draws using cluster-level score perturbations and a
#' one-step Newton approximation. Also returns patient-level averaging weights
#' (shared with score perturbations) for `avg_sops` marginalization.
#'
#' @param model A `robcov_vglm` object with stored `scores`, `bread`, and
#'   `cluster` components, or an `orm` object fitted with `x = TRUE, y = TRUE`.
#'   For `orm`, `cluster` must be supplied.
#' @param baseline_data Baseline data (one row per patient).
#' @param id_var Name of patient ID variable in `baseline_data`.
#' @param n_sim Number of simulation draws.
#' @param score_weight_dist Cluster weight distribution. Currently only
#'   `"exponential"` is supported.
#' @param cluster Optional row-level cluster vector for `orm` models.
#'
#' @return A list with:
#'   \itemize{
#'     \item `beta_draws`: Matrix `[n_sim x p]` of perturbed coefficients.
#'     \item `baseline_weights`: Matrix `[n_sim x n_patients]` of normalized
#'       patient averaging weights.
#'   }
#'
#' @keywords internal
generate_score_bootstrap_draws <- function(
  model,
  baseline_data,
  id_var,
  n_sim,
  score_weight_dist = "exponential",
  cluster = NULL
) {
  if (!inherits(model, "robcov_vglm") && !inherits(model, "orm")) {
    stop(
      "`engine = \"score_bootstrap\"` requires a 'robcov_vglm' model ",
      "or an 'orm' model with `cluster` supplied. Wrap vglm fits with ",
      "robcov_vglm(fit, cluster = <id>), or call inferences(..., ",
      "engine = \"score_bootstrap\", cluster = <id>) for orm fits."
    )
  }

  score_weight_dist <- match.arg(score_weight_dist, choices = "exponential")

  if (!id_var %in% names(baseline_data)) {
    stop("ID variable '", id_var, "' not found in baseline data.")
  }

  components <- score_bootstrap_components(model, cluster = cluster)
  beta_hat <- components$coefficients
  bread <- components$bread
  scores <- components$scores
  cluster <- components$cluster

  if (!is.matrix(scores)) {
    stop("`model$scores` must be a matrix for score bootstrap.")
  }
  if (!is.matrix(bread)) {
    stop("`model$bread` must be a matrix for score bootstrap.")
  }
  if (length(cluster) != nrow(scores)) {
    stop(
      "Length of `model$cluster` (",
      length(cluster),
      ") does not match ",
      "the number of score rows (",
      nrow(scores),
      ")."
    )
  }
  if (length(beta_hat) != ncol(scores)) {
    stop(
      "Coefficient length (",
      length(beta_hat),
      ") does not match ",
      "score columns (",
      ncol(scores),
      ")."
    )
  }

  # Fix cluster ordering to first occurrence and aggregate score contributions.
  cluster_chr <- as.character(cluster)
  cluster_ids <- unique(cluster_chr)
  cluster_fac <- factor(cluster_chr, levels = cluster_ids)
  scores_clustered <- rowsum(scores, cluster_fac, reorder = FALSE)

  baseline_ids <- as.character(baseline_data[[id_var]])
  cluster_idx <- match(baseline_ids, cluster_ids)
  if (anyNA(cluster_idx)) {
    missing_ids <- unique(baseline_ids[is.na(cluster_idx)])
    stop(
      "Some baseline IDs used in avg_sops() were not found in the cluster ",
      "variable stored in robcov_vglm: ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) " ..." else "",
      ". Ensure avg_sops(id_var = ...) matches the clustering used in ",
      "robcov_vglm(cluster = ...)."
    )
  }

  n_clusters <- length(cluster_ids)
  n_pat <- length(baseline_ids)
  p <- length(beta_hat)

  beta_draws <- matrix(NA_real_, nrow = n_sim, ncol = p)
  colnames(beta_draws) <- names(beta_hat)
  baseline_weights <- matrix(NA_real_, nrow = n_sim, ncol = n_pat)

  for (i in seq_len(n_sim)) {
    # Exponential multipliers are centered at 1; centered form drives score perturbation.
    w_cluster <- stats::rexp(n_clusters, rate = 1)
    u_cluster <- w_cluster - 1

    # multiple each centered patient weight by the patient-level score
    # contribution, then sum across patients to get the overall score
    # perturbation for the draw.

    # Before applying the weights, colSums of the scores would be (by
    # definition) zero (because we're looking at the first derivative at the MLE
    #).

    # U_star is the perturbed total score vector for this draw. It
    # represents the total score perturbation under resampled patient weights.
    U_star <- as.vector(crossprod(u_cluster, scores_clustered))

    # One-Step Newton-Raphson Update:
    # If we ignore covariance: if var(beta_x) very low, then even if U_star is
    # large, then the update will be small.
    #
    # Reminder to self bread %*% U_star is equivalent to sum(bread[i, ] *
    # U_star) for each row of bread.

    # We can use Newton-Raphson to find root of the score equation, which gives
    # us the new MLE estimates for this draw.
    # We want to find where $f(x) = 0$ (our score)
    # $x_{new} = x_{old} - \frac{f(x_{old})}{f'(x_{old})}$
    # $x$ is is beta
    # $f(x)$ is the score equation, which is zero at the MLE
    # $f'(x)$ is the Hessian
    # $$\beta_{new} = \beta_{old} - \frac{U(\beta_{old})}{H(\beta_{old})}$$

    # Dividing by a matrix is the same as multiplying by the inverse, so we can
    # write:
    # $\beta_{new} = \beta_{old} - H^{-1} \times U(\beta_{old})$
    # vcov (V) is defined as the inverse of the negative Hessian ($V = (-H)^{-1}$), so we can rewrite:
    # $\beta_{new} = \beta_{old} + V \times U(\beta_{old})$
    # The new beta (x1) is where the tangent line at the old beta (x0) intersects the x-axis (where score = 0).

    # The bread component is vcov(fit), so this directly applies the
    # one-step update from the total score perturbation.
    beta_draws[i, ] <- beta_hat + as.vector(bread %*% U_star)

    w_patient <- w_cluster[cluster_idx]
    w_sum <- sum(w_patient)
    if (!is.finite(w_sum) || w_sum <= 0) {
      stop("Invalid baseline weights generated for score bootstrap.")
    }
    baseline_weights[i, ] <- w_patient / w_sum
  }

  list(beta_draws = beta_draws, baseline_weights = baseline_weights)
}

score_bootstrap_components <- function(model, cluster = NULL) {
  if (inherits(model, "robcov_vglm")) {
    required <- c("coefficients", "bread", "scores", "cluster")
    missing_required <- required[vapply(
      required,
      function(x) is.null(model[[x]]),
      logical(1)
    )]
    if (length(missing_required) > 0) {
      stop(
        "The supplied robcov_vglm model is missing required components for ",
        "score bootstrap: ",
        paste(missing_required, collapse = ", "),
        ". Refit with robcov_vglm()."
      )
    }

    return(list(
      coefficients = model$coefficients,
      bread = model$bread,
      scores = model$scores,
      cluster = model$cluster
    ))
  }

  if (!inherits(model, "orm")) {
    stop("score_bootstrap_components() supports robcov_vglm and orm models.")
  }

  if (is.null(cluster)) {
    stop(
      "`cluster` is required for orm score-bootstrap inference. ",
      "Supply the row-level patient ID vector used to fit the orm model, ",
      "for example `cluster = data$id`."
    )
  }

  scores <- compute_scores_orm(model)
  cluster <- align_cluster_orm(model, cluster, nrow(scores))
  bread <- get_orm_model_vcov(model)

  list(
    coefficients = stats::coef(model),
    bread = bread,
    scores = scores,
    cluster = cluster
  )
}

compute_scores_orm <- function(model) {
  if (!inherits(model, "orm")) {
    stop("'model' must be an orm object.")
  }
  if (is.null(model$x) || is.null(model$y)) {
    stop(
      "orm score bootstrap requires the model to be fitted with ",
      "`x = TRUE, y = TRUE`."
    )
  }

  beta <- stats::coef(model)
  Gamma <- get_effective_coefs(model, beta = beta)
  X <- cbind("(Intercept)" = 1, as.matrix(model$x))
  X <- X[, colnames(Gamma), drop = FALSE]

  eta <- X %*% t(Gamma)
  cum_probs <- stats::plogis(eta)
  probs <- lp_to_probs(eta, nrow(Gamma))
  probs <- pmax(probs, .Machine$double.eps)

  y_levels <- as_state_labels(model$yunique)
  y_idx <- match(as_state_labels(model$y), y_levels)
  if (anyNA(y_idx)) {
    stop("Could not match orm response values to fitted outcome levels.")
  }

  n <- nrow(X)
  M <- nrow(Gamma)
  d_eta <- matrix(0, nrow = n, ncol = M)

  for (i in seq_len(n)) {
    k <- y_idx[i]
    if (k == 1L) {
      d_eta[i, 1L] <- -cum_probs[i, 1L] * (1 - cum_probs[i, 1L]) /
        probs[i, 1L]
    } else if (k == M + 1L) {
      d_eta[i, M] <- cum_probs[i, M] * (1 - cum_probs[i, M]) /
        probs[i, M + 1L]
    } else {
      d_eta[i, k - 1L] <- cum_probs[i, k - 1L] *
        (1 - cum_probs[i, k - 1L]) / probs[i, k]
      d_eta[i, k] <- -cum_probs[i, k] *
        (1 - cum_probs[i, k]) / probs[i, k]
    }
  }

  intercept_scores <- d_eta
  slope_scores <- as.matrix(model$x) * rowSums(d_eta)
  scores <- cbind(intercept_scores, slope_scores)
  colnames(scores) <- names(beta)
  scores
}

align_cluster_orm <- function(model, cluster, n_scores) {
  if (inherits(cluster, "formula")) {
    stop("`cluster` for orm score bootstrap must be a row-level vector.")
  }

  if (length(cluster) == n_scores) {
    return(cluster)
  }

  na_action <- model$na.action
  if (!is.null(na_action) && length(cluster) > n_scores) {
    keep <- seq_along(cluster)
    keep <- keep[-as.integer(na_action)]
    if (length(keep) == n_scores) {
      return(cluster[keep])
    }
  }

  stop(
    "Length of `cluster` (",
    length(cluster),
    ") does not match the number of orm score rows (",
    n_scores,
    "). Supply the row-level patient ID vector used to fit the model."
  )
}

get_orm_model_vcov <- function(model) {
  beta <- stats::coef(model)

  if (!is.null(model$orig.var)) {
    return(validate_coef_vcov(beta, model$orig.var, arg = "orm model vcov"))
  }

  if (requireNamespace("rms", quietly = TRUE)) {
    robust <- rms::robcov(model)
    if (!is.null(robust$orig.var)) {
      return(validate_coef_vcov(beta, robust$orig.var, arg = "orm model vcov"))
    }
  }

  stop(
    "Could not obtain a full model-based covariance matrix for the orm fit. ",
    "Refit with rms::orm(..., x = TRUE, y = TRUE)."
  )
}


# =============================================================================
# BOOTSTRAP-BASED INFERENCE
# =============================================================================

#' Bootstrap-Based Inference for SOPs
#'
#' Internal function that implements bootstrap inference for SOPs. Resamples
#' patients with replacement, refits the model, and computes SOPs for each
#' bootstrap sample.
#'
#' @param object A `markov_avg_sops` object.
#' @param n_boot Number of bootstrap iterations.
#' @param workers Number of parallel workers. If NULL or 1, uses sequential
#'   processing. If > 1, uses parallel processing.
#' @param conf_level Confidence level.
#' @param return_draws Store individual bootstrap draws?
#' @param update_datadist Update datadist for rms models?
#' @param use_coefstart Use original coefficients as starting values?
#'
#' @return The input object with added confidence intervals.
#'
#' @keywords internal
inferences_bootstrap <- function(
  object,
  n_boot,
  workers,
  conf_level,
  return_draws,
  update_datadist,
  use_coefstart
) {
  # Bootstrap only supports avg_sops for now
  if (!inherits(object, "markov_avg_sops")) {
    stop(
      "Bootstrap inference currently only supports 'markov_avg_sops' objects. ",
      "For individual-level SOPs, use method = 'simulation'."
    )
  }

  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  call_args <- attr(object, "call_args")
  avg_args <- attr(object, "avg_args")

  tvarname <- attr(object, "tvarname")
  pvarname <- attr(object, "pvarname")
  p2varname <- attr(object, "p2varname")
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
  gap <- attr(object, "gap")
  t_covs <- attr(object, "t_covs")

  variables <- avg_args$variables
  by <- avg_args$by
  times <- avg_args$times
  id_var <- avg_args$id_var

  if (is.null(newdata_orig)) {
    stop("Original newdata not stored. Cannot perform bootstrap.")
  }

  if (!id_var %in% names(newdata_orig)) {
    stop("ID variable '", id_var, "' not found in data.")
  }

  # --- 2. Validate Data is Longitudinal (Not Just Baseline) ---
  # Check if tvarname exists and has multiple time points per patient
  if (tvarname %in% names(newdata_orig)) {
    rows_per_patient <- tabulate(match(newdata_orig[[id_var]], unique(newdata_orig[[id_var]])))

    if (all(rows_per_patient == 1)) {
      stop(
        "Bootstrap inference requires full longitudinal data (all time points), ",
        "but the data passed to avg_sops() appears to be baseline only ",
        "(one row per patient).\\n\\n",
        "For bootstrap: avg_sops(model, newdata = data, ...)\\n",
        "For simulation: avg_sops(model, newdata = data |> filter(time == 1), ...)"
      )
    }
  }

  # Dimensions for result expansion
  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 2. Identify Factor Columns ---
  factor_cols <- c("y", pvarname, p2varname)
  factor_cols <- intersect(factor_cols, names(newdata_orig))

  # --- 3. Generate Bootstrap ID Samples ---
  boot_ids <- fast_group_bootstrap(
    data = newdata_orig,
    id_var = id_var,
    n_boot = n_boot
  )

  # --- 4. Define Analysis Function ---
  analysis_fn <- function(boot_data) {
    # A. Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = newdata_orig,
      ylevels = ylevels,
      absorb = absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )

    m_boot <- boot_result$model
    boot_data <- boot_result$data
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }

    # Get states present in original numbering (for mapping back later)
    ylevel_names <- as_state_labels(ylevels)
    states_present <- ylevel_names[
      !ylevel_names %in% as_state_labels(missing_states)
    ]

    # B. Compute standardized SOPs using G-computation on bootstrap data
    # Extract baseline data (one row per patient)
    baseline_boot <- boot_data[!duplicated(boot_data[["new_id"]]), ]

    # Create counterfactual datasets for each variable value
    grid <- do.call(expand.grid, variables)
    cf_data_list <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
      dt_copy <- baseline_boot
      for (v in names(grid)) {
        dt_copy[[v]] <- grid[i, v]
      }
      cf_data_list[[i]] <- dt_copy
    }
    newdata_cf <- do.call(rbind, cf_data_list)

    # Compute individual SOPs for counterfactual data
    sops_array <- tryCatch(
      soprob_markov(
        object = m_boot,
        data = newdata_cf,
        times = times,
        ylevels = factor(boot_ylevels),
        absorb = boot_absorb,
        tvarname = tvarname,
        pvarname = pvarname,
        p2varname = p2varname,
        gap = gap,
        t_covs = t_covs
      ),
      error = function(e) {
        warning("soprob_markov failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(sops_array)) {
      return(NULL)
    }

    # C. Marginalize (average) across patients for each counterfactual
    # sops_array is [n_pat, n_times, n_boot_states]
    n_pat <- dim(sops_array)[1]
    n_boot_states <- dim(sops_array)[3]
    n_cf <- nrow(grid)
    n_each <- n_pat %/% n_cf

    # Average within each counterfactual group
    avg_sops_list <- vector("list", n_cf)
    for (cf_i in seq_len(n_cf)) {
      start_idx <- (cf_i - 1) * n_each + 1
      end_idx <- cf_i * n_each

      # Subset and average [n_each x n_times x n_boot_states]
      sops_cf <- sops_array[start_idx:end_idx, , , drop = FALSE]
      avg_sops_mat <- apply(sops_cf, c(2, 3), mean) # [n_times x n_boot_states]

      avg_sops_list[[cf_i]] <- avg_sops_mat
    }

    # D. Expand back to original state space with zero-padding
    if (length(missing_states) > 0) {
      for (cf_i in seq_len(n_cf)) {
        avg_sops_mat <- avg_sops_list[[cf_i]]
        full_mat <- matrix(0, nrow = n_times, ncol = n_states)

        # Fill in states that were present
        for (i in seq_along(states_present)) {
          original_state <- states_present[i]
          original_idx <- which(ylevel_names == as_state_labels(original_state))
          full_mat[, original_idx] <- avg_sops_mat[, i]
        }

        avg_sops_list[[cf_i]] <- full_mat
      }
    }

    # E. Format as data frame
    result_list <- vector("list", n_cf)
    for (cf_i in seq_len(n_cf)) {
      avg_sops_mat <- avg_sops_list[[cf_i]]

      df <- expand.grid(time = times, state = ylevels)
      df$estimate <- as.vector(avg_sops_mat)

      # Add variable values for this counterfactual
      for (v in names(grid)) {
        df[[v]] <- grid[cf_i, v]
      }

      result_list[[cf_i]] <- df
    }

    bind_rows_fill(result_list)
  }

  # --- 5. Apply to Bootstrap Samples ---
  boot_results <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = newdata_orig,
    id_var = id_var,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "stats"),
    globals = c(
      "model",
      "variables",
      "times",
      "ylevels",
      "absorb",
      "pvarname",
      "p2varname",
      "tvarname",
      "gap",
      "t_covs",
      "n_times",
      "n_states",
      "update_datadist",
      "factor_cols",
      "use_coefstart"
    )
  )

  # --- 6. Combine and Compute Summary Statistics ---
  boot_results <- Filter(Negate(is.null), boot_results)

  if (length(boot_results) == 0) {
    stop("All bootstrap iterations failed.")
  }

  # Add draw_id and combine
  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- bind_rows_fill(boot_results)

  # Compute confidence intervals
  group_cols <- c("time", "state", names(variables))
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = boot_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = "perc" # Bootstrap always uses percentile
  )

  # Merge with original object
  final_result <- merge(object, summary_stats, by = group_cols, all.x = TRUE)

  # Restore attributes
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(final_result, a) <- attr(object, a)
    }
  }
  class(final_result) <- class(object)

  # Add bootstrap metadata
  attr(final_result, "n_boot") <- n_boot
  attr(final_result, "n_successful") <- length(boot_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "method") <- "bootstrap"

  # Store full bootstrap draws if requested
  if (return_draws) {
    attr(final_result, "bootstrap_draws") <- boot_df
  }

  final_result
}


# =============================================================================
# HELPER FUNCTIONS FOR INFERENCE
# =============================================================================

#' Create Counterfactual Datasets for G-Computation
#'
#' Creates copies of baseline data with treatment variable set to each level.
#'
#' @param baseline_data Data frame with one row per patient.
#' @param grid Data frame of variable combinations.
#' @param variables Named list of variable values.
#'
#' @return Data frame with counterfactual data stacked.
#'
#' @keywords internal
create_counterfactual_data <- function(baseline_data, grid, variables) {
  cf_data_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt_copy <- baseline_data
    for (v in names(grid)) {
      dt_copy[[v]] <- grid[i, v]
    }
    cf_data_list[[i]] <- dt_copy
  }
  do.call(rbind, cf_data_list)
}


#' Marginalize SOPs Array Over Patients
#'
#' Averages individual-level SOPs to get population-average (marginal) SOPs.
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param grid Data frame of variable combinations.
#' @param times Vector of time points.
#' @param ylevels Vector of state levels.
#' @param variables Named list of variables.
#' @param n_cf Number of counterfactual scenarios.
#' @param n_each Number of patients per scenario.
#'
#' @return Data frame with marginalized SOPs.
#'
#' @keywords internal
marginalize_sops_array <- function(
  sops_array,
  grid,
  times,
  ylevels,
  variables,
  n_cf,
  n_each,
  weights = NULL
) {
  n_times <- length(times)
  n_states <- length(ylevels)

  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != n_each || any(weights < 0)) {
      stop("`weights` must be a non-negative numeric vector of length n_each.")
    }
    weight_sum <- sum(weights)
    if (!is.finite(weight_sum) || weight_sum <= 0) {
      stop("`weights` must have a positive finite sum.")
    }
    weights <- weights / weight_sum
  }

  # Average within each counterfactual group
  avg_sops_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    start_idx <- (cf_i - 1) * n_each + 1
    end_idx <- cf_i * n_each

    sops_cf <- sops_array[start_idx:end_idx, , , drop = FALSE]
    if (is.null(weights)) {
      avg_sops_mat <- apply(sops_cf, c(2, 3), mean)
    } else {
      sops_cf_mat <- matrix(
        aperm(sops_cf, c(2, 3, 1)),
        nrow = n_times * n_states,
        ncol = n_each
      )
      avg_sops_mat <- matrix(
        as.vector(sops_cf_mat %*% weights),
        nrow = n_times,
        ncol = n_states
      )
    }

    avg_sops_list[[cf_i]] <- avg_sops_mat
  }

  # Format as data frame
  result_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    avg_sops_mat <- avg_sops_list[[cf_i]]

    df <- expand.grid(time = times, state = ylevels)
    df$estimate <- as.vector(avg_sops_mat)

    for (v in names(grid)) {
      df[[v]] <- grid[cf_i, v]
    }

    result_list[[cf_i]] <- df
  }

  bind_rows_fill(result_list)
}


#' Convert SOPs Array to Individual-Level Data Frame
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param times Vector of time points.
#' @param ylevels Vector of state levels.
#' @param newdata Original data with rowid.
#' @param by Optional character vector of variables to aggregate by.
#'
#' @return Data frame with individual SOPs (or aggregated if by is specified).
#'
#' @keywords internal
array_to_df_individual <- function(
  sops_array,
  times,
  ylevels,
  newdata,
  by = NULL
) {
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]

  # Flatten array
  probs_flat <- as.vector(sops_array)

  # Construct indices
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)

  # Build result by repeating newdata rows
  result <- newdata[idx_pat, , drop = FALSE]

  # Add rowid if it exists
  if ("rowid" %in% names(newdata)) {
    result$rowid <- newdata$rowid[idx_pat]
  } else {
    result$rowid <- idx_pat
  }

  result$time <- idx_time
  result$state <- idx_state
  result$estimate <- probs_flat
  rownames(result) <- NULL

  # Apply stratified aggregation if 'by' is specified
  if (!is.null(by)) {
    # Validate that by variables exist
    missing_vars <- setdiff(by, names(result))
    if (length(missing_vars) > 0) {
      warning(
        "Variables in 'by' not found in result: ",
        paste(missing_vars, collapse = ", "),
        ". Skipping aggregation."
      )
      return(result)
    }

    # Group by time, state, and stratification variables, then average
    group_cols <- unique(c("time", "state", by))

    agg_formula <- stats::as.formula(
      paste("estimate ~", paste(group_cols, collapse = " + "))
    )
    result <- stats::aggregate(
      agg_formula,
      data = result,
      FUN = mean,
      na.rm = TRUE
    )
  }

  result
}


#' Compute Confidence Intervals from Draws
#'
#' Computes confidence intervals from simulation or bootstrap draws.
#'
#' @param draws_df Data frame with draws (must have `estimate` and `draw_id` columns).
#' @param group_cols Character vector of grouping columns.
#' @param conf_level Confidence level (default 0.95).
#' @param conf_type Type of CI: "perc" (percentile) or "wald".
#'
#' @return Data frame with summary statistics.
#'
#' @keywords internal
compute_ci_from_draws <- function(
  draws_df,
  group_cols,
  conf_level = 0.95,
  conf_type = "perc"
) {
  alpha <- 1 - conf_level

  # Aggregate to get summary statistics
  agg_formula <- stats::as.formula(
    paste("estimate ~", paste(group_cols, collapse = " + "))
  )

  if (conf_type == "perc") {
    # Percentile confidence intervals
    summary_stats <- stats::aggregate(
      agg_formula,
      data = draws_df,
      FUN = function(x) {
        c(
          conf.low = stats::quantile(x, alpha / 2, na.rm = TRUE),
          conf.high = stats::quantile(x, 1 - alpha / 2, na.rm = TRUE),
          std.error = stats::sd(x, na.rm = TRUE)
        )
      }
    )
  } else if (conf_type == "wald") {
    # Wald confidence intervals
    summary_stats <- stats::aggregate(
      agg_formula,
      data = draws_df,
      FUN = function(x) {
        se <- stats::sd(x, na.rm = TRUE)
        critical <- abs(stats::qnorm(alpha / 2))
        mean_est <- mean(x, na.rm = TRUE)
        c(
          conf.low = mean_est - critical * se,
          conf.high = mean_est + critical * se,
          std.error = se
        )
      }
    )
  } else {
    stop("conf_type must be 'perc' or 'wald'")
  }

  # Fix aggregate's matrix column output
  mat <- summary_stats$estimate
  summary_stats$estimate <- NULL
  summary_stats$conf.low <- mat[, 1]
  summary_stats$conf.high <- mat[, 2]
  summary_stats$std.error <- mat[, 3]

  summary_stats
}


#' Extract Individual Draws from Inference Objects
#'
#' Extracts the individual draws (bootstrap samples, simulated values from MVN, etc.)
#' from an object returned by `inferences()` with `return_draws = TRUE`.
#' This function joins the draws back to the original point estimate object,
#' preserving all covariates, grouping variables, and summary statistics.
#'
#' @param object An object returned by `inferences()` with `return_draws = TRUE`.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item draw_id: Simulation or bootstrap iteration number
#'     \item time: Time point
#'     \item state: State number
#'     \item draw: Draw-specific estimate of state occupation probability
#'     \item estimate: The original point estimate from the model
#'     \item conf.low, conf.high, std.error: Summary statistics from the point estimate object
#'     \item Additional columns from the original object (covariates, etc.)
#'   }
#'   Each row represents one draw for a specific time-state combination.
#'
#' @details
#' This function retrieves the draws from objects created by `inferences()`.
#' It performs a join between the draws and the point estimate object to ensure
#' that all metadata is preserved. To avoid conflict, the `estimate` column from
#' the draws is renamed to `draw`.
#'
#' @examples
#' \dontrun{
#' # Create object with bootstrap draws
#' result <- avg_sops(
#'   model = fit,
#'   newdata = data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:30,
#'   ylevels = 1:6,
#'   absorb = 6,
#'   id_var = "id"
#' ) |>
#'   inferences(n_boot = 1000, return_draws = TRUE)
#'
#' # Extract draws
#' draws <- get_draws(result)
#'
#' # Plot distribution for state 1 at time 10 under treatment
#' library(ggplot2)
#' draws |>
#'   filter(time == 10, state == 1, tx == 1) |>
#'   ggplot(aes(x = draw)) +
#'   geom_histogram(bins = 30) +
#'   labs(title = "Distribution: P(State 1 | Time 10, Treatment)")
#'
#' # Compute mean time in state 1 with bootstrap CI
#' library(dplyr)
#' time_in_state_boot <- draws |>
#'   filter(state == 1, tx == 1) |>
#'   group_by(draw_id) |>
#'   summarise(total_time = sum(draw))
#'
#' quantile(time_in_state_boot$total_time, c(0.025, 0.5, 0.975))
#'
#' # Compare treatment effect on time in state
#' treatment_effect <- draws |>
#'   filter(state == 1) |>
#'   group_by(draw_id, tx) |>
#'   summarise(total_time = sum(draw), .groups = "drop") |>
#'   pivot_wider(names_from = tx, values_from = total_time,
#'               names_prefix = "tx") |>
#'   mutate(effect = tx1 - tx0)
#'
#' quantile(treatment_effect$effect, c(0.025, 0.5, 0.975))
#' }
#'
#' @seealso [inferences()], [avg_sops()]
#'
#' @export
get_draws <- function(object) {
  if (!inherits(object, c("markov_avg_sops", "markov_sops"))) {
    stop(
      "get_draws() requires an object from inferences(). ",
      "Got: ",
      paste(class(object), collapse = ", ")
    )
  }

  # 1. Extract draws attribute
  draws <- attr(object, "bootstrap_draws")
  if (is.null(draws)) {
    draws <- attr(object, "simulation_draws")
  }
  if (is.null(draws)) {
    draws <- attr(object, "draws")
  }

  if (is.null(draws)) {
    method <- attr(object, "method")
    msg <- "No draws found. Run inferences() with return_draws = TRUE."
    if (!is.null(method)) {
      msg <- paste0(msg, " (Method used: '", method, "')")
    }
    stop(msg)
  }

  # 2. Prepare metadata from the original object
  meta <- as.data.frame(object)

  # Rename 'estimate' in draws to 'draw' to avoid conflict with point estimates
  if ("estimate" %in% names(draws)) {
    names(draws)[names(draws) == "estimate"] <- "draw"
  }

  # 3. Identify join keys
  # Keys are columns present in both draws and meta (excluding estimate/draw)
  common_cols <- intersect(names(draws), names(meta))
  keys <- setdiff(common_cols, c("draw", "estimate"))

  if (length(keys) == 0) {
    return(draws)
  }

  # 4. Join metadata back to draws
  draws <- merge(draws, meta, by = keys, all.x = TRUE, sort = FALSE)

  return(draws)
}


# =============================================================================
# FAST PATH HELPER FUNCTIONS
# =============================================================================
# These internal functions provide the fast path for simulation-based inference
# on VGLM models by pre-computing design matrices and bypassing VGAM's
# predict() function.

#' Build Fast Markov Simulation Components (Internal)
#'
#' Pre-calculates design matrices for fast Markov simulation. This is the
#' "heavy lifting" step that should be done ONCE before running thousands of
#' simulation draws.
#'
#' @param model Fitted vglm model
#' @param data Baseline data frame (one row per patient)
#' @param t_covs Time-dependent covariate lookup (optional)
#' @param times Time points to predict
#' @param ylevels State labels
#' @param absorb Absorbing state(s). These transition matrices are not needed
#'   because absorbing states keep their prior probability mass.
#' @param tvarname Name of time variable
#' @param pvarname Name of previous state variable
#' @param gap Optional gap variable. With factor visit time, numeric gap values
#'   must be supplied in `t_covs`.
#' @param ... Ignored
#'
#' @return A list containing pre-calculated matrices and metadata.
#'
#' @keywords internal
markov_msm_build <- function(
  model,
  data,
  t_covs = NULL,
  times = NULL,
  ylevels = 1:6,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  gap = NULL,
  ...
) {
  n_pat <- nrow(data)
  ylevel_names <- as_state_labels(ylevels)
  absorb_names <- as_state_labels(absorb)
  n_states <- length(ylevel_names)
  M <- n_states - 1
  absorb_idx <- if (length(absorb_names)) {
    which(ylevel_names %in% absorb_names)
  } else {
    integer(0)
  }

  time_res <- resolve_sop_times(
    model,
    data,
    times,
    tvarname,
    t_covs = t_covs,
    default = "fast"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap, t_covs, time_info)

  # Need Gamma structure to know which columns to keep
  Gamma_template <- get_effective_coefs(model)
  common_cols <- colnames(Gamma_template)

  # Prepare data
  if (!pvarname %in% names(data)) {
    data[[pvarname]] <- factor(ylevel_names[1], levels = ylevel_names)
  }
  if (!tvarname %in% names(data)) {
    data <- assign_sop_time(data, tvarname, times[1], time_info)
  }

  is_vglm <- inherits(model, "vglm")
  is_orm <- inherits(model, "orm")

  # Terms object
  if (is_vglm) {
    tt <- stats::terms(model)
    tt <- stats::delete.response(tt)
    contrasts_arg <- if (length(model@contrasts)) {
      model@contrasts
    } else {
      NULL
    }
    xlev_arg <- model@xlevels
  } else if (is_orm) {
    tt <- NULL
    contrasts_arg <- NULL
    xlev_arg <- NULL
  } else {
    stop("markov_msm_build() supports only vglm and orm models.")
  }

  # Helper to get design matrix columns used by the effective coefficient
  # matrix. For inline transforms such as rms::rcs(), model.matrix() uses the
  # fitted terms object's predvars, including stored knot locations.
  get_X <- function(d) {
    if (is_orm) {
      X <- orm_model_matrix(model, d, include_intercept = TRUE)
    } else {
      X <- stats::model.matrix(
        tt,
        data = d,
        contrasts.arg = contrasts_arg,
        xlev = xlev_arg
      )
    }
    common <- intersect(colnames(X), common_cols)
    if (length(common) == 0) {
      return(NULL)
    }
    X[, common, drop = FALSE]
  }

  set_time <- function(d, time_idx) {
    assign_sop_visit(
      d,
      tvarname = tvarname,
      times = times,
      index = time_idx,
      t_covs = t_covs,
      gap = gap,
      time_info = time_info
    )
  }

  # A. Baseline matrix for T=1 uses each patient's observed previous state.
  d_init <- set_time(data, 1)
  d_init[[pvarname]] <- normalize_previous_state_column(
    d_init[[pvarname]],
    data[[pvarname]],
    pvarname
  )
  X_init <- get_X(d_init)
  if (is.null(X_init)) {
    stop("Could not construct a design matrix for the fitted model.")
  }
  col_names <- colnames(X_init)

  align_X <- function(X) {
    missing_cols <- setdiff(col_names, colnames(X))
    if (length(missing_cols) > 0) {
      stop(
        "Design matrix columns missing during fast path build: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    X[, col_names, drop = FALSE]
  }

  # B. Transition matrices: one matrix for each time and previous-state value.
  X_transition <- vector("list", length(times))
  for (time_idx in seq_along(times)) {
    X_transition[[time_idx]] <- vector("list", n_states)
    d_time <- set_time(data, time_idx)

    for (k in seq_len(n_states)) {
      if (k %in% absorb_idx) {
        X_transition[[time_idx]][[k]] <- NULL
        next
      }

      d_k <- d_time
      d_k[[pvarname]] <- make_previous_state_column(
        states = ylevel_names[k],
        prototype = data[[pvarname]],
        n = nrow(d_k),
        pvarname = pvarname
      )
      X_transition[[time_idx]][[k]] <- align_X(get_X(d_k))
    }
  }

  list(
    X_init = align_X(X_init),
    X_transition = X_transition,
    n_pat = n_pat,
    n_states = n_states,
    M = M,
    ylevels = ylevel_names,
    col_names = col_names
  )
}


#' Run Fast Markov Simulation (Internal)
#'
#' Runs the Markov loop using pre-calculated components and a specific
#' coefficient vector.
#'
#' @param components List returned by `markov_msm_build`
#' @param Gamma Effective coefficient matrix with dimensions M x P (M thresholds by P predictors)
#' @param times Vector of time points
#' @param absorb Absorbing state(s). Can be NULL for no absorbing states.
#'
#' @return Array of state probabilities with dimensions n_pat x n_times x n_states
#'
#' @keywords internal
markov_msm_run <- function(components, Gamma, times, absorb = NULL) {
  # Unpack
  X_init <- components$X_init
  X_transition <- components$X_transition
  n_pat <- components$n_pat
  n_states <- components$n_states
  M <- components$M
  ylevels <- components$ylevels
  col_names <- components$col_names
  n_times <- length(times)

  # Align Gamma to components
  Gamma <- Gamma[, col_names, drop = FALSE]
  Gamma_t <- t(Gamma) # [P x M]

  # Helper: X * Gamma_t -> LP [N x M]
  calc_lp <- function(X) {
    X %*% Gamma_t
  }

  # Simulation
  P_out <- array(0, dim = c(n_pat, n_times, n_states))
  dimnames(P_out)[[3]] <- ylevels

  absorb_idx <- if (!is.null(absorb)) {
    which(as.character(ylevels) %in% as.character(absorb))
  } else {
    integer(0)
  }
  non_absorb_idx <- setdiff(1:n_states, absorb_idx)

  # T=1
  P_out[, 1, ] <- lp_to_probs(calc_lp(X_init), M)

  # Loop 2..T
  if (n_times >= 2) {
    for (t_idx in 2:n_times) {
      p_current <- matrix(0, nrow = n_pat, ncol = n_states)

      for (k in non_absorb_idx) {
        p_prev_k <- P_out[, t_idx - 1, k]
        if (max(p_prev_k) < 1e-12) {
          next
        }

        LP_k <- calc_lp(X_transition[[t_idx]][[k]])
        trans_k <- lp_to_probs(LP_k, M)
        p_current <- p_current + (trans_k * p_prev_k)
      }

      for (a in absorb_idx) {
        p_current[, a] <- p_current[, a] + P_out[, t_idx - 1, a]
      }
      P_out[, t_idx, ] <- p_current
    }
  }

  P_out
}


#' Convert Cumulative Log-Odds to Probabilities (Internal)
#'
#' Converts a matrix of cumulative log-odds to category probabilities.
#'
#' @param eta Matrix of linear predictors with dimensions N x M (N observations by M thresholds)
#' @param M Number of thresholds (one less than the number of categories)
#' @return Matrix of probabilities with dimensions N x (M+1)
#'
#' @keywords internal
lp_to_probs <- function(eta, M) {
  cum_probs <- stats::plogis(eta) # [N x M]
  n_rows <- nrow(eta)
  probs <- matrix(0, nrow = n_rows, ncol = M + 1)

  probs[, 1] <- 1 - cum_probs[, 1]

  if (M > 1) {
    probs[, 2:M] <- cum_probs[, 1:(M - 1)] - cum_probs[, 2:M]
  }

  probs[, M + 1] <- cum_probs[, M]
  probs[probs < 0] <- 0

  probs
}


#' Compute Effective Coefficients from Beta Vector (Internal)
#'
#' Transforms a vector of VGLM coefficients into an effective coefficient matrix
#' by applying the constraint matrices.
#'
#' @param beta Coefficient vector from a VGLM fit
#' @param C_list Constraint matrices list from VGAM::constraints()
#' @return Effective coefficient matrix with dimensions M x P (M linear predictors by P terms)
#'
#' @keywords internal
compute_Gamma <- function(beta, C_list) {
  M <- nrow(C_list[[1]])
  P <- length(C_list)
  term_names <- names(C_list)

  G <- matrix(0, nrow = M, ncol = P)
  colnames(G) <- term_names
  rownames(G) <- paste0("eta", 1:M)

  curr <- 1
  for (j in seq_along(C_list)) {
    k <- ncol(C_list[[j]])
    chunk <- beta[curr:(curr + k - 1)]
    G[, j] <- C_list[[j]] %*% chunk
    curr <- curr + k
  }
  G
}
