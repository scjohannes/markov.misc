# SOP model adapters and time helpers.

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

  if (inherits(model_chk, "blrm")) {
    return(invisible(NULL))
  }

  # Check if it's a vglm model
  if (inherits(model_chk, "vglm")) {
    # Extract family object
    fam <- model_chk@family

    # Check if it's cumulative family
    if (!grepl("cumulative", fam@vfamily[1], ignore.case = TRUE)) {
      stop(
        "For vglm models, only the cumulative family is supported.\n",
        "Current family: ",
        fam@vfamily[1],
        "\n",
        "Please refit your model with: family = cumulative(reverse = TRUE, ...)"
      )
    }

    link <- unique(as.character(model_chk@misc$link))
    if (length(link) != 1L || !identical(tolower(link), "logitlink")) {
      stop_unsupported_markov_link(link)
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
        "For vglm models, the cumulative family must use reverse = TRUE.\n",
        "Current setting: reverse = ",
        reverse_param,
        "\n",
        "Please refit your model with: family = cumulative(reverse = TRUE, ...)"
      )
    }
  }

  if (inherits(model_chk, "orm")) {
    family <- model_chk$family
    if (
      is.null(family) ||
        length(family) != 1L ||
        !identical(tolower(as.character(family)), "logistic")
    ) {
      stop_unsupported_markov_link(family)
    }
  }

  # Other model types will be caught by the existing class check in soprob_markov

  invisible(NULL)
}

stop_unsupported_markov_link <- function(link = NULL) {
  message <- paste0(
    "Only cumulative-logit models are supported; ",
    "refit the model with a logit link."
  )
  condition <- structure(
    list(
      message = message,
      call = NULL,
      link = if (length(link)) as.character(link) else NA_character_
    ),
    class = c("markov_misc_unsupported_link", "error", "condition")
  )
  stop(condition)
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
    if (
      !is.numeric(n_draws) ||
        length(n_draws) != 1 ||
        is.na(n_draws) ||
        n_draws < 1
    ) {
      stop("`n_draws` must be a positive integer or NULL.")
    }
    n <- min(as.integer(n_draws), n_available)
  }

  if (n == n_available) {
    return(seq_len(n_available))
  }

  if (!is.null(seed)) {
    old_seed_exists <- exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
    old_seed <- if (old_seed_exists) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit(
      {
        if (old_seed_exists) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
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
  id_var <- markov_model_id_var(model, id_var)
  if (!is.null(id_var)) {
    return(id_var)
  }

  cluster_name <- model$clusterInfo$name
  if (
    !is.null(cluster_name) && length(cluster_name) == 1 && nzchar(cluster_name)
  ) {
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

resolve_blrm_cppo <- function(cppo, envir = parent.frame()) {
  if (is.function(cppo)) {
    return(cppo)
  }

  if (
    !is.character(cppo) || length(cppo) != 1L || is.na(cppo) || !nzchar(cppo)
  ) {
    stop(
      "Only constrained partial proportional odds `blrm` models are supported."
    )
  }

  is_syntactic_name <- make.names(cppo) == cppo && !grepl("^\\.[0-9]", cppo)
  if (!is_syntactic_name) {
    stop(
      "String `cppo` values must name an available function. ",
      "Inline expressions are not supported."
    )
  }

  tryCatch(
    get(cppo, mode = "function", envir = envir, inherits = TRUE),
    error = function(e) {
      stop(
        "String `cppo` values must name an available function. ",
        "Could not find function `",
        cppo,
        "`.",
        call. = FALSE
      )
    }
  )
}

predict_blrm_response_markov <- function(
  object,
  newdata,
  include_re = FALSE,
  id_var = NULL,
  draw_indices = NULL,
  gamma_draws = NULL,
  prediction_cache = NULL
) {
  if (is.null(draw_indices)) {
    draw_indices <- seq_len(nrow(object$draws))
  }
  manual_error <- NULL

  manual <- tryCatch(
    {
      cache_index <- NULL
      if (!is.null(prediction_cache)) {
        prediction_cache$cursor <- prediction_cache$cursor + 1L
        cache_index <- prediction_cache$cursor
      }
      X <- if (
        !is.null(cache_index) && length(prediction_cache$X) >= cache_index
      ) {
        prediction_cache$X[[cache_index]]
      } else {
        value <- blrm_design_matrix(object, newdata = newdata, second = FALSE)
        if (!is.null(cache_index)) {
          prediction_cache$X[[cache_index]] <- value
        }
        value
      }
      draws <- object$draws[draw_indices, , drop = FALSE]
      ndraws <- nrow(draws)
      ns <- object$non.slopes
      K <- ns + 1L
      cn <- colnames(draws)
      tauinfo <- object$tauInfo
      tau_names <- if (!is.null(tauinfo) && length(tauinfo$name)) {
        tauinfo$name
      } else {
        character()
      }

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
        cppo <- resolve_blrm_cppo(object$cppo, envir = parent.frame())
        Z <- if (
          !is.null(cache_index) && length(prediction_cache$Z) >= cache_index
        ) {
          prediction_cache$Z[[cache_index]]
        } else {
          value <- blrm_design_matrix(object, newdata = newdata, second = TRUE)
          if (!is.null(cache_index)) {
            prediction_cache$Z[[cache_index]] <- value
          }
          value
        }
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
        cum_probs[,, k] <- stats::plogis(ep)
      }

      out[,, 1] <- 1 - cum_probs[,, 1]
      if (K > 2) {
        out[,, 2:(K - 1)] <- cum_probs[,, 1:(K - 2), drop = FALSE] -
          cum_probs[,, 2:(K - 1), drop = FALSE]
      }
      out[,, K] <- cum_probs[,, K - 1]
      out[out < 0] <- 0

      normalize_probability_array(out)
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
  on.exit(
    {
      options(contrasts = oldopts$contrasts)
      options(Design.attr = oldopts$Design.attr)
    },
    add = TRUE
  )

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
        paste(
          utils::head(unique(as.character(values[bad])), 5),
          collapse = ", "
        ),
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

coerce_previous_state_values <- function(values, prototype, p_var) {
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
        "`y_levels` must be integer-compatible when `",
        p_var,
        "` is an integer previous-state variable."
      )
    }
    return(out)
  }

  if (is.numeric(prototype)) {
    out <- suppressWarnings(as.numeric(labels))
    if (anyNA(out) && any(!is.na(labels))) {
      stop(
        "`y_levels` must be numeric-compatible when `",
        p_var,
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

make_previous_state_column <- function(states, prototype, n, p_var) {
  values <- coerce_previous_state_values(states, prototype, p_var)
  rep(values, each = n)
}

normalize_previous_state_column <- function(values, prototype, p_var) {
  coerce_previous_state_values(values, prototype, p_var)
}

get_model_factor_levels <- function(model, varname) {
  model_chk <- if (inherits(model, "robcov_vglm")) {
    model$vglm_fit
  } else {
    model
  }

  if (inherits(model_chk, "vglm")) {
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

get_time_info <- function(model, data, time_var) {
  values <- if (!is.null(time_var) && time_var %in% names(data)) {
    data[[time_var]]
  } else {
    NULL
  }
  model_levels <- if (!is.null(model)) {
    get_model_factor_levels(model, time_var)
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

default_sop_times <- function(data, time_var, time_covariates, default) {
  if (default == "max_sequence") {
    return(seq_len(max(data[[time_var]], na.rm = TRUE)))
  }

  if (default == "fast") {
    if (!is.null(time_covariates)) {
      return(seq_len(nrow(time_covariates)))
    }
    if (!is.null(time_var) && time_var %in% names(data)) {
      return(sort(unique(data[[time_var]])))
    }
    return(1)
  }

  if (is.null(time_var) || !time_var %in% names(data)) {
    stop("`times` must be specified if `time_var` is not in data.")
  }

  sort(unique(data[[time_var]]))
}

coerce_visit_times <- function(times, time_info, time_var) {
  labels <- as.character(times)
  bad <- is.na(match(labels, time_info$levels))
  if (any(bad)) {
    bad_values <- unique(labels[bad])
    stop(
      "Requested visit values are not among fitted factor levels for `",
      time_var,
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
  time_var,
  time_covariates = NULL,
  default = c("unique", "max_sequence", "fast")
) {
  default <- match.arg(default)
  time_info <- get_time_info(model, data, time_var)

  if (is.null(times)) {
    times <- if (time_info$is_factor) {
      time_info$levels
    } else {
      default_sop_times(data, time_var, time_covariates, default)
    }
  }

  if (time_info$is_factor) {
    times <- coerce_visit_times(times, time_info, time_var)
  }

  if (!is.null(time_covariates) && nrow(time_covariates) != length(times)) {
    stop("`time_covariates` must have one row per requested time point.")
  }

  list(times = times, time_info = time_info)
}

assign_sop_time <- function(data, time_var, value, time_info) {
  if (isTRUE(time_info$is_factor)) {
    data[[time_var]] <- factor(
      as.character(value),
      levels = time_info$levels,
      ordered = isTRUE(time_info$ordered)
    )
  } else {
    data[[time_var]] <- value
  }
  data
}

validate_factor_gap <- function(gap_var, time_covariates, time_info) {
  if (is.null(gap_var) || !isTRUE(time_info$is_factor)) {
    return(invisible(NULL))
  }

  if (
    is.null(time_covariates) ||
      !gap_var %in% names(time_covariates) ||
      !is.numeric(time_covariates[[gap_var]])
  ) {
    stop(
      "Factor visit time with `gap_var` requires numeric gap_var values in ",
      "`time_covariates[[\"",
      gap_var,
      "\"]]`; elapsed gaps are not inferred from ",
      "factor visit labels."
    )
  }

  invisible(NULL)
}

assign_sop_gap <- function(
  data,
  gap_var,
  times,
  index,
  time_covariates,
  time_info
) {
  if (is.null(gap_var)) {
    return(data)
  }

  validate_factor_gap(gap_var, time_covariates, time_info)
  if (
    !isTRUE(time_info$is_factor) &&
      (is.null(time_covariates) || !gap_var %in% names(time_covariates))
  ) {
    data[[gap_var]] <- if (index == 1L) {
      times[index]
    } else {
      times[index] - times[index - 1L]
    }
  }

  data
}

assign_sop_t_covs <- function(data, time_covariates, index) {
  if (is.null(time_covariates)) {
    return(data)
  }

  for (nm in names(time_covariates)) {
    data[[nm]] <- time_covariates[index, nm]
  }
  data
}

assign_sop_visit <- function(
  data,
  time_var,
  times,
  index,
  time_covariates,
  gap_var,
  time_info
) {
  data <- assign_sop_time(data, time_var, times[index], time_info)
  data <- assign_sop_gap(
    data,
    gap_var,
    times,
    index,
    time_covariates,
    time_info
  )
  assign_sop_t_covs(data, time_covariates, index)
}
