# Markov model data metadata helpers.

#' Fit an rms Markov Model
#'
#' `orm_markov()` is a thin package-aware wrapper around [rms::orm()] for
#' Markov SOP workflows. It fits with `x = TRUE` and `y = TRUE` by default,
#' rejects offsets, stores the row-aligned fitting data on the returned model,
#' and computes cluster-robust standard errors automatically when `id_var` is
#' supplied.
#'
#' @inheritParams rms::orm
#' @param id_var Optional character scalar naming the patient or cluster ID
#'   column in `data`. When supplied, the returned object is wrapped with
#'   [rms::robcov()] using this column as the cluster. When omitted, the fit is
#'   returned without cluster-robust covariance and a warning is issued.
#'
#' @return A fitted `orm` object. If `id_var` is supplied, the object contains
#'   the robust covariance computed by [rms::robcov()].
#'
#' @examples
#' \dontrun{
#' dd <- rms::datadist(data)
#' options(datadist = "dd")
#'
#' fit <- orm_markov(
#'   ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
#'   data = data,
#'   id_var = "id"
#' )
#'
#' sops(fit, times = 1:60, ylevels = factor(1:6), absorb = "6")
#' }
#'
#' @seealso [blrm_markov()], [vglm_markov()], [sops()], [avg_sops()]
#' @export
orm_markov <- function(
  formula,
  data,
  ...,
  id_var = NULL,
  x = TRUE,
  y = TRUE
) {
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be supplied as a data frame.")
  }

  mt <- stats::terms(formula, specials = "offset", data = data)
  if (terms_has_offset(mt)) {
    stop_unsupported_offset()
  }

  orm_call <- match.call(expand.dots = FALSE)
  dots <- orm_call$...
  if (length(dots) > 0 && "offset" %in% names(dots)) {
    stop_unsupported_offset()
  }

  id_var <- validate_markov_id_var(id_var, data, "data")
  if (is.null(id_var)) {
    warn_missing_markov_id_var("orm_markov()")
  }

  fit <- rms::orm(
    formula = formula,
    data = data,
    ...,
    x = x,
    y = y
  )
  fit <- markov_normalize_rms_call(
    fit,
    fun = quote(rms::orm),
    formula = formula,
    data = orm_call$data,
    dots = dots,
    x = x,
    y = y
  )
  fit_frame <- markov_rms_model_frame(orm_call, dots, parent.frame())
  fit_data <- markov_align_model_data(data, fit_frame)
  fit <- markov_attach_model_data(fit, fit_data, id_var)

  if (!is.null(id_var)) {
    robust <- rms::robcov(fit, cluster = fit_data[[id_var]])
    robust <- markov_normalize_rms_call(
      robust,
      fun = quote(rms::orm),
      formula = formula,
      data = orm_call$data,
      dots = dots,
      x = x,
      y = y
    )
    robust <- markov_attach_model_data(robust, fit_data, id_var)
    return(robust)
  }

  fit
}

#' Fit an rmsb Markov Model
#'
#' `blrm_markov()` is a thin package-aware wrapper around [rmsb::blrm()] for
#' Markov SOP workflows. It fits with `x = TRUE` and `y = TRUE` by default,
#' rejects offsets, and stores the row-aligned fitting data and optional patient
#' ID column on the returned model. Bayesian fits do not receive
#' cluster-robust covariance; posterior uncertainty comes from the fitted
#' `blrm` draws.
#'
#' @inheritParams rmsb::blrm
#' @param id_var Optional character scalar naming the patient or cluster ID
#'   column in `data`. When supplied, SOP calls can reuse the stored full data
#'   and infer the ID column for baseline extraction and random-effect
#'   prediction.
#'
#' @return A fitted `blrm` object with `markov.misc` stored-data metadata.
#'
#' @examples
#' \dontrun{
#' dd <- rms::datadist(data)
#' options(datadist = "dd")
#'
#' fit <- blrm_markov(
#'   y ~ yprev + rms::rcs(time, 4) * tx,
#'   ppo = ~ time,
#'   cppo = function(y) y,
#'   data = data,
#'   id_var = "id"
#' )
#'
#' avg_sops(fit, variables = list(tx = c(0, 1)), times = 1:60)
#' }
#'
#' @seealso [orm_markov()], [vglm_markov()], [sops()], [avg_sops()]
#' @export
blrm_markov <- function(
  formula,
  ppo = NULL,
  cppo = NULL,
  data,
  ...,
  id_var = NULL,
  x = TRUE,
  y = TRUE
) {
  if (!requireNamespace("rmsb", quietly = TRUE)) {
    stop("Package 'rmsb' is required for `blrm_markov()`.")
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be supplied as a data frame.")
  }

  mt <- stats::terms(formula, specials = "offset", data = data)
  if (terms_has_offset(mt)) {
    stop_unsupported_offset()
  }

  blrm_call <- match.call(expand.dots = FALSE)
  dots <- blrm_call$...
  if (length(dots) > 0 && "offset" %in% names(dots)) {
    stop_unsupported_offset()
  }

  id_var <- validate_markov_id_var(id_var, data, "data")
  if (is.null(id_var)) {
    warn_missing_markov_id_var("blrm_markov()")
  }

  fit <- rmsb::blrm(
    formula = formula,
    ppo = ppo,
    cppo = cppo,
    data = data,
    ...,
    x = x,
    y = y
  )
  fit <- markov_normalize_rms_call(
    fit,
    fun = quote(rmsb::blrm),
    formula = formula,
    data = blrm_call$data,
    dots = dots,
    ppo = ppo,
    cppo = cppo,
    x = x,
    y = y
  )
  fit_frame <- markov_rms_model_frame(blrm_call, dots, parent.frame())
  fit_data <- markov_align_model_data(data, fit_frame)
  markov_attach_model_data(fit, fit_data, id_var)
}

markov_normalize_rms_call <- function(
  fit,
  fun,
  formula,
  data,
  dots,
  ...,
  ppo = NULL,
  cppo = NULL
) {
  call_args <- c(
    list(formula = formula),
    if (!is.null(ppo)) list(ppo = ppo) else NULL,
    if (!is.null(cppo)) list(cppo = cppo) else NULL,
    list(data = data),
    as.list(dots),
    list(...)
  )
  fit$call <- as.call(c(list(fun), call_args))
  fit
}

markov_rms_model_frame <- function(rms_call, dots, eval_env) {
  mf <- rms_call[c(1L, match(c("formula", "data"), names(rms_call), 0L))]
  for (nm in intersect(c("subset", "weights", "na.action"), names(dots))) {
    mf[[nm]] <- dots[[nm]]
  }
  if (is.null(mf$na.action)) {
    mf$na.action <- utils::getFromNamespace("na.delete", "Hmisc")
  }
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  eval(mf, eval_env)
}

validate_markov_id_var <- function(id_var, data, data_arg = "data") {
  if (is.null(id_var)) {
    return(NULL)
  }
  if (!is.character(id_var) || length(id_var) != 1 || is.na(id_var)) {
    stop("`id_var` must be a single column name.")
  }
  if (!id_var %in% names(data)) {
    stop("ID variable '", id_var, "' not found in ", data_arg, ".")
  }
  id_var
}

warn_missing_markov_id_var <- function(wrapper) {
  warning(
    "`id_var` was not supplied to ",
    wrapper,
    "; stored patient-ID metadata, cluster-robust standard errors for ",
    "frequentist wrappers, and automatic FWB refits will not be available ",
    "from the fitted object.",
    call. = FALSE
  )
}

markov_attach_model_data <- function(model, data, id_var = NULL) {
  if (!is.null(data)) {
    attr(model, "markov_data") <- data
  }
  if (!is.null(id_var)) {
    attr(model, "markov_id_var") <- id_var
  }
  model
}

markov_model_data <- function(model) {
  data <- attr(model, "markov_data", exact = TRUE)
  if (!is.null(data)) {
    return(data)
  }

  if (inherits(model, "robcov_vglm") && !is.null(model$vglm_fit)) {
    data <- attr(model$vglm_fit, "markov_data", exact = TRUE)
    if (!is.null(data)) {
      return(data)
    }
  }

  NULL
}

markov_model_id_var <- function(model, id_var = NULL) {
  if (!is.null(id_var)) {
    return(id_var)
  }

  out <- attr(model, "markov_id_var", exact = TRUE)
  if (!is.null(out)) {
    return(out)
  }

  if (inherits(model, "robcov_vglm") && !is.null(model$vglm_fit)) {
    out <- attr(model$vglm_fit, "markov_id_var", exact = TRUE)
    if (!is.null(out)) {
      return(out)
    }
  }

  NULL
}

markov_align_model_data <- function(data, model_or_frame) {
  if (!is.data.frame(data)) {
    return(NULL)
  }

  rn <- NULL
  if (is.data.frame(model_or_frame)) {
    rn <- rownames(model_or_frame)
  } else {
    mf <- tryCatch(stats::model.frame(model_or_frame), error = function(e) NULL)
    if (!is.null(mf)) {
      rn <- rownames(mf)
    }
  }

  if (!is.null(rn) && length(rn) > 0 && all(rn %in% rownames(data))) {
    return(data[rn, , drop = FALSE])
  }

  model_n <- tryCatch(stats::nobs(model_or_frame), error = function(e) NULL)
  if (!is.null(model_n) && nrow(data) != model_n) {
    model_vars <- tryCatch(
      all.vars(stats::formula(model_or_frame)),
      error = function(e) character()
    )
    model_vars <- intersect(model_vars, names(data))
    if (length(model_vars) > 0) {
      keep <- stats::complete.cases(data[, model_vars, drop = FALSE])
      if (sum(keep) == model_n) {
        return(data[keep, , drop = FALSE])
      }
    }
  }

  data
}

resolve_markov_source_data <- function(model, newdata, refit_data = NULL) {
  refit_data <- refit_data %||% markov_model_data(model)
  source_data <- newdata %||% refit_data

  if (is.null(source_data)) {
    stop(
      "Provide newdata, or fit with `orm_markov()`, `blrm_markov()`, or ",
      "`vglm_markov(id_var = ...)` so the model stores its original data."
    )
  }
  if (!is.data.frame(source_data)) {
    stop("`newdata` must be a data frame.")
  }

  list(
    source_data = source_data,
    refit_data = refit_data
  )
}

resolve_markov_prediction_data <- function(
  data,
  id_var,
  tvarname,
  data_label = "newdata"
) {
  if (is.null(id_var) || !id_var %in% names(data)) {
    return(data)
  }

  id <- data[[id_var]]
  if (!anyDuplicated(id)) {
    return(data)
  }

  if (tvarname %in% names(data)) {
    split_idx <- split(seq_len(nrow(data)), id)
    has_longitudinal_repeats <- any(vapply(
      split_idx,
      function(idx) {
        length(unique(data[[tvarname]][idx])) > 1
      },
      logical(1)
    ))
    if (!has_longitudinal_repeats) {
      return(data)
    }

    baseline_idx <- vapply(
      split_idx,
      function(idx) {
        time_values <- data[[tvarname]][idx]
        idx[order(time_values, na.last = TRUE)[1]]
      },
      integer(1)
    )
    return(data[baseline_idx, , drop = FALSE])
  }

  warning(
    data_label,
    " contains repeated IDs but no `",
    tvarname,
    "` column; using the first row per ID as the prediction row.",
    call. = FALSE
  )
  data[!duplicated(id), , drop = FALSE]
}

ensure_markov_rowid <- function(data) {
  if (!"rowid" %in% names(data)) {
    data$rowid <- seq_len(nrow(data))
  }
  data
}
