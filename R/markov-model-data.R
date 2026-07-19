# Markov model data metadata helpers.

#' Fit an rms Markov Model
#'
#' `orm_markov()` is a thin package-aware wrapper around [rms::orm()] for
#' Markov SOP workflows. It fits with `x = TRUE` and `y = TRUE` by default,
#' rejects offsets, stores separate likelihood, refit, and designated-start
#' profile data on the returned model, and computes cluster-robust standard
#' errors automatically when `id_var` is supplied.
#'
#' @inheritParams rms::orm
#' @param id_var Optional character scalar naming the patient or cluster ID
#'   column in `data`. When supplied, the returned object is wrapped with
#'   [rms::robcov()] using this column as the cluster. When omitted, the fit is
#'   returned without cluster-robust covariance and a warning is issued.
#' @param time_var Character scalar naming the modeled time column used to
#'   identify the designated starting-profile row.
#' @param start_time Optional starting-profile visit. `NULL` uses the first
#'   scheduled numeric time or the first factor level across the pre-NA fitting
#'   data. This is a cohort-wide visit, never each patient's first retained
#'   likelihood row.
#' @param origin_time Origin assigned to the previous-state value carried by the
#'   starting profile. The default is 0. This metadata may differ from
#'   `start_time`, for example when the row predicting day 1 carries a day-0
#'   state in `yprev`.
#' @param time_map Optional visit-to-real-time mapping retained with the fitted
#'   model. For factor time with `start_time = NULL`, the visit with the smallest
#'   mapped real time is designated when the mapping is complete.
#'
#' @details Starting profiles are retained before response-driven model-frame
#'   omission. Every fitted patient must have exactly one complete profile at
#'   the cohort-wide `start_time`, including ID, model predictors, time, and the
#'   previous state; the transition response on that row may be missing. A
#'   patient is included only when at least one usable likelihood transition is
#'   fitted somewhere. Patients with no fitted transition are excluded, and a
#'   fitted patient without a complete designated profile is an error. Later
#'   likelihood rows are not substituted for that profile.
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
#' sops(fit, times = 1:60, y_levels = factor(1:6), absorb = "6")
#' }
#'
#' @seealso [blrm_markov()], [vglm_markov()], [sops()], [avg_sops()]
#' @export
orm_markov <- function(
  formula,
  data,
  ...,
  id_var = NULL,
  time_var = "time",
  start_time = NULL,
  origin_time = 0,
  time_map = NULL,
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
  } else {
    warn_duplicate_markov_id_time(
      data,
      id_var,
      time_var = time_var,
      wrapper = "orm_markov()"
    )
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
  stored <- markov_prepare_stored_data(
    data = data,
    formula = formula,
    subset = dots$subset,
    eval_env = parent.frame(),
    fit_data = fit_data,
    id_var = id_var,
    time_var = time_var,
    start_time = start_time,
    origin_time = origin_time,
    time_map = time_map
  )
  fit <- markov_attach_model_data(
    fit,
    data = fit_data,
    id_var = id_var,
    refit_data = stored$refit_data,
    origin_data = stored$origin_data,
    origin_metadata = stored$origin_metadata
  )

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
    robust <- markov_attach_model_data(
      robust,
      data = fit_data,
      id_var = id_var,
      refit_data = stored$refit_data,
      origin_data = stored$origin_data,
      origin_metadata = stored$origin_metadata
    )
    return(robust)
  }

  fit
}

#' Fit an rmsb Markov Model
#'
#' `blrm_markov()` is a thin package-aware wrapper around [rmsb::blrm()] for
#' Markov SOP workflows. It fits with `x = TRUE` and `y = TRUE` by default,
#' rejects offsets, and stores separate likelihood, refit, and designated-start
#' profile data plus the optional patient ID on the returned model. Bayesian
#' fits do not receive
#' cluster-robust covariance; posterior uncertainty comes from the fitted
#' `blrm` draws.
#'
#' @inheritParams rmsb::blrm
#' @param id_var Optional character scalar naming the patient or cluster ID
#'   column in `data`. When supplied, SOP calls can reuse the stored full data
#'   and infer the ID column for baseline extraction and random-effect
#'   prediction.
#' @inheritParams orm_markov
#'
#' @details The designated-start profile and fitted-patient eligibility contract
#'   is the same as for [orm_markov()]. In particular, the transition response
#'   at the starting row may be missing when that patient contributes another
#'   usable likelihood transition.
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
  time_var = "time",
  start_time = NULL,
  origin_time = 0,
  time_map = NULL,
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
  } else {
    warn_duplicate_markov_id_time(
      data,
      id_var,
      time_var = time_var,
      wrapper = "blrm_markov()"
    )
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
  stored <- markov_prepare_stored_data(
    data = data,
    formula = formula,
    subset = dots$subset,
    eval_env = parent.frame(),
    fit_data = fit_data,
    id_var = id_var,
    time_var = time_var,
    start_time = start_time,
    origin_time = origin_time,
    time_map = time_map
  )
  markov_attach_model_data(
    fit,
    data = fit_data,
    id_var = id_var,
    refit_data = stored$refit_data,
    origin_data = stored$origin_data,
    origin_metadata = stored$origin_metadata
  )
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

warn_duplicate_markov_id_time <- function(
  data,
  id_var,
  time_var = "time",
  wrapper = "the Markov wrapper"
) {
  if (
    is.null(id_var) ||
      !id_var %in% names(data) ||
      !time_var %in% names(data)
  ) {
    return(invisible(NULL))
  }

  key <- data[c(id_var, time_var)]
  if (anyDuplicated(key) == 0L) {
    return(invisible(NULL))
  }

  warning(
    wrapper,
    " received data with duplicate `",
    id_var,
    "` x `",
    time_var,
    "` combinations. This is usually accidental in Markov transition data; ",
    "automatic SOP prediction from stored data will use one prediction row ",
    "per `",
    id_var,
    "`.",
    call. = FALSE
  )

  invisible(NULL)
}

markov_attach_model_data <- function(
  model,
  data,
  id_var = NULL,
  refit_data = NULL,
  origin_data = NULL,
  origin_metadata = NULL
) {
  if (!is.null(data)) {
    attr(model, "markov_data") <- data
  }
  if (!is.null(id_var)) {
    attr(model, "markov_id_var") <- id_var
  }
  if (!is.null(refit_data)) {
    attr(model, "markov_refit_data") <- refit_data
  }
  if (!is.null(origin_data)) {
    attr(model, "markov_origin_data") <- origin_data
  }
  if (!is.null(origin_metadata)) {
    attr(model, "markov_origin_metadata") <- origin_metadata
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

markov_model_metadata_attr <- function(model, name) {
  value <- attr(model, name, exact = TRUE)
  if (!is.null(value)) {
    return(value)
  }

  if (inherits(model, "robcov_vglm") && !is.null(model$vglm_fit)) {
    return(attr(model$vglm_fit, name, exact = TRUE))
  }

  NULL
}

markov_model_refit_data <- function(model) {
  markov_model_metadata_attr(model, "markov_refit_data") %||%
    markov_model_data(model)
}

markov_model_origin_metadata <- function(model) {
  markov_model_metadata_attr(model, "markov_origin_metadata")
}

markov_format_patient_ids <- function(ids, limit = 5L) {
  ids <- unique(as.character(ids))
  paste0(
    paste(utils::head(ids, limit), collapse = ", "),
    if (length(ids) > limit) " ..." else ""
  )
}

markov_validate_origin_profiles <- function(model, required_variables = NULL) {
  profiles <- markov_model_metadata_attr(model, "markov_origin_data")
  metadata <- markov_model_origin_metadata(model)
  if (is.null(profiles) || is.null(metadata)) {
    stop(
      "Provide newdata, or fit with `orm_markov()`, `blrm_markov()`, or ",
      "`vglm_markov()` so ",
      "automatic SOP prediction has stored starting-profile metadata.",
      call. = FALSE
    )
  }

  id_var <- metadata$id_var
  fitted_ids <- as.character(metadata$fitted_ids)
  if (
    is.null(id_var) ||
      !id_var %in% names(profiles) ||
      anyNA(fitted_ids) ||
      any(fitted_ids == "") ||
      anyDuplicated(fitted_ids)
  ) {
    stop("Stored starting-profile patient metadata is invalid.", call. = FALSE)
  }

  profile_ids <- as.character(profiles[[id_var]])
  missing_ids <- setdiff(fitted_ids, profile_ids)
  if (length(missing_ids)) {
    stop(
      "Fitted patient(s) lack a row at the designated starting visit `",
      as.character(metadata$start_time),
      "`: ",
      markov_format_patient_ids(missing_ids),
      ". A later likelihood row is not substituted for the starting profile.",
      call. = FALSE
    )
  }

  duplicated_ids <- unique(profile_ids[
    duplicated(profile_ids) |
      duplicated(profile_ids, fromLast = TRUE)
  ])
  if (length(duplicated_ids)) {
    stop(
      "Fitted patient(s) have multiple rows at the designated starting visit ",
      "`",
      as.character(metadata$start_time),
      "`: ",
      markov_format_patient_ids(duplicated_ids),
      ". Exactly one starting profile is required per fitted patient.",
      call. = FALSE
    )
  }

  required <- unique(c(metadata$required_variables, required_variables))
  missing_variables <- setdiff(required, names(profiles))
  if (length(missing_variables)) {
    stop(
      "Stored starting profiles are missing required prediction variables: ",
      paste(missing_variables, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  incomplete <- !stats::complete.cases(profiles[, required, drop = FALSE])
  if (any(incomplete)) {
    stop(
      "Fitted patient(s) have incomplete starting-profile predictors: ",
      markov_format_patient_ids(profile_ids[incomplete]),
      ". The transition response is not required, but ID, modeled predictors, ",
      "time/origin metadata, and the previous state must be complete.",
      call. = FALSE
    )
  }

  fitting_data <- markov_model_data(model)
  factor_predictors <- intersect(
    metadata$predictor_variables %||% character(),
    names(profiles)[vapply(profiles, is.factor, logical(1))]
  )
  for (variable in factor_predictors) {
    if (is.null(fitting_data) || !variable %in% names(fitting_data)) {
      next
    }
    supported <- unique(as.character(fitting_data[[variable]]))
    unsupported <- !as.character(profiles[[variable]]) %in% supported
    if (any(unsupported)) {
      stop(
        "Fitted patient(s) have unsupported factor value(s) for starting ",
        "predictor `",
        variable,
        "`: ",
        markov_format_patient_ids(profile_ids[unsupported]),
        ".",
        call. = FALSE
      )
    }
  }

  unexpected_ids <- setdiff(profile_ids, fitted_ids)
  if (length(unexpected_ids)) {
    stop(
      "Stored starting profiles contain patient(s) with no usable likelihood ",
      "transition: ",
      markov_format_patient_ids(unexpected_ids),
      ".",
      call. = FALSE
    )
  }

  profiles <- profiles[match(fitted_ids, profile_ids), , drop = FALSE]
  rownames(profiles) <- NULL
  profiles
}

markov_subset_source_data <- function(data, subset, eval_env) {
  if (is.null(subset)) {
    return(list(data = data, rows = seq_len(nrow(data))))
  }

  selection <- eval(subset, data, eval_env)
  rows <- seq_len(nrow(data))[selection]
  rows <- rows[!is.na(rows)]
  list(data = data[rows, , drop = FALSE], rows = rows)
}

markov_start_time <- function(time, start_time, time_map) {
  if (!is.null(start_time)) {
    if (length(start_time) != 1L || is.na(start_time)) {
      stop("`start_time` must contain one non-missing value.", call. = FALSE)
    }
    return(start_time)
  }

  if (is.factor(time)) {
    levels <- levels(time)
    if (!is.null(time_map)) {
      map <- standardize_time_map(time_map)
      mapped <- map$real_time[match(levels, map$visit)]
      if (!anyNA(mapped)) {
        return(levels[which.min(mapped)])
      }
    }
    return(levels[[1L]])
  }

  finite <- time[!is.na(time)]
  if (!length(finite)) {
    stop("The stored time variable has no non-missing scheduled values.")
  }
  min(finite)
}

markov_prepare_stored_data <- function(
  data,
  formula,
  subset,
  eval_env,
  fit_data,
  id_var,
  time_var,
  start_time,
  origin_time,
  time_map
) {
  if (!is.data.frame(data)) {
    return(list(
      refit_data = NULL,
      origin_data = NULL,
      origin_metadata = NULL
    ))
  }
  source <- markov_subset_source_data(data, subset, eval_env)
  refit_data <- source$data
  if (is.null(id_var)) {
    return(list(
      refit_data = refit_data,
      origin_data = NULL,
      origin_metadata = NULL
    ))
  }
  validate_markov_id_var(time_var, refit_data, "data")
  validate_markov_id_var(id_var, fit_data, "fitted data")
  if (length(origin_time) != 1L || is.na(origin_time)) {
    stop("`origin_time` must contain one non-missing value.", call. = FALSE)
  }

  fitted_ids <- unique(as.character(fit_data[[id_var]]))
  if (anyNA(fitted_ids) || any(fitted_ids == "")) {
    stop("Fitted patient IDs must be non-missing and non-empty.", call. = FALSE)
  }
  scheduled_start <- markov_start_time(
    refit_data[[time_var]],
    start_time,
    time_map
  )
  at_start <- as.character(refit_data[[time_var]]) ==
    as.character(scheduled_start)
  in_fit <- as.character(refit_data[[id_var]]) %in% fitted_ids
  origin_data <- refit_data[at_start & in_fit, , drop = FALSE]
  origin_data$.markov_source_row <- source$rows[at_start & in_fit]

  predictor_terms <- stats::delete.response(stats::terms(formula, data = data))
  predictor_variables <- all.vars(predictor_terms)
  required <- unique(c(
    id_var,
    predictor_variables,
    time_var
  ))

  list(
    refit_data = refit_data,
    origin_data = origin_data,
    origin_metadata = list(
      id_var = id_var,
      time_var = time_var,
      start_time = scheduled_start,
      origin_time = origin_time,
      time_map = time_map,
      required_variables = required,
      predictor_variables = predictor_variables,
      fitted_ids = fitted_ids,
      factor_levels = lapply(
        refit_data[vapply(refit_data, is.factor, logical(1))],
        levels
      )
    )
  )
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

resolve_markov_source_data <- function(
  model,
  newdata,
  refit_data = NULL,
  time_var = NULL,
  p_var = NULL
) {
  newdata_supplied <- !is.null(newdata)
  refit_data <- refit_data %||% markov_model_refit_data(model)
  source_data <- if (newdata_supplied) {
    newdata
  } else {
    metadata <- markov_model_origin_metadata(model)
    if (
      !is.null(time_var) &&
        !is.null(metadata$time_var) &&
        !identical(time_var, metadata$time_var)
    ) {
      stop(
        "Automatic SOP prediction requested `time_var = \"",
        time_var,
        "\"`, but the fitted wrapper stored starting profiles using `",
        metadata$time_var,
        "`. Refit with the intended `time_var` or supply explicit `newdata`.",
        call. = FALSE
      )
    }
    markov_validate_origin_profiles(
      model,
      required_variables = c(time_var, p_var)
    )
  }

  if (is.null(source_data)) {
    stop(
      "Provide `newdata`, or fit with `orm_markov()`, `blrm_markov()`, or ",
      "`vglm_markov()` so ",
      "the model stores its designated starting profiles."
    )
  }
  if (!is.data.frame(source_data)) {
    stop("`newdata` must be a data frame.")
  }

  list(
    source_data = source_data,
    refit_data = refit_data,
    newdata_supplied = newdata_supplied
  )
}

resolve_markov_prediction_data <- function(
  data,
  id_var,
  time_var,
  data_label = "newdata"
) {
  if (is.null(id_var) || !id_var %in% names(data)) {
    return(data)
  }

  id <- data[[id_var]]
  if (!anyDuplicated(id)) {
    return(data)
  }

  if (time_var %in% names(data)) {
    split_idx <- split(seq_len(nrow(data)), id)
    baseline_idx <- vapply(
      split_idx,
      function(idx) {
        time_values <- data[[time_var]][idx]
        idx[order(time_values, na.last = TRUE)[1]]
      },
      integer(1)
    )
    return(data[baseline_idx, , drop = FALSE])
  }

  warning(
    data_label,
    " contains repeated IDs but no `",
    time_var,
    "` column; using the first row per ID as the prediction row.",
    call. = FALSE
  )
  data[!duplicated(id), , drop = FALSE]
}

ensure_markov_rowid <- function(data) {
  data$rowid <- seq_len(nrow(data))
  data
}
