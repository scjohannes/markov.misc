#' Calculate State Occupancy Probabilities
#'
#' Computes patient-level state occupancy probabilities from a supported
#' Markov transition model. Frequentist first-order models use the optimized
#' native engine; other supported workflows use the reference engine.
#'
#' @param model A fitted `vglm`, `robcov_vglm`, `orm`, or `blrm` Markov
#'   transition model.
#' @param newdata A data frame containing one baseline prediction row per
#'   patient or profile.
#' @param times Required visit-scale time points.
#' @param y_levels Ordered outcome-state levels.
#' @param absorb Optional absorbing state levels.
#' @param time_var Name of the model time variable.
#' @param p_var Name of the previous-state variable.
#' @param p2_var Optional second previous-state variable. A non-`NULL` value
#'   selects second-order recursion.
#' @param gap_var Optional elapsed-time-gap variable.
#' @param time_covariates Optional data frame with one row per requested visit
#'   containing precomputed time-varying covariates.
#' @param include_re For `blrm` fits, include fitted random-effect draws for
#'   known IDs.
#' @param id_var ID column used for `blrm` random effects.
#' @param n_draws Number of `blrm` posterior draws, or `NULL` for all draws.
#' @param seed Optional seed used when selecting posterior draws.
#' @param ... Reserved for internal posterior draw batching.
#' @return An array of state probabilities. Frequentist fits return
#'   `[patient, time, state]`; Bayesian fits return
#'   `[draw, patient, time, state]`.
#' @export
soprob_markov <- function(
  model,
  newdata,
  times = NULL,
  y_levels,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  include_re = FALSE,
  id_var = NULL,
  n_draws = 100L,
  seed = NULL,
  ...
) {
  dots <- list(...)
  model_fast <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  vglm_fast <- inherits(model_fast, "vglm") &&
    methods::is(model_fast, "vglm") &&
    isTRUE(tryCatch(
      !is.null(methods::slot(model_fast, "constraints")),
      error = function(e) FALSE
    ))
  orm_fast <- inherits(model_fast, "orm") && !inherits(model_fast, "blrm")
  use_fast <- (vglm_fast || orm_fast) &&
    !inherits(model_fast, "blrm") &&
    is.null(p2_var) &&
    length(dots) == 0L

  if (!use_fast) {
    return(soprob_markov_reference(
      model = model,
      newdata = newdata,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      include_re = include_re,
      id_var = id_var,
      n_draws = n_draws,
      seed = seed,
      ...
    ))
  }

  validate_markov_model(model_fast)
  validate_fast_previous_state_support(
    model_fast,
    p_var = p_var,
    y_levels = y_levels,
    absorb = absorb
  )
  plan <- compile_sop_execution_plan(
    model = model_fast,
    newdata = newdata,
    time_covariates = time_covariates,
    times = times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    gap_var = gap_var,
    builder = "batched",
    output = "array"
  )
  resolved_times <- plan$times
  out <- run_sop_execution_plan(plan, get_effective_coefs(model_fast))
  dimnames(out) <- list(
    rownames(newdata),
    resolved_times,
    as_state_labels(y_levels)
  )
  out
}

validate_fast_previous_state_support <- function(
  model,
  p_var,
  y_levels,
  absorb
) {
  model_data <- markov_model_data(model)
  if (
    is.null(model_data) ||
      !p_var %in% names(model_data) ||
      !is.factor(model_data[[p_var]])
  ) {
    return(invisible(NULL))
  }
  observed <- unique(as.character(model_data[[p_var]]))
  unsupported <- setdiff(as_state_labels(y_levels), observed)
  unsupported <- setdiff(unsupported, as_state_labels(absorb))
  if (length(unsupported) > 0L) {
    stop(
      "Model prediction requires transitions from state level(s) not ",
      "represented in the fitted transition data: ",
      paste(unsupported, collapse = ", "),
      ". If a level is absorbing, pass it via `absorb`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
