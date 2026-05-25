# Legacy SOP standardization wrapper.

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
