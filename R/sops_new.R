#' Calculate Individual State Occupation Probabilities
#'
#' Computes individual-level state occupation probabilities (SOPs) for each row
#' in a dataset. This function serves as the foundation for `avg_sops()` and
#' follows a `marginaleffects`-style workflow.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`).
#' @param newdata Optional. A data frame of new data for prediction. If NULL,
#'   uses the data used to fit the model.
#' @param times A numeric vector of time points to estimate.
#' @param ylevels A vector of state levels. If NULL, attempts to infer from model.
#' @param absorb Character or integer. The absorbing state.
#' @param tvarname Name of the time variable in the model.
#' @param pvarname Name of the previous state variable in the model.
#' @param gap Name of the time gap variable (if used).
#' @param t_covs Optional time-varying covariate lookup table.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class `markov_sops` containing:
#'   \item{rowid}{Row identifier from newdata}
#'   \item{time}{Time point}
#'   \item{state}{State name}
#'   \item{estimate}{Probability of being in the state}
#'   Plus all columns from `newdata`.
#'
#' @details
#' This function wraps `soprob_markov()` and converts its array output to a tidy
#' data frame. The output contains one row per patient-time-state combination.
#'
#' For computing marginal/standardized SOPs (G-computation), use `avg_sops()`
#' instead, which creates counterfactual datasets and averages over individuals.
#'
#' @seealso [avg_sops()] for marginal SOPs, [soprob_markov()] for the underlying
#'   computation.
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
  gap = NULL,
  t_covs = NULL,
  ...
) {
  # --- 1. Setup & Defaults ---
  if (is.null(newdata)) {
    newdata <- model$x
    if (is.null(newdata)) {
      newdata <- tryCatch(model.frame(model), error = function(e) NULL)
    }
    if (is.null(newdata)) {
      stop("Could not extract data from model. Please provide `newdata`.")
    }
  }

  # Ensure rowid exists for tracking
  if (!"rowid" %in% names(newdata)) {
    newdata$rowid <- seq_len(nrow(newdata))
  }

  if (is.null(times)) {
    if (is.null(tvarname) || !tvarname %in% names(newdata)) {
      stop("`times` must be specified if `tvarname` is not in data.")
    }
    times <- sort(unique(newdata[[tvarname]]))
  }

  if (is.null(ylevels)) {
    # Try to infer from model
    if (inherits(model, "vglm")) {
      # VGAM stores response levels
      ylevels <- model@extra$colnames.y
      if (is.null(ylevels)) ylevels <- 1:6
    } else if (inherits(model, "robcov_vglm")) {
      # robcov_vglm stores extra slot in list
      ylevels <- model$extra$colnames.y
      if (is.null(ylevels)) ylevels <- 1:6
    } else if (inherits(model, "orm")) {
      # rms orm stores levels
      ylevels <- model$yunique
      if (is.null(ylevels)) ylevels <- 1:6
    } else {
      ylevels <- 1:6
    }
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
    gap = gap,
    t_covs = t_covs
  )

  # --- 3. Tidy the Output ---
  # Array dims: [n_pat, n_times, n_states]
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]

  # Flatten array (R is column-major: iterates dim1 fastest)
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

  # Store attributes for downstream use
  attr(result, "model") <- model
  attr(result, "tvarname") <- tvarname
  attr(result, "pvarname") <- pvarname
  attr(result, "ylevels") <- ylevels
  attr(result, "absorb") <- absorb
  attr(result, "gap") <- gap
  attr(result, "t_covs") <- t_covs
  attr(result, "call_args") <- list(
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    gap = gap,
    t_covs = t_covs
  )
  # Store original newdata for inferences() (needed for simulation-based inference)
  attr(result, "newdata_orig") <- newdata

  class(result) <- c("markov_sops", class(result))
  return(result)
}


#' Calculate Averaged State Occupation Probabilities (Marginal Effects)
#'
#' Computes standardized (marginal) state occupation probabilities using
#' G-computation. Creates counterfactual cohorts by setting all individuals
#' to each level of the treatment variable and averaging over the covariate
#' distribution.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`).
#' @param newdata Data frame for prediction. For **simulation inference**, pass
#'   baseline data only (one row per patient). For **bootstrap inference**, you
#'   must pass the full longitudinal dataset (all time points) since the model
#'   needs to be refit on bootstrap samples. If NULL, extracts from model.
#' @param variables A named list specifying the variable(s) to standardize over.
#'   E.g., `list(tx = c(0, 1))` creates counterfactual datasets for treatment
#'   and control. **Required** - G-computation needs a treatment variable.
#' @param by Optional character vector of additional variables to group by.
#' @param times Numeric vector of time points. If NULL, inferred from data.
#' @param id_var Name of the patient ID variable (default "id"). Required for
#'   bootstrap inference.
#' @param ... Additional arguments passed to `sops()` (e.g., `ylevels`, `absorb`,
#'   `tvarname`, `pvarname`, `t_covs`).
#'
#' @return A data frame of class `markov_avg_sops` with columns:
#'   \item{time}{Time point}
#'   \item{state}{State level}
#'   \item{(variables)}{Value of standardization variable (e.g., tx)}
#'   \item{estimate}{Average probability across individuals}
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
#' ) |> inferences(method = "simulation", n_sim = 1000)
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
  id_var = "id",
  ...
) {
  # --- 1. Input Validation ---
  if (is.null(variables)) {
    stop(
      "`variables` is required for G-computation. ",
      "Specify the treatment variable, e.g., `variables = list(tx = c(0, 1))`."
    )
  }

  if (is.null(newdata)) {
    newdata <- model$x
    if (is.null(newdata)) {
      newdata <- tryCatch(model.frame(model), error = function(e) NULL)
    }
    if (is.null(newdata)) {
      stop("Provide newdata or ensure model stores data (x = TRUE).")
    }
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

  # Ensure rowid exists
  if (!"rowid" %in% names(baseline_data)) {
    baseline_data$rowid <- seq_len(nrow(baseline_data))
  }

  # --- 3. Create Counterfactual Datasets ---
  # For each combination in variables, create a copy of baseline_data with
  # the variable(s) set to that value
  grid <- do.call(expand.grid, variables)

  expanded_data_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt_copy <- baseline_data
    for (v in names(grid)) {
      dt_copy[[v]] <- grid[i, v]
    }
    expanded_data_list[[i]] <- dt_copy
  }
  newdata_expanded <- do.call(rbind, expanded_data_list)

  # --- 4. Compute Individual SOPs ---
  sops_ind <- sops(model, newdata = newdata_expanded, times = times, ...)

  # --- 5. Aggregate (Marginalize) ---
  # Group by time, state, and the variables used for standardization
  group_cols <- c("time", "state", names(variables))

  # Add optional 'by' variables
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }

  # Validate grouping columns exist
  missing_groups <- setdiff(group_cols, names(sops_ind))
  if (length(missing_groups) > 0) {
    stop("Grouping variables missing: ", paste(missing_groups, collapse = ", "))
  }

  # Aggregate using formula interface
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
  attr(result, "ylevels") <- attr(sops_ind, "ylevels")
  attr(result, "absorb") <- attr(sops_ind, "absorb")
  attr(result, "gap") <- attr(sops_ind, "gap")
  attr(result, "t_covs") <- attr(sops_ind, "t_covs")

  # Specific attributes for avg_sops/inferences
  attr(result, "avg_args") <- list(
    variables = variables,
    by = by,
    times = times,
    id_var = id_var
  )
  # Store ORIGINAL newdata for bootstrap (not the expanded counterfactual)
  attr(result, "newdata_orig") <- newdata

  class(result) <- c("markov_avg_sops", class(result))
  return(result)
}


#' Inference for State Occupation Probabilities
#'
#' Adds confidence intervals to SOP objects using simulation-based (MVN) or
#' bootstrap methods. The default method is simulation, which draws coefficient
#' vectors from a multivariate normal distribution centered at the point
#' estimates. This is typically much faster than bootstrap as it does not
#' require refitting the model.
#'
#' @param object A `markov_avg_sops` object from `avg_sops()` or a
#'   `markov_sops` object from `sops()`.
#' @param method Character. Inference method:
#'   \itemize{
#'     \item `"simulation"` (default): Draws from MVN(beta_hat, Sigma). Fast,
#'       does not refit models. Works for both individual and averaged SOPs.
#'     \item `"bootstrap"`: Resamples patients with replacement, refits model.
#'       More robust but slower. Only works for `markov_avg_sops` objects.
#'   }
#' @param n_sim Number of simulation draws (for simulation) or bootstrap
#'   iterations (for bootstrap). Default is 1000.
#' @param n_boot Deprecated alias for `n_sim` when using bootstrap method.
#'   For backward compatibility.
#' @param vcov Optional custom variance-covariance matrix. If provided,
#'   overrides the vcov extracted from the model. This is useful for sensitivity
#'   analyses with inflated variance.
#' @param parallel Logical. Use parallel processing? Default is FALSE.
#' @param workers Number of parallel workers. If NULL, uses detectCores() - 1.
#' @param conf_level Confidence level for intervals (default 0.95).
#' @param conf_type Type of confidence interval (simulation method only):
#'   \itemize{
#'     \item `"perc"` (default): Percentile-based intervals from the simulation
#'       distribution.
#'     \item `"wald"`: Uses simulation standard errors with normal quantiles.
#'   }
#' @param return_draws Logical. If TRUE, stores all individual simulation/bootstrap
#'   draws as an attribute. Extract with `get_draws()`. Default is FALSE.
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
#' ## Simulation Method (Default)
#'
#' The simulation method (Krinsky & Robb, 1986) proceeds as follows:
#' 1. Extract coefficient vector beta_hat and (robust) variance-covariance Sigma
#' 2. Draw n_sim coefficient vectors from MVN(beta_hat, Sigma)
#' 3. For each draw, replace model coefficients and compute SOPs
#' 4. Compute confidence intervals from the empirical distribution
#'
#' **Advantages:**
#' - Much faster than bootstrap (no model refitting)
#' - Works for both individual-level (`sops()`) and averaged (`avg_sops()`) SOPs
#' - Supports robust (cluster) variance estimation via `cluster` argument
#'
#' **When to use:**
#' - Default choice for most applications
#' - When speed is important
#' - When the normal approximation for coefficients is reasonable
#'
#' ## Bootstrap Method
#'
#' The bootstrap method resamples patients and refits the model:
#' 1. Sample patient IDs with replacement
#' 2. Refit model on bootstrap sample (handles missing states)
#' 3. Compute SOPs using G-computation
#' 4. Compute percentile-based confidence intervals
#'
#' **Advantages:**
#' - Makes no distributional assumptions about coefficients
#' - Properly handles rare states that may be missing from samples
#'
#' **When to use:**
#' - When the normal approximation may be poor (small samples, near boundaries)
#' - When you need to verify simulation-based results
#'
#' **Important:** Bootstrap requires the full longitudinal dataset (all time
#' points) in the original `newdata` passed to `avg_sops()`, not just baseline
#' data. This is because the model must be refit on each bootstrap sample.
#'
#' ## Cluster-Robust Variance
#'
#' For longitudinal data with repeated measurements per patient, you should
#' use cluster-robust standard errors. To do this, wrap your model with
#' `robcov_vglm()` (for vglm) or `rms::robcov()` (for orm) **before** passing
#' it to `sops()` or `avg_sops()`:
#'
#' ```r
#' # Fit model
#' fit <- vglm(y ~ time + tx + yprev, family = cumulative(...), data = data)
#'
#' # Wrap with cluster-robust vcov
#' fit_robust <- robcov_vglm(fit, cluster = data$id)
#'
#' # Use in workflow - vcov is automatically extracted
#' avg_sops(model = fit_robust, ...) |> inferences()
#' ```
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
#' fit <- vglm(y ~ time + tx + yprev, family = cumulative(...), data = data)
#' fit_robust <- robcov_vglm(fit, cluster = data$id)
#'
#' # Step 2a: Simulation inference (use baseline data only)
#' baseline_data <- data |> filter(time == 1)
#' result_sim <- avg_sops(
#'   model = fit_robust,
#'   newdata = baseline_data,  # Baseline only for simulation
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6,
#'   id_var = "id"
#' ) |>
#'   inferences(method = "simulation", n_sim = 1000)
#'
#' # Step 2b: Bootstrap inference (must use full data)
#' result_boot <- avg_sops(
#'   model = fit_robust,
#'   newdata = data,  # Full longitudinal data for bootstrap
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6,
#'   id_var = "id"
#' ) |>
#'   inferences(method = "bootstrap", n_sim = 500, parallel = TRUE)
#'
#' # Individual-level SOPs with simulation inference
#' ind_result <- sops(model = fit_robust, newdata = test_data, times = 1:30) |>
#'   inferences(n_sim = 500)
#'
#' # Extract draws for custom analyses
#' result_draws <- avg_sops(model = fit_robust, ...) |>
#'   inferences(return_draws = TRUE)
#' draws <- get_draws(result_draws)
#'
#' # Compute treatment effect on time in state
#' library(dplyr)
#' library(tidyr)
#' treatment_effect <- draws |>
#'   filter(state == 1) |>
#'   group_by(draw_id, tx) |>
#'   summarise(total_time = sum(estimate), .groups = "drop") |>
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
  n_sim = 1000,
  n_boot = NULL,
  vcov = NULL,
  parallel = FALSE,
  workers = NULL,
  conf_level = 0.95,
  conf_type = "perc",
  return_draws = FALSE,
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

  method <- match.arg(method, choices = c("simulation", "bootstrap"))

  # Handle backward compatibility: n_boot overrides n_sim for bootstrap
  if (!is.null(n_boot) && method == "bootstrap") {
    n_sim <- n_boot
  }

  # --- Dispatch to Method-Specific Implementation ---
  if (method == "simulation") {
    inferences_simulation(
      object = object,
      n_sim = n_sim,
      vcov = vcov,
      conf_level = conf_level,
      conf_type = conf_type,
      parallel = parallel,
      workers = workers,
      return_draws = return_draws
    )
  } else if (method == "bootstrap") {
    inferences_bootstrap(
      object = object,
      n_boot = n_sim,
      parallel = parallel,
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
#' @param n_sim Number of simulation draws.
#' @param vcov Custom variance-covariance matrix (optional, overrides model vcov).
#' @param conf_level Confidence level.
#' @param conf_type Type of confidence interval ("perc" or "wald").
#' @param parallel Use parallel processing?
#' @param workers Number of parallel workers.
#' @param return_draws Store individual draws?
#'
#' @return The input object with added confidence intervals.
#'
#' @keywords internal
inferences_simulation <- function(
  object,
  n_sim,
  vcov,
  conf_level,
  conf_type,
  parallel,
  workers,
  return_draws
) {
  # Check for mvtnorm package
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop(
      "Package 'mvtnorm' is required for simulation-based inference.\n",
      "Install with: install.packages('mvtnorm')"
    )
  }

  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  call_args <- attr(object, "call_args")
  avg_args <- attr(object, "avg_args")

  tvarname <- attr(object, "tvarname")
  pvarname <- attr(object, "pvarname")
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
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
    by <- NULL
    id_var <- NULL
  }

  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform simulation inference.")
  }

  # --- 2. Get Coefficients and VCov from Model ---
  # The vcov is extracted directly from the model object:
  # - For robcov_vglm: uses the stored cluster-robust vcov
  # - For orm with rms::robcov(): uses the stored robust vcov
  # - For plain vglm/orm: uses model-based vcov
  # Users who want cluster-robust SEs should wrap their model with
  # robcov_vglm() or rms::robcov() BEFORE passing to sops()/avg_sops().
  beta_hat <- get_coef(model)

  if (!is.null(vcov) && is.matrix(vcov)) {
    # User-provided custom vcov matrix (for sensitivity analyses)
    Sigma <- vcov
  } else {
    # Extract vcov from model (robust if model was wrapped with robcov_vglm)
    Sigma <- get_vcov_robust(model, cluster = NULL)
  }

  # Validate dimensions
  if (length(beta_hat) != nrow(Sigma) || length(beta_hat) != ncol(Sigma)) {
    stop(
      "Dimension mismatch: coefficients (",
      length(beta_hat),
      ") vs ",
      "vcov matrix (",
      nrow(Sigma),
      " x ",
      ncol(Sigma),
      ")"
    )
  }

  # --- 3. Draw n_sim Coefficient Vectors from MVN ---
  beta_draws <- mvtnorm::rmvnorm(n_sim, mean = beta_hat, sigma = Sigma)

  # --- 4. Prepare Prediction Data (COMPUTED ONCE!) ---
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

  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 5. Define Analysis Function for Each Draw ---
  analysis_fn <- function(i) {
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
        n_each = n_each
      )
    } else {
      # Individual-level: convert array to data frame
      result <- array_to_df_individual(sops_array, times, ylevels, newdata_pred)
    }

    result
  }

  # --- 6. Apply Across All Draws ---
  if (parallel) {
    # Setup parallel
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("Package 'furrr' is required for parallel processing")
    }

    if (is.null(workers)) {
      workers <- max(1, parallel::detectCores() - 1)
    }

    future::plan(future.callr::callr, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    sim_results <- furrr::future_map(
      seq_len(n_sim),
      analysis_fn,
      .options = furrr::furrr_options(
        seed = TRUE,
        globals = c(
          "model",
          "beta_draws",
          "newdata_pred",
          "times",
          "ylevels",
          "absorb",
          "tvarname",
          "pvarname",
          "t_covs",
          "is_avg",
          "grid",
          "variables",
          "n_cf",
          "n_each"
        ),
        packages = c("rms", "VGAM", "dplyr", "stats")
      ),
      .progress = FALSE
    )
  } else {
    sim_results <- lapply(seq_len(n_sim), analysis_fn)
  }

  # --- 7. Combine Results ---
  sim_results <- Filter(Negate(is.null), sim_results)

  if (length(sim_results) == 0) {
    stop("All simulation draws failed.")
  }

  # Add draw_id
  for (i in seq_along(sim_results)) {
    sim_results[[i]]$draw_id <- i
  }
  draws_df <- dplyr::bind_rows(sim_results)

  # --- 8. Compute Confidence Intervals ---
  if (is_avg) {
    group_cols <- c("time", "state", names(variables))
    if (!is.null(by)) {
      group_cols <- unique(c(group_cols, by))
    }
  } else {
    group_cols <- c("rowid", "time", "state")
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = draws_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )

  # --- 9. Merge with Original Object ---
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

  if (return_draws) {
    attr(final_result, "simulation_draws") <- draws_df
  }

  final_result
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
#' @param parallel Use parallel processing?
#' @param workers Number of parallel workers.
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
  parallel,
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
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
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
    rows_per_patient <- newdata_orig |>
      dplyr::group_by(!!rlang::sym(id_var)) |>
      dplyr::summarise(n_rows = dplyr::n(), .groups = "drop") |>
      dplyr::pull(n_rows)
    
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
  factor_cols <- c("y", pvarname)
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
    states_present <- setdiff(ylevels, as.numeric(missing_states))

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
          original_idx <- which(ylevels == original_state)
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

    dplyr::bind_rows(result_list)
  }

  # --- 5. Apply to Bootstrap Samples ---
  boot_results <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = newdata_orig,
    id_var = id_var,
    parallel = parallel,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "dplyr", "stats"),
    globals = c(
      "model",
      "variables",
      "times",
      "ylevels",
      "absorb",
      "pvarname",
      "tvarname",
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
  boot_df <- dplyr::bind_rows(boot_results)

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
  n_each
) {
  n_times <- length(times)
  n_states <- length(ylevels)

  # Average within each counterfactual group
  avg_sops_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    start_idx <- (cf_i - 1) * n_each + 1
    end_idx <- cf_i * n_each

    sops_cf <- sops_array[start_idx:end_idx, , , drop = FALSE]
    avg_sops_mat <- apply(sops_cf, c(2, 3), mean)

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

  dplyr::bind_rows(result_list)
}


#' Convert SOPs Array to Individual-Level Data Frame
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param times Vector of time points.
#' @param ylevels Vector of state levels.
#' @param newdata Original data with rowid.
#'
#' @return Data frame with individual SOPs.
#'
#' @keywords internal
array_to_df_individual <- function(sops_array, times, ylevels, newdata) {
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]

  # Flatten array
  probs_flat <- as.vector(sops_array)

  # Construct indices
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)

  # Build result
  result <- data.frame(
    rowid = if ("rowid" %in% names(newdata)) {
      newdata$rowid[idx_pat]
    } else {
      idx_pat
    },
    time = idx_time,
    state = idx_state,
    estimate = probs_flat
  )

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
  if (is.null(draws)) draws <- attr(object, "simulation_draws")
  if (is.null(draws)) draws <- attr(object, "draws")

  if (is.null(draws)) {
    method <- attr(object, "method")
    msg <- "No draws found. Run inferences() with return_draws = TRUE."
    if (!is.null(method)) msg <- paste0(msg, " (Method used: '", method, "')")
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
  if (requireNamespace("dplyr", quietly = TRUE)) {
    draws <- dplyr::left_join(draws, meta, by = keys)
  } else {
    draws <- merge(draws, meta, by = keys, all.x = TRUE, sort = FALSE)
  }

  return(draws)
}
