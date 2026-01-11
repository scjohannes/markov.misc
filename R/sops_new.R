#' Calculate Individual State Occupation Probabilities (Tidy)
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
#' @param newdata Optional data frame. If NULL, extracts from model.
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
#' # Compute marginal SOPs comparing tx=1 vs tx=0
#' avg_sops(
#'   model = fit,
#'   newdata = data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6
#' )
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


#' Bootstrap Inference for State Occupation Probabilities
#'
#' Adds bootstrap confidence intervals to `avg_sops` objects. Implements the
#' same robust approach as `bootstrap_standardized_sops()` for handling rare
#' states in bootstrap samples.
#'
#' @param object A `markov_avg_sops` object from `avg_sops()`.
#' @param method Character. Currently only "bootstrap" is supported.
#' @param n_boot Number of bootstrap iterations.
#' @param parallel Logical. Use parallel processing?
#' @param workers Number of parallel workers. If NULL, uses detectCores() - 1.
#' @param update_datadist Logical. Whether to update datadist for rms models.
#' @param conf_level Confidence level for intervals (default 0.95).
#' @param return_draws Logical. If TRUE, stores all individual bootstrap
#'   estimates as an attribute that can be extracted with `draw_bootstrap()`.
#'   This allows users to plot distributions, compute custom statistics, or
#'   aggregate across time. Default is FALSE to save memory.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object with added columns:
#'   \item{conf.low}{Lower confidence bound}
#'   \item{conf.high}{Upper confidence bound}
#'   \item{std.error}{Bootstrap standard error}
#'
#'   If `return_draws = TRUE`, the object also has a "bootstrap_draws"
#'   attribute containing a data frame with all individual bootstrap estimates.
#'   Extract with `get_draws()`.
#'
#' @details
#' This function implements group bootstrap resampling (by patient ID) with
#' proper handling of rare states that may be missing from bootstrap samples:
#'
#' **Bootstrap Procedure:**
#' 1. Sample patient IDs with replacement using `fast_group_bootstrap()`
#' 2. For each bootstrap sample:
#'    - Relevel factors to consecutive integers if states are missing
#'    - Refit the model on the bootstrap data
#'    - Predict on the ORIGINAL counterfactual grid (not bootstrap sample)
#'    - Map predictions back to original state space
#'    - Zero-pad missing states (they have 0 probability in that bootstrap)
#' 3. Compute quantile-based confidence intervals across bootstrap iterations
#'
#' **Missing State Handling:**
#' When a state is absent from a bootstrap sample:
#' - The model is refit without that state level
#' - The `pvarname` values in the prediction data are mapped to the new levels
#' - Results are expanded back to the original state space
#' - Missing states receive probability 0 (not NA)
#'
#' This approach maintains valid probability distributions (sum to 1) and
#' avoids the need for imputation or jittering strategies.
#'
#' @seealso [avg_sops()], [bootstrap_standardized_sops()],
#'   [fast_group_bootstrap()], [relevel_factors_consecutive()]
#'
#' @examples
#' \dontrun{
#' # Compute marginal SOPs with bootstrap CIs
#' result <- avg_sops(
#'   model = fit,
#'   newdata = data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6,
#'   id_var = "id"
#' ) |>
#'   inferences(n_boot = 500, parallel = TRUE, workers = 4)
#'
#' # Extract individual bootstrap draws for custom analyses
#' result_with_draws <- avg_sops(...) |>
#'   inferences(n_boot = 500, return_draws = TRUE)
#'
#' # Get all bootstrap estimates
#' boot_draws <- draw_bootstrap(result_with_draws)
#'
#' # Example: Compute mean time in state 1 for each bootstrap iteration
#' library(dplyr)
#' time_in_state <- boot_draws |>
#'   filter(state == 1, tx == 1) |>
#'   group_by(boot_id) |>
#'   summarise(mean_time = sum(estimate))
#' }
#'
#' @seealso [get_draws()] to extract all individual draws from any sampling procedure
#'
#' @export
inferences <- function(
  object,
  method = "bootstrap",
  n_boot = 1000,
  parallel = FALSE,
  workers = NULL,
  update_datadist = TRUE,
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  if (!inherits(object, "markov_avg_sops")) {
    stop("inferences() currently only supports 'markov_avg_sops' objects.")
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

  # Dimensions for result expansion
  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 2. Extract Baseline Data and Create Counterfactual Grid ---
  # For SOP prediction, we need one row per patient (baseline)
  baseline_data <- newdata_orig[!duplicated(newdata_orig[[id_var]]), ]

  # Create counterfactual grid for prediction (FIXED across bootstraps)
  grid <- do.call(expand.grid, variables)
  cf_data_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt_copy <- baseline_data
    for (v in names(grid)) {
      dt_copy[[v]] <- grid[i, v]
    }
    cf_data_list[[i]] <- dt_copy
  }
  newdata_cf <- do.call(rbind, cf_data_list)

  # --- 3. Generate Bootstrap ID Samples ---
  # Bootstrap from the FULL longitudinal data (for model refitting)
  boot_ids <- fast_group_bootstrap(
    data = newdata_orig,
    id_var = id_var,
    n_boot = n_boot
  )

  # --- 4. Define Analysis Function ---
  # This runs on each bootstrap sample (refit model, predict on CF data)
  analysis_fn <- function(boot_data) {
    # A. Relevel factors to handle missing states
    factor_cols <- c("y", pvarname)
    factor_cols <- intersect(factor_cols, names(boot_data))

    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = newdata_orig,
      ylevels = ylevels,
      absorb = absorb,
      update_datadist = update_datadist
    )

    m_boot <- boot_result$model
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }

    # B. Prepare prediction data with mapped pvarname levels
    # Get states present in bootstrap (in original numbering)
    states_present <- setdiff(ylevels, as.numeric(missing_states))
    state_mapping <- stats::setNames(seq_along(states_present), states_present)

    # Copy counterfactual data and map pvarname
    newdata_pred <- newdata_cf

    # Get pvarname values
    pvar_vals <- newdata_pred[[pvarname]]
    if (is.factor(pvar_vals)) {
      pvar_vals <- as.numeric(as.character(pvar_vals))
    }

    # Check which values are present in the model
    is_present <- pvar_vals %in% states_present

    if (!all(is_present)) {
      # Some pvarname values are not in the model
      # For G-computation prediction, we need valid levels
      # Map to nearest present level (this is for PREDICTION, not training)
      for (idx in which(!is_present)) {
        old_val <- pvar_vals[idx]
        # Find nearest present level
        diffs <- abs(states_present - old_val)
        new_val <- states_present[which.min(diffs)]
        pvar_vals[idx] <- new_val
      }
    }

    # Apply mapping to consecutive integers for the bootstrap model
    mapped_vals <- state_mapping[as.character(pvar_vals)]
    newdata_pred[[pvarname]] <- factor(
      mapped_vals,
      levels = seq_along(states_present)
    )

    # C. Compute SOPs on counterfactual data
    sops_result <- tryCatch(
      soprob_markov(
        object = m_boot,
        data = newdata_pred,
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

    if (is.null(sops_result)) {
      return(NULL)
    }

    # D. Expand results back to original state space
    # sops_result is array [n_pat, n_times, n_boot_states]
    n_pat <- dim(sops_result)[1]
    n_boot_states <- dim(sops_result)[3]

    # Initialize full array with zeros for all original states
    sops_full <- array(0, dim = c(n_pat, n_times, n_states))

    # Fill in present states
    for (i in seq_along(states_present)) {
      original_idx <- which(ylevels == states_present[i])
      sops_full[,, original_idx] <- sops_result[,, i]
    }

    # E. Marginalize and format as avg_sops output
    # For each treatment value, average across patients
    n_cf <- nrow(grid)
    # newdata_cf has n_baseline_patients * n_cf rows total
    n_each <- n_pat %/% n_cf # Number of patients per counterfactual group

    result_list <- vector("list", n_cf)
    for (cf_i in seq_len(n_cf)) {
      # Indices for this counterfactual group
      start_idx <- (cf_i - 1) * n_each + 1
      end_idx <- cf_i * n_each

      # Subset array and average across patients
      sops_cf <- sops_full[start_idx:end_idx, , , drop = FALSE]
      avg_sops_mat <- apply(sops_cf, c(2, 3), mean) # [n_times x n_states]

      # Convert to data frame
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
      "newdata_cf",
      "variables",
      "grid",
      "times",
      "ylevels",
      "absorb",
      "pvarname",
      "tvarname",
      "t_covs",
      "n_times",
      "n_states",
      "update_datadist"
    )
  )

  # --- 6. Combine and Compute Summary Statistics ---
  # Filter out NULL results
  boot_results <- Filter(Negate(is.null), boot_results)

  if (length(boot_results) == 0) {
    stop("All bootstrap iterations failed.")
  }

  # Add boot_id and combine
  for (i in seq_along(boot_results)) {
    boot_results[[i]]$boot_id <- i
  }
  boot_df <- dplyr::bind_rows(boot_results)

  # Compute confidence intervals
  alpha <- 1 - conf_level
  group_cols <- c("time", "state", names(variables))
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }

  # Aggregate to get quantiles
  agg_formula <- stats::as.formula(
    paste("estimate ~", paste(group_cols, collapse = " + "))
  )
  summary_stats <- stats::aggregate(
    agg_formula,
    data = boot_df,
    FUN = function(x) {
      c(
        conf.low = stats::quantile(x, alpha / 2, na.rm = TRUE),
        conf.high = stats::quantile(x, 1 - alpha / 2, na.rm = TRUE),
        std.error = stats::sd(x, na.rm = TRUE)
      )
    }
  )

  # Fix aggregate's matrix column output
  mat <- summary_stats$estimate
  summary_stats$estimate <- NULL
  summary_stats$conf.low <- mat[, 1]
  summary_stats$conf.high <- mat[, 2]
  summary_stats$std.error <- mat[, 3]

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
  attr(final_result, "method") <- method

  # Store full bootstrap draws if requested
  if (return_draws) {
    attr(final_result, "bootstrap_draws") <- boot_df
  }

  return(final_result)
}


#' Extract Individual Draws from Inference Objects
#'
#' Extracts the individual draws (bootstrap samples, simulated values from MVN, etc.)
#' from an object returned by `inferences()` with `return_draws = TRUE`.
#' This allows users to analyze the full distribution, compute custom statistics,
#' or create visualizations, regardless of the sampling method used.
#'
#' @param object An object returned by `inferences()` with `return_draws = TRUE`
#'   (or similar inference functions that support storing draws).
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item boot_id: Bootstrap iteration number
#'     \item time: Time point
#'     \item state: State number
#'     \item estimate: Bootstrap estimate of state occupation probability
#'     \item Additional columns for standardization variables (e.g., tx)
#'   }
#'   Each row represents one bootstrap estimate for a specific time-state-treatment
#'   combination.
#'
#' @details
#' This function retrieves the draws from objects created by `inferences()`
#' (or similar inference functions). The exact name of the draws attribute
#' depends on the sampling method used (e.g., "bootstrap_draws" for bootstrap,
#' potentially "mvn_draws" for MVN simulation). If the attribute is not present
#' (i.e., the object was created without `return_draws = TRUE`), an error is raised.
#'
#' **Common Use Cases:**
#' \itemize{
#'   \item Plot distributions with histograms or density plots
#'   \item Compute custom summary statistics (e.g., median, quantiles)
#'   \item Aggregate across time to calculate total time in state
#'   \item Compute treatment effects with CIs
#'   \item Perform sensitivity analyses on the sampled values
#'   \item Combine draws across multiple inference results
#' }
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
#'   ggplot(aes(x = estimate)) +
#'   geom_histogram(bins = 30) +
#'   labs(title = "Distribution: P(State 1 | Time 10, Treatment)")
#'
#' # Compute mean time in state 1 with bootstrap CI
#' library(dplyr)
#' time_in_state_boot <- draws |>
#'   filter(state == 1, tx == 1) |>
#'   group_by(boot_id) |>
#'   summarise(total_time = sum(estimate))
#'
#' quantile(time_in_state_boot$total_time, c(0.025, 0.5, 0.975))
#'
#' # Compare treatment effect on time in state
#' treatment_effect <- draws |>
#'   filter(state == 1) |>
#'   group_by(boot_id, tx) |>
#'   summarise(total_time = sum(estimate), .groups = "drop") |>
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
  if (!inherits(object, "markov_avg_sops")) {
    stop(
      "get_draws() requires an object from inferences(). ",
      "Did you forget to call inferences() first?"
    )
  }

  # Try to find draws attribute (could be bootstrap_draws, mvn_draws, etc.)
  draws <- attr(object, "bootstrap_draws")

  # Not yet implemented
  # if (is.null(draws)) {
  #   draws <- attr(object, "mvn_draws")
  # }

  # if (is.null(draws)) {
  #   draws <- attr(object, "draws")
  # }

  if (is.null(draws)) {
    stop(
      "No draws found. ",
      "Run inferences() with return_draws = TRUE (or similar) to store draws."
    )
  }

  return(draws)
}
