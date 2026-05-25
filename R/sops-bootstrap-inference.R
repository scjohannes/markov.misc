# SOP refit-bootstrap inference.

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
    newdata_cf <- create_counterfactual_data(baseline_boot, grid, variables)

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

  final_result <- restore_sops_attrs(final_result, object)

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
