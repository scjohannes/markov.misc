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
#' @param object A `markov_avg_sops` object, or a `markov_sops` object when
#'   `engine = "fwb"`.
#' @param engine Bootstrap engine, `"standard"` or `"fwb"`.
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
  engine,
  n_boot,
  workers,
  conf_level,
  conf_type = "perc",
  return_draws,
  update_datadist,
  use_coefstart
) {
  engine <- match.arg(engine, choices = c("standard", "fwb"))

  if (!inherits(object, "markov_avg_sops")) {
    if (engine != "fwb") {
      stop(
        "Standard refit bootstrap is not supported for individual-level ",
        "SOPs because ordinary resampling can drop outcome-state support ",
        "needed by fixed prediction rows. Use `engine = \"fwb\"`, or call ",
        "`avg_sops()` for marginal bootstrap inference."
      )
    }
    return(inferences_bootstrap_sops_fwb(
      object = object,
      n_boot = n_boot,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    ))
  }

  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  prediction_data <- attr(object, "newdata_pred")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))
  refit_data <- attr(object, "refit_data") %||% newdata_orig
  call_args <- attr(object, "call_args")
  avg_args <- attr(object, "avg_args")

  time_var <- attr(object, "time_var")
  p_var <- attr(object, "p_var")
  p2_var <- attr(object, "p2_var")
  y_levels <- attr(object, "y_levels")
  absorb <- attr(object, "absorb")
  gap_var <- attr(object, "gap_var")
  time_covariates <- attr(object, "time_covariates")

  variables <- avg_args$variables
  by <- avg_args$by
  times <- avg_args$times
  id_var <- avg_args$id_var

  if (is.null(refit_data)) {
    stop("Full refit data not stored. Cannot perform bootstrap.")
  }

  validate_refit_bootstrap_data(refit_data, id_var, time_var)

  if (newdata_supplied && engine == "fwb") {
    warn_fixed_profile_bootstrap_weights("FWB")
  }

  # Dimensions for result expansion
  n_times <- length(times)
  n_states <- length(y_levels)

  # --- 2. Identify Factor Columns ---
  factor_cols <- c("y", p_var, p2_var)
  factor_cols <- intersect(factor_cols, names(refit_data))

  # --- 4. Define Analysis Function ---
  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    fit_weights <- NULL
    if (engine == "fwb") {
      fit_weights <- boot_data$fwb_weight
    }

    # A. Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      y_levels = y_levels,
      absorb = absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart,
      fit_weights = fit_weights
    )

    m_boot <- boot_result$model
    boot_data <- boot_result$data
    boot_ylevels <- boot_result$ylevels %||% boot_result$y_levels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }

    # Get states present in original numbering (for mapping back later)
    ylevel_names <- as_state_labels(y_levels)
    states_present <- ylevel_names[
      !ylevel_names %in% as_state_labels(missing_states)
    ]

    grid <- do.call(expand.grid, variables)
    if (newdata_supplied) {
      newdata_cf <- prediction_data
      n_cf <- nrow(grid)
      if (is.null(newdata_cf) || nrow(newdata_cf) %% n_cf != 0L) {
        stop("Stored prediction data is not aligned with counterfactual grid.")
      }
      baseline_weights <- NULL
    } else {
      # B. Compute standardized SOPs using G-computation on bootstrap data
      # Extract baseline data (one row per patient)
      if (engine == "standard") {
        baseline_boot <- resolve_markov_prediction_data(
          boot_data,
          id_var = "new_id",
          time_var = time_var,
          data_label = "bootstrap data"
        )
        baseline_weights <- NULL
      } else {
        baseline_boot <- resolve_markov_prediction_data(
          boot_data,
          id_var = id_var,
          time_var = time_var,
          data_label = "bootstrap data"
        )
        baseline_weights <- fwb_baseline_weights(
          fwb_weights = fwb_weights,
          baseline_data = baseline_boot,
          id_var = id_var
        )
      }

      # Create counterfactual datasets for each variable value
      newdata_cf <- create_counterfactual_data(baseline_boot, grid, variables)
    }

    # Compute individual SOPs for counterfactual data
    sops_array <- tryCatch(
      soprob_markov(
        model = m_boot,
        newdata = newdata_cf,
        times = times,
        y_levels = factor(boot_ylevels),
        absorb = boot_absorb,
        time_var = time_var,
        p_var = p_var,
        p2_var = p2_var,
        gap_var = gap_var,
        time_covariates = time_covariates
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
    n_cf <- nrow(grid)
    n_each <- nrow(newdata_cf) / n_cf
    boot_avg <- marginalize_sops_array(
      sops_array = sops_array,
      grid = grid,
      times = times,
      y_levels = factor(boot_ylevels),
      variables = variables,
      n_cf = n_cf,
      n_each = n_each,
      weights = baseline_weights,
      by = by,
      newdata = newdata_cf
    )

    if (!is.null(by)) {
      if (length(missing_states) > 0) {
        boot_avg <- complete_bootstrap_sop_states(
          boot_avg,
          times = times,
          y_levels = y_levels,
          variables = variables,
          by = by
        )
      }
      return(boot_avg)
    }

    avg_sops_list <- bootstrap_avg_df_to_matrices(
      boot_avg = boot_avg,
      grid = grid,
      times = times,
      states = factor(boot_ylevels),
      variables = variables
    )

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

      df <- expand.grid(time = times, state = y_levels)
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
  if (engine == "standard") {
    boot_ids <- fast_group_bootstrap(
      data = refit_data,
      id_var = id_var,
      n_boot = n_boot
    )

    boot_results <- apply_to_bootstrap(
      boot_samples = boot_ids,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats"),
      globals = c(
        "model",
        "prediction_data",
        "newdata_supplied",
        "variables",
        "times",
        "y_levels",
        "absorb",
        "p_var",
        "p2_var",
        "time_var",
        "gap_var",
        "time_covariates",
        "n_times",
        "n_states",
        "update_datadist",
        "factor_cols",
        "use_coefstart",
        "engine",
        "id_var",
        "by",
        "complete_bootstrap_sop_states"
      )
    )
  } else {
    fwb_samples <- generate_fwb_bootstrap_weights(
      data = refit_data,
      id_var = id_var,
      n_boot = n_boot
    )

    boot_results <- apply_to_fwb_bootstrap(
      fwb_samples = fwb_samples,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats"),
      globals = c(
        "model",
        "prediction_data",
        "newdata_supplied",
        "variables",
        "times",
        "y_levels",
        "absorb",
        "p_var",
        "p2_var",
        "time_var",
        "gap_var",
        "time_covariates",
        "n_times",
        "n_states",
        "update_datadist",
        "factor_cols",
        "use_coefstart",
        "engine",
        "id_var",
        "by",
        "complete_bootstrap_sop_states"
      )
    )
  }

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
    conf_type = conf_type,
    point_estimates = as.data.frame(object)
  )

  # Merge with original object
  object$.sop_order <- seq_len(nrow(object))
  final_result <- merge(
    object,
    summary_stats,
    by = group_cols,
    all.x = TRUE,
    sort = FALSE
  )
  final_result <- final_result[order(final_result$.sop_order), , drop = FALSE]
  final_result$.sop_order <- NULL
  rownames(final_result) <- NULL

  final_result <- restore_sops_attrs(final_result, object)

  # Add bootstrap metadata
  attr(final_result, "n_boot") <- n_boot
  attr(final_result, "n_successful") <- length(boot_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "method") <- if (engine == "fwb") "fwb" else "bootstrap"
  attr(final_result, "engine") <- engine
  if (engine == "fwb") {
    attr(final_result, "fwb_weight_type") <- "exponential"
    attr(final_result, "fwb_weight_scale") <- "cluster_mean_1"
    attr(final_result, "draw_weights_attached") <- FALSE
    attr(final_result, "draw_weight_col") <- NULL
    attr(final_result, "draw_weight_omission_reason") <- "averaged_sops"
  }

  # Store full bootstrap draws if requested
  if (return_draws) {
    attr(final_result, "draws") <- boot_df
  }

  final_result
}

inferences_bootstrap_sops_fwb <- function(
  object,
  n_boot,
  workers,
  conf_level,
  conf_type = "perc",
  return_draws,
  update_datadist,
  use_coefstart
) {
  model <- attr(object, "model")
  prediction_data <- attr(object, "newdata_pred") %||%
    attr(object, "newdata_orig")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))
  refit_data <- attr(object, "refit_data") %||% markov_model_data(model)
  call_args <- attr(object, "call_args")

  time_var <- attr(object, "time_var")
  p_var <- attr(object, "p_var")
  p2_var <- attr(object, "p2_var")
  y_levels <- attr(object, "y_levels")
  absorb <- attr(object, "absorb")
  gap_var <- attr(object, "gap_var")
  time_covariates <- attr(object, "time_covariates")
  by <- call_args$by %||% attr(object, "by")
  times <- call_args$times
  id_var <- attr(object, "id_var") %||% markov_model_id_var(model)

  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform bootstrap.")
  }
  if (is.null(prediction_data)) {
    stop("Prediction data not stored. Cannot perform bootstrap.")
  }
  if (is.null(refit_data)) {
    stop(
      "Full refit data not stored. Fit with `orm_markov()` or ",
      "`vglm_markov(id_var = ...)`, or pass full data to `sops()`."
    )
  }
  if (is.null(id_var)) {
    stop(
      "`id_var` is required for fractional weighted bootstrap. Fit with ",
      "`orm_markov(id_var = ...)`, `vglm_markov(id_var = ...)`, or pass ",
      "`id_var` to `sops()`."
    )
  }
  validate_refit_bootstrap_data(refit_data, id_var, time_var)

  prediction_data <- ensure_markov_rowid(prediction_data)
  if (!newdata_supplied) {
    validate_prediction_weight_ids(prediction_data, id_var)
  } else if (!is.null(by)) {
    warn_fixed_profile_bootstrap_weights("FWB")
  } else if (isTRUE(return_draws)) {
    warn_fixed_profile_draw_weights("FWB")
  }

  draw_weight_col <- if (!newdata_supplied && is.null(by)) {
    "fwb_weight"
  } else {
    NULL
  }
  validate_draw_weight_column_available(prediction_data, draw_weight_col)
  draw_weights_attached <- !is.null(draw_weight_col) && isTRUE(return_draws)
  draw_weight_omission_reason <- if (draw_weights_attached) {
    NULL
  } else if (newdata_supplied) {
    "user_supplied_newdata"
  } else if (!is.null(by)) {
    "grouped_sops"
  } else {
    NULL
  }

  factor_cols <- c("y", p_var, p2_var)
  factor_cols <- intersect(factor_cols, names(refit_data))

  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    fit_weights <- boot_data$fwb_weight

    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      y_levels = y_levels,
      absorb = absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart,
      fit_weights = fit_weights
    )

    m_boot <- boot_result$model
    boot_ylevels <- boot_result$ylevels %||% boot_result$y_levels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }
    if (length(missing_states) > 0) {
      warning(
        "Fractional weighted bootstrap unexpectedly dropped outcome-state ",
        "support; skipping this draw.",
        call. = FALSE
      )
      return(NULL)
    }

    sops_array <- tryCatch(
      soprob_markov(
        model = m_boot,
        newdata = prediction_data,
        times = times,
        y_levels = factor(boot_ylevels),
        absorb = boot_absorb,
        time_var = time_var,
        p_var = p_var,
        p2_var = p2_var,
        gap_var = gap_var,
        time_covariates = time_covariates
      ),
      error = function(e) {
        warning("soprob_markov failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(sops_array)) {
      return(NULL)
    }

    prediction_weights <- if (newdata_supplied) {
      NULL
    } else {
      fwb_baseline_weights(
        fwb_weights = fwb_weights,
        baseline_data = prediction_data,
        id_var = id_var
      )
    }

    array_to_df_individual(
      sops_array,
      times,
      factor(boot_ylevels),
      prediction_data,
      by = by,
      weights = prediction_weights,
      weight_col = draw_weight_col
    )
  }

  fwb_samples <- generate_fwb_bootstrap_weights(
    data = refit_data,
    id_var = id_var,
    n_boot = n_boot
  )

  boot_results <- apply_to_fwb_bootstrap(
    fwb_samples = fwb_samples,
    analysis_fn = analysis_fn,
    data = refit_data,
    id_var = id_var,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "stats"),
    globals = c(
      "model",
      "prediction_data",
      "refit_data",
      "times",
      "y_levels",
      "absorb",
      "p_var",
      "p2_var",
      "time_var",
      "gap_var",
      "time_covariates",
      "update_datadist",
      "factor_cols",
      "use_coefstart",
      "id_var",
      "by",
      "newdata_supplied",
      "draw_weight_col"
    )
  )

  boot_results <- Filter(Negate(is.null), boot_results)
  if (length(boot_results) == 0) {
    stop("All bootstrap iterations failed.")
  }

  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- bind_rows_fill(boot_results)

  if (is.null(by)) {
    group_cols <- c("rowid", "time", "state")
  } else {
    group_cols <- unique(c("time", "state", by))
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = boot_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type,
    point_estimates = as.data.frame(object)
  )

  object$.sop_order <- seq_len(nrow(object))
  final_result <- merge(
    object,
    summary_stats,
    by = group_cols,
    all.x = TRUE,
    sort = FALSE
  )
  final_result <- final_result[order(final_result$.sop_order), , drop = FALSE]
  final_result$.sop_order <- NULL
  rownames(final_result) <- NULL
  final_result <- restore_sops_attrs(final_result, object)

  attr(final_result, "n_boot") <- n_boot
  attr(final_result, "n_successful") <- length(boot_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "method") <- "fwb"
  attr(final_result, "engine") <- "fwb"
  attr(final_result, "fwb_weight_type") <- "exponential"
  attr(final_result, "fwb_weight_scale") <- "cluster_mean_1"
  attr(final_result, "draw_weights_attached") <- draw_weights_attached
  attr(final_result, "draw_weight_col") <- draw_weight_col
  if (!is.null(draw_weight_omission_reason)) {
    attr(final_result, "draw_weight_omission_reason") <-
      draw_weight_omission_reason
  }

  if (return_draws) {
    attr(final_result, "draws") <- boot_df
  }

  final_result
}


# =============================================================================
# HELPER FUNCTIONS FOR INFERENCE
# =============================================================================

validate_refit_bootstrap_data <- function(
  refit_data,
  id_var,
  time_var,
  data_arg = "refit_data"
) {
  validate_markov_id_var(id_var, refit_data, data_arg)

  if (!is.null(time_var) && time_var %in% names(refit_data)) {
    rows_per_patient <- tabulate(match(
      refit_data[[id_var]],
      unique(refit_data[[id_var]])
    ))

    if (all(rows_per_patient == 1)) {
      stop(
        "Bootstrap inference requires full longitudinal data (all time points), ",
        "but `",
        data_arg,
        "` appears to be baseline only (one row per patient).\n\n",
        "For bootstrap: pass the full longitudinal data as `newdata` or ",
        "`refit_data`.\n",
        "For simulation: baseline-only prediction profiles are supported.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}

bootstrap_avg_df_to_matrices <- function(
  boot_avg,
  grid,
  times,
  states,
  variables
) {
  out <- vector("list", nrow(grid))

  for (cf_i in seq_len(nrow(grid))) {
    keep <- rep(TRUE, nrow(boot_avg))
    for (v in names(variables)) {
      keep <- keep & boot_avg[[v]] == grid[[v]][cf_i]
    }

    df <- boot_avg[keep, , drop = FALSE]
    mat <- matrix(
      df$estimate,
      nrow = length(times),
      ncol = length(states),
      byrow = FALSE
    )
    out[[cf_i]] <- mat
  }

  out
}

complete_bootstrap_sop_states <- function(
  boot_avg,
  times,
  y_levels,
  variables,
  by = NULL
) {
  group_cols <- unique(c(names(variables), by))
  group_values <- unique(boot_avg[, group_cols, drop = FALSE])
  base <- expand.grid(
    time = times,
    state = y_levels,
    KEEP.OUT.ATTRS = FALSE
  )

  pieces <- vector("list", nrow(group_values))
  for (i in seq_len(nrow(group_values))) {
    piece <- base
    for (col in group_cols) {
      piece[[col]] <- group_values[[col]][i]
    }
    pieces[[i]] <- piece
  }

  full <- bind_rows_fill(pieces)
  keys <- unique(c("time", "state", group_cols))
  out <- merge(full, boot_avg, by = keys, all.x = TRUE, sort = FALSE)
  out$estimate[is.na(out$estimate)] <- 0
  out
}
