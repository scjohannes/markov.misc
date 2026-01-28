# =============================================================================
# ARCHIVED: Old Slow Implementation for Benchmarking
# =============================================================================
#
# This file contains the original slow implementation of SOPs computation
# and inference. It's kept for benchmarking purposes only.
#
# The slow path calls VGAM::predict() for each simulation draw, which has
# significant overhead from formula parsing and model matrix construction.
#
# The fast path (now in R/sops.R) pre-computes design matrix decompositions
# once and uses matrix multiplication for each draw.
# =============================================================================

#' Slow Simulation-Based Inference (FOR BENCHMARKING ONLY)
#'
#' This is the original slow implementation that calls predict() for each
#' simulation draw. It's kept for benchmarking against the fast path.
#'
#' @keywords internal
inferences_simulation_slow <- function(
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
    stop("Package 'mvtnorm' is required for simulation-based inference.")
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

  is_avg <- inherits(object, "markov_avg_sops")
  if (is_avg) {
    variables <- avg_args$variables
    by <- avg_args$by
    times <- avg_args$times
    id_var <- avg_args$id_var
  } else {
    times <- call_args$times
    variables <- NULL
    by <- NULL
    id_var <- NULL
  }

  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform simulation inference.")
  }

  # --- 2. Get Coefficients and VCov from Model ---
  beta_hat <- markov.misc::get_coef(model)

  if (!is.null(vcov) && is.matrix(vcov)) {
    Sigma <- vcov
  } else {
    Sigma <- markov.misc::get_vcov_robust(model, cluster = NULL)
  }

  # Validate dimensions
  if (length(beta_hat) != nrow(Sigma) || length(beta_hat) != ncol(Sigma)) {
    stop("Dimension mismatch: coefficients vs vcov matrix")
  }

  # --- 3. Draw n_sim Coefficient Vectors from MVN ---
  beta_draws <- mvtnorm::rmvnorm(n_sim, mean = beta_hat, sigma = Sigma)

  # --- 4. Prepare Prediction Data ---
  if (is_avg) {
    baseline_data <- newdata_orig[!duplicated(newdata_orig[[id_var]]), ]
    grid <- do.call(expand.grid, variables)
    newdata_pred <- create_counterfactual_data_slow(baseline_data, grid, variables)
    n_cf <- nrow(grid)
    n_each <- nrow(baseline_data)
  } else {
    newdata_pred <- newdata_orig
    grid <- NULL
    n_cf <- 1
    n_each <- nrow(newdata_pred)
  }

  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 5. SLOW PATH: Use soprob_markov with coefficient replacement ---
  # This is the slow part - calling predict() for each draw
  analysis_fn <- function(i) {
    # Replace coefficients
    model_i <- markov.misc::set_coef(model, beta_draws[i, ])

    # Compute SOPs using soprob_markov (which calls predict internally)
    sops_array <- tryCatch(
      markov.misc::soprob_markov(
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
      result <- marginalize_sops_array_slow(
        sops_array = sops_array,
        grid = grid,
        times = times,
        ylevels = ylevels,
        variables = variables,
        n_cf = n_cf,
        n_each = n_each
      )
    } else {
      result <- array_to_df_individual_slow(sops_array, times, ylevels, newdata_pred)
    }

    result
  }

  # --- 6. Apply Across All Draws ---
  if (parallel) {
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
        packages = c("rms", "VGAM", "dplyr", "stats", "markov.misc")
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

  summary_stats <- compute_ci_from_draws_slow(
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
  attr(final_result, "method") <- "simulation_slow"

  if (return_draws) {
    attr(final_result, "simulation_draws") <- draws_df
  }

  final_result
}


# Helper functions (copies for the slow path)

create_counterfactual_data_slow <- function(baseline_data, grid, variables) {
  cf_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt <- baseline_data
    for (v in names(grid)) dt[[v]] <- grid[i, v]
    cf_list[[i]] <- dt
  }
  do.call(rbind, cf_list)
}

marginalize_sops_array_slow <- function(sops_array, grid, times, ylevels, variables, n_cf, n_each) {
  out_list <- vector("list", n_cf)
  n_states <- dim(sops_array)[3]
  use_ylevels <- if(is.null(ylevels)) 1:n_states else ylevels

  for (i in seq_len(n_cf)) {
    idx_start <- (i - 1) * n_each + 1
    idx_end <- i * n_each
    sub_arr <- sops_array[idx_start:idx_end, , , drop = FALSE]
    avg_mat <- apply(sub_arr, c(2, 3), mean)

    df <- expand.grid(time = times, state = use_ylevels)
    df$estimate <- as.vector(avg_mat)
    for (v in names(grid)) df[[v]] <- grid[i, v]
    out_list[[i]] <- df
  }
  do.call(rbind, out_list)
}

array_to_df_individual_slow <- function(sops_array, times, ylevels, newdata) {
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]

  probs_flat <- as.vector(sops_array)

  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)

  result <- newdata[idx_pat, , drop = FALSE]
  result$time <- idx_time
  result$state <- idx_state
  result$estimate <- probs_flat
  rownames(result) <- NULL

  result
}

compute_ci_from_draws_slow <- function(draws_df, group_cols, conf_level, conf_type) {
  alpha <- (1 - conf_level) / 2

  summary_df <- draws_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      std.error = stats::sd(estimate, na.rm = TRUE),
      conf.low = if (conf_type == "perc") {
        stats::quantile(estimate, alpha, na.rm = TRUE)
      } else {
        mean(estimate, na.rm = TRUE) - stats::qnorm(1 - alpha) * stats::sd(estimate, na.rm = TRUE)
      },
      conf.high = if (conf_type == "perc") {
        stats::quantile(estimate, 1 - alpha, na.rm = TRUE)
      } else {
        mean(estimate, na.rm = TRUE) + stats::qnorm(1 - alpha) * stats::sd(estimate, na.rm = TRUE)
      },
      .groups = "drop"
    )

  summary_df
}
