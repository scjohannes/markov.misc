#' Prepare Markov Data for Proportional Odds Modeling
#'
#' Prepares longitudinal Markov trajectory data for proportional odds regression
#' by removing absorbing states and converting variables to appropriate types.
#' Works with both VGAM::vglm and rms::orm models.
#'
#' @param data A data frame containing Markov trajectory data with columns
#'   `y` (current state) and `yprev` (previous state).
#' @param absorbing_state Integer. The absorbing state to filter out (default: 6
#'   for death). Observations where `yprev` equals this value are removed.
#' @param ordered_response Logical. Should `y` be converted to ordered factor?
#'   (default: TRUE). Required for clm and vglm, not required for orm.
#' @param factor_previous Logical. Should `yprev` be converted to factor?
#'   (default: TRUE).
#'
#' @return A data frame with modified `y` and `yprev` columns and absorbing
#'   states removed.
#'
#' @details
#' This is a convenience function that performs standard data preparation for
#' proportional odds models on Markov data:
#' - Removes rows where patients were in an absorbing state at the previous time
#' - Converts current state to ordered factor (for cumulative models)
#' - Converts previous state to factor (for including as predictor)
#'
#' **Package compatibility**: This function prepares data for:
#' - **VGAM::vglm**: Use `cumulative(parallel = TRUE, reverse = TRUE)` default
#'   datasets.
#' - **rms::orm**: Automatically treats integers as ordered factors.
#' - **ordinal::clm**: Only accepts ordered factors.
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' model_data <- prepare_markov_data(trajectories)
#'
#' # Fit with VGAM (no robust SEs; use bootstrap for inference)
#' library(VGAM)
#' fit <- vglm(y ~ tx + rcs(time, 4) + yprev,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = model_data)
#'
#' # Fit with rms (supports robust SEs via robcov)
#' library(rms)
#' dd <- datadist(model_data)
#' options(datadist = "dd")
#' fit <- orm(y ~ tx + rcs(time, 4) + yprev, data = model_data, x = TRUE, y = TRUE)
#' fit_robust <- robcov(fit, cluster = model_data$id)
#'
#' # Fit with ordinal
#' library(ordinal)
#' fit <- clm(y ~ tx + time + yprev, data = model_data)
#' }
#'
#' @importFrom dplyr filter mutate
#'
#' @export
prepare_markov_data <- function(
  data,
  absorbing_state = 6,
  ordered_response = TRUE,
  factor_previous = TRUE
) {
  if (!is.null(absorbing_state)) {
    data <- data |>
      dplyr::filter(yprev != absorbing_state)
  }

  if (factor_previous) {
    data <- data |>
      dplyr::mutate(yprev = factor(yprev))
  }

  if (ordered_response) {
    data <- data |>
      dplyr::mutate(y = ordered(y))
  }

  data
}


#' Relevel Factors to Consecutive Integers
#'
#' Handles missing factor levels in data by releveling ordered factors to
#' consecutive integers. This is useful when bootstrap samples or subsets
#' are missing certain levels, which can cause model fitting failures.
#'
#' @param data A data frame containing the data to relevel.
#' @param factor_cols Character vector of column names to relevel. Only columns
#'   that are factors with numeric-coercible levels will be releveled.
#' @param original_data Optional data frame containing the original data before
#'   subsetting. Used to determine which levels are present in the full dataset.
#' @param ylevels Optional integer vector of original state levels (e.g., 1:6).
#'   If provided along with absorb, these will be updated to match the new levels.
#' @param absorb Optional integer specifying the absorbing state in the original
#'   levels. Will be updated to the new position if provided.
#'
#' @return A list with components:
#'   - data: The data frame with releveled factors
#'   - ylevels: Updated ylevels (if provided as input), or NULL
#'   - absorb: Updated absorb (if provided as input), or NULL
#'   - missing_levels: Character vector of levels that were missing
#'
#' @details
#' When certain factor levels are absent from a dataset (e.g., in bootstrap
#' samples), this function:
#' 1. Identifies which levels are present
#' 2. Creates a mapping from old to new consecutive integers
#' 3. Relevels the specified factor columns
#' 4. Updates ylevels and absorb parameters if provided
#'
#' This is particularly useful for ordered factors representing health states
#' in Markov models, where missing states would otherwise cause model fitting
#' failures.
#'
#' @examples
#' \dontrun{
#' # Relevel y and yprev in bootstrap sample
#' releveled <- relevel_factors_consecutive(
#'   data = boot_data,
#'   factor_cols = c("y", "yprev"),
#'   original_data = full_data,
#'   ylevels = 1:6,
#'   absorb = 6
#' )
#'
#' boot_data <- releveled$data
#' boot_ylevels <- releveled$ylevels
#' boot_absorb <- releveled$absorb
#' }
#'
#' @export
relevel_factors_consecutive <- function(
  data,
  factor_cols = c("y", "yprev"),
  original_data = NULL,
  ylevels = NULL,
  absorb = NULL
) {
  # Determine original levels
  if (!is.null(original_data)) {
    original_levels <- sort(unique(unlist(lapply(
      factor_cols,
      function(col) {
        if (col %in% names(original_data)) {
          if (is.factor(original_data[[col]])) {
            as.numeric(levels(original_data[[col]]))
          } else {
            unique(original_data[[col]])
          }
        }
      }
    ))))
  } else if (!is.null(ylevels)) {
    original_levels <- ylevels
  } else {
    # Infer from data
    original_levels <- sort(unique(unlist(lapply(
      factor_cols,
      function(col) {
        if (col %in% names(data)) {
          unique(data[[col]])
        }
      }
    ))))
  }

  # Get levels present in current data
  states_present <- sort(unique(unlist(lapply(
    factor_cols,
    function(col) {
      if (col %in% names(data)) {
        if (is.factor(data[[col]])) {
          as.numeric(as.character(data[[col]]))
        } else {
          unique(data[[col]])
        }
      }
    }
  ))))

  missing_levels <- setdiff(original_levels, states_present)

  # If no levels are missing, return data as-is
  if (length(missing_levels) == 0) {
    return(list(
      data = data,
      ylevels = ylevels,
      absorb = absorb,
      missing_levels = character(0)
    ))
  }

  # Create mapping from old to new state numbers
  state_mapping <- setNames(seq_along(states_present), states_present)

  # Relevel each factor column
  for (col in factor_cols) {
    if (col %in% names(data)) {
      # Convert to numeric if factor
      if (is.factor(data[[col]])) {
        data[[col]] <- as.numeric(as.character(data[[col]]))
      }

      # Apply mapping
      data[[col]] <- state_mapping[as.character(data[[col]])]

      # Convert back to factor with new levels
      data[[col]] <- factor(data[[col]], levels = seq_along(states_present))
    }
  }

  # Update ylevels if provided
  new_ylevels <- if (!is.null(ylevels)) {
    seq_along(states_present)
  } else {
    NULL
  }

  # Update absorbing state if provided and present
  new_absorb <- if (!is.null(absorb)) {
    if (absorb %in% states_present) {
      state_mapping[as.character(absorb)]
    } else {
      warning(
        "Absorbing state ",
        absorb,
        " not present in data. Set to NULL."
      )
      NULL
    }
  } else {
    NULL
  }

  return(list(
    data = data,
    ylevels = new_ylevels,
    absorb = new_absorb,
    missing_levels = as.character(missing_levels)
  ))
}


#' Calculate Jackknife Monte Carlo Standard Error
#'
#' Computes the Monte Carlo standard error (MCSE) of a statistic using the
#' jackknife method. This provides an estimate of the variability in the
#' statistic across simulation repetitions.
#'
#' @param estimates Numeric vector of estimates from simulation runs.
#' @param statistic A function to apply to the estimates (default: mean).
#'   Should accept a numeric vector and return a single numeric value.
#'
#' @return A numeric value representing the jackknife MCSE.
#'
#' @details
#' The jackknife MCSE is calculated by:
#' 1. Creating leave-one-out estimates by removing each simulation run in turn
#' 2. Calculating the statistic on the remaining runs
#' 3. Computing the variance of these leave-one-out estimates
#' 4. Applying the jackknife variance formula: sqrt(((n-1)/n) * sum((x_i - mean)^2))
#'
#' This method is useful for assessing the Monte Carlo uncertainty in simulation
#' studies and can be applied to various statistics beyond the mean.
#'
#' @examples
#' # Calculate MCSE for mean of simulation estimates
#' sim_results <- rnorm(1000, mean = 5, sd = 2)
#' mcse_mean <- jackknife_mcse(sim_results)
#'
#' # Calculate MCSE for median
#' mcse_median <- jackknife_mcse(sim_results, statistic = median)
#'
#' @export
jackknife_mcse <- function(estimates, statistic = mean) {
  # Number of simulation repetitions
  nsim <- length(estimates)

  # Vector to store the "leave-one-out" estimates
  leave_one_out_estimates <- numeric(nsim)

  # Step 1: Create the "leave-one-out" estimates
  # Loop through each simulation run, remove it, and recalculate the statistic
  for (i in 1:nsim) {
    leave_one_out_estimates[i] <- statistic(estimates[-i])
  }

  # Step 2: Calculate the average of all the leave-one-out estimates
  mean_of_estimates <- mean(leave_one_out_estimates)

  # Step 3: Calculate the sum of squared differences
  sum_sq_diff <- sum((leave_one_out_estimates - mean_of_estimates)^2)

  # Step 4: Apply the final jackknife formula
  mcse <- sqrt(((nsim - 1) / nsim) * sum_sq_diff)

  return(mcse)
}


#' Convert State Trajectory Data to T-Test Format
#'
#' This function summarizes patient trajectories by counting the number of time
#' periods each patient spent in a specific state (typically state 1 = Home/Discharged).
#' The resulting data is suitable for simple t-tests or Wilcoxon tests comparing
#' treatment groups.
#'
#' @param data A data frame containing trajectory data, typically the output from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`. Must contain
#'   columns: `id`, `y` (state), and `tx` (treatment).
#' @param target_state Integer. The state to count (default: 1, representing Home/Discharged).
#'
#' @return A data frame with columns:
#'   - id: patient identifier
#'   - y: total number of time periods spent in `target_state`
#'   - tx: treatment assignment
#'
#' @details
#' This function aggregates longitudinal trajectory data into a single summary
#' value per patient: the total number of days (or time periods) spent in the
#' target state. This outcome can be analyzed with standard two-sample tests.
#'
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, lp_function = my_lp)
#' t_data <- states_to_ttest(trajectories, target_state = 1)
#' }
#'
#' @importFrom dplyr group_by summarise first
#'
#' @export
states_to_ttest <- function(data, target_state = 1) {
  # Input validation
  required_cols <- c("id", "y", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Base summary
  result <- data |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      y = sum(y == target_state),
      tx = dplyr::first(tx),
      .groups = "drop"
    )

  return(result)
}


#' Convert State Trajectory Data to Days Returned to Baseline (DRS) Format
#'
#' This function calculates "Days Returned to Baseline" (or "Days at Home") for
#' each patient. This is defined as the number of days spent in the target state
#' (typically Home) after the last discharge from a worse state. Patients who die
#' receive a score of -1.
#'
#' @param data A data frame containing trajectory data, typically the output from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`. Must contain
#'   columns: `id`, `time`, `y` (state), and `tx` (treatment).
#' @param follow_up_time Integer. Total follow-up time used in the simulation.
#'   This is needed to calculate days at home from the last discharge.
#' @param target_state Integer. The state representing "home" or "baseline" (default: 1).
#' @param death_state Integer. The state representing death (default: 6).
#' @param covariates Character vector of additional covariate names to include in
#'   the output (default: c("age", "sofa")). The first value of each covariate
#'   per patient will be retained.
#'
#' @return A data frame with columns:
#'   - id: patient identifier
#'   - tx: treatment assignment
#'   - drs: days returned to baseline (or -1 if died)
#'   - any additional covariates specified
#'
#' @details
#' The DRS (Days Returned to Baseline) outcome focuses on sustained recovery by
#' counting only the days at home after the final discharge. The calculation:
#' - Find the last time point where the patient was NOT at home
#' - Count the remaining days from that point to end of follow-up
#' - If the patient died at any point, assign -1
#'
#' This outcome is more sensitive to sustained recovery than simply counting
#' total days at home, as it emphasizes staying home rather than bouncing
#' in and out of the hospital.
#'
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, follow_up_time = 60,
#'                                         lp_function = my_lp)
#' drs_data <- states_to_drs(trajectories, follow_up_time = 60)
#'
#' # With different covariates
#' drs_data <- states_to_drs(trajectories, follow_up_time = 60,
#'                           covariates = c("age", "sofa", "baseline_severity"))
#' }
#'
#' @importFrom dplyr group_by arrange summarise if_else select across all_of left_join mutate
#'
#' @export
states_to_drs <- function(
  data,
  follow_up_time,
  target_state = 1,
  death_state = 6,
  covariates = c("age", "sofa")
) {
  # Input validation
  required_cols <- c("id", "time", "y", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate DRS
  result <- data |>
    dplyr::group_by(id) |>
    dplyr::arrange(time) |>
    dplyr::summarise(
      last_not_home = max(0, which(y != target_state)),
      dead = any(y == death_state),
      tx = dplyr::first(tx),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      drs = dplyr::if_else(dead, -1, follow_up_time - last_not_home)
    ) |>
    dplyr::select(id, tx, drs)

  # Add covariates if specified
  if (!is.null(covariates)) {
    available_covs <- intersect(covariates, names(data))
    if (length(available_covs) > 0) {
      cov_data <- data |>
        dplyr::group_by(id) |>
        dplyr::summarise(
          dplyr::across(dplyr::all_of(available_covs), dplyr::first),
          .groups = "drop"
        )
      result <- result |>
        dplyr::left_join(cov_data, by = "id")
    }
  }

  return(result)
}


#' Convert State Trajectory Data to Time-to-event data
#'
#' This function calculates time-to-event for usage in Cox or recurrent event
#' models.
#'
#' @param data A data frame containing trajectory data, typically the output from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`. Must contain
#'   columns: `id`, `time`, `y` (state), and `tx` (treatment).
#' @param covariates Character vector of additional covariate names to include in
#'   the output (default: c("age", "sofa")). The first value of each covariate
#'   per patient will be retained.
#'
#' @return A data frame with columns:
#'   - id: patient identifier
#'   - tx: treatment assignment
#'   - start: start of the interval, first value per participant is zero
#'   - stop: end of the interval
#'   - y: state at the end of the interval
#'   - any additional covariates specified
#'
#' @details
#' Each event change within a participant concludes a time interval indicated by
#' start and stop and the state at the end of the interval. Columns given in
#' covariates will be kept in the output.
#'
#' @importFrom dplyr group_by arrange summarise if_else select across all_of left_join mutate
#'
#' @export
states_to_tte_old <- function(
  data,
  covariates = c("age", "sofa")
) {
  # Input validation
  required_cols <- c("id", "time", "y", "tx", "yprev")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate Andersen-Gill count data format ----------------------------------

  # Indicate rows with state changes
  data$state_change <- data[["y"]] - data[["yprev"]]

  # Keep state changes and ...
  data <- data[
    data$state_change != 0 |
      # ... the last entry of each individual, unless they died
      (data[["time"]] == ave(data[["time"]], data[["id"]], FUN = max) &
        !data[["y"]] == 6),
  ]

  data[["start"]] <- ave(data$time, data$id, FUN = function(t) {
    c(0, t[-length(t)])
  })
  data[["stop"]] <- data[["time"]]
  data <- data[, c("id", "tx", "yprev", "start", "stop", "y")]

  # Add covariates if specified ------------------------------------------------
  if (!is.null(covariates)) {
    available_covs <- intersect(covariates, colnames(data))
    data <- data[,
      colnames(data) %in%
        c(
          c("id", "tx", "yprev", "start", "stop", "y"),
          available_covs
        )
    ]
  }

  return(result = data)
}


#' Convert Markov Data to Time-to-Event (Start-Stop) Format
#'
#' Converts a longitudinal dataset of daily states into a collapsed
#' start-stop format suitable for survival analysis
#'
#' @param data A data frame containing trajectory data, typically the output from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`. Must contain
#'   columns: `id`, `time`, `y` (state), and `tx` (treatment).
#' @param covariates Character vector of additional covariate names to include in
#'   the output (default: c("age", "sofa")). The first value of each covariate
#'   per patient will be retained.
#'
#' @return A tibble in start-stop format with columns:
#'   \code{id}, \code{tx}, \code{start}, \code{stop}, \code{y}, \code{yprev},
#'   and specified covariates.
#'
#' @details
#' Each event change within a participant concludes a time interval indicated by
#' start and stop and the state at the end of the interval. Columns given in
#' covariates will be kept in the output.
#'
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, follow_up_time = 60,
#'                                         lp_function = my_lp)
#' tte_data <- states_to_tte(trajectories)
#'
#' # With different covariates
#' tte_data <- states_to_tte(trajectories, covariates = c("age", "sofa", "baseline_severity"))
#' }
#' @importFrom dplyr arrange group_by mutate lag summarize first select any_of everything
#' @export
states_to_tte <- function(
  data,
  covariates = c("age", "sofa"),
  absorbing_state = 6
) {
  # 1. Input Validation
  required_cols <- c("id", "time", "y", "yprev", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # 2. Logic: Collapse consecutive runs of identical states
  # We cannot just filter; we must aggregate.
  result <- data |>
    dplyr::filter(yprev != absorbing_state) |>
    dplyr::arrange(id, time) |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      # Define the raw interval for every single row
      # If time is 1, interval is 0 -> 1.
      raw_start = dplyr::lag(time, default = 0),
      raw_stop = time,

      # Create a unique ID for every "run" of identical states
      # If y changes, we increment the ID.
      state_change = y != dplyr::lag(y, default = -999),
      run_id = cumsum(state_change)
    ) |>
    # Collapse: Squash all rows in the same run into one interval
    dplyr::group_by(id, run_id) |>
    dplyr::summarize(
      # Variables that are constant for the run
      y = dplyr::first(y),
      tx = dplyr::first(tx),
      yprev = dplyr::first(yprev),

      # The interval becomes the Min Start to the Max Stop of the run
      start = min(raw_start),
      stop = max(raw_stop),

      # Capture covariates (taking the value at the start of the run)
      dplyr::across(dplyr::any_of(covariates), dplyr::first),

      .groups = "drop"
    ) |>
    # Clean up and reorder
    dplyr::select(-run_id) |>
    dplyr::select(id, tx, start, stop, y, yprev, dplyr::everything())

  return(result)
}


#' Calculate Difference in Time Alive and Out Of Hospital
#'
#' Computes the true treatment effect as the difference in time spent in a
#' target state (typically "home" or "alive and out of hospital") between
#' treatment and control groups in simulated Markov trajectory data.
#'
#' @param data A data frame containing trajectory data with columns: `id`,
#'   `y` (state), `tx` (treatment), and `time`.
#' @param target_state Integer. The state representing the outcome of interest
#'   (default: 1, typically representing Home/Discharged/Alive and Out of Hospital).
#'
#' @return A tibble with the following columns:
#'   - `true_effect`: Difference in state occupation probability (tx=1 minus tx=0)
#'     summed across all time points
#'   - `tx0_mean_time`: Mean time spent in target state for control group
#'   - `tx0_sd_time`: Standard deviation of time in target state for control group
#'   - `tx1_mean_time`: Mean time spent in target state for treatment group
#'   - `tx1_sd_time`: Standard deviation of time in target state for treatment group
#'
#' @details
#' This function calculates the true treatment effect in two ways:
#'
#' 1. **State Occupation Probability Method**: Calculates the proportion of
#'    patients in the target state at each time point for each treatment group,
#'    then sums the difference across all time points.
#'
#' 2. **Individual Time in State**: For each patient, counts the total number
#'    of time periods spent in the target state, then computes group-level
#'    means and standard deviations.
#'
#' The `true_effect` represents the area between the state occupation curves
#' and can be interpreted as the expected difference in total days spent in
#' the target state over the follow-up period.
#'
#' @examples
#' \dontrun{
#' # Calculate true effect from simulated data
#' trajectories <- sim_trajectories_markov(baseline_data, lp_function = my_lp)
#' effect_summary <- calc_time_in_state_diff(trajectories, target_state = 1)
#'
#' # View results
#' print(effect_summary)
#' }
#'
#' @importFrom dplyr group_by summarise filter select mutate ungroup pull
#' @importFrom tidyr pivot_wider
#' @importFrom stats sd
#'
#' @export
calc_time_in_state_diff <- function(data, target_state = 1) {
  # Input validation
  required_cols <- c("id", "y", "tx", "time")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate state occupation probabilities
  sops <- data |>
    dplyr::group_by(y, tx, time) |>
    dplyr::summarise(count = n(), .groups = "drop") |>
    dplyr::group_by(tx, time) |>
    dplyr::mutate(sop = count / sum(count)) |>
    dplyr::ungroup()

  # Calculate estimand (difference in probability of being in target state)
  estimand <- sops |>
    dplyr::filter(y %in% target_state) |>
    dplyr::group_by(tx, time) |>
    dplyr::summarise(sop = sum(sop), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = tx,
      values_from = sop,
      names_prefix = "tx_"
    ) |>
    dplyr::summarise(true_effect = sum(tx_1) - sum(tx_0)) |>
    dplyr::pull(true_effect)

  # Calculate time spent in target state for each group
  time_in_state <- data |>
    dplyr::group_by(id, tx) |>
    dplyr::summarise(
      time_in_state = sum(y %in% target_state),
      .groups = "drop"
    ) |>
    dplyr::group_by(tx) |>
    dplyr::summarise(
      mean_time = mean(time_in_state),
      sd_time = sd(time_in_state),
      .groups = "drop"
    )

  # Extract values for each treatment group
  tx0_mean <- time_in_state |> dplyr::filter(tx == 0) |> dplyr::pull(mean_time)
  tx0_sd <- time_in_state |> dplyr::filter(tx == 0) |> dplyr::pull(sd_time)
  tx1_mean <- time_in_state |> dplyr::filter(tx == 1) |> dplyr::pull(mean_time)
  tx1_sd <- time_in_state |> dplyr::filter(tx == 1) |> dplyr::pull(sd_time)

  return(tibble::tibble(
    true_effect = estimand,
    tx0_mean_time = tx0_mean,
    tx0_sd_time = tx0_sd,
    tx1_mean_time = tx1_mean,
    tx1_sd_time = tx1_sd
  ))
}


#' Tidy Bootstrap Coefficient Estimates
#'
#' Summarizes bootstrap coefficient estimates by computing confidence intervals
#' and point estimates. Typically used with output from
#' \code{\link{bootstrap_model_coefs}}.
#'
#' @param boot_coefs A data frame or tibble containing bootstrap coefficient
#'   estimates, typically the output from \code{\link{bootstrap_model_coefs}}.
#'   Should have a column for bootstrap iteration ID (default: "boot_id") and
#'   one column per model coefficient.
#' @param id_col Name of the column containing bootstrap iteration IDs
#'   (default: "boot_id"). This column will be excluded from the summary.
#' @param probs Numeric vector of length 2 specifying the lower and upper
#'   quantiles for confidence intervals. Default is c(0.025, 0.975) for 95% CI.
#'   Can be changed to other levels, e.g., c(0.05, 0.95) for 90% CI.
#' @param estimate Character string specifying the point estimate to use:
#'   "median" (default), "mean", or NULL to omit point estimate. This is
#'   computed separately from the quantiles in \code{probs}.
#' @param na.rm Logical indicating whether to remove NA values before computing
#'   statistics. Default is FALSE, which will result in NA output if any
#'   bootstrap iteration failed. Set to TRUE to compute statistics on
#'   non-missing values only.
#'
#' @return A tibble with one row per coefficient containing:
#'   \itemize{
#'     \item term: Coefficient name
#'     \item estimate: Point estimate (median or mean, if requested)
#'     \item lower: Lower confidence limit (first value in \code{probs})
#'     \item upper: Upper confidence limit (second value in \code{probs})
#'   }
#'
#' @details
#' This function computes quantile-based confidence intervals from bootstrap
#' coefficient estimates. The default settings provide:
#' \itemize{
#'   \item Median as point estimate
#'   \item 2.5th percentile as lower bound
#'   \item 97.5th percentile as upper bound
#' }
#'
#' This corresponds to the percentile bootstrap confidence interval method,
#' which is appropriate when the bootstrap distribution is roughly symmetric
#' and unbiased.
#'
#' For highly skewed or biased bootstrap distributions, consider using
#' bias-corrected and accelerated (BCa) intervals instead (not currently
#' implemented in this function).
#'
#' **Usage with assess_operating_characteristics()**: When using this function
#' within fitting functions for \code{\link{assess_operating_characteristics}},
#' you'll need to add additional columns and rename \code{lower}/\code{upper}
#' to \code{conf_low}/\code{conf_high}. See examples below.
#'
#' @examples
#' \dontrun{
#' # After running bootstrap
#' boot_coefs <- bootstrap_model_coefs(
#'   model = m,
#'   data = my_data,
#'   n_boot = 1000
#' )
#'
#' # Get 95% CI with median
#' tidy_boot <- tidy_bootstrap_coefs(boot_coefs)
#'
#' # Get 90% CI with mean
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   probs = c(0.05, 0.95),
#'   estimate = "mean"
#' )
#'
#' # Get 99% CI with median
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   probs = c(0.005, 0.995)
#' )
#'
#' # Remove NA values if some bootstrap iterations failed
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   na.rm = TRUE
#' )
#'
#' # Format for assess_operating_characteristics()
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   probs = c(0.025, 0.975),
#'   estimate = "median"
#' ) |>
#'   mutate(
#'     iter = 1,
#'     analysis = "markov",
#'     se_type = "boot",
#'     std_error = NULL,
#'     statistic = NULL,
#'     p_value = NULL,
#'     conf_low = lower,
#'     conf_high = upper,
#'     .before = 1
#'   ) |>
#'   select(iter, analysis, se_type, term, estimate, std_error,
#'          statistic, p_value, conf_low, conf_high)
#' }
#'
#' @importFrom dplyr select across summarise everything
#' @importFrom tidyr pivot_longer
#' @importFrom stats quantile median
#'
#' @export
tidy_bootstrap_coefs <- function(
  boot_coefs,
  id_col = "boot_id",
  probs = c(0.025, 0.975),
  estimate = "median",
  na.rm = TRUE
) {
  # Input validation
  if (!is.data.frame(boot_coefs)) {
    stop("boot_coefs must be a data frame or tibble")
  }

  if (!id_col %in% names(boot_coefs)) {
    stop("id_col '", id_col, "' not found in boot_coefs")
  }

  if (
    !is.numeric(probs) || length(probs) != 2 || any(probs < 0) || any(probs > 1)
  ) {
    stop(
      "probs must be a numeric vector of length 2 with values between 0 and 1"
    )
  }

  if (!is.null(estimate) && !estimate %in% c("median", "mean")) {
    stop("estimate must be 'median', 'mean', or NULL")
  }

  if (!is.logical(na.rm) || length(na.rm) != 1) {
    stop("na.rm must be a single logical value (TRUE or FALSE)")
  }

  # Ensure probs are sorted
  probs <- sort(probs)

  # Remove the ID column and summarize
  result <- boot_coefs |>
    dplyr::select(-dplyr::all_of(id_col)) |>
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        list(
          estimate = ~ if (!is.null(estimate)) {
            if (estimate == "median") {
              median(.x, na.rm = na.rm)
            } else {
              mean(.x, na.rm = na.rm)
            }
          } else {
            NA_real_
          },
          lower = ~ quantile(.x, probs[1], na.rm = na.rm),
          upper = ~ quantile(.x, probs[2], na.rm = na.rm),
          n_iter = ~ n()
        ),
        .names = "{.col}___{.fn}"
      )
    )

  # Reshape to long format
  result_long <- result |>
    tidyr::pivot_longer(
      cols = dplyr::everything(),
      names_to = c("term", "statistic"),
      names_sep = "___",
      values_to = "value"
    ) |>
    tidyr::pivot_wider(
      names_from = statistic,
      values_from = value
    )

  # Remove estimate column if it's all NA (when estimate = NULL)
  if ("estimate" %in% names(result_long) && all(is.na(result_long$estimate))) {
    result_long <- result_long |>
      dplyr::select(-estimate)
  }

  # Reorder columns: term, estimate (if present), lower, upper
  if ("estimate" %in% names(result_long)) {
    result_long <- result_long |>
      dplyr::select(term, estimate, lower, upper, n_iter)
  } else {
    result_long <- result_long |>
      dplyr::select(term, lower, upper, n_iter)
  }

  return(result_long)
}


#' Format Time-to-Event Output for Competing Risks
#'
#' Prepares the output from `states_to_tte()` for competing-risks analysis
#' (e.g., Fine-Gray or Cause-Specific Cox) by identifying transitions from
#' a risk state (e.g., hospital) to an absorbing state (e.g., discharge or death).
#'
#' @description
#' This function transforms "state occupancy" data into "time-to-event" data.
#' It looks ahead to the subsequent interval to detect when a subject transitions
#' from a generic state to an event state. The interval ending at that transition
#' is assigned the corresponding event status. Rows representing the absorbing
#' states themselves are excluded, as the subject is no longer at risk.
#'
#' @param data A data frame produced by `states_to_tte()` containing at least
#'   columns `id`, `start`, `stop`, `y` (state), and `tx`.
#' @param event_status Integer code identifying the event of interest (e.g., discharge)
#'   in the `y` column (default: 1).
#' @param death_status Integer code identifying the competing event (e.g., death)
#'   in the `y` column (default: 6).
#'
#' @return A data frame where each row represents an interval at risk. The output
#'   includes a `status` column:
#'   \itemize{
#'     \item \code{0}: Censored (no event observed or lost to follow-up)
#'     \item \code{1}: Event of interest occurred at the end of this interval
#'     \item \code{2}: Competing event occurred at the end of this interval
#'   }
#'   Rows where `y` matches `event_status` or `death_status` are removed.
#'
#' @details
#' The function performs the following steps for each subject (`id`):
#' \enumerate{
#'   \item **Sorts** the data by time (`start`) to ensure chronological order.
#'   \item **Looks ahead** to the state (`y`) of the *next* interval.
#'   \item **Assigns status**:
#'     \itemize{
#'       \item If the *next* state is `event_status`, the current interval gets `status = 1`.
#'       \item If the *next* state is `death_status`, the current interval gets `status = 2`.
#'       \item Otherwise, `status = 0`.
#'     }
#'   \item **Filters**: Removes rows where the subject is already in an absorbing
#'     state (`event_status` or `death_status`), as they are no longer at risk.
#'   \item **Truncates**: Ensures only the trajectory up to the first observed event
#'     is retained per subject.
#' }
#'
#' @examples
#' \dontrun{
#' tte <- states_to_tte(trajectories)
#' # Prepare for Fine-Gray model (1=Discharge vs 2=Death)
#' fg_ready <- format_competing_risks(tte, event_status = 1, death_status = 6)
#' }
#'
#' @importFrom dplyr arrange group_by mutate lead filter ungroup select everything
#' @export
format_competing_risks <- function(
  data,
  event_status = 1,
  death_status = 6,
  version_old = TRUE
) {
  # Basic validation
  required_cols <- c("id", "start", "stop", "y", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!version_old) {
    data <- data |>
      dplyr::arrange(id, start) |>
      dplyr::group_by(id) |>
      dplyr::mutate(
        # Look at the NEXT state to see if an event happens at the end of THIS interval
        next_y = dplyr::lead(y),

        # Assign Status to the interval PRECEDING the event
        # E.g., 1 = Discharge, 2 = Death, 0 = Censored (no event or loss to follow up)
        status = dplyr::case_when(
          next_y == event_status ~ 1L,
          next_y == death_status ~ 2L,
          TRUE ~ 0L
        ),
        # 2. Drop all rows appearing AFTER the first event
        # Calculate cumulative events.
        # Once this hits 1, the current row is the event.
        cum_events = cumsum(status > 0),

        # We want to remove rows where the event occurred PREVIOUSLY.
        # lag(cum_events, default=0) will be 0 for the event row,
        # but 1 for all subsequent rows.
        prev_events = dplyr::lag(cum_events, default = 0)
      ) |>
      # Keep rows where no event has happened in the past
      dplyr::filter(prev_events == 0) |>
      dplyr::ungroup() |>
      dplyr::select(
        id,
        tx,
        start,
        stop,
        status,
        y,
        dplyr::everything(),
        -next_y,
        -cum_events,
        -prev_events
      )
  }

  if (version_old) {
    # Split data into id groups
    data_split <- split(data, data[["id"]])

    # For each individual ...
    data_split <- lapply(data_split, FUN = function(d) {
      # ...if any event occurred in this group, select the according rows...
      event_or_death <- d[["y"]] == event_status | d[["y"]] == death_status
      if (any(event_or_death)) {
        d <- d[event_or_death, ]
        # ...else select the last row
      } else {
        d <- d[nrow(d), ]
      }
      # introduce status variable, which indicates event (1) or death (2)
      d[["status"]] <- factor(ifelse(
        d[["y"]] == event_status,
        1,
        ifelse(d[["y"]] == death_status, 2, 0)
      ))
      d <- d[1, ]
      d[["etime"]] <- as.numeric(d[["stop"]])
      return(d)
    })
    return(do.call(rbind, data_split))
  }
}
