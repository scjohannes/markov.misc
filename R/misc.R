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
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, follow_up_time = 60,
#'                                         lp_function = my_lp)
#' tte_data <- states_to_tte(trajectories, follow_up_time = 60)
#'
#' # With different covariates
#' tte_data <- states_to_tte(trajectories, follow_up_time = 60,
#'                           covariates = c("age", "sofa", "baseline_severity"))
#' }
#'
#' @importFrom dplyr group_by arrange summarise if_else select across all_of left_join mutate
#'
#' @export
states_to_tte <- function(
    data,
    covariates = c("age", "sofa")
) {
  # Input validation
  required_cols <- c("id", "time", "y", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate Andersen-Gill count data format
  data$state_change <- ave(data[["y"]], data[["id"]],
                           FUN = function(x) c(0, diff(x)))

  data <- data[data$state_change != 0, ]

  data[["start"]] <- lag(data[["time"]], default = 0)
  data[["stop"]] <- data[["time"]]
  data <- data[, !(colnames(data) %in% c("state_change", "time"))]

  # Add covariates if specified
  if (!is.null(covariates)) {
    available_covs <- intersect(covariates, colnames(data))
    data <- data[, colnames(data) %in% c(
      c("id", "y", "tx", "start", "stop", "yprev"),
      available_covs)
    ]
  }

  return(result = data)
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
    dplyr::filter(y == target_state) |>
    dplyr::select(tx, time, sop) |>
    tidyr::pivot_wider(
      names_from = tx,
      values_from = sop,
      names_prefix = "tx_"
    ) |>
    dplyr::summarise(true_effect = sum(tx_1) - sum(tx_0)) |>
    dplyr::pull(true_effect)

  # Calculate time spent in target state for each group
  time_in_state <- data |>
    dplyr::filter(y == target_state) |>
    dplyr::group_by(id, tx) |>
    dplyr::summarise(time_in_state = n(), .groups = "drop") |>
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
