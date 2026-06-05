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


first_by_id <- function(x, id) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  tapply(x, id, function(values) values[1])
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
#' @param target_state Integer or integer vector. The state to count (default: 1, representing Home/Discharged).
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
#'
#' # Compare performance to bind_rows_fill + group_by + summarise
#' microbenchmark::microbenchmark(
#'   original = {
#'  ids <- unique(trajectories$id)
#'  result <- bind_rows_fill(lapply(ids, function(id) {
#'    rows <- trajectories$id == id
#'    data.frame(
#'      id = id,
#'      y = sum(trajectories$y[rows] == 1),
#'      tx = trajectories$tx[rows][1],
#'      check.names = FALSE
#'    )
#' return(result)
#'  }))
#'  },
#'  vectorized = {
#'  y_count <- tapply(trajectories$y == 1, trajectories$id, sum)
#'
#'  # get tx per id (first observed)
#'  tx_first <- tapply(trajectories$tx, trajectories$id, function(x) x[1])
#'
#'  return(data.frame(
#'    id = as.integer(names(y_count)),
#'    y = as.integer(y_count),
#'    tx = tx_first,
#'    check.names = FALSE
#'  ))}
#'  )
#' }
#'
#' @export
states_to_ttest <- function(data, target_state = 1) {
  # Input validation
  required_cols <- c("id", "y", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # count days in target_state per id
  y_count <- tapply(data$y %in% target_state, data$id, sum)

  # get tx per id (first observed)
  tx_first <- first_by_id(data$tx, data$id)

  data.frame(
    id = as.integer(names(y_count)),
    y = as.integer(y_count),
    tx = tx_first,
    check.names = FALSE
  )
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
#' @param target_state Integer or integer vector. The state representing "home" or "baseline" (default: 1).
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
#' trajectories <- sim_actt2_brownian()
#' drs_data <- states_to_drs(trajectories, follow_up_time = 60)
#'
#' # With different covariates
#' drs_data <- states_to_drs(trajectories, follow_up_time = 60,
#'                           covariates = c("age", "sofa", "baseline_severity"))
#'
#' # Compare performance to bind_rows_fill + group_by + summarise
#' follow_up_time <- 28
#' microbenchmark::microbenchmark(
#'   original = {
#'  ids <- unique(trajectories$id)
#'  result <- bind_rows_fill(lapply(ids, function(id) {
#'    rows <- trajectories$id == id
#'    last_not_home <- max(trajectories$time[rows][trajectories$y[rows] != 1], na.rm = TRUE)
#'    dead <- any(trajectories$y[rows] == 8)
#'    data.frame(
#'      id = id,
#'      tx = trajectories$tx[rows][1],
#'      drs = ifelse(
#'        dead,
#'        follow_up_time + 2,
#'         ifelse(
#'            follow_up_time == last_not_home,
#'            follow_up_time + 1,
#'            last_not_home
#'      )
#'    ),
#'      check.names = FALSE
#'    )
#'  }))
#'  },
#'  vectorized = {
#'  last_not_home <- tapply(
#'    ifelse(trajectories$y != 1, trajectories$time, NA_real_),
#'    trajectories$id,
#'    function(x) if (all(is.na(x))) 0 else max(x, na.rm = TRUE))
#'  dead <- tapply(trajectories$y == 8, trajectories$id, any)
#'  tx_first <- tapply(trajectories$tx, trajectories$id, function(x) x[1])
#'
#'  result <- data.frame(
#'    id = as.integer(names(last_not_home)),
#'    tx = tx_first,
#'    drs = ifelse(
#'      dead,
#'      follow_up_time + 2,
#'      ifelse(
#'        follow_up_time == last_not_home,
#'        follow_up_time + 1,
#'        last_not_home
#'      )
#'    ),
#'    check.names = FALSE
#'  )
#' result
#' })
#' }
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

  # ensure sorted so "first" is consistent
  data <- data[order(data$id, data$time), , drop = FALSE]

  # last time not in target_state (0 if never left target)
  last_not_home <- tapply(
    ifelse(!(data$y %in% target_state), data$time, NA_real_),
    data$id,
    function(x) if (all(is.na(x))) 0 else max(x, na.rm = TRUE)
  )

  dead <- tapply(data$y == death_state, data$id, any)
  tx_first <- first_by_id(data$tx, data$id)

  result <- data.frame(
    id = as.integer(names(last_not_home)),
    tx = tx_first,
    drs = ifelse(
      dead,
      follow_up_time + 2,
      ifelse(
        follow_up_time == last_not_home,
        follow_up_time + 1,
        last_not_home
      )
    ),
    check.names = FALSE
  )

  # add covariates (first value per id)
  if (!is.null(covariates)) {
    available_covs <- intersect(covariates, names(data))
    for (cov in available_covs) {
      result[[cov]] <- first_by_id(data[[cov]], data$id)
    }
  }

  return(result)
}
