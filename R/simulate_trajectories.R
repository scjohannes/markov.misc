#' Simulate Individual Patient Trajectories Using Markov Model
#'
#' This function generates individual patient trajectories over time based on
#' a proportional odds model with customizable intercepts, linear predictor
#' function, and baseline patient characteristics.
#'
#' @param baseline_data A data frame containing baseline patient characteristics.
#'   Must include columns: `id`, `yprev` (initial state), and any covariates
#'   needed by `lp_function`. Typically includes `tx`, `age`, `sofa`, etc.
#' @param follow_up_time Integer. Number of time periods to simulate (default: 60).
#' @param intercepts Numeric vector of intercepts for the proportional odds model.
#'   Should have length = (number of states - 1). Default values are from VIOLET
#'   study with 6-state expansion:
#'   c(-9.353, -4.294, -1.389, -0.556, 3.127)
#' @param lp_function A function that calculates the linear predictor for each
#'   patient at each time point. Should accept parameters:
#'   - yprev: previous state
#'   - t: current time
#'   - tx: treatment indicator
#'   - Additional named arguments matching columns in `baseline_data`
#'   - parameter: treatment effect (log odds ratio)
#'   - extra_params: named vector of coefficients
#'   Must return a numeric value representing the linear predictor.
#' @param extra_params Named numeric vector of model coefficients used by
#'   `lp_function`. Default values from VIOLET study include time, spline,
#'   age, sofa, and previous state effects with interactions.
#' @param parameter Numeric. Treatment effect on log odds ratio scale (default: 0,
#'   meaning odds ratio = 1).
#' @param absorbing_states Integer vector of states that are absorbing (once
#'   entered, patients cannot leave). Default is 6 (death).
#' @param seed Integer. Random seed for reproducibility (default: NULL, no seed set).
#' @param covariate_names Character vector of covariate names from `baseline_data`
#'   to pass to `lp_function`. Default is c("age", "sofa", "tx"). The function
#'   will automatically detect and pass these to `lp_function`.
#'
#' @return A data frame with columns:
#'   - id: patient identifier
#'   - time: time point (1 to follow_up_time)
#'   - y: observed state at this time
#'   - yprev: state at previous time point
#'   - all columns from baseline_data
#'
#' @details
#' The function implements a discrete-time Markov model where transition
#' probabilities are determined by a proportional odds model. At each time step:
#' 1. Calculate linear predictor for each active patient
#' 2. Apply intercepts to get cumulative probabilities
#' 3. Convert to state probabilities
#' 4. Sample next state for each patient
#' 5. Patients in absorbing states remain there
#'
#' States are assumed to be ordered integers (e.g., 1=Home, 2=Hospital mild,
#' 3=Hospital oxygen, 4=Hospital NIV, 5=Ventilator, 6=Death). The model
#' uses a proportional odds structure where higher intercepts represent
#' more severe health states.
#'
#' @examples
#' \dontrun{
#' # Using default VIOLET parameters
#' baseline <- data.frame(
#'   id = 1:100,
#'   yprev = sample(2:5, 100, replace = TRUE),
#'   tx = rbinom(100, 1, 0.5),
#'   age = rnorm(100, 60, 15),
#'   sofa = rpois(100, 5)
#' )
#'
#' # Define a simple linear predictor function
#' my_lp <- function(yprev, t, age, sofa, tx, parameter = 0, extra_params) {
#'   tx_effect <- parameter * tx
#'   sofa_effect <- extra_params["sofa"] * sofa
#'   age_effect <- extra_params["age"] * age
#'   time_effect <- extra_params["time"] * t
#'   time_2_effect <- extra_params["time'"] * pmax(t - 2, 0)
#'
#'   yprev_effect <- 0
#'   yprev_time_effect <- 0
#'   if (yprev != 2) {
#'     yprev_coef_name <- paste0('yprev=', yprev)
#'     if (yprev_coef_name %in% names(extra_params)) {
#'       yprev_effect <- extra_params[[yprev_coef_name]]
#'     }
#'     yprev_time_name <- paste0('yprev=', yprev, ' * time')
#'     if (yprev_time_name %in% names(extra_params)) {
#'       yprev_time_effect <- extra_params[[yprev_time_name]] * t
#'     }
#'   }
#'
#'   lp <- as.numeric(tx_effect + sofa_effect + age_effect +
#'                    time_effect + time_2_effect +
#'                    yprev_effect + yprev_time_effect)
#'   return(lp)
#' }
#'
#' trajectories <- simulate_trajectories(
#'   baseline_data = baseline,
#'   follow_up_time = 30,
#'   lp_function = my_lp,
#'   parameter = log(0.8),  # OR = 0.8
#'   seed = 12345
#' )
#' }
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate arrange group_by lag ungroup filter left_join
#' @importFrom stats plogis
#'
#' @export
simulate_trajectories <- function(
  baseline_data,
  follow_up_time = 60,
  intercepts = c(-9.353417, -4.294121, -1.389221, -0.555688, 3.127056),
  lp_function,
  extra_params = c(
    "time" = -0.738194,
    "time'" = 0.7464006,
    "age" = 0.010321,
    "sofa" = 0.046901,
    "yprev=1" = -9.464827,
    "yprev=3" = 0.438444,
    "yprev=4" = 0.657666,
    "yprev=5" = 6.576662,
    "yprev=1 * time" = 0,
    "yprev=3 * time" = 0,
    "yprev=4 * time" = 0,
    "yprev=5 * time" = 0
  ),
  parameter = 0,
  absorbing_states = 6,
  seed = NULL,
  covariate_names = c("age", "sofa", "tx")
) {
  # Input validation
  if (!is.data.frame(baseline_data)) {
    stop("baseline_data must be a data frame")
  }

  required_cols <- c("id", "yprev")
  missing_cols <- setdiff(required_cols, names(baseline_data))
  if (length(missing_cols) > 0) {
    stop(
      "baseline_data must contain columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.numeric(intercepts) || length(intercepts) < 1) {
    stop("intercepts must be a numeric vector with at least one element")
  }

  if (!is.function(lp_function)) {
    stop("lp_function must be a function")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Setup
  N <- nrow(baseline_data)
  times <- 1:follow_up_time
  n_states <- length(intercepts) + 1

  # Initialize state matrix
  # Rows = patients, Columns = time points (0 to follow_up_time)
  state_matrix <- matrix(
    as.numeric(baseline_data$yprev),
    nrow = N,
    ncol = length(times) + 1,
    dimnames = list(baseline_data$id, 0:max(times))
  )

  # Detect which covariates from baseline_data need to be passed to lp_function
  available_covariates <- intersect(covariate_names, names(baseline_data))
  if (length(available_covariates) == 0) {
    warning("None of the specified covariate_names found in baseline_data")
  }

  # Main simulation loop
  for (t in times) {
    # Current states
    yprev <- state_matrix[, t]

    # Identify patients not in absorbing states
    active_idx <- which(!yprev %in% absorbing_states)

    # If all patients are in absorbing states, fill remaining time and exit
    if (length(active_idx) == 0) {
      if (t <= max(times)) {
        for (i in (t + 1):ncol(state_matrix)) {
          state_matrix[, i] <- state_matrix[, t]
        }
      }
      break
    }

    # Extract data for active patients
    yprev_active <- yprev[active_idx]
    X_active <- baseline_data[active_idx, , drop = FALSE]

    # Build argument list for lp_function
    # This allows dynamic passing of covariates
    lp_args <- list(
      yprev = yprev_active
    )

    # MoreArgs for scalar values that don't vary across patients
    more_args <- list(
      t = t,
      parameter = parameter,
      extra_params = extra_params
    )

    # Add covariates dynamically
    for (cov_name in available_covariates) {
      lp_args[[cov_name]] <- X_active[[cov_name]]
    }

    # Special handling for treatment variable (might be tx_val in function)
    if ("tx" %in% names(X_active) && !"tx" %in% names(lp_args)) {
      lp_args[["tx"]] <- X_active[["tx"]]
    }

    # Calculate linear predictor for each active patient
    lp <- do.call(
      mapply,
      c(
        list(FUN = lp_function, SIMPLIFY = TRUE, USE.NAMES = FALSE),
        lp_args,
        list(MoreArgs = more_args)
      )
    )

    # Calculate transition probabilities using proportional odds model
    # For each intercept threshold, calculate cumulative probability
    thresholds_matrix <- outer(intercepts, lp, "+")
    cum_probs <- t(plogis(thresholds_matrix)) # (n_active x n_states-1)

    # Convert cumulative probabilities to individual state probabilities
    prob_matrix <- cbind(cum_probs, 1) - cbind(0, cum_probs)

    # States in descending order (highest to lowest)
    states_vec <- n_states:1

    # Sample next state for each active patient
    y_new_active <- apply(prob_matrix, 1, function(p) {
      sample(states_vec, size = 1, prob = p)
    })

    # Update state matrix
    # First, carry forward previous state for all patients
    state_matrix[, t + 1] <- state_matrix[, t]
    # Then update only active patients with new states
    state_matrix[active_idx, t + 1] <- y_new_active
  }

  # Convert matrix to long format data frame
  result <- as.data.frame(state_matrix) |>
    tibble::rownames_to_column(var = "id") |>
    tidyr::pivot_longer(
      cols = -id,
      names_to = "time",
      values_to = "y"
    ) |>
    dplyr::mutate(
      id = as.integer(id),
      time = as.integer(time)
    ) |>
    dplyr::left_join(
      baseline_data,
      by = "id"
    ) |>
    dplyr::arrange(id, time) |>
    dplyr::group_by(id) |>
    dplyr::mutate(yprev = dplyr::lag(y)) |>
    dplyr::ungroup() |>
    dplyr::filter(time > 0) # Remove time 0 (initial state)

  return(result)
}
