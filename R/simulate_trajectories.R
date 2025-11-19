#' Simulate Individual Patient Trajectories Using Markov Model
#'
#' This function generates individual patient trajectories over time based on
#' a proportional odds model with customizable intercepts, linear predictor
#' function, and baseline patient characteristics.
#'
#' @param baseline_data A data frame containing baseline patient characteristics.
#'   Must include columns: `id`, `yprev` (initial state), and any covariates
#'   needed by `lp_function`. Typically includes `tx`, `age`, `sofa`, etc.
#'   Default is `violet_baseline`, a dataset of 10,000 patients derived from
#'   the VIOLET trial (see `?violet_baseline` for details).
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
#'   Default is `lp_violet`, which implements the VIOLET study model
#'   (see `?lp_violet` for details).
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
#' # Using all defaults (violet_baseline data, lp_violet function, VIOLET params)
#' trajectories <- sim_trajectories_markov(
#'   parameter = log(0.8),  # OR = 0.8 for treatment effect
#'   seed = 12345
#' )
#'
#' # Using default function with custom baseline data
#' custom_baseline <- data.frame(
#'   id = 1:100,
#'   yprev = sample(2:5, 100, replace = TRUE),
#'   tx = rbinom(100, 1, 0.5),
#'   age = rnorm(100, 60, 15),
#'   sofa = rpois(100, 5)
#' )
#'
#' trajectories_custom <- sim_trajectories_markov(
#'   baseline_data = custom_baseline,
#'   follow_up_time = 30,
#'   seed = 12345
#' )
#'
#' # Define a custom linear predictor function
#' my_custom_lp <- function(yprev, t, age, sofa, tx, parameter = 0, extra_params) {
#'   # Custom implementation
#'   tx_effect <- parameter * tx
#'   time_effect <- extra_params["time"] * t
#'   # ... your custom logic ...
#'   return(tx_effect + time_effect)
#' }
#'
#' trajectories_custom_lp <- sim_trajectories_markov(
#'   lp_function = my_custom_lp,
#'   parameter = log(0.75),
#'   seed = 12345
#' )
#' }
#'
#' @seealso
#' \code{\link{lp_violet}} for the default linear predictor function
#'
#' \code{\link{violet_baseline}} for the default baseline dataset
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate arrange group_by lag ungroup filter left_join
#' @importFrom stats plogis
#'
#' @export
sim_trajectories_markov <- function(
  baseline_data = violet_baseline,
  follow_up_time = 60,
  intercepts = c(-9.353417, -4.294121, -1.389221, -0.555688, 3.127056),
  lp_function = lp_violet,
  extra_params = c(
    "time" = -0.738194,
    "time'" = 0.7464006,
    "age" = 0.010321,
    "sofa" = 0.046901,
    "yprev=1" = -8.518344,
    "yprev=3" = 0,
    "yprev=4" = 1.315332,
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


#' Simulate Patient Trajectories Using Brownian Motion Latent Ordinal Model
#'
#' Generates individual patient trajectories based on a latent continuous severity
#' variable that evolves as a Gaussian random walk (Brownian motion) with drift.
#' Observed states are determined by thresholding the latent variable using an
#' ordered logistic or ordered probit model.
#'
#' @param n_patients Integer. Number of patients to simulate (default: 1000).
#' @param follow_up_time Integer. Number of days to simulate per patient (default: 60).
#' @param n_states Integer. Number of ordered states/categories (default: 6).
#' @param mu_drift Numeric. Baseline drift per day for control group (default: -0.18).
#'   Negative values indicate improvement (movement toward lower severity).
#' @param mu_treatment_effect Numeric. Additional drift for treatment group
#'   (default: 0). Negative values indicate faster improvement with treatment.
#' @param sigma_rw Numeric. Random walk innovation standard deviation per day
#'   (default: 0.12). Lower values create stronger autocorrelation.
#' @param x0_sd Numeric. Standard deviation for baseline latent severity
#'   (default: 1.0).
#' @param thresholds Numeric vector of K-1 cutpoints for the ordered model
#'   (default: c(-3.5, 0, 1, 3, 4.5) for 6 states). Higher values
#'   represent more severe states.
#' @param treatment_prob Numeric. Probability of assignment to treatment group
#'   (default: 0.5).
#' @param absorbing_state Integer. State number that represents death or other
#'   absorbing state (default: 6). Once entered, patients remain there.
#' @param drift_change_times Numeric vector of length 2 specifying when drift
#'   begins to decline and when it reaches zero (default: c(25, 50)).
#'   Drift is constant until first time, then linearly declines to 0 by second time.
#' @param latent_dist Character. The distribution used to link the latent variable
#'   to observed states. Options are "logistic" (default) or "normal".
#'   "logistic" uses the standard logistic CDF (plogis), while "normal" uses
#'   the standard normal CDF (pnorm).
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#'
#' @return A data frame (tibble) with columns:
#'   - id: patient identifier
#'   - time: time/day (1 to follow_up_time)
#'   - tx: treatment assignment (0 = control, 1 = treatment)
#'   - y: observed ordinal state at time
#'   - x: latent continuous severity (NA after entering absorbing state)
#'
#' @details
#' This function implements a latent variable model where:
#' 1. Each patient has an unobserved continuous "severity" X(t) that evolves
#'    as a Gaussian random walk: X(t) = X(t-1) + μ(t) + ε, where ε ~ N(0, σ²)
#' 2. The drift μ(t) depends on treatment assignment and time
#' 3. Observed states Y(t) are generated by thresholding X(t).
#'    If `latent_dist = "logistic"`, P(Y ≤ k) = logit⁻¹(c_k - X).
#'    If `latent_dist = "normal"`, P(Y ≤ k) = Φ(c_k - X).
#' 4. Once a patient enters the absorbing state (death), they remain there.
#'
#' @examples
#' \dontrun{
#' # Standard Logistic model
#' traj_logit <- sim_trajectories_brownian(
#'   n_patients = 1000,
#'   latent_dist = "logistic"
#' )
#'
#' # Probit model which approximates viral respiratory trajectories okay-ishly
#' # Because standard normal is narrower than logistic, thresholds need to be closer together, drift smaller, varaince lower etc.
#' traj_norm <- sim_trajectories_brownian(n_patients = 1000,
#'   latent_dist = "normal",
#'   thresholds = c(-3, -1.0, 0.5, 2, 3),
#'   mu_drift = -0.12,
#'   drift_change_times = c(40, 50),
#'   sigma_rw = 0.02)
#' }
#'
#' @importFrom stats rnorm plogis pnorm
#' @importFrom tibble tibble
#' @importFrom dplyr mutate left_join select arrange group_by ungroup filter
#' @importFrom tidyr pivot_longer
#'
#' @export
sim_trajectories_brownian <- function(
  n_patients = 1000,
  follow_up_time = 60,
  n_states = 6,
  mu_drift = -0.18,
  mu_treatment_effect = 0,
  sigma_rw = 0.12,
  x0_sd = 1.0,
  thresholds = c(-3.5, 0, 1, 3, 4.5),
  treatment_prob = 0.5,
  absorbing_state = 6,
  drift_change_times = c(25, 50),
  latent_dist = c("logistic", "normal"),
  seed = NULL
) {
  # Input validation
  if (n_patients < 1 || !is.numeric(n_patients)) {
    stop("n_patients must be a positive integer")
  }

  if (follow_up_time < 1 || !is.numeric(follow_up_time)) {
    stop("follow_up_time must be a positive integer")
  }

  if (n_states < 2 || !is.numeric(n_states)) {
    stop("n_states must be at least 2")
  }

  if (length(thresholds) != (n_states - 1)) {
    stop(
      "thresholds must have length n_states - 1 (",
      n_states - 1,
      "), but has length ",
      length(thresholds)
    )
  }

  if (!is.numeric(thresholds) || !all(diff(thresholds) > 0)) {
    stop("thresholds must be a strictly increasing numeric vector")
  }

  if (treatment_prob < 0 || treatment_prob > 1) {
    stop("treatment_prob must be between 0 and 1")
  }

  if (absorbing_state < 1 || absorbing_state > n_states) {
    stop("absorbing_state must be between 1 and n_states")
  }

  if (
    length(drift_change_times) != 2 ||
      drift_change_times[1] >= drift_change_times[2]
  ) {
    stop(
      "drift_change_times must be a vector of length 2 with increasing values"
    )
  }

  # Validate and select distribution
  latent_dist <- match.arg(latent_dist)

  # Select the appropriate Cumulative Distribution Function
  # Assigning this before the loop improves performance
  cdf_fun <- if (latent_dist == "logistic") {
    stats::plogis
  } else {
    stats::pnorm
  }

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Randomly assign treatment
  treatment <- sample(
    c(0, 1),
    n_patients,
    replace = TRUE,
    prob = c(1 - treatment_prob, treatment_prob)
  )

  # Initialize matrices for latent X and observed Y
  # Need follow_up_time + 1 columns to include time 0
  X <- matrix(NA_real_, n_patients, follow_up_time + 1)
  Y <- matrix(NA_integer_, n_patients, follow_up_time + 1)

  # Simulate trajectories for each patient
  for (i in 1:n_patients) {
    # Initialize latent severity at day 0 (baseline)
    X[i, 1] <- rnorm(1, 0, x0_sd)

    # Generate observation at day 0 using selected CDF
    pcat <- diff(c(0, cdf_fun(thresholds - X[i, 1]), 1))
    Y[i, 1] <- sample.int(n_states, 1, prob = pcat)

    # Simulate remaining days (1 to follow_up_time)
    for (t in 2:(follow_up_time + 1)) {
      # Check if patient is in absorbing state
      if (Y[i, t - 1] == absorbing_state) {
        Y[i, t] <- absorbing_state
        X[i, t] <- NA_real_ # Latent variable no longer evolves
      } else {
        # Calculate individual-specific drift
        mu_i <- mu_drift + treatment[i] * mu_treatment_effect

        # Apply time-varying drift
        if (t <= drift_change_times[1]) {
          mu_t <- mu_i
        } else if (t <= drift_change_times[2]) {
          # Linear decline from mu_i to 0
          mu_t <- mu_i *
            (drift_change_times[2] - t) /
            (drift_change_times[2] - drift_change_times[1])
        } else {
          mu_t <- 0
        }

        # Update latent severity (random walk with drift)
        X[i, t] <- rnorm(1, mean = X[i, t - 1] + mu_t, sd = sigma_rw)

        # Generate observation using selected CDF
        pcat <- diff(c(0, cdf_fun(thresholds - X[i, t]), 1))
        Y[i, t] <- sample.int(n_states, 1, prob = pcat)
      }
    }
  }

  # Convert matrices to long format data frame
  # Create data frame for latent X
  dat_latent <- as.data.frame(X)
  colnames(dat_latent) <- as.character(0:follow_up_time)
  dat_latent <- dat_latent |>
    dplyr::mutate(
      id = seq_len(n_patients),
      tx = treatment
    ) |>
    tidyr::pivot_longer(
      cols = -c(id, tx),
      names_to = "time",
      values_to = "x"
    ) |>
    dplyr::mutate(time = as.integer(time))

  # Create data frame for observed Y
  dat_observed <- as.data.frame(Y)
  colnames(dat_observed) <- as.character(0:follow_up_time)
  dat_observed <- dat_observed |>
    dplyr::mutate(
      id = seq_len(n_patients),
      tx = treatment
    ) |>
    tidyr::pivot_longer(
      cols = -c(id, tx),
      names_to = "time",
      values_to = "y"
    ) |>
    dplyr::mutate(time = as.integer(time)) |>
    dplyr::group_by(id) |>
    dplyr::mutate(yprev = dplyr::lag(y, 1)) |>
    dplyr::ungroup() |>
    dplyr::filter(time > 0)

  # Combine latent and observed data
  result <- dat_observed |>
    dplyr::left_join(
      dat_latent |> dplyr::select(id, time, x),
      by = c("id", "time")
    ) |>
    dplyr::arrange(id, time)

  return(tibble::tibble(result))
}
