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
    # We pass yprev as a factor to ensure compatibility with lp functions that expect it
    lp_args <- list(
      yprev = factor(yprev_active, levels = 1:n_states)
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
#' @param mu_drift Numeric. Baseline (constant) drift per day for control group (default: -0.18).
#'   Negative values indicate improvement (movement toward lower severity).
#' @param mu_quad Numeric. Quadratic coefficient for time-varying drift (default: 0).
#'   The drift at time t is: `mu_drift + mu_quad * t + mu_quad2 * t^2`.
#'   Positive values cause the drift to become less negative (slower improvement) over time.
#' @param mu_quad2 Numeric. Second-order quadratic coefficient for drift (default: 0).
#'   Allows for smooth non-linear drift dynamics (e.g., initial increase then decrease).
#'   For example, setting `mu_quad > 0` and `mu_quad2 < 0` creates a drift that
#'   initially increases then decreases.
#' @param mu_treatment_effect Numeric. Additional constant drift for treatment group
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
#' @param allowed_start_state Integer vector or NULL. Defines one or more states allowed at
#'   time 0. If NULL, participants can start in all states (default: 2:5).
#' @param absorbing_state Integer or NULL. State number that represents death or other
#'   absorbing state (default: 6). Once entered, patients remain there.
#'   Set to NULL if no absorbing state is desired.
#' @param drift_change_times Numeric vector of length 2 specifying when drift
#'   begins to decline and when it reaches zero (default: c(25, 50)).
#'   Drift is constant until first time, then linearly declines to 0 by second time.
#'   Note: This piecewise adjustment is applied AFTER the quadratic time effects
#'   (`mu_quad`, `mu_quad2`), allowing both smooth and piecewise dynamics.
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
#' 2. The drift μ(t) depends on treatment assignment and time:
#'    μ(t) = (mu_drift + mu_quad * t + mu_quad2 * t² + tx * mu_treatment_effect) * piecewise_factor
#'    where piecewise_factor applies the drift_change_times logic (1 before first time,
#'    linearly declining to 0 between first and second time, 0 after).
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
#'
#' # Smooth non-linear drift: initial worsening that reverses to improvement
#' # The quadratic terms create a U-shaped or inverted-U drift curve
#' traj_nonlinear <- sim_trajectories_brownian(
#'   n_patients = 1000,
#'   mu_drift = 0.05,       # slight initial worsening
#'   mu_quad = -0.002,      # drift becomes more negative over time
#'   mu_quad2 = 0,          # no second-order curvature
#'   drift_change_times = c(60, 61)  # disable piecewise decline
#' )
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
  mu_quad = 0,
  mu_quad2 = 0,
  mu_treatment_effect = 0,
  sigma_rw = 0.12,
  x0_sd = 1.0,
  thresholds = c(-3.5, 0, 1, 3, 4.5),
  treatment_prob = 0.5,
  allowed_start_state = 2:5,
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

  if (!is.null(allowed_start_state)) {
    if (!is.integer(allowed_start_state)) {
      stop("allowed_start_state must be NULL or an integer")
    }
    if (!all(allowed_start_state %in% 1:n_states)) {
      stop("All allowed_start_state must be between 1 and n_states")
    }
  }

  if (!is.null(absorbing_state)) {
    if (absorbing_state < 1 || absorbing_state > n_states) {
      stop("absorbing_state must be between 1 and n_states")
    }
  }

  if (
    !is.null(drift_change_times) &
      (length(drift_change_times) != 2 ||
        drift_change_times[1] >= drift_change_times[2])
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
    if (!is.null(allowed_start_state)) {
      pcat[-allowed_start_state] <- 0 # states nobody should start in
    }
    Y[i, 1] <- sample.int(n_states, 1, prob = pcat)

    # Simulate remaining days (1 to follow_up_time)
    for (t in 2:(follow_up_time + 1)) {
      # Check if patient is in absorbing state
      if (!is.null(absorbing_state) && Y[i, t - 1] == absorbing_state) {
        Y[i, t] <- absorbing_state
        X[i, t] <- NA_real_ # Latent variable no longer evolves
      } else {
        # Calculate time-varying base drift with quadratic terms
        # Note: t here is 1-indexed (t=2 means day 1), so we use (t-1) for actual day
        day <- t - 1
        mu_base <- mu_drift + mu_quad * day + mu_quad2 * day^2

        # Add treatment effect to base drift (so it also gets scaled by piecewise factor)
        mu_i <- mu_base + treatment[i] * mu_treatment_effect

        # Apply piecewise time-varying factor (for backward compatibility)
        if (is.null(drift_change_times)) {
          mu_t <- mu_i
        } else if (day <= drift_change_times[1]) {
          mu_t <- mu_i
        } else if (day <= drift_change_times[2]) {
          # Linear decline from mu_i to 0
          mu_t <- mu_i *
            (drift_change_times[2] - day) /
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


#' Simulate Patient Trajectories with Non-Linear (Quadratic) Dynamics
#'
#' Generates deterministic patient trajectories where the underlying latent
#' variable follows a quadratic curve (time + time^2) with random walk noise.
#'
#' Default values are for a **mild to moderate respiratory viral disease population**,
#' the default parameters simulate a general recovery trend (negative linear
#' slope).
#'
#' @param n_patients Integer. Number of patients (default: 1000).
#' @param follow_up_time Integer. Days to simulate (default: 60).
#' @param n_states Integer. Number of states (default: 6).
#' @param mu_linear Numeric. Initial linear slope at t=0 (default: -0.1).
#'   Negative values indicate improvement.
#' @param linear_sd Numeric. SD of the linear slope between patients (default: 0.02).
#' @param mu_quad Numeric. Quadratic coefficient (acceleration) (default: 0).
#'   The default is 0 (linear recovery). Positive values combined with a negative
#'   linear slope will cause patients to improve, bottom out, and then worsen (relapse).
#' @param quad_sd Numeric. SD of the quadratic term between patients (default: 0.0006).
#' @param mu_treatment_effect Numeric. Effect of treatment on the *linear* slope
#'   (default: 0). If set to a negative value (e.g., -0.05), treatment accelerates recovery.
#' @param sigma_rw Numeric. Daily random walk noise (default: 0.25).
#' @param x0_mean Numeric. Baseline severity mean (default: 0.1).
#'   Represents mild/moderate starting severity.
#' @param x0_sd Numeric. Baseline severity SD (default: 1.0).
#' @param thresholds Numeric vector. Cutpoints for states (default: c(-1.8, -0.5, 0.5, 1.2, 2.3)).
#' @param treatment_prob Numeric. Probability of treatment (default: 0.5).
#' @param seed Integer. Seed for reproducibility.
#'
#' @return A tibble with id, time, tx, y, x.
#'
#' @importFrom stats rnorm rbinom
#' @importFrom tibble tibble
#' @importFrom dplyr mutate left_join arrange
#' @importFrom tidyr pivot_longer
#'
#' @export
sim_trajectories_deterministic <- function(
  n_patients = 1000,
  follow_up_time = 60,
  n_states = 6,
  mu_linear = -0.1,
  linear_sd = 0.02,
  mu_quad = 0,
  quad_sd = 0.0006,
  mu_treatment_effect = 0,
  sigma_rw = 0.25,
  x0_mean = 0.1,
  x0_sd = 1.0,
  thresholds = c(-1.8, -0.5, 0.5, 1.2, 2.3),
  treatment_prob = 0.5,
  seed = NULL
) {
  # --- Validation ---
  if (length(thresholds) != n_states - 1) {
    stop("Thresholds length error")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Setup ---
  X <- matrix(NA_real_, n_patients, follow_up_time + 1)
  Y <- matrix(NA_integer_, n_patients, follow_up_time + 1)
  treatment <- rbinom(n_patients, 1, treatment_prob)

  # --- Assign Individual Parameters ---
  # 1. Linear term (Velocity): Varies by patient + Treatment effect
  # Note: We apply treatment only to the linear term (initial response)
  beta1_i <- rnorm(n_patients, mu_linear, linear_sd) +
    (treatment * mu_treatment_effect)

  # 2. Quadratic term (Acceleration/Curvature): Varies by patient
  beta2_i <- rnorm(n_patients, mu_quad, quad_sd)

  # --- Simulation Loop ---
  for (i in 1:n_patients) {
    # Baseline
    X[i, 1] <- rnorm(1, x0_mean, x0_sd)
    Y[i, 1] <- findInterval(X[i, 1], thresholds) + 1

    absorbed <- FALSE
    if (Y[i, 1] == n_states) {
      absorbed <- TRUE
    }

    for (t in 2:(follow_up_time + 1)) {
      if (absorbed) {
        Y[i, t] <- n_states
        X[i, t] <- NA
      } else {
        # Calculate the deterministic drift for this specific timepoint
        # The derivative of (b1*t + b2*t^2) is (b1 + 2*b2*t)
        # Note: t here represents the step index.
        # Since we are adding to the *previous* X, we add the instantaneous slope.
        current_drift <- beta1_i[i] + (2 * beta2_i[i] * (t - 1))

        # Update Latent X
        X[i, t] <- X[i, t - 1] + current_drift + rnorm(1, 0, sigma_rw)

        # Map to Observed Y
        state <- findInterval(X[i, t], thresholds) + 1

        if (state >= n_states) {
          Y[i, t] <- n_states
          absorbed <- TRUE
        } else {
          Y[i, t] <- state
        }
      }
    }
  }

  # --- Formatting ---
  dat_x <- as.data.frame(X)
  colnames(dat_x) <- as.character(0:follow_up_time)
  dat_x <- dat_x |>
    dplyr::mutate(id = 1:n_patients, tx = treatment) |>
    tidyr::pivot_longer(-c(id, tx), names_to = "time", values_to = "x")

  dat_y <- as.data.frame(Y)
  colnames(dat_y) <- as.character(0:follow_up_time)
  dat_y <- dat_y |>
    dplyr::mutate(id = 1:n_patients, tx = treatment) |>
    tidyr::pivot_longer(-c(id, tx), names_to = "time", values_to = "y")

  result <- dplyr::left_join(dat_y, dat_x, by = c("id", "tx", "time")) |>
    dplyr::mutate(time = as.integer(time)) |>
    dplyr::arrange(id, time)

  return(tibble::tibble(result))
}


#' Generate recurrent event times for individuals (helper for latent TTE DGM)
#'
#' Simulate recurrent event times in long format for a set of individuals.
#' For each individual the waiting times between events are drawn from
#' exponential distributions with state-specific rates and then cumulated to
#' event times. Event times that exceed the administrative follow-up are
#' censored (removed).
#'
#' @param id Optional numeric vector of subject identifiers. If supplied, these
#'   identifiers are repeated `max_events` times internally (one block of
#'   rates per subject). If missing, identifiers 1:n are generated where
#'   `n` is provided or inferred.
#' @param n Optional integer. Number of participants to simulate. If missing,
#'   `n` is set to `length(id)`. If both `id` and `n` are provided their
#'   lengths must match.
#' @param dist Character. Name of the waiting-time distribution. Currently only
#'   "Exponential" is supported (default).
#' @param param Numeric. Parameter vector for the waiting-time distribution.
#'   The first element is taken as the baseline rate \u03bb (lambda) used to
#'   construct the per-event rates. Only `param[1]` is used by the current
#'   implementation.
#' @param b Numeric scalar. Increment added to the rate for every subsequent
#'   event (autoregressive Poisson / accelerating rates). Rate for event j is
#'   lambda + b * (j - 1). Default `0` (constant rates).
#' @param follow_up Numeric scalar. Administrative follow-up time. Event times
#'   \code{>= follow_up} are censored and not returned. Default `60`.
#' @param max_events Optional integer. Maximum number of events to generate per
#'   individual. If NULL (default) the function chooses a suitable `max_events`
#'   between 3 and 20 so that the probability of observing all events within
#'   `follow_up` is very small (threshold 1e-4). When provided, the function
#'   uses the given value.
#'
#' @return A two-column numeric matrix / data.frame with columns
#'   - \code{id}: subject identifier (numeric)
#'   - \code{event_time}: event time (numeric) for events strictly less than \code{follow_up}
#'
#' @details
#' Algorithm summary:
#' 1. Construct a vector of per-event rates of length \code{max_events}:
#'    \eqn{rate_j = lambda + b * (j - 1)}.
#' 2. If \code{max_events} is not supplied, choose the smallest value in
#'    3:20 for which the probability of experiencing all \code{max_events}
#'    within \code{follow_up} is < 1e-4. When rates are constant the gamma
#'    distribution (pgamma) is used; otherwise \code{sdprisk::phypoexp} is used.
#' 3. Simulate waiting times with \code{rexp(n * max_events, rate = rates)}.
#'    Rates are recycled to length \code{n * max_events} so that each subject
#'    receives the block of per-event rates in order.
#' 4. For each subject compute the cumulative sums of waiting times to obtain
#'    event times, and retain only those event times < \code{follow_up}.
#'
#' Notes and caveats:
#' - Event times exactly equal to \code{follow_up} are excluded (strict <).
#'
#' @examples
#' # small example
#' recurr_event(id = 1:3, param = 0.1, b = 0, follow_up = 10, max_events = 5)
#'
#' # let function choose max_events automatically
#' recurr_event(id = 1:10, param = 0.05, b = 0.01, follow_up = 30)
#'
#' @seealso \code{\link[sdprisk]{phypoexp}} for hypoexponential CDF used internally
#' @keywords recurr_event
#' @importFrom sdprisk phypoexp
#' @export

recurr_event <- function(
  id,
  dist = "Exponential",
  param,
  b = 0,
  follow_up = 60,
  max_events = NULL
) {
  # prepare variables
  n <- length(id)

  lambda_i <- param

  #______________________________________________________________________________#
  #____Find max_events and corresponding rates___________________________________#
  #______________________________________________________________________________#
  # Find max_events so that the probability of experiencing all events within
  # the study interval is very small (here set to 0.0001, can be changed below)
  if (is.null(max_events)) {
    for (max_events in 3:20) {
      rates <- numeric(max_events)

      if (length(param) == 1) {
        for (j in 1:max_events) {
          rates[j] <- lambda_i + b * (j - 1)
        }
      } else {
        if (length(param) == max_events) {
          rates <- param
        } else {
          stop("Length of param must be one or equal to max_events.")
        }
      }

      if (length(unique(rates)) == 1) {
        p <- pgamma(follow_up, shape = length(rates), rate = unique(rates))
      } else {
        p <- phypoexp(follow_up, rate = rates)
      }

      if (p < 0.0001) {
        break()
      }

      if (max_events == 20) {
        warning(paste0(
          "max_events = 20, but propability of experiencing all events within time ",
          follow_up,
          " is still ",
          p,
          ". To lower this probability, try different follow-up-time or lambda_i or override max_events."
        ))
      }
    }
  } else {
    # define rates for waiting time prior to event i
    rates <- numeric(max_events)
    for (i in 1:max_events) {
      rates[i] <- lambda_i + b * (i - 1)
    }
  }

  #______________________________________________________________________________#
  #____ Model waiting times between events and store as vector   ________________#
  #______________________________________________________________________________#
  # Assign patient ID for later dataframe
  if (missing(id)) {
    id <- rep(1:n, each = max_events)
  } else {
    id <- rep(id, each = max_events)
  }

  # generate the waiting times between events for each individual
  rates_rep <- rep(rates, times = n)
  waiting_times_orig <- rexp(n * max_events, rate = rates_rep)

  # vector of cumulative sums of individuals, i.e., time points of state change
  event_times <- ave(waiting_times_orig, id, FUN = cumsum)
  names(event_times) <- id

  #   - censor intervals which exceed the follow-up (administrative censoring)
  censored_event_times <- event_times[event_times < follow_up]

  return(cbind(
    id = as.numeric(names(censored_event_times)),
    event_time = censored_event_times
  ))
}


#' Simulate Patient Trajectories Using Latent Time-to-Event Model
#'
#' Generates individual patient trajectories using a latent time-to-event mechanism
#' with recurrent events. Each patient transitions through ordered states based on
#' state-specific event rates, with support for treatment effects via hazard ratios.
#'
#' @param baseline_data A data frame containing baseline patient characteristics.
#'   Must include columns: `id`, `tx` (treatment indicator), and `state` (initial state).
#'   Default is `NULL` (user must provide). Each row represents one patient.
#' @param baseline_states Integer vector. States that patients can start in at time 0.
#'   Default is NULL, which allows starting in any state. If provided, only these
#'   states will be assigned at baseline.
#' @param prob Probabilities for each baseline_state.
#' @param states Integer vector. Ordered states in the model (default: 1:6).
#'   States should be numbered consecutively.
#' @param absorbing_states Integer vector. States that are absorbing (once entered,
#'   patients cannot leave). Default is 6 (typically death).
#' @param follow_up_time Integer. Number of days to simulate per patient (default: 60).
#' @param param Numeric vector. Baseline event rates (hazards) for each state.
#'   Length must equal `length(states)`. These represent the rate parameter λ for
#'   exponential waiting times in the control group.
#' @param hazard_ratios List of numeric vectors. Treatment effects as multiplicative
#'   factors on the baseline rates for each treatment arm.
#'   Each list element corresponds to a treatment arm (excluding control, which is 1.0).
#'   Length of each vector must equal `length(states)`.
#'   Default is `list(c(1, 1, 1, 1, 1, 1))` (no treatment effect).
#' @param b Numeric. Autoregressive coefficient for recurrent events within a state
#'   (default: 0). The rate for event j in a given state is:
#'   rate_j = param + b * (j - 1). Positive b means events become more likely over time
#'   (Poisson process acceleration); b = 0 means independent events.
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#'
#' @return A data frame (tibble) with columns:
#'   - id: patient identifier
#'   - tx: treatment assignment (0 = control, 1 = treatment)
#'   - time: time/day (1 to follow_up_time)
#'   - y: observed state at time
#'   - yprev: state at previous time point
#'
#' @details
#' This function implements a latent time-to-event data generating mechanism:
#' 1. For each patient, state, and treatment arm, waiting times to recurrent events
#'    are generated from exponential distributions with state-specific rates.
#' 2. Rates are adjusted by the corresponding hazard ratio (treatment effect).
#' 3. Event times are cumulated to determine when state transitions occur.
#' 4. Events occurring after follow-up are censored (administrative censoring).
#' 5. Each patient's trajectory is expanded to include all days 1 to follow_up_time,
#'    with states carried forward between event times (last observation carried forward, LOCF).
#' 6. Once an absorbing state is entered, no further transitions occur.
#'
#' Patients in absorbing states at baseline are handled by advancing their first
#' absorbing event to day 1.
#'
#' @examples
#' \dontrun{
#' # Basic example with control parameters
#' baseline <- data.frame(
#'   id = 1:100,
#'   tx = rbinom(100, 1, 0.5),
#'   state = sample(2:5, 100, replace = TRUE)
#' )
#'
#' trajectories <- sim_trajectories_tte(
#'   baseline_data = baseline,
#'   param = c(0.05, 0.0035, 0.0025, 0.002, 0.002, 0.005),
#'   hazard_ratios = list(c(0.8, 1, 1, 1, 1, 1)),
#'   follow_up_time = 30,
#'   seed = 12345
#' )
#'
#' # Generate baseline_data within the function
#' test_traj <- sim_trajectories_tte(
#'   baseline_data = NULL,
#' states = 1:6,
#'   baseline_states = c(2:5),
#'   prob = c(0.55, 0.2, 0.15, 0.1),
#'   n = 10000,
#'   absorbing_states = 6,
#'   follow_up_time = 60,
#'   param = c(0.05, 0.003, 0.001, 0.001, 0.001, 0.0015),
#'   hazard_ratios = list(c(1.145, 1, 1, 1, 1, 1)),
#'   b = 0,
#'   seed = NULL
#' )
#' plot_sops(test_traj)
#' calc_time_in_state_diff(test_traj)
#'
#' # With multiple treatment arms
#' trajectories_multi <- sim_trajectories_tte(
#'   baseline_data = baseline,
#'   param = c(0.05, 0.0035, 0.0025, 0.002, 0.002, 0.005),
#'   hazard_ratios = list(
#'     c(0.8, 1, 1, 1, 1, 1),  # Treatment arm 1
#'     c(0.9, 1, 1, 1, 1, 1)   # Treatment arm 2
#'   ),
#'   follow_up_time = 30,
#'   seed = 12345
#' )
#' }
#'
#' @seealso
#' \code{\link{recurr_event}} for the underlying recurrent event generation
#'
#' @importFrom tidyr tibble
#' @importFrom stats rnorm
#'
#' @export

sim_trajectories_tte <- function(
  baseline_data,
  baseline_states = NULL,
  prob = c(0.55, 0.2, 0.15, 0.1),
  n = 1000,
  states = 1:6,
  absorbing_states = 6,
  follow_up_time = 60,
  param = c(0.05, 0.0035, 0.0025, 0.002, 0.002, 0.005),
  hazard_ratios = list(c(1, 1, 1, 1, 1, 1)),
  b = 0,
  seed = NULL
) {
  # Input validation
  if (!is.null(baseline_data) && !is.data.frame(baseline_data)) {
    stop("baseline_data must be a data frame or NULL")
  }

  if (is.null(baseline_data)) {
    baseline_data <- data.frame(
      id = 1:n,
      event_time = 0,
      tx = sample(c(0, 1), n, replace = TRUE),
      state = sample(
        if (!is.null(baseline_states)) baseline_states else states,
        n,
        replace = TRUE,
        prob = prob
      )
    )
  }

  required_cols <- c("id", "tx", "state", "event_time")
  missing_cols <- setdiff(required_cols, names(baseline_data))
  if (length(missing_cols) > 0) {
    stop(
      "baseline_data must contain columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.numeric(param) || length(param) != length(states)) {
    stop("param length must equal length(states)")
  }

  if (!is.list(hazard_ratios)) {
    stop("hazard_ratios must be a list of numeric vectors")
  }

  if (
    !all(vapply(
      hazard_ratios,
      function(x) {
        is.numeric(x) && length(x) == length(states)
      },
      logical(1)
    ))
  ) {
    stop(
      "Each element of hazard_ratios must have length equal to length(states)"
    )
  }

  if (any(baseline_data$state %in% absorbing_states)) {
    warning(
      "Some individuals are in absorbing state at baseline. ",
      "Consider changing the baseline distribution."
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Expand param to include control arm
  param_list <- list(ctrl = param)
  for (j in seq_along(hazard_ratios)) {
    param_list[[j + 1]] <- param * hazard_ratios[[j]]
  }

  # Generate recurrent events for each treatment arm × state combination
  state_changes <- list()
  m <- 1
  tx_levels <- sort(unique(baseline_data$tx))

  for (j in seq_along(tx_levels)) {
    tx_val <- tx_levels[j]
    ids_tx <- baseline_data$id[baseline_data$tx == tx_val]

    for (i in states) {
      state_changes[[m]] <- recurr_event(
        id = ids_tx,
        param = param_list[[j]][i],
        b = b,
        follow_up = follow_up_time,
        max_events = NULL
      )
      state_changes[[m]] <- cbind(state_changes[[m]], state = i, tx = tx_val)
      m <- m + 1
    }
  }

  # Combine event times, add baseline, and order
  state_changes_long <- rbind(
    do.call(rbind, state_changes),
    baseline_data
  )

  state_changes_long <- state_changes_long[
    order(state_changes_long[, "id"], state_changes_long[, "event_time"]),
  ]
  state_changes_long[, "time"] <- floor(state_changes_long[, "event_time"])

  # Split by patient ID
  split_state_changes_long <- split(
    as.data.frame(state_changes_long),
    state_changes_long[, "id"]
  )

  # LOCF helper function
  locf <- function(v) {
    ind <- which(!is.na(v))
    if (length(ind) == 0L) {
      return(v)
    }
    fi <- findInterval(seq_along(v), ind)
    res <- v[ind][pmax(fi, 1L)]
    res[fi == 0L] <- NA
    res
  }

  # Process each patient: remove post-absorbing events, handle duplicates, expand days, and LOCF
  state_changes_list <- lapply(split_state_changes_long, function(x) {
    this_id <- unique(x$id[!is.na(x$id)])
    if (length(this_id) != 1L) {
      this_id <- this_id[1L]
    }

    this_tx <- unique(x$tx[!is.na(x$tx)])
    if (length(this_tx) != 1L) {
      this_tx <- this_tx[1L]
    }

    # Absorbing state at baseline: set to day 1 (observed next day)
    if (any(x$time == 0 & x$state %in% absorbing_states)) {
      x$time[x$time == 0 & x$state %in% absorbing_states] <- 1
    }

    # Remove events after absorbing state
    absorbing_selector <- ave(
      x[, "state"] %in% absorbing_states,
      x[, "id"],
      FUN = cumsum
    )
    absorbing_selector <- c(0, absorbing_selector[-length(absorbing_selector)])
    x <- x[absorbing_selector == 0, , drop = FALSE]

    # Keep last event of each day
    x <- x[!duplicated(x$time, fromLast = TRUE), , drop = FALSE]

    # Expand to all days (0:follow_up_time)
    complete_times <- 0:follow_up_time
    x <- merge(
      data.frame(time = complete_times),
      x,
      by = "time",
      all.x = TRUE
    )

    # Restore id and tx
    x$id <- this_id
    x$tx <- this_tx

    # Convert state to numeric and apply LOCF
    if (is.factor(x$state)) {
      x$state <- as.integer(as.character(x$state))
    }
    x$y <- locf(x$state)

    # Generate yprev
    x$yprev <- c(NA, x$y[-nrow(x)])

    # Drop day 0 (keep only days 1 to follow_up_time)
    x <- x[x$time > 0, , drop = FALSE]

    # Select relevant columns
    x[, c("id", "tx", "time", "y", "yprev"), drop = FALSE]
  })

  state_changes_long <- do.call(rbind, state_changes_list)
  rownames(state_changes_long) <- NULL

  tibble::tibble(state_changes_long)
}
