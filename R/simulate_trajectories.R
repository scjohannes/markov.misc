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

  calc_piecewise_factor <- function(day) {
    if (is.null(drift_change_times) || day <= drift_change_times[1]) {
      return(1)
    }

    if (day <= drift_change_times[2]) {
      return(
        (drift_change_times[2] - day) /
          (drift_change_times[2] - drift_change_times[1])
      )
    }

    0
  }

  calc_shared_drift <- function(day, tx) {
    mu_base <- mu_drift + mu_quad * day + mu_quad2 * day^2
    (mu_base + tx * mu_treatment_effect) * calc_piecewise_factor(day)
  }

  threshold_matrix <- matrix(
    thresholds,
    nrow = follow_up_time + 1,
    ncol = n_states - 1,
    byrow = TRUE
  )

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
    pcat <- diff(c(0, cdf_fun(threshold_matrix[1, ] - X[i, 1]), 1))
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
        # Note: t here is 1-indexed (t = 2 means day 1), so use (t - 1).
        day <- t - 1
        mu_t <- calc_shared_drift(day, treatment[i])

        # Update latent severity (random walk with drift)
        X[i, t] <- rnorm(1, mean = X[i, t - 1] + mu_t, sd = sigma_rw)

        # Generate observation using selected CDF
        pcat <- diff(c(0, cdf_fun(threshold_matrix[t, ] - X[i, t]), 1))
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


#' Simulate Patient Trajectories Using Brownian Motion with Refresh Gaps
#'
#' Generates individual patient trajectories based on a latent continuous severity
#' variable that evolves daily as a Gaussian random walk with drift, but only
#' refreshes non-absorbing observed states on scheduled refresh days. This
#' reduces excessive observed switching while preserving daily latent evolution.
#'
#' @param n_patients Integer. Number of patients to simulate (default: 1000).
#' @param follow_up_time Integer. Number of days to simulate per patient (default: 60).
#' @param n_states Integer. Number of ordered states/categories (default: 6).
#' @param mu_drift Numeric. Baseline (constant) drift per day for control group
#'   (default: -0.26). Negative values indicate improvement.
#' @param mu_drift_sd Numeric. Standard deviation of a patient-specific random
#'   drift component added to `mu_drift` (default: 0). This induces persistent
#'   heterogeneity in individual improvement or worsening rates over follow-up.
#' @param mu_quad Numeric. Quadratic coefficient for time-varying drift
#'   (default: 0.006).
#' @param mu_quad2 Numeric. Second-order quadratic coefficient for drift (default: 0).
#' @param mu_treatment_effect Numeric. Additional constant drift for treatment group
#'   (default: -0.04). Negative values indicate faster improvement with treatment.
#' @param sigma_rw Numeric. Random walk innovation standard deviation per day
#'   (default: 0.12).
#' @param x0_sd Numeric. Standard deviation for baseline latent severity
#'   (default: 1.0).
#' @param thresholds Numeric vector of K-1 cutpoints for the ordered model
#'   (default: c(-0.9, 0.7, 1.95, 3, 4.25) for 6 states). Higher values
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
#' @param latent_dist Character. The distribution used to link the latent variable
#'   to observed states. Options are "logistic" (default) or "normal".
#' @param refresh_rate Numeric scalar. Exponential rate per day for the next
#'   non-absorbing state refresh (default: 0.048). Lower values create longer
#'   dwell times between observed state updates.
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#'
#' @return A data frame (tibble) with columns:
#'   - id: patient identifier
#'   - time: time/day (1 to follow_up_time)
#'   - tx: treatment assignment (0 = control, 1 = treatment)
#'   - y: observed ordinal state at time
#'   - x: latent continuous severity
#'
#' @details
#' This function modifies the daily Brownian simulator by separating latent
#' evolution from observed-state refresh:
#' 1. Each patient has an unobserved continuous "severity" X(t) that evolves
#'    daily as a Gaussian random walk with drift.
#' 2. A patient-specific waiting time to the next refresh is drawn from an
#'    exponential distribution with rate `refresh_rate`.
#' 3. Between refresh dates, non-absorbing observed states are carried forward.
#' 4. On refresh dates, the observed state is resampled from the latent-variable
#'    category probabilities excluding the absorbing state.
#' 5. Entry into the absorbing state can occur on any day according to the
#'    daily latent death probability; once entered, it is carried forward.
#'
#' @examples
#' \dontrun{
#' traj_gap <- sim_trajectories_brownian_gap(
#'   n_patients = 1000,
#'   seed = 12345
#' )
#'
#' traj_gap_low_switch <- sim_trajectories_brownian_gap(
#'   n_patients = 1000,
#'   refresh_rate = 0.02,
#'   seed = 12345
#' )
#' }
#'
#' @importFrom stats rnorm plogis pnorm rexp runif
#' @importFrom tibble tibble
#' @importFrom dplyr mutate left_join select arrange group_by ungroup filter
#' @importFrom tidyr pivot_longer
#'
#' @export
sim_trajectories_brownian_gap <- function(
  n_patients = 1000,
  follow_up_time = 60,
  n_states = 6,
  mu_drift = -0.26,
  mu_drift_sd = 0,
  mu_quad = 0.006,
  mu_quad2 = 0,
  mu_treatment_effect = -0.04,
  sigma_rw = 0.12,
  x0_sd = 1.0,
  thresholds = c(-0.9, 0.7, 1.95, 3, 4.25),
  treatment_prob = 0.5,
  allowed_start_state = 2:5,
  absorbing_state = 6,
  drift_change_times = c(25, 50),
  latent_dist = c("logistic", "normal"),
  refresh_rate = 0.048,
  seed = NULL
) {
  if (n_patients < 1 || !is.numeric(n_patients)) {
    stop("n_patients must be a positive integer")
  }

  if (follow_up_time < 1 || !is.numeric(follow_up_time)) {
    stop("follow_up_time must be a positive integer")
  }

  if (n_states < 2 || !is.numeric(n_states)) {
    stop("n_states must be at least 2")
  }

  if (!is.numeric(mu_drift_sd) || length(mu_drift_sd) != 1 || mu_drift_sd < 0) {
    stop("mu_drift_sd must be a non-negative numeric scalar")
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

  if (!is.numeric(refresh_rate) || length(refresh_rate) != 1 || refresh_rate <= 0) {
    stop("refresh_rate must be a positive numeric scalar")
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

  latent_dist <- match.arg(latent_dist)

  cdf_fun <- if (latent_dist == "logistic") {
    stats::plogis
  } else {
    stats::pnorm
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  treatment <- sample(
    c(0, 1),
    n_patients,
    replace = TRUE,
    prob = c(1 - treatment_prob, treatment_prob)
  )
  drift_random_effect <- stats::rnorm(n_patients, mean = 0, sd = mu_drift_sd)

  X <- matrix(NA_real_, n_patients, follow_up_time + 1)
  Y <- matrix(NA_integer_, n_patients, follow_up_time + 1)

  draw_wait <- function() {
    max(1L, ceiling(stats::rexp(1, rate = refresh_rate)))
  }

  calc_piecewise_factor <- function(day) {
    if (is.null(drift_change_times) || day <= drift_change_times[1]) {
      return(1)
    }

    if (day <= drift_change_times[2]) {
      return(
        (drift_change_times[2] - day) /
          (drift_change_times[2] - drift_change_times[1])
      )
    }

    0
  }

  calc_drift <- function(day, tx, drift_i) {
    mu_base <- mu_drift + drift_i + mu_quad * day + mu_quad2 * day^2
    mu_i <- mu_base + tx * mu_treatment_effect

    mu_i * calc_piecewise_factor(day)
  }

  threshold_matrix <- matrix(
    thresholds,
    nrow = follow_up_time + 1,
    ncol = n_states - 1,
    byrow = TRUE
  )

  for (i in seq_len(n_patients)) {
    X[i, 1] <- stats::rnorm(1, 0, x0_sd)

    # sample initial state occupancy probabilities
    pcat0 <- diff(c(0, cdf_fun(threshold_matrix[1, ] - X[i, 1]), 1))

    # set probabilities to 0 for states that are not allowed at baseline
    if (!is.null(allowed_start_state)) {
      pcat0[-allowed_start_state] <- 0
    }

    # sample initial state (baseline observation)
    Y[i, 1] <- sample.int(n_states, 1, prob = pcat0)

    # Initialize next refresh time for this patient
    next_refresh <- 1L + draw_wait()

    # still calculate latent X daily (could maybe be optimized to only calculate on refresh days)
    for (t in 2:(follow_up_time + 1)) {
      prev_y <- Y[i, t - 1]

      if (!is.null(absorbing_state) && prev_y == absorbing_state) {
        Y[i, t] <- absorbing_state
        X[i, t] <- NA_real_
        next
      }

      day <- t - 1
      mu_t <- calc_drift(day, treatment[i], drift_random_effect[i])
      X[i, t] <- stats::rnorm(1, mean = X[i, t - 1] + mu_t, sd = sigma_rw)

      pcat <- diff(c(0, cdf_fun(threshold_matrix[t, ] - X[i, t]), 1))

      # patients can enter absorbing state on any day, not just refresh days
      died_today <- !is.null(absorbing_state) &&
        stats::runif(1) < pcat[absorbing_state]

      if (died_today) {
        Y[i, t] <- absorbing_state
      } else if (t >= next_refresh) {
        # if did not enter absorbing state and a state transition is allowed, sample from latent variable probabilities
        if (!is.null(absorbing_state)) {
          pcat[absorbing_state] <- 0
        }

        if (sum(pcat) <= 0) {
          Y[i, t] <- prev_y
        } else {
          Y[i, t] <- sample.int(n_states, 1, prob = pcat)
        }

        # draw next refresh time
        next_refresh <- t + draw_wait()
      } else {
        # if before next time that allows transition, carry previous state forward
        Y[i, t] <- prev_y
      }
    }
  }

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

  result <- dat_observed |>
    dplyr::left_join(
      dat_latent |> dplyr::select(id, time, x),
      by = c("id", "time")
    ) |>
    dplyr::arrange(id, time)

  return(tibble::tibble(result))
}


#' Simulate ACTT-2-Inspired Ordinal Trajectories
#'
#' Generates synthetic ACTT-2-like ordinal outcome data using the Brownian gap
#' simulator with defaults tuned to roughly match the control-arm state
#' occupancy, cumulative death, and low rehospitalization patterns described in
#' the ACTT-2 trial. The defaults are null-aligned, so there is no built-in
#' treatment difference unless `mu_treatment_effect` is changed.
#'
#' This is intended as a convenient starting point for simulation studies, not
#' as an exact recreation of the original trial.
#'
#' @param n_patients Integer. Number of patients to simulate (default: 1000).
#' @param treatment_prob Numeric. Probability of assignment to treatment group
#'   (default: 0.5).
#' @param mu_treatment_effect Numeric. Additional constant drift for treatment
#'   group (default: 0, meaning no arm difference).
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#' @param ... Additional arguments passed to
#'   [sim_trajectories_brownian_gap()]. These can be used to override any of the
#'   ACTT-2-inspired defaults such as thresholds, drift terms, or refresh rate.
#'
#' @return A tibble with the same columns returned by
#'   [sim_trajectories_brownian_gap()].
#'
#' @examples
#' \dontrun{
#' actt2_null <- sim_actt2_brownian(n_patients = 1000, seed = 123)
#'
#' actt2_alt <- sim_actt2_brownian(
#'   n_patients = 1000,
#'   mu_treatment_effect = -0.03,
#'   seed = 123
#' )
#' }
#'
#' @export
sim_actt2_brownian <- function(
  n_patients = 1000,
  treatment_prob = 0.5,
  mu_treatment_effect = 0,
  seed = NULL,
  ...
) {
  args <- modifyList(
    list(
      n_patients = n_patients,
      follow_up_time = 28,
      n_states = 8,
      mu_drift = -2.65,
      mu_drift_sd = 1, #0.28,
      mu_quad = 0.1,
      mu_quad2 = -0.001,
      mu_treatment_effect = mu_treatment_effect,
      sigma_rw = 0.055,
      x0_sd = 0.35,
      thresholds = c(-6.42, -4.88, -3.98, -2.16, 0.5, 1.16, 1.8),
      treatment_prob = treatment_prob,
      allowed_start_state = 4:7,
      absorbing_state = 8L,
      drift_change_times = NULL,
      latent_dist = "logistic",
      refresh_rate = 0.067,
      seed = seed
    ),
    list(...)
  )

  do.call(sim_trajectories_brownian_gap, args)
}


#' Simulate Patient Trajectories Using "Line of Destiny" Model
#'
#' Generates deterministic patient trajectories where each patient follows their own
#' random line: θ₀ᵢ + θ₁ᵢ × log(d), where d is the day since randomization.
#' The observed ordinal score for day d is given as floor of θ₀ᵢ + θ₁ᵢ × log(d).
#' Death (highest state) and recovery (state 1) are absorbing states.
#'
#' This implements the data generating mechanism from Dodd et al. (2020) for
#' COVID-19 ordinal outcome power analysis, where subjects deterministically slide
#' up or down their own "line of destiny" and report the integer value each day.
#'
#' @param n_patients Integer. Number of patients (default: 1000).
#' @param follow_up_time Integer. Days to simulate (default: 28).
#' @param n_states Integer. Number of ordinal states (default: 6).
#'   State 1 = recovery (absorbing), State 6 = death (absorbing).
#' @param B0 Numeric. Fixed intercept for the model (default: 4).
#'   This should be set so that baseline states are distributed across
#'   the middle range (e.g., states 2-5 for hospitalized patients).
#' @param B1 Numeric. Fixed slope for log(time) (default: -0.4.
#'   Negative values indicate improvement over log-time.
#' @param B2 Numeric. Fixed treatment effect on slope (default: 0).
#'   The treatment group's slope is B1 + B2×Z where Z=1 for treatment.
#' @param b0_sd Numeric. SD of random intercept b₀ᵢ (default: 1).
#' @param mortality_prob_control Numeric. Probability of death trajectory for control
#'   group (default: 0.08). These patients get large positive b₁ᵢ values.
#' @param mortality_prob_treatment Numeric. Probability of death trajectory for treatment
#'   group (default: 0.08). Should be ≤ mortality_prob_control.
#' @param b1_mortality_mean Numeric. Mean of b₁ᵢ for patients destined to die (default: 1.5).
#'   Large positive values ensure progression to death.
#' @param b1_mortality_sd Numeric. SD of b₁ᵢ for mortality trajectories (default: 0.2).
#' @param b1_recovery_mean Numeric. Mean of b₁ᵢ for patients destined to recover (default: -0.35).
#'   Large negative values ensure fast recovery.
#' @param b1_recovery_sd Numeric. SD of b₁ᵢ for recovery trajectories (default: 0.12).
#' @param treatment_prob Numeric. Probability of treatment assignment (default: 0.5).
#' @param absorbing_states Integer vector. States that are absorbing (default: c(1, 6)).
#'   Typically c(1, n_states) for recovery and death. Note: patients cannot start
#'   in absorbing states at baseline (day 0); they are constrained to transient states.
#' @param seed Integer. Seed for reproducibility (default: NULL).
#'
#' @return A tibble with columns:
#'   - id: patient identifier
#'   - time: day since randomization (0 to follow_up_time, where 0 = baseline/randomization)
#'   - tx: treatment assignment (0 = control, 1 = treatment)
#'   - y: observed ordinal state
#'   - yprev: previous state (lagged y, non-missing for time >= 1)
#'   - x: continuous latent value (θ₀ᵢ + θ₁ᵢ × log(d) for d >= 1, θ₀ᵢ for d = 0)
#'   - theta0: patient-specific total intercept (B0 + b₀ᵢ)
#'   - theta1: patient-specific total slope (B1 + B2×Z + b₁ᵢ)
#'   - b0: random intercept b₀ᵢ
#'   - b1: random slope b₁ᵢ
#'
#' @details
#' The "line of destiny" model generates trajectories according to:
#'
#' **Y_id = B0 + B1×log(d) + B2×Z×log(d) + b₀ᵢ + b₁ᵢ×log(d)**
#'
#' where:
#' - B0, B1, B2 are fixed effects (population-level parameters)
#' - b₀ᵢ ~ N(0, b0_sd²) is the random intercept
#' - b₁ᵢ is a random slope with mixture distribution:
#'   * With probability p_death: b₁ᵢ ~ N(b1_mortality_mean, b1_mortality_sd²) → death
#'   * With probability (1 - p_death): b₁ᵢ ~ N(b1_recovery_mean, b1_recovery_sd²) → recovery
#' - Z is treatment indicator (0 or 1)
#' - d is day since randomization
#'
#' **Implementation steps:**
#'
#' 1. Assign treatment Z ~ Bernoulli(treatment_prob)
#'
#' 2. Sample random effects for each patient:
#'    - b₀ᵢ ~ N(0, b0_sd²)
#'    - Mortality indicator I ~ Bernoulli(p_death) where p_death depends on treatment
#'    - If I=1: b₁ᵢ ~ N(b1_mortality_mean, b1_mortality_sd²)
#'    - If I=0: b₁ᵢ ~ N(b1_recovery_mean, b1_recovery_sd²)
#'
#' 3. Calculate total slope for each patient:
#'    θ₁ᵢ = B1 + B2×Z + b₁ᵢ
#'
#' 4. Generate trajectories:
#'    - Day 0: xᵢ(0) = B0 + b₀ᵢ, yᵢ(0) = floor(xᵢ(0))
#'    - Days 1+: xᵢ(d) = B0 + b₀ᵢ + θ₁ᵢ×log(d), yᵢ(d) = floor(xᵢ(d))
#'    - Constrain y to 1, n_states
#'
#' 5. Once a patient reaches an absorbing state (1 or 7), they remain there.
#'
#' The log-time transformation log(d) creates realistic non-linear progression:
#' rapid initial change that slows over time. The floor() operation discretizes
#' the continuous trajectory into ordinal states.
#'
#' @examples
#' \dontrun{
#' # Default: 10% mortality in control, 5% in treatment
#' traj <- sim_trajectories_deterministic(
#'   n_patients = 1000,
#'   follow_up_time = 28,
#'   seed = 123
#' )
#'
#' # Higher mortality scenario (20% control, 10% treatment)
#' traj_high_mort <- sim_trajectories_deterministic(
#'   n_patients = 1000,
#'   mortality_prob_control = 0.20,
#'   mortality_prob_treatment = 0.10,
#'   seed = 456
#' )
#'
#' # Treatment effect on recovery slope (faster recovery)
#' traj_tx_effect <- sim_trajectories_deterministic(
#'   n_patients = 1000,
#'   theta1_treatment_effect = -0.05,  # More negative = faster improvement
#'   seed = 789
#' )
#'
#' # Visualize trajectories
#' library(ggplot2)
#' traj |>
#'   dplyr::filter(id <= 20) |>
#'   ggplot(aes(x = time, y = y, group = id, color = factor(tx))) +
#'   geom_line() +
#'   labs(title = "Lines of Destiny", color = "Treatment")
#' }
#'
#' @references
#' Dodd, Lori E., Dean Follmann, Jing Wang, Franz Koenig, Lisa L. Korn,
#' Christian Schoergenhofer, Michael Proschan, et al. 2020. “Endpoints for
#' Randomized Controlled Clinical Trials for COVID-19 Treatments.” Clinical
#' Trials (London, England) 17 (5): 472–82.

#'
#' @return A tibble with id, time, tx, y, yprev, x, theta0, theta1.
#'
#' @importFrom stats rbinom rnorm
#' @importFrom tibble tibble
#' @importFrom dplyr mutate left_join arrange group_by ungroup
#' @importFrom tidyr pivot_longer
#'
#' @export
sim_trajectories_deterministic <- function(
  n_patients = 1000,
  follow_up_time = 28,
  n_states = 6,
  B0 = 4,
  B1 = -0.40,
  B2 = 0,
  b0_sd = 1.0,
  mortality_prob_control = 0.08,
  mortality_prob_treatment = 0.08,
  b1_mortality_mean = 1.5,
  b1_mortality_sd = 0.2,
  b1_recovery_mean = -0.35,
  b1_recovery_sd = 0.12,
  treatment_prob = 0.5,
  absorbing_states = c(1, 7),
  seed = NULL
) {
  # --- Input Validation ---
  if (n_patients < 1 || !is.numeric(n_patients)) {
    stop("n_patients must be a positive integer")
  }

  if (follow_up_time < 1 || !is.numeric(follow_up_time)) {
    stop("follow_up_time must be a positive integer")
  }

  if (n_states < 2 || !is.numeric(n_states)) {
    stop("n_states must be at least 2")
  }

  if (mortality_prob_control < 0 || mortality_prob_control > 1) {
    stop("mortality_prob_control must be between 0 and 1")
  }

  if (mortality_prob_treatment < 0 || mortality_prob_treatment > 1) {
    stop("mortality_prob_treatment must be between 0 and 1")
  }

  if (treatment_prob < 0 || treatment_prob > 1) {
    stop("treatment_prob must be between 0 and 1")
  }

  if (!all(absorbing_states %in% 1:n_states)) {
    stop("All absorbing_states must be between 1 and n_states")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # --- Treatment Assignment ---
  treatment <- rbinom(n_patients, 1, treatment_prob)

  # --- Generate Patient-Specific Random Effects ---
  # Random intercepts: b₀ᵢ ~ N(0, b0_sd²)
  b0 <- rnorm(n_patients, mean = 0, sd = b0_sd)

  # Random slopes: b₁ᵢ from mixture distribution
  b1 <- numeric(n_patients)

  for (i in 1:n_patients) {
    # Determine mortality probability based on treatment
    mort_prob <- if (treatment[i] == 1) {
      mortality_prob_treatment
    } else {
      mortality_prob_control
    }

    # Decide if patient is destined to die or recover
    is_mortality_trajectory <- runif(1) < mort_prob

    if (is_mortality_trajectory) {
      # Patient destined to die: b₁ᵢ ~ N(b1_mortality_mean, b1_mortality_sd²)
      b1[i] <- rnorm(1, b1_mortality_mean, b1_mortality_sd)
    } else {
      # Patient destined to recover: b₁ᵢ ~ N(b1_recovery_mean, b1_recovery_sd²)
      b1[i] <- rnorm(1, b1_recovery_mean, b1_recovery_sd)
    }
  }

  # --- Compute Total Intercepts and Slopes ---
  # Total intercept: θ₀ᵢ = B0 + b₀ᵢ
  theta0 <- B0 + b0

  # Total slope: θ₁ᵢ = B1 + B2×Z + b₁ᵢ
  theta1 <- B1 + B2 * treatment + b1 # --- Simulate Trajectories ---
  # Initialize matrices for latent X and observed Y
  # Need follow_up_time + 1 columns to include day 0 (baseline/randomization)
  X <- matrix(NA_real_, n_patients, follow_up_time + 1)
  Y <- matrix(NA_integer_, n_patients, follow_up_time + 1)

  for (i in 1:n_patients) {
    absorbed <- FALSE
    absorbed_state <- NA_integer_

    # Day 0 (baseline/randomization): use only intercept θ₀ᵢ
    X[i, 1] <- theta0[i]
    y_raw <- floor(X[i, 1])
    y_constrained <- max(1, min(n_states, y_raw))

    # Ensure patient doesn't start in absorbing state at baseline
    # If they would start in an absorbing state, move them to nearest transient state
    if (y_constrained %in% absorbing_states) {
      if (y_constrained == 1) {
        y_constrained <- 2 # Move up from home to hospital
      } else if (y_constrained == n_states) {
        y_constrained <- n_states - 1 # Move down from death to ICU
      }
    }

    Y[i, 1] <- y_constrained # Days 1 to follow_up_time: apply log-time model
    for (d in 1:follow_up_time) {
      if (absorbed) {
        # Remain in absorbing state
        Y[i, d + 1] <- absorbed_state
        X[i, d + 1] <- NA_real_
      } else {
        # Calculate latent value: θ₀ᵢ + θ₁ᵢ × log(d)
        X[i, d + 1] <- theta0[i] + theta1[i] * log(d)

        # Observed state: floor of latent value
        y_raw <- floor(X[i, d + 1])

        # Constrain to valid state range [1, n_states]
        y_constrained <- max(1, min(n_states, y_raw))

        Y[i, d + 1] <- y_constrained

        # Check if entered absorbing state
        if (y_constrained %in% absorbing_states) {
          absorbed <- TRUE
          absorbed_state <- y_constrained
        }
      }
    }
  }

  # --- Format Output as Long Data Frame ---
  # Convert Y matrix to data frame
  dat_y <- as.data.frame(Y)
  colnames(dat_y) <- as.character(0:follow_up_time)
  dat_y <- dat_y |>
    dplyr::mutate(
      id = 1:n_patients,
      tx = treatment,
      theta0 = theta0,
      theta1 = theta1,
      b0 = b0,
      b1 = b1
    ) |>
    tidyr::pivot_longer(
      cols = -c(id, tx, theta0, theta1, b0, b1),
      names_to = "time",
      values_to = "y"
    ) |>
    dplyr::mutate(time = as.integer(time))

  # Convert X matrix to data frame
  dat_x <- as.data.frame(X)
  colnames(dat_x) <- as.character(0:follow_up_time)
  dat_x <- dat_x |>
    dplyr::mutate(id = 1:n_patients) |>
    tidyr::pivot_longer(
      cols = -id,
      names_to = "time",
      values_to = "x"
    ) |>
    dplyr::mutate(time = as.integer(time))

  # Combine and add lagged y (yprev)
  result <- dplyr::left_join(dat_y, dat_x, by = c("id", "time")) |>
    dplyr::arrange(id, time) |>
    dplyr::group_by(id) |>
    dplyr::mutate(yprev = dplyr::lag(y, 1)) |>
    dplyr::ungroup() |>
    dplyr::filter_out(time == 0)

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
#'   (Poisson process acceleration); b = 0 means independent events. If a vector is given, 
#'   it must be the same length as the states vector (event-type-specific acceleration). 
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

  # Expand b vector for each event type 
  if (length(b) == 1) {
    # Fixed acceleration
    b <- rep(b, times = length(states))
    } else {
    # Event-type-specific acceleration: check input
    if (length(b) != length(states)) {
      stop(
      "Autoregressive parameter b must have length one or equal to length(states)"
    )
      }
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
        b = b[i],
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
