#' Simulate Patient Trajectories Using "Line of Destiny" Model
#'
#' Generates deterministic patient trajectories where each patient follows their own
#' random line: theta0_i + theta1_i * log(d), where d is the day since
#' randomization.
#' The observed ordinal score for day d is given as
#' floor(theta0_i + theta1_i * log(d)).
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
#'   The treatment group's slope is B1 + B2 * Z where Z = 1 for treatment.
#' @param b0_sd Numeric. SD of random intercept b0_i (default: 1).
#' @param mortality_prob_control Numeric. Probability of death trajectory for control
#'   group (default: 0.08). These patients get large positive b1_i values.
#' @param mortality_prob_treatment Numeric. Probability of death trajectory for treatment
#'   group (default: 0.08). Should be <= mortality_prob_control.
#' @param b1_mortality_mean Numeric. Mean of b1_i for patients destined to die
#'   (default: 1.5).
#'   Large positive values ensure progression to death.
#' @param b1_mortality_sd Numeric. SD of b1_i for mortality trajectories
#'   (default: 0.2).
#' @param b1_recovery_mean Numeric. Mean of b1_i for patients destined to
#'   recover (default: -0.35).
#'   Large negative values ensure fast recovery.
#' @param b1_recovery_sd Numeric. SD of b1_i for recovery trajectories
#'   (default: 0.12).
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
#'   - x: continuous latent value (theta0_i + theta1_i * log(d) for d >= 1,
#'     theta0_i for d = 0)
#'   - theta0: patient-specific total intercept (B0 + b0_i)
#'   - theta1: patient-specific total slope (B1 + B2 * Z + b1_i)
#'   - b0: random intercept b0_i
#'   - b1: random slope b1_i
#'
#' @details
#' The "line of destiny" model generates trajectories according to:
#'
#' **Y_id = B0 + B1 * log(d) + B2 * Z * log(d) + b0_i + b1_i * log(d)**
#'
#' where:
#' - B0, B1, B2 are fixed effects (population-level parameters)
#' - b0_i ~ N(0, b0_sd^2) is the random intercept
#' - b1_i is a random slope with mixture distribution:
#'   * With probability p_death: b1_i ~
#'     N(b1_mortality_mean, b1_mortality_sd^2) -> death
#'   * With probability (1 - p_death): b1_i ~
#'     N(b1_recovery_mean, b1_recovery_sd^2) -> recovery
#' - Z is treatment indicator (0 or 1)
#' - d is day since randomization
#'
#' **Implementation steps:**
#'
#' 1. Assign treatment Z ~ Bernoulli(treatment_prob)
#'
#' 2. Sample random effects for each patient:
#'    - b0_i ~ N(0, b0_sd^2)
#'    - Mortality indicator I ~ Bernoulli(p_death) where p_death depends on treatment
#'    - If I = 1: b1_i ~ N(b1_mortality_mean, b1_mortality_sd^2)
#'    - If I = 0: b1_i ~ N(b1_recovery_mean, b1_recovery_sd^2)
#'
#' 3. Calculate total slope for each patient:
#'    theta1_i = B1 + B2 * Z + b1_i
#'
#' 4. Generate trajectories:
#'    - Day 0: x_i(0) = B0 + b0_i, y_i(0) = floor(x_i(0))
#'    - Days 1+: x_i(d) = B0 + b0_i + theta1_i * log(d),
#'      y_i(d) = floor(x_i(d))
#'    - Constrain y to 1, n_states
#'
#' 5. Once a patient reaches an absorbing state (typically 1 or `n_states`),
#'    they remain there.
#'
#' The log-time transformation log(d) creates realistic non-linear progression:
#' rapid initial change that slows over time. The floor() operation discretizes
#' the continuous trajectory into ordinal states.
#'
#' @examples
#' # Default mortality settings
#' traj <- sim_trajectories_deterministic(
#'   n_patients = 20,
#'   follow_up_time = 7,
#'   seed = 123
#' )
#'
#' # Higher mortality scenario
#' traj_high_mort <- sim_trajectories_deterministic(
#'   n_patients = 20,
#'   follow_up_time = 7,
#'   mortality_prob_control = 0.20,
#'   mortality_prob_treatment = 0.10,
#'   seed = 456
#' )
#'
#' # Treatment effect on recovery slope
#' traj_tx_effect <- sim_trajectories_deterministic(
#'   n_patients = 20,
#'   follow_up_time = 7,
#'   B2 = -0.05,
#'   seed = 789
#' )
#'
#' @references
#' Dodd, Lori E., Dean Follmann, Jing Wang, Franz Koenig, Lisa L. Korn,
#' Christian Schoergenhofer, Michael Proschan, et al. 2020. "Endpoints for
#' Randomized Controlled Clinical Trials for COVID-19 Treatments." Clinical
#' Trials (London, England) 17 (5): 472-82.
#'
#' @return A tibble with id, time, tx, y, yprev, x, theta0, theta1.
#'
#' @importFrom stats rbinom rnorm
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
  absorbing_states = c(1, 6),
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
  # Random intercepts: b0_i ~ N(0, b0_sd^2)
  b0 <- rnorm(n_patients, mean = 0, sd = b0_sd)

  # Random slopes: b1_i from mixture distribution
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
      # Patient destined to die: b1_i ~ N(b1_mortality_mean, b1_mortality_sd^2)
      b1[i] <- rnorm(1, b1_mortality_mean, b1_mortality_sd)
    } else {
      # Patient destined to recover: b1_i ~ N(b1_recovery_mean, b1_recovery_sd^2)
      b1[i] <- rnorm(1, b1_recovery_mean, b1_recovery_sd)
    }
  }

  # --- Compute Total Intercepts and Slopes ---
  # Total intercept: theta0_i = B0 + b0_i
  theta0 <- B0 + b0

  # Total slope: theta1_i = B1 + B2 * Z + b1_i
  theta1 <- B1 + B2 * treatment + b1 # --- Simulate Trajectories ---
  # Initialize matrices for latent X and observed Y
  # Need follow_up_time + 1 columns to include day 0 (baseline/randomization)
  X <- matrix(NA_real_, n_patients, follow_up_time + 1)
  Y <- matrix(NA_integer_, n_patients, follow_up_time + 1)

  for (i in 1:n_patients) {
    absorbed <- FALSE
    absorbed_state <- NA_integer_

    # Day 0 (baseline/randomization): use only intercept theta0_i
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
        # Calculate latent value: theta0_i + theta1_i * log(d)
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
  colnames(Y) <- as.character(0:follow_up_time)
  colnames(X) <- as.character(0:follow_up_time)

  dat_y <- matrix_to_long(Y, value_name = "y")
  dat_y$id <- as.integer(dat_y$id)
  dat_y$time <- as.integer(dat_y$time)
  dat_y$tx <- treatment[dat_y$id]
  dat_y$theta0 <- theta0[dat_y$id]
  dat_y$theta1 <- theta1[dat_y$id]
  dat_y$b0 <- b0[dat_y$id]
  dat_y$b1 <- b1[dat_y$id]
  dat_y <- reorder_columns(dat_y, c("id", "tx", "theta0", "theta1", "b0", "b1", "time", "y"))

  dat_x <- matrix_to_long(X, value_name = "x")
  dat_x$id <- as.integer(dat_x$id)
  dat_x$time <- as.integer(dat_x$time)

  result <- left_join_preserve_order(dat_y, dat_x, by = c("id", "time"))
  result <- result[order(result$id, result$time), , drop = FALSE]
  result$yprev <- ave(result$y, result$id, FUN = function(x) c(NA, utils::head(x, -1)))
  result <- result[result$time != 0, , drop = FALSE]
  rownames(result) <- NULL

  return(result)
}
