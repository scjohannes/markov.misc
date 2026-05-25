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

  result <- data.frame(
    id = rep(seq_len(n_patients), each = follow_up_time),
    tx = rep(treatment, each = follow_up_time),
    time = rep(seq_len(follow_up_time), times = n_patients),
    y = as.integer(as.vector(t(Y[, -1, drop = FALSE]))),
    yprev = as.integer(as.vector(t(Y[, -ncol(Y), drop = FALSE]))),
    x = as.vector(t(X[, -1, drop = FALSE]))
  )

  return(result)
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
      thresholds = c(-14.42, -4.88, -3.98, -2.16, 0.5, 1.16, 1.8),
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
