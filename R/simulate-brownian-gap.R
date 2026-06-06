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
#' @param mu_drift Numeric scalar or vector of length `n_states - 1`. Baseline
#'   threshold-specific drift per day for the control group (default: -0.26).
#'   Negative values indicate improvement. Vector element `j` affects threshold
#'   `j`, the boundary between states `j` and `j + 1`.
#' @param mu_drift_sd Numeric. Standard deviation of a patient-specific random
#'   drift component added to `mu_drift` (default: 0). This induces persistent
#'   heterogeneity in individual improvement or worsening rates over follow-up.
#' @param mu_quad Numeric scalar or vector of length `n_states - 1`. Quadratic
#'   coefficient for threshold-specific time-varying drift (default: 0.006).
#' @param mu_quad2 Numeric scalar or vector of length `n_states - 1`.
#'   Second-order quadratic coefficient for threshold-specific drift (default:
#'   0).
#' @param mu_treatment_effect Numeric scalar or vector of length `n_states - 1`.
#'   Additional threshold-specific drift for the treatment group (default:
#'   -0.04). Negative values indicate faster improvement with treatment. To
#'   target a single adjacent-state boundary, set only that vector element to a
#'   nonzero value.
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
#'   begins to decline and when it reaches zero after drift has started
#'   (default: c(25, 50)).
#' @param drift_start Numeric scalar, numeric vector of length `n_patients`, or
#'   function. Day when deterministic drift starts for each patient (default:
#'   0). A function is called as `drift_start(n_patients)` and must return a
#'   non-negative numeric scalar or vector of length `n_patients`.
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
#'   - x: reference latent continuous severity
#'
#' @details
#' This function modifies the daily Brownian simulator by separating latent
#' evolution from observed-state refresh:
#' 1. Each patient has an unobserved continuous "severity" X(t) that evolves
#'    daily as a Gaussian random walk with drift.
#' 2. When drift parameters are threshold-specific vectors, `x` follows a
#'    reference latent drift while observed-state probabilities are generated
#'    from threshold-specific cutpoint trajectories.
#' 3. A patient-specific waiting time to the next refresh is drawn from an
#'    exponential distribution with rate `refresh_rate`.
#' 4. Between refresh dates, non-absorbing observed states are carried forward.
#' 5. On refresh dates, the observed state is resampled from the latent-variable
#'    category probabilities excluding the absorbing state.
#' 6. Entry into the absorbing state can occur on any day according to the
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
#'
#' # Treatment effect only at the state 5/6 boundary
#' traj_gap_targeted <- sim_trajectories_brownian_gap(
#'   n_patients = 1000,
#'   mu_treatment_effect = c(0, 0, 0, 0, -0.04),
#'   seed = 12345
#' )
#' }
#'
#' @importFrom stats rnorm plogis pnorm rexp runif rpois
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
  drift_start = 0,
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

  normalize_threshold_effect <- function(x, name) {
    if (
      !is.numeric(x) ||
        !(length(x) %in% c(1, n_states - 1)) ||
        any(!is.finite(x))
    ) {
      stop(
        name,
        " must be a finite numeric scalar or vector with length n_states - 1"
      )
    }

    rep(x, length.out = n_states - 1)
  }

  mu_drift_threshold <- normalize_threshold_effect(mu_drift, "mu_drift")
  mu_quad_threshold <- normalize_threshold_effect(mu_quad, "mu_quad")
  mu_quad2_threshold <- normalize_threshold_effect(mu_quad2, "mu_quad2")
  mu_treatment_threshold <- normalize_threshold_effect(
    mu_treatment_effect,
    "mu_treatment_effect"
  )

  mu_drift_reference <- mean(mu_drift_threshold)
  mu_quad_reference <- mean(mu_quad_threshold)
  mu_quad2_reference <- mean(mu_quad2_threshold)
  mu_treatment_reference <- mean(mu_treatment_threshold)

  if (treatment_prob < 0 || treatment_prob > 1) {
    stop("treatment_prob must be between 0 and 1")
  }

  if (
    !is.numeric(refresh_rate) ||
      length(refresh_rate) != 1 ||
      any(!is.finite(refresh_rate)) ||
      any(refresh_rate <= 0)
  ) {
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

  if (is.function(drift_start)) {
    drift_start <- drift_start(n_patients)
  }

  if (
    !is.numeric(drift_start) ||
      !(length(drift_start) %in% c(1, n_patients)) ||
      any(!is.finite(drift_start)) ||
      any(drift_start < 0)
  ) {
    stop(
      "drift_start must be a non-negative numeric scalar, a numeric vector ",
      "of length n_patients, or a function returning one of those"
    )
  }

  drift_start <- rep(drift_start, length.out = n_patients)

  X <- matrix(NA_real_, n_patients, follow_up_time + 1)
  Y <- matrix(NA_integer_, n_patients, follow_up_time + 1)

  draw_wait <- function() {
    max(1L, ceiling(stats::rexp(1, rate = refresh_rate)))
  }

  calc_piecewise_factor <- function(elapsed_day) {
    if (elapsed_day < 0) {
      return(0)
    }

    if (is.null(drift_change_times) || elapsed_day <= drift_change_times[1]) {
      return(1)
    }

    if (elapsed_day <= drift_change_times[2]) {
      return(
        (drift_change_times[2] - elapsed_day) /
          (drift_change_times[2] - drift_change_times[1])
      )
    }

    0
  }

  calc_drift <- function(day, tx, drift_i, drift_start_i) {
    elapsed_day <- day - drift_start_i
    if (elapsed_day < 0) {
      return(0)
    }

    mu_base <- mu_drift_reference +
      drift_i +
      mu_quad_reference * elapsed_day +
      mu_quad2_reference * elapsed_day^2
    mu_i <- mu_base + tx * mu_treatment_reference

    mu_i * calc_piecewise_factor(elapsed_day)
  }

  threshold_matrix_cache <- new.env(parent = emptyenv())
  make_threshold_matrix <- function(drift_start_i, tx) {
    cache_key <- paste(drift_start_i, tx, sep = ":")
    if (exists(cache_key, envir = threshold_matrix_cache, inherits = FALSE)) {
      return(get(cache_key, envir = threshold_matrix_cache, inherits = FALSE))
    }

    reference_drift <- numeric(follow_up_time + 1)
    threshold_drift <- matrix(
      0,
      nrow = follow_up_time + 1,
      ncol = n_states - 1
    )
    for (t in 2:(follow_up_time + 1)) {
      day <- t - 1
      elapsed_day <- day - drift_start_i
      if (elapsed_day < 0) {
        reference_drift[t] <- reference_drift[t - 1]
        threshold_drift[t, ] <- threshold_drift[t - 1, ]
      } else {
        piecewise_factor <- calc_piecewise_factor(elapsed_day)
        reference_step <- (mu_drift_reference +
          mu_quad_reference * elapsed_day +
          mu_quad2_reference * elapsed_day^2 +
          tx * mu_treatment_reference) *
          piecewise_factor
        threshold_step <- (mu_drift_threshold +
          mu_quad_threshold * elapsed_day +
          mu_quad2_threshold * elapsed_day^2 +
          tx * mu_treatment_threshold) *
          piecewise_factor

        reference_drift[t] <- reference_drift[t - 1] + reference_step
        threshold_drift[t, ] <- threshold_drift[t - 1, ] + threshold_step
      }
    }

    threshold_matrix <- matrix(
      thresholds,
      nrow = follow_up_time + 1,
      ncol = n_states - 1,
      byrow = TRUE
    ) +
      reference_drift -
      threshold_drift

    if (!all(apply(threshold_matrix, 1, function(x) all(diff(x) > 0)))) {
      stop(
        "threshold-specific drift terms induce crossing thresholds; ",
        "use effects that keep thresholds strictly increasing over follow-up"
      )
    }

    assign(cache_key, threshold_matrix, envir = threshold_matrix_cache)
    threshold_matrix
  }

  for (i in seq_len(n_patients)) {
    threshold_matrix_i <- make_threshold_matrix(drift_start[i], treatment[i])
    X[i, 1] <- stats::rnorm(1, 0, x0_sd)

    # sample initial state occupancy probabilities
    pcat0 <- diff(c(0, cdf_fun(threshold_matrix_i[1, ] - X[i, 1]), 1))

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
      mu_t <- calc_drift(
        day,
        treatment[i],
        drift_random_effect[i],
        drift_start[i]
      )
      X[i, t] <- stats::rnorm(1, mean = X[i, t - 1] + mu_t, sd = sigma_rw)

      pcat <- diff(c(0, cdf_fun(threshold_matrix_i[t, ] - X[i, t]), 1))

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
#' simulator with defaults tuned to roughly preserve control-arm state
#' occupancy while reducing rehospitalization-like movement from recovered
#' states 1/2 into states 3:7. The defaults are null-aligned, so there is no
#' built-in treatment difference unless `mu_treatment_effect` is changed.
#' `mu_treatment_effect` may be a scalar, which affects every threshold equally,
#' or a vector with one value per threshold.
#'
#' This is intended as a convenient starting point for simulation studies, not
#' as an exact recreation of the original trial.
#'
#' @param n_patients Integer. Number of patients to simulate (default: 1000).
#' @param treatment_prob Numeric. Probability of assignment to treatment group
#'   (default: 0.5).
#' @param mu_treatment_effect Numeric scalar or vector of length 7. Additional
#'   threshold-specific drift for the treatment group (default: 0, meaning no
#'   arm difference). Vector element `j` affects the boundary between states `j`
#'   and `j + 1`.
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
#'
#' actt2_state5_to_6 <- sim_actt2_brownian(
#'   n_patients = 1000,
#'   mu_treatment_effect = c(0, 0, 0, 0, 0.03, 0, 0),
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
      mu_drift = c(-0.99, -1.1, -1.1, -1.1, -1.1, -1.078, -1.1),
      mu_drift_sd = 1, #0.28,
      mu_quad = c(0, 0, 0, 0, 0, 0, 0),
      mu_quad2 = c(
        -0.00045,
        -0.0005,
        -0.0005,
        -0.0005,
        -0.0005,
        -0.00049,
        -0.0005
      ),
      mu_treatment_effect = mu_treatment_effect,
      sigma_rw = 0.055,
      x0_sd = 0.35,
      thresholds = c(-7.0, -4.88, -3.98, -2.16, 0.3, 1.16, 9),
      treatment_prob = treatment_prob,
      allowed_start_state = 4:7,
      absorbing_state = 8L,
      drift_change_times = c(10, 28),
      drift_start = function(n) rpois(n, 3.2),
      latent_dist = "logistic",
      refresh_rate = 0.32,
      seed = seed
    ),
    list(...)
  )

  do.call(sim_trajectories_brownian_gap, args)
}
