test_that("sim_trajectories_brownian preserves default behavior", {
  traj_default <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 12,
    seed = 123
  )

  expect_trajectory_contract(
    traj_default,
    expected_cols = c("id", "tx", "time", "y", "yprev", "x"),
    n_patients = 40,
    follow_up_time = 12,
    states = 1:6
  )
})

test_that("sim_trajectories_brownian keeps the absorbing state", {
  traj <- sim_trajectories_brownian(
    n_patients = 4,
    follow_up_time = 5,
    allowed_start_state = 6L,
    seed = 321,
    absorbing_state = 6
  )

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "tx", "time", "y", "yprev", "x"),
    n_patients = 4,
    follow_up_time = 5,
    states = 1:6
  )
  expect_absorbing_state_sticky(traj, absorbing_state = 6)
  expect_equal(unique(traj$y), 6L)
})

test_that("sim_trajectories_brownian_gap preserves default behavior", {
  traj_default <- sim_trajectories_brownian_gap(
    n_patients = 40,
    follow_up_time = 12,
    seed = 123
  )

  expect_trajectory_contract(
    traj_default,
    expected_cols = c("id", "tx", "time", "y", "yprev", "x"),
    n_patients = 40,
    follow_up_time = 12,
    states = 1:6
  )
})

test_that("sim_trajectories_brownian_gap supports patient-specific drift heterogeneity", {
  traj_zero <- sim_trajectories_brownian_gap(
    n_patients = 80,
    follow_up_time = 15,
    mu_drift_sd = 0,
    seed = 123
  )

  traj_random <- sim_trajectories_brownian_gap(
    n_patients = 80,
    follow_up_time = 15,
    mu_drift_sd = 0.08,
    seed = 123
  )

  expect_identical(names(traj_zero), names(traj_random))
  expect_identical(dim(traj_zero), dim(traj_random))
  expect_gt(sum(traj_zero$y != traj_random$y), 0)

  expect_error(
    sim_trajectories_brownian_gap(
      n_patients = 5,
      follow_up_time = 10,
      mu_drift_sd = -0.1,
      seed = 1
    ),
    "mu_drift_sd must be a non-negative numeric scalar",
    fixed = TRUE
  )
})

test_that("sim_trajectories_brownian_gap supports delayed drift starts", {
  traj_immediate <- sim_trajectories_brownian_gap(
    n_patients = 40,
    follow_up_time = 8,
    drift_start = 0,
    seed = 123
  )

  traj_delayed <- sim_trajectories_brownian_gap(
    n_patients = 40,
    follow_up_time = 8,
    drift_start = function(n) sample(0:4, n, replace = TRUE),
    seed = 123
  )

  expect_identical(names(traj_immediate), names(traj_delayed))
  expect_identical(dim(traj_immediate), dim(traj_delayed))
  expect_gt(sum(traj_immediate$y != traj_delayed$y), 0)
})

test_that("sim_trajectories_brownian_gap validates refresh rate", {
  traj <- sim_trajectories_brownian_gap(
    n_patients = 30,
    follow_up_time = 10,
    refresh_rate = 0.03,
    seed = 7
  )

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "tx", "time", "y", "yprev", "x"),
    n_patients = 30,
    follow_up_time = 10,
    states = 1:6
  )
  expect_error(
    sim_trajectories_brownian_gap(
      n_patients = 5,
      follow_up_time = 10,
      refresh_rate = c(0.1, 0.2),
      seed = 1
    ),
    "refresh_rate must be a positive numeric scalar",
    fixed = TRUE
  )
})

test_that("sim_trajectories_brownian_gap keeps the absorbing state", {
  traj <- sim_trajectories_brownian_gap(
    n_patients = 4,
    follow_up_time = 5,
    allowed_start_state = 6L,
    seed = 321,
    absorbing_state = 6
  )

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "tx", "time", "y", "yprev", "x"),
    n_patients = 4,
    follow_up_time = 5,
    states = 1:6
  )
  expect_absorbing_state_sticky(traj, absorbing_state = 6)
  expect_equal(unique(traj$y), 6L)
})

test_that("sim_actt2_brownian returns ACTT-2-shaped output", {
  traj <- sim_actt2_brownian(n_patients = 60, seed = 1987)

  expect_trajectory_contract(
    traj,
    expected_cols = c("id", "tx", "time", "y", "yprev", "x"),
    n_patients = 60,
    follow_up_time = 28,
    states = 1:8
  )
  expect_true(all(traj$yprev[traj$time == 1] %in% 4:7))
})

test_that("sim_actt2_brownian uses a null treatment effect by default", {
  traj_default <- sim_actt2_brownian(n_patients = 40, seed = 99)
  traj_null <- sim_actt2_brownian(
    n_patients = 40,
    mu_treatment_effect = 0,
    seed = 99
  )

  expect_identical(traj_default, traj_null)
})

test_that("sim_trajectories_markov validates inputs and absent covariates", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1))

  expect_error(
    sim_trajectories_markov(baseline_data = list()),
    "baseline_data must be a data frame",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_markov(baseline_data = data.frame(id = 1)),
    "baseline_data must contain columns: yprev",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_markov(baseline_data = baseline, intercepts = character()),
    "intercepts must be a numeric vector",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_markov(baseline_data = baseline, lp_function = "nope"),
    "lp_function must be a function",
    fixed = TRUE
  )

  expect_warning(
    result <- sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 2,
      intercepts = 0,
      lp_function = function(yprev, t, parameter, extra_params) 0,
      covariate_names = "missing",
      seed = 1
    ),
    "None of the specified covariate_names found"
  )

  expect_equal(nrow(result), 4)
  expect_equal(names(result), c("id", "time", "y", "yprev"))
})

test_that("sim_trajectories_markov carries forward once all patients are absorbing", {
  baseline <- data.frame(id = 1:3, yprev = c(2, 2, 2), tx = 0)

  result <- sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = 4,
    intercepts = -Inf,
    lp_function = function(yprev, t, tx, parameter, extra_params) 0,
    absorbing_states = 2,
    seed = 1,
    covariate_names = "tx"
  )

  expect_equal(unique(result$y), 2)
  expect_equal(nrow(result), 12)
})

test_that("sim_trajectories_markov passes tx even when not requested as a covariate", {
  baseline <- data.frame(id = 1:2, yprev = c(1, 1), tx = c(0, 1))

  expect_warning(
    result <- sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 2,
      intercepts = 0,
      lp_function = function(yprev, t, tx, parameter, extra_params) {
        ifelse(tx == 1, 10, -10)
      },
      covariate_names = character(0),
      seed = 1
    ),
    "None of the specified covariate_names found",
    fixed = TRUE
  )

  expect_equal(nrow(result), 4)
  expect_true(all(result$tx %in% c(0, 1)))
})

test_that("sim_trajectories_markov preserves character patient IDs", {
  baseline <- data.frame(id = c("b", "a"), yprev = c(1, 1), tx = c(0, 1))

  result <- sim_trajectories_markov(
    baseline_data = baseline,
    follow_up_time = 2,
    intercepts = 0,
    lp_function = function(yprev, t, tx, parameter, extra_params) 0,
    covariate_names = "tx",
    seed = 1
  )

  expect_identical(result$id, rep(c("b", "a"), each = 2))
  expect_equal(result$time, rep(1:2, 2))
  expect_equal(result$tx, rep(c(0, 1), each = 2))
})

test_that("sim_trajectories_brownian validates inputs and supports normal link options", {
  expect_error(
    sim_trajectories_brownian(n_patients = 0),
    "n_patients must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(follow_up_time = 0),
    "follow_up_time must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(n_states = 1),
    "n_states must be at least 2",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(n_states = 3, thresholds = 0),
    "thresholds must have length n_states - 1",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(thresholds = c(1, 1, 2, 3, 4)),
    "thresholds must be a strictly increasing numeric vector",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(treatment_prob = 2),
    "treatment_prob must be between 0 and 1",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(allowed_start_state = c(1, 2)),
    "allowed_start_state must be NULL or an integer",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(n_states = 3, allowed_start_state = 4L, thresholds = c(0, 1)),
    "All allowed_start_state must be between 1 and n_states",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(
      n_states = 3,
      absorbing_state = 4,
      thresholds = c(0, 1),
      allowed_start_state = NULL
    ),
    "absorbing_state must be between 1 and n_states",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_brownian(drift_change_times = c(2, 1)),
    "drift_change_times must be a vector of length 2",
    fixed = TRUE
  )

  result <- sim_trajectories_brownian(
    n_patients = 4,
    follow_up_time = 3,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = NULL,
    absorbing_state = NULL,
    drift_change_times = NULL,
    latent_dist = "normal",
    seed = 10
  )

  expect_equal(nrow(result), 12)
  expect_setequal(unique(result$y), unique(result$y[result$y %in% 1:3]))
})

test_that("sim_trajectories_brownian applies late piecewise drift", {
  result <- sim_trajectories_brownian(
    n_patients = 4,
    follow_up_time = 4,
    drift_change_times = c(1, 2),
    seed = 101
  )

  expect_equal(max(result$time), 4)
  expect_named(result, c("id", "tx", "time", "y", "yprev", "x"))
})

test_that("sim_trajectories_brownian_gap validates additional inputs", {
  expect_error(sim_trajectories_brownian_gap(n_patients = 0), "n_patients must be")
  expect_error(sim_trajectories_brownian_gap(follow_up_time = 0), "follow_up_time must be")
  expect_error(sim_trajectories_brownian_gap(n_states = 1), "n_states must be")
  expect_error(
    sim_trajectories_brownian_gap(n_states = 3, thresholds = 0),
    "thresholds must have length n_states - 1"
  )
  expect_error(
    sim_trajectories_brownian_gap(thresholds = c(1, 1, 2, 3, 4)),
    "thresholds must be a strictly increasing"
  )
  expect_error(
    sim_trajectories_brownian_gap(threshold_time_effect_factors = c(0, 1)),
    "threshold_time_effect_factors must be NULL"
  )
  expect_error(
    sim_trajectories_brownian_gap(drift_start = -1),
    "drift_start must be a non-negative numeric scalar"
  )
  expect_error(
    sim_trajectories_brownian_gap(n_patients = 3, drift_start = c(0, 1)),
    "drift_start must be a non-negative numeric scalar"
  )
  expect_error(
    sim_trajectories_brownian_gap(
      thresholds = c(-1, 0, 1, 2, 3),
      threshold_time_effect_factors = c(0, 10, 0, 0, 0)
    ),
    "threshold_time_effect_factors induce crossing thresholds"
  )
  expect_error(sim_trajectories_brownian_gap(treatment_prob = -1), "treatment_prob must")
  expect_error(
    sim_trajectories_brownian_gap(allowed_start_state = c(1, 2)),
    "allowed_start_state must"
  )
  expect_error(
    sim_trajectories_brownian_gap(n_states = 3, allowed_start_state = 4L, thresholds = c(0, 1)),
    "All allowed_start_state"
  )
  expect_error(
    sim_trajectories_brownian_gap(
      n_states = 3,
      absorbing_state = 4,
      thresholds = c(0, 1),
      allowed_start_state = NULL
    ),
    "absorbing_state must"
  )
  expect_error(
    sim_trajectories_brownian_gap(drift_change_times = c(2, 1)),
    "drift_change_times must"
  )

  result <- sim_trajectories_brownian_gap(
    n_patients = 4,
    follow_up_time = 3,
    n_states = 3,
    thresholds = c(-1, 1),
    allowed_start_state = NULL,
    absorbing_state = NULL,
    drift_change_times = NULL,
    latent_dist = "normal",
    refresh_rate = 10,
    seed = 11
  )

  expect_equal(nrow(result), 12)
  expect_true(all(result$y %in% 1:3))
})

test_that("sim_trajectories_brownian_gap handles piecewise drift and zero transition mass", {
  result <- sim_trajectories_brownian_gap(
    n_patients = 3,
    follow_up_time = 4,
    n_states = 2,
    thresholds = 10,
    allowed_start_state = 2L,
    absorbing_state = 1L,
    refresh_rate = 100,
    drift_change_times = c(1, 2),
    seed = 102
  )

  expect_equal(max(result$time), 4)
  expect_true(all(result$y %in% 1:2))
})

test_that("sim_trajectories_brownian_gap reaches both piecewise drift tails", {
  result <- sim_trajectories_brownian_gap(
    n_patients = 3,
    follow_up_time = 5,
    absorbing_state = NULL,
    drift_change_times = c(2, 4),
    refresh_rate = 100,
    seed = 108
  )

  expect_equal(max(result$time), 5)
  expect_named(result, c("id", "tx", "time", "y", "yprev", "x"))
})

test_that("sim_trajectories_deterministic validates inputs and returns trajectories", {
  expect_error(sim_trajectories_deterministic(n_patients = 0), "n_patients must")
  expect_error(sim_trajectories_deterministic(follow_up_time = 0), "follow_up_time must")
  expect_error(sim_trajectories_deterministic(n_states = 1), "n_states must")
  expect_error(
    sim_trajectories_deterministic(mortality_prob_control = -0.1),
    "mortality_prob_control must"
  )
  expect_error(
    sim_trajectories_deterministic(mortality_prob_treatment = 1.1),
    "mortality_prob_treatment must"
  )
  expect_error(sim_trajectories_deterministic(treatment_prob = 2), "treatment_prob must")
  expect_error(
    sim_trajectories_deterministic(n_states = 5, absorbing_states = 7),
    "All absorbing_states must"
  )

  result <- sim_trajectories_deterministic(
    n_patients = 5,
    n_states = 7,
    follow_up_time = 4,
    mortality_prob_control = 1,
    mortality_prob_treatment = 0,
    treatment_prob = 0,
    seed = 12
  )

  expect_named(result, c("id", "tx", "theta0", "theta1", "b0", "b1", "time", "y", "x", "yprev"))
  expect_equal(nrow(result), 20)
  expect_equal(min(result$time), 1)
  expect_true(any(result$y %in% c(1, 6)))
})

test_that("sim_trajectories_deterministic handles baseline and later absorbing states", {
  high_start <- sim_trajectories_deterministic(
    n_patients = 2,
    n_states = 6,
    follow_up_time = 3,
    B0 = 100,
    b0_sd = 0,
    B1 = 0,
    b1_recovery_mean = 0,
    b1_recovery_sd = 0,
    mortality_prob_control = 0,
    mortality_prob_treatment = 0,
    treatment_prob = 0,
    absorbing_states = c(1, 6),
    seed = 103
  )
  low_start <- sim_trajectories_deterministic(
    n_patients = 2,
    n_states = 6,
    follow_up_time = 3,
    B0 = -100,
    b0_sd = 0,
    B1 = 0,
    b1_recovery_mean = 0,
    b1_recovery_sd = 0,
    mortality_prob_control = 0,
    mortality_prob_treatment = 0,
    treatment_prob = 0,
    absorbing_states = c(1, 6),
    seed = 104
  )
  later_absorb <- sim_trajectories_deterministic(
    n_patients = 2,
    n_states = 6,
    follow_up_time = 3,
    B0 = 3,
    b0_sd = 0,
    B1 = 10,
    b1_recovery_mean = 0,
    b1_recovery_sd = 0,
    mortality_prob_control = 0,
    mortality_prob_treatment = 0,
    treatment_prob = 0,
    absorbing_states = c(1, 6),
    seed = 105
  )

  expect_true(all(high_start$yprev[high_start$time == 1] == 5))
  expect_true(all(low_start$yprev[low_start$time == 1] == 2))
  expect_true(any(later_absorb$y == 6))

  treated <- sim_trajectories_deterministic(
    n_patients = 2,
    n_states = 6,
    follow_up_time = 2,
    treatment_prob = 1,
    mortality_prob_treatment = 1,
    b0_sd = 0,
    b1_mortality_sd = 0,
    absorbing_states = c(1, 6),
    seed = 109
  )
  expect_true(all(treated$tx == 1))
})

test_that("recurr_event simulates events and validates vector parameter lengths", {
  withr::local_seed(1)
  result <- recurr_event(id = 1:2, param = 0.5, b = 0, follow_up = 10, max_events = 3)

  expect_equal(ncol(result), 2)
  expect_true(all(result[, "id"] %in% 1:2))
  expect_true(all(result[, "event_time"] < 10))

  expect_error(
    recurr_event(id = 1:2, param = c(0.1, 0.2), follow_up = 10, max_events = NULL),
    "Length of param must be one or equal to max_events",
    fixed = TRUE
  )
})

test_that("recurr_event auto-selects event counts for scalar and vector rates", {
  withr::local_seed(2)
  scalar <- recurr_event(param = 0.01, b = 0.001, follow_up = 1, max_events = NULL)
  vector <- recurr_event(
    id = 1,
    param = c(0.01, 0.02, 0.03),
    b = 0,
    follow_up = 1,
    max_events = NULL
  )

  expect_equal(ncol(scalar), 2)
  expect_equal(ncol(vector), 2)

  expect_warning(
    recurr_event(id = 1, param = 100, b = 0, follow_up = 100, max_events = NULL),
    "max_events = 20",
    fixed = TRUE
  )
})

test_that("sim_trajectories_tte validates inputs", {
  expect_error(
    sim_trajectories_tte(baseline_data = list()),
    "baseline_data must be a data frame or NULL",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_tte(baseline_data = data.frame(id = 1)),
    "baseline_data must contain columns",
    fixed = TRUE
  )
  baseline <- data.frame(id = 1, tx = 0, state = 2, event_time = 0)
  expect_error(
    sim_trajectories_tte(baseline, states = 1:3, param = c(0.1, 0.2)),
    "param length must equal length\\(states\\)"
  )
  expect_error(
    sim_trajectories_tte(baseline, hazard_ratios = c(1, 1, 1, 1, 1, 1)),
    "hazard_ratios must be a list",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_tte(baseline, states = 1:3, param = rep(0.1, 3), hazard_ratios = list(c(1, 1))),
    "Each element of hazard_ratios must have length",
    fixed = TRUE
  )
  expect_error(
    sim_trajectories_tte(
      baseline,
      states = 1:3,
      param = rep(0.1, 3),
      hazard_ratios = list(rep(1, 3)),
      b = c(0, 0)
    ),
    "Autoregressive parameter b must have length one",
    fixed = TRUE
  )
})

test_that("sim_trajectories_tte expands event histories into daily states", {
  baseline <- data.frame(
    id = c(1, 2),
    tx = c(0, 1),
    state = c(2, 2),
    event_time = 0
  )

  with_mocked_bindings(
    recurr_event = function(id, param, b, follow_up, max_events) {
      if (length(id) == 0) {
        return(cbind(id = numeric(0), event_time = numeric(0)))
      }
      cbind(id = id, event_time = seq_along(id))
    },
    {
      result <- sim_trajectories_tte(
        baseline,
        states = 1:3,
        absorbing_states = 3,
        follow_up_time = 3,
        param = c(0.1, 0.2, 0.3),
        hazard_ratios = list(c(1, 2, 3)),
        b = c(0, 0.01, 0.02),
        seed = 1
      )
    }
  )

  expect_named(result, c("id", "tx", "time", "y", "yprev"))
  expect_equal(nrow(result), 6)
  expect_equal(result$time, rep(1:3, 2))
  expect_true(all(result$y %in% 1:3))
})

test_that("sim_trajectories_tte maps positive tx codes to matching hazard ratios", {
  baseline <- data.frame(
    id = 1,
    tx = 1,
    state = 1,
    event_time = 0
  )
  observed_param <- numeric()

  with_mocked_bindings(
    recurr_event = function(id, param, b, follow_up, max_events) {
      observed_param <<- c(observed_param, param)
      data.frame(id = numeric(0), event_time = numeric(0))
    },
    {
      sim_trajectories_tte(
        baseline,
        states = 1:2,
        absorbing_states = 2,
        follow_up_time = 1,
        param = c(10, 20),
        hazard_ratios = list(c(0.5, 2)),
        b = 0,
        seed = 1
      )
    }
  )

  expect_equal(observed_param, c(5, 40))
})

test_that("sim_trajectories_tte requires hazard ratios for observed treatment arms", {
  baseline <- data.frame(
    id = 1,
    tx = 2,
    state = 1,
    event_time = 0
  )

  expect_error(
    sim_trajectories_tte(
      baseline,
      states = 1:2,
      param = c(10, 20),
      hazard_ratios = list(c(1, 1)),
      b = 0
    ),
    "Missing arm\\(s\\): 2"
  )
})

test_that("sim_trajectories_tte generates baseline data and handles absorbing baselines", {
  with_mocked_bindings(
    recurr_event = function(id, param, b, follow_up, max_events) {
      cbind(id = numeric(0), event_time = numeric(0))
    },
    {
      generated <- sim_trajectories_tte(
        baseline_data = NULL,
        baseline_states = 1:3,
        prob = c(0.2, 0.3, 0.5),
        n = 2,
        states = 1:3,
        absorbing_states = 3,
        follow_up_time = 2,
        param = c(0.01, 0.02, 0.03),
        hazard_ratios = list(rep(1, 3)),
        b = 0,
        seed = 106
      )
    }
  )

  expect_equal(nrow(generated), 4)
  expect_named(generated, c("id", "tx", "time", "y", "yprev"))

  baseline <- data.frame(
    id = c(1, 2),
    tx = c(0, 1),
    state = factor(c(NA, 3), levels = 1:3),
    event_time = 0
  )
  with_mocked_bindings(
    recurr_event = function(id, param, b, follow_up, max_events) {
      cbind(id = numeric(0), event_time = numeric(0))
    },
    {
      expect_warning(
        absorbing <- sim_trajectories_tte(
          baseline,
          states = 1:3,
          absorbing_states = 3,
          follow_up_time = 2,
          param = c(0.01, 0.02, 0.03),
          hazard_ratios = list(rep(1, 3)),
          b = c(0, 0, 0),
          seed = 107
        ),
        "absorbing state at baseline",
        fixed = TRUE
      )
    }
  )

  expect_equal(nrow(absorbing), 4)
  expect_true(any(is.na(absorbing$y)))

  numeric_absorbing <- data.frame(
    id = 1,
    tx = 0,
    state = 3,
    event_time = 0
  )
  with_mocked_bindings(
    recurr_event = function(id, param, b, follow_up, max_events) {
      data.frame(id = numeric(0), event_time = numeric(0))
    },
    {
      expect_warning(
        numeric_out <- sim_trajectories_tte(
          numeric_absorbing,
          states = 1:3,
          absorbing_states = 3,
          follow_up_time = 2,
          param = c(0.01, 0.02, 0.03),
          hazard_ratios = list(rep(1, 3)),
          b = 0,
          seed = 110
        ),
        "absorbing state at baseline",
        fixed = TRUE
      )
    }
  )
  expect_equal(numeric_out$y, c(3, 3))
})
