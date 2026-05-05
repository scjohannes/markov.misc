error_message <- function(expr) {
  tryCatch(
    {
      eval.parent(substitute(expr))
      NULL
    },
    error = \(e) conditionMessage(e)
  )
}

test_that("sim_trajectories_brownian preserves default behavior", {
  traj_default <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 12,
    seed = 123
  )

  expect_named(traj_default, c("id", "tx", "time", "y", "yprev", "x"))
  expect_equal(max(traj_default$time), 12)
})

test_that("sim_trajectories_brownian keeps the absorbing state", {
  traj <- sim_trajectories_brownian(
    n_patients = 80,
    follow_up_time = 20,
    seed = 321
  )

  deaths <- traj[traj$y == 6, ]
  if (nrow(deaths) > 0) {
    for (patient_id in unique(deaths$id)) {
      patient_data <- traj[traj$id == patient_id, ]
      first_death <- min(patient_data$time[patient_data$y == 6])
      subsequent_states <- patient_data$y[patient_data$time >= first_death]
      expect_setequal(unique(subsequent_states), 6)
    }
  }
})

test_that("sim_trajectories_brownian_gap preserves default behavior", {
  traj_default <- sim_trajectories_brownian_gap(
    n_patients = 40,
    follow_up_time = 12,
    seed = 123
  )

  expect_named(traj_default, c("id", "tx", "time", "y", "yprev", "x"))
  expect_equal(max(traj_default$time), 12)
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

  expect_match(
    error_message(
      sim_trajectories_brownian_gap(
        n_patients = 5,
        follow_up_time = 10,
        mu_drift_sd = -0.1,
        seed = 1
      )
    ),
    "mu_drift_sd must be a non-negative numeric scalar"
  )
})

test_that("sim_trajectories_brownian_gap validates refresh rate", {
  traj <- sim_trajectories_brownian_gap(
    n_patients = 30,
    follow_up_time = 10,
    refresh_rate = 0.03,
    seed = 7
  )

  expect_named(traj, c("id", "tx", "time", "y", "yprev", "x"))
  expect_equal(max(traj$time), 10)

  expect_match(
    error_message(
      sim_trajectories_brownian_gap(
        n_patients = 5,
        follow_up_time = 10,
        refresh_rate = c(0.1, 0.2),
        seed = 1
      )
    ),
    "refresh_rate must be a positive numeric scalar"
  )
})

test_that("sim_trajectories_brownian_gap keeps the absorbing state", {
  traj <- sim_trajectories_brownian_gap(
    n_patients = 80,
    follow_up_time = 20,
    seed = 321
  )

  deaths <- traj[traj$y == 6, ]
  if (nrow(deaths) > 0) {
    for (patient_id in unique(deaths$id)) {
      patient_data <- traj[traj$id == patient_id, ]
      first_death <- min(patient_data$time[patient_data$y == 6])
      subsequent_states <- patient_data$y[patient_data$time >= first_death]
      expect_setequal(unique(subsequent_states), 6)
    }
  }
})

test_that("sim_actt2_brownian returns ACTT-2-shaped output", {
  traj <- sim_actt2_brownian(n_patients = 60, seed = 42)

  expect_named(traj, c("id", "tx", "time", "y", "yprev", "x"))
  expect_equal(min(traj$time), 1)
  expect_equal(max(traj$time), 28)
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
