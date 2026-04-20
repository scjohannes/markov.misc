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

test_that("sim_trajectories_actt2 returns ACTT-2-shaped output", {
  traj <- sim_trajectories_actt2(n_patients = 60, seed = 42)

  expect_named(traj, c("id", "tx", "time", "y", "yprev", "x"))
  expect_equal(min(traj$time), 1)
  expect_equal(max(traj$time), 28)
  expect_true(all(traj$yprev[traj$time == 1] %in% 4:7))
})

test_that("sim_trajectories_actt2 uses a null treatment effect by default", {
  traj_default <- sim_trajectories_actt2(n_patients = 40, seed = 99)
  traj_null <- sim_trajectories_actt2(
    n_patients = 40,
    mu_treatment_effect = 0,
    seed = 99
  )

  expect_identical(traj_default, traj_null)
})
