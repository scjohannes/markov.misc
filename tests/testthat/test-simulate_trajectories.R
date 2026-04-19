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

  traj_null <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 12,
    threshold_drift_deviations = NULL,
    seed = 123
  )

  traj_zero <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 12,
    threshold_drift_deviations = rep(0, 5),
    seed = 123
  )

  expect_identical(traj_default, traj_null)
  expect_identical(traj_default, traj_zero)
})

test_that("sim_trajectories_brownian validates threshold drift deviations", {
  expect_match(
    error_message(
      sim_trajectories_brownian(
        n_patients = 5,
        follow_up_time = 3,
        threshold_drift_deviations = "bad",
        seed = 1
      )
    ),
    "threshold_drift_deviations must be a numeric vector without NA values"
  )

  expect_match(
    error_message(
      sim_trajectories_brownian(
        n_patients = 5,
        follow_up_time = 3,
        threshold_drift_deviations = c(0.1, 0.2),
        seed = 1
      )
    ),
    "threshold_drift_deviations must have length n_states - 1"
  )

  expect_match(
    error_message(
      sim_trajectories_brownian(
        n_patients = 5,
        follow_up_time = 20,
        thresholds = c(0, 1, 2, 3, 4),
        threshold_drift_deviations = c(0, 0, -0.2, -0.4, -0.6),
        drift_change_times = NULL,
        seed = 1
      )
    ),
    "threshold_drift_deviations produce non-increasing thresholds over follow-up"
  )
})

test_that("sim_trajectories_brownian supports threshold-specific time effects", {
  traj_default <- sim_trajectories_brownian(
    n_patients = 300,
    follow_up_time = 15,
    drift_change_times = NULL,
    seed = 999
  )

  traj_nonparallel <- sim_trajectories_brownian(
    n_patients = 300,
    follow_up_time = 15,
    threshold_drift_deviations = c(0.04, 0.02, 0, -0.02, -0.04),
    drift_change_times = NULL,
    seed = 999
  )

  expect_identical(names(traj_default), names(traj_nonparallel))
  expect_identical(dim(traj_default), dim(traj_nonparallel))
  expect_gt(sum(traj_default$y != traj_nonparallel$y), 0)

  mean_default <- mean(traj_default$y[traj_default$time == 15])
  mean_nonparallel <- mean(traj_nonparallel$y[traj_nonparallel$time == 15])
  expect_gt(abs(mean_default - mean_nonparallel), 0.01)
})

test_that("sim_trajectories_brownian keeps the absorbing state with deviations", {
  traj <- sim_trajectories_brownian(
    n_patients = 80,
    follow_up_time = 20,
    threshold_drift_deviations = c(0.02, 0.01, 0, -0.01, -0.02),
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
