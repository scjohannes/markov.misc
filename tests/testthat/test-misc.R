test_that("calc_time_in_state_diff returns all observed states by default", {
  data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
    time = rep(1:3, 4),
    y = c(
      0, 1, 1,
      0, 0, 1,
      1, 1, 1,
      0, 1, 1
    ),
    tx = c(rep(0, 6), rep(1, 6))
  )

  result <- calc_time_in_state_diff(data)

  expect_equal(result$state, c(0, 1))
  expect_equal(result$reference_level, c(0, 0))
  expect_equal(result$treatment_level, c(1, 1))
  expect_equal(result$true_effect, c(-1, 1))
  expect_equal(result$reference_mean_time, c(1.5, 1.5))
  expect_equal(result$treatment_mean_time, c(0.5, 2.5))
  expect_equal(result$reference_sd_time, c(stats::sd(c(1, 2)), stats::sd(c(2, 1))))
  expect_equal(result$treatment_sd_time, c(stats::sd(c(0, 1)), stats::sd(c(3, 2))))
})

test_that("calc_time_in_state_diff works for an explicit single target state", {
  data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4),
    time = rep(1:3, 4),
    y = c(
      0, 1, 1,
      0, 0, 1,
      1, 1, 1,
      0, 1, 1
    ),
    tx = c(rep(0, 6), rep(1, 6))
  )

  result <- calc_time_in_state_diff(data, target_state = 1)

  expect_equal(nrow(result), 1)
  expect_equal(result$state, 1)
  expect_equal(result$reference_level, 0)
  expect_equal(result$treatment_level, 1)
  expect_equal(result$true_effect, 1)
  expect_equal(result$reference_mean_time, 1.5)
  expect_equal(result$treatment_mean_time, 2.5)
  expect_equal(result$reference_sd_time, stats::sd(c(2, 1)))
  expect_equal(result$treatment_sd_time, stats::sd(c(3, 2)))
})

test_that("calc_time_in_state_diff supports custom columns and multi-level treatment", {
  data <- data.frame(
    id = rep(1:6, each = 3),
    day = rep(1:3, 6),
    y = c(
      0, 1, 1,
      0, 0, 1,
      1, 1, 1,
      0, 1, 1,
      0, 0, 0,
      0, 1, 0
    ),
    arm = factor(
      rep(c("A", "A", "B", "B", "C", "C"), each = 3),
      levels = c("A", "B", "C")
    )
  )

  result <- calc_time_in_state_diff(
    data,
    tvarname = "day",
    txvarname = "arm"
  )

  expect_equal(as.character(result$state), c("0", "0", "1", "1"))
  expect_equal(as.character(result$reference_level), c("A", "A", "A", "A"))
  expect_equal(as.character(result$treatment_level), c("B", "C", "B", "C"))
  expect_equal(result$true_effect, c(-1, 1, 1, -1))
  expect_equal(result$reference_mean_time, c(1.5, 1.5, 1.5, 1.5))
  expect_equal(result$reference_sd_time, c(
    stats::sd(c(1, 2)),
    stats::sd(c(1, 2)),
    stats::sd(c(2, 1)),
    stats::sd(c(2, 1))
  ))
  expect_equal(result$treatment_mean_time, c(0.5, 2.5, 2.5, 0.5))
  expect_equal(
    result$treatment_sd_time,
    c(
      stats::sd(c(0, 1)),
      stats::sd(c(3, 2)),
      stats::sd(c(3, 2)),
      stats::sd(c(0, 1))
    )
  )
})

test_that("calc_time_in_state_diff allows an explicit reference level", {
  data <- data.frame(
    id = rep(1:6, each = 3),
    day = rep(1:3, 6),
    y = c(
      0, 1, 1,
      0, 0, 1,
      1, 1, 1,
      0, 1, 1,
      0, 0, 0,
      0, 1, 0
    ),
    arm = factor(
      rep(c("A", "A", "B", "B", "C", "C"), each = 3),
      levels = c("A", "B", "C")
    )
  )

  result <- calc_time_in_state_diff(
    data,
    target_state = 1,
    tvarname = "day",
    txvarname = "arm",
    reference_level = "B"
  )

  expect_equal(result$state, c(1, 1))
  expect_equal(as.character(result$reference_level), c("B", "B"))
  expect_equal(as.character(result$treatment_level), c("A", "C"))
  expect_equal(result$true_effect, c(-1, -2))
  expect_equal(result$reference_mean_time, c(2.5, 2.5))
  expect_equal(result$treatment_mean_time, c(1.5, 0.5))
})

test_that("state conversion helpers validate required columns", {
  expect_error(states_to_ttest(data.frame(id = 1)), "data must contain columns: y, tx")
  expect_error(states_to_drs(data.frame(id = 1), follow_up_time = 3), "data must contain columns")
  expect_error(states_to_tte_old(data.frame(id = 1)), "data must contain columns")
  expect_error(states_to_tte(data.frame(id = 1)), "data must contain columns")
  expect_error(format_competing_risks(data.frame(id = 1)), "data must contain columns")
})

test_that("prepare_markov_data can keep previous state numeric", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(2, 3, 4, 5),
    yprev = c(1, 2, 4, 6)
  )

  prepared <- prepare_markov_data(
    data,
    absorbing_state = 6,
    factor_previous = FALSE
  )

  expect_type(prepared$yprev, "double")
  expect_s3_class(prepared$y, "ordered")
  expect_equal(prepared$yprev, c(1, 2, 4))
})

test_that("states_to_ttest counts target-state time per patient", {
  data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    y = c(1, 2, 1, 2, 2, 1),
    tx = c(0, 0, 0, 1, 1, 1)
  )

  result <- states_to_ttest(data, target_state = 1)

  expect_equal(result$id, c(1, 2))
  expect_equal(result$y, c(2, 1))
  expect_equal(result$tx, c(0, 1))
})

test_that("states_to_drs handles death, sustained home time, and covariates", {
  data <- data.frame(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    y = c(2, 1, 1, 1, 2, 2, 1, 1, 2, 6, 6, 6),
    tx = rep(c(0, 1, 1), each = 4),
    age = rep(c(50, 60, 70), each = 4)
  )

  result <- states_to_drs(data, follow_up_time = 4, covariates = c("age", "missing"))

  expect_equal(result$id, 1:3)
  expect_equal(result$drs, c(3, 2, -1))
  expect_equal(result$age, c(50, 60, 70))

  no_covs <- states_to_drs(data, follow_up_time = 4, covariates = NULL)
  expect_equal(names(no_covs), c("id", "tx", "drs"))
})

test_that("states_to_tte_old and states_to_tte collapse trajectories", {
  data <- data.frame(
    id = rep(1:2, each = 5),
    time = rep(1:5, 2),
    y = c(2, 2, 1, 1, 1, 3, 3, 4, 6, 6),
    yprev = c(2, 2, 2, 1, 1, 3, 3, 3, 4, 6),
    tx = rep(c(0, 1), each = 5),
    age = rep(c(50, 60), each = 5)
  )

  old <- states_to_tte_old(data, covariates = "age")
  current <- states_to_tte(data, covariates = "age", absorbing_state = 6)

  expect_equal(names(old), c("id", "tx", "yprev", "start", "stop", "y"))
  expect_equal(current$id, c(1, 1, 2, 2, 2))
  expect_equal(current$start, c(0, 2, 0, 2, 3))
  expect_equal(current$stop, c(2, 5, 2, 3, 4))
  expect_equal(current$age, c(50, 50, 60, 60, 60))
})

test_that("format_competing_risks supports old and new status logic", {
  data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2, 3, 3),
    tx = c(0, 0, 0, 1, 1, 1, 0, 0),
    start = c(0, 1, 2, 0, 1, 2, 0, 1),
    stop = c(1, 2, 3, 1, 2, 3, 1, 2),
    y = c(3, 2, 1, 3, 2, 6, 3, 2)
  )

  old <- format_competing_risks(data, event_status = 1, death_status = 6)
  new <- format_competing_risks(
    data,
    event_status = 1,
    death_status = 6,
    version_old = FALSE
  )

  expect_equal(as.character(old$status), c("1", "2", "0"))
  expect_equal(old$etime, c(3, 3, 2))
  expect_equal(new$status, c(0, 1, 0, 2, 0, 0))
  expect_false("next_y" %in% names(new))
})

test_that("calc_time_in_state_diff validates inputs", {
  data <- data.frame(
    id = c(1, 2),
    time = c(1, 1),
    y = c(0, 1),
    tx = c(0, 0)
  )

  expect_error(calc_time_in_state_diff(data.frame(id = 1)), "data must contain columns")
  expect_error(
    calc_time_in_state_diff(data, target_state = 3),
    "target_state must be drawn from observed y values"
  )
  expect_error(
    calc_time_in_state_diff(data),
    "txvarname must contain at least 2 observed treatment levels"
  )
  expect_error(
    calc_time_in_state_diff(
      transform(data, tx = c(0, 1)),
      reference_level = 2
    ),
    "reference_level must be one of"
  )
})
