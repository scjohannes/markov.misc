test_that("calc_time_in_state_diff returns all observed states by default", {
  data <- tibble::tibble(
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
  data <- tibble::tibble(
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
  data <- tibble::tibble(
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
  data <- tibble::tibble(
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
