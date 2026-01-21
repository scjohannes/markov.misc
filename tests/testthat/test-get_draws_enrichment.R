
test_that("get_draws() preserves covariates for sops objects", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("dplyr")

  # Simulate test data
  data <- make_test_data(n_patients = 20, seed = 123, follow_up_time = 30)

  # Add dummy age covariate manually since brownian sim doesn't provide it
  data$age <- rnorm(nrow(data), 60, 10)

  # Ensure we have a covariate to check
  expect_true("tx" %in% names(data))
  expect_true("id" %in% names(data))
  expect_true("age" %in% names(data))

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev + age,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  m_vglm_rob <- robcov_vglm(m_vglm, cluster = data$id)

  # Compute individual sops
  sops_result <- sops(
    model = m_vglm_rob,
    newdata = data |> filter(time == 1),
    times = 1:5,
    ylevels = 1:6,
    absorb = 6,
    tvarname = "time",
    pvarname = "yprev"
  )

  # Add simulation-based inference
  result_ci <- inferences(
    sops_result,
    method = "simulation",
    n_sim = 10,
    return_draws = TRUE
  )

  # Extract draws
  draws <- get_draws(result_ci)

  # Check structure
  expect_s3_class(draws, "data.frame")

  # Core columns
  expect_true("draw_id" %in% names(draws))
  expect_true("rowid" %in% names(draws))
  expect_true("time" %in% names(draws))
  expect_true("state" %in% names(draws))
  expect_true("estimate" %in% names(draws))

  # Covariates
  expect_true("tx" %in% names(draws))
  expect_true("age" %in% names(draws))

  # Check content
  # Pick a random row from draws and check if tx matches the original data for that rowid
  sample_row <- draws[1, ]
  orig_row <- data[sample_row$rowid, ]

  expect_equal(sample_row$tx, orig_row$tx)
  expect_equal(sample_row$age, orig_row$age)
})

