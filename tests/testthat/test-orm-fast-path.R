local_orm_datadist <- function(data) {
  env <- parent.frame()
  dd_name <- "dd_orm_test"
  assign(dd_name, rms::datadist(data), envir = globalenv())
  withr::defer(
    if (exists(dd_name, envir = globalenv())) {
      rm(list = dd_name, envir = globalenv())
    },
    envir = env
  )
  withr::local_options(list(datadist = dd_name), .local_envir = env)
}

add_test_age <- function(data) {
  data$age <- 45 + (data$id %% 25)
  data
}

test_that("orm fast path supports categorical previous state", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 60, follow_up_time = 7, seed = 6101)
  data <- add_test_age(data)
  local_orm_datadist(data)

  fit <- rms::orm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    x = TRUE,
    y = TRUE
  )
  baseline <- data[!duplicated(data$id), ][1:12, , drop = FALSE]

  direct <- markov.misc:::predict_orm_response_markov(fit, baseline)
  expect_equal(dim(direct), c(nrow(baseline), length(fit$yunique)))
  expect_equal(rowSums(direct), rep(1, nrow(baseline)), tolerance = 1e-12)

  slow <- markov.misc::soprob_markov(
    fit,
    baseline,
    times = 1:5,
    ylevels = fit$yunique,
    absorb = "6",
    pvarname = "yprev"
  )
  components <- markov.misc:::markov_msm_build(
    model = fit,
    data = baseline,
    times = 1:5,
    ylevels = fit$yunique,
    absorb = "6",
    pvarname = "yprev"
  )
  fast <- markov.misc:::markov_msm_run(
    components,
    markov.misc::get_effective_coefs(fit),
    times = 1:5,
    absorb = "6"
  )

  expect_equal(unname(fast), unname(slow), tolerance = 1e-12)
})

test_that("orm fast path agrees with predict.orm for numeric previous state", {
  skip_if_not_installed("rms")

  data <- markov.misc::sim_trajectories_brownian(
    n_patients = 70,
    follow_up_time = 7,
    seed = 6102,
    absorbing_state = 6
  ) |>
    markov.misc::prepare_markov_data(
      absorbing_state = 6,
      factor_previous = FALSE
    )
  local_orm_datadist(data)

  fit <- suppressWarnings(
    rms::orm(
      y ~ rms::rcs(time, 3) * tx + rms::rcs(yprev, 3),
      data = data,
      x = TRUE,
      y = TRUE
    )
  )
  newdata <- data[1:15, , drop = FALSE]

  direct <- markov.misc:::predict_orm_response_markov(fit, newdata)
  predicted <- stats::predict(fit, newdata = newdata, type = "fitted.ind")

  expect_equal(unname(direct), unname(predicted), tolerance = 1e-12)
})

test_that("orm supports MVN and score-bootstrap inference", {
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  data <- make_test_data(n_patients = 55, follow_up_time = 6, seed = 6103)
  data <- add_test_age(data)
  local_orm_datadist(data)

  fit <- rms::orm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    x = TRUE,
    y = TRUE
  )
  robust_fit <- rms::robcov(fit, cluster = data$id)
  baseline <- data[!duplicated(data$id), ]

  avg <- markov.misc::avg_sops(
    robust_fit,
    newdata = baseline,
    variables = "tx",
    times = 1:4,
    ylevels = fit$yunique,
    absorb = "6",
    id_var = "id"
  )

  expect_error(
    markov.misc::inferences(
      avg,
      method = "simulation",
      vcov = stats::vcov(fit),
      n_sim = 2
    ),
    "Dimension mismatch"
  )

  withr::local_seed(61031)
  mvn <- markov.misc::inferences(
    avg,
    method = "simulation",
    engine = "mvn",
    n_sim = 2,
    return_draws = TRUE
  )
  expect_equal(attr(mvn, "n_successful"), 2L)
  expect_false(anyNA(mvn$conf.low))
  expect_s3_class(markov.misc::get_draws(mvn), "data.frame")

  expect_error(
    markov.misc::inferences(
      avg,
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 2
    ),
    "`cluster` is required"
  )

  withr::local_seed(61032)
  score <- markov.misc::inferences(
    avg,
    method = "simulation",
    engine = "score_bootstrap",
    cluster = data$id,
    n_sim = 2,
    return_draws = TRUE
  )
  expect_equal(attr(score, "engine"), "score_bootstrap")
  expect_equal(attr(score, "n_successful"), 2L)
  expect_true(any(score$std.error > 0))
})

test_that("orm supports refit bootstrap inference", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 45, follow_up_time = 6, seed = 6104)
  data <- add_test_age(data)
  local_orm_datadist(data)

  fit <- rms::orm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    x = TRUE,
    y = TRUE
  )

  avg <- markov.misc::avg_sops(
    fit,
    newdata = data,
    variables = "tx",
    times = 1:4,
    ylevels = fit$yunique,
    absorb = "6",
    id_var = "id"
  )

  withr::local_seed(61041)
  boot <- markov.misc::inferences(
    avg,
    method = "bootstrap",
    n_sim = 1,
    update_datadist = FALSE,
    return_draws = TRUE
  )

  expect_equal(attr(boot, "method"), "bootstrap")
  expect_equal(attr(boot, "n_successful"), 1L)
  expect_false(anyNA(boot$conf.low))
})
