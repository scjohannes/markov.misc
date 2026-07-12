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
    y_levels = fit$yunique,
    absorb = "6",
    p_var = "yprev"
  )
  components <- markov.misc:::markov_msm_build(
    model = fit,
    data = baseline,
    times = 1:5,
    y_levels = fit$yunique,
    absorb = "6",
    p_var = "yprev"
  )
  fast <- markov.misc:::markov_msm_run(
    components,
    markov.misc::get_effective_coefs(fit),
    times = 1:5,
    absorb = "6"
  )

  expect_equal(unname(fast), unname(slow), tolerance = 1e-12)
})

test_that("orm fast path validates state support and initial probabilities", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 45, follow_up_time = 6, seed = 6111)
  data <- add_test_age(data)
  local_orm_datadist(data)
  fit <- rms::orm(
    y ~ time + tx + yprev + age,
    data = data,
    x = TRUE,
    y = TRUE
  )
  baseline <- data[!duplicated(data$id), ][1:3, , drop = FALSE]

  expect_snapshot(
    markov.misc::soprob_markov(
      fit,
      baseline,
      times = 1,
      y_levels = fit$yunique[-length(fit$yunique)]
    ),
    error = TRUE
  )

  baseline$age[1L] <- NA_real_
  expect_snapshot(
    markov.misc::soprob_markov(
      fit,
      baseline[1L, , drop = FALSE],
      times = 1,
      y_levels = fit$yunique
    ),
    error = TRUE
  )
})

test_that("execution-plan budgets bound all retained first-order designs", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 40, follow_up_time = 6, seed = 6112)
  data <- add_test_age(data)
  local_orm_datadist(data)
  fit <- rms::orm(
    y ~ time + tx + yprev + age,
    data = data,
    x = TRUE,
    y = TRUE
  )
  baseline <- data[!duplicated(data$id), ][1:4, , drop = FALSE]
  args <- list(
    model = fit,
    newdata = baseline,
    times = 1:5,
    y_levels = fit$yunique,
    absorb = "6"
  )
  expected <- do.call(markov.misc:::soprob_markov_reference, args)

  withr::local_options(markov.misc.execution_plan_max_bytes = 1)
  for (builder in c("batched", "streamed")) {
    condition <- tryCatch(
      do.call(
        markov.misc:::compile_sop_execution_plan,
        c(args, list(builder = builder))
      ),
      error = function(e) e
    )
    expect_s3_class(condition, "markov_misc_execution_plan_too_large")
    expect_gt(condition$required_bytes, condition$max_bytes)
  }
  expect_message(
    actual <- do.call(markov.misc::soprob_markov, args),
    "Compiled C\\+\\+ SOP calculations were not used"
  )
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("fallback notifications are emitted once per inference call", {
  messages <- testthat::capture_messages({
    markov.misc:::with_sop_fallback_notification_scope({
      markov.misc:::notify_sop_reference_fallback("first draw failed.")
      markov.misc:::with_sop_fallback_notification_scope({
        markov.misc:::notify_sop_reference_fallback("second draw failed.")
      })
    })
  })
  expect_length(messages, 1L)
  expect_match(messages, "Compiled C\\+\\+ SOP calculations were not used")
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

test_that("orm scale=TRUE agrees with default scale in SOP workflows", {
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  data <- make_test_data(n_patients = 45, follow_up_time = 6, seed = 6105)
  data <- add_test_age(data)
  local_orm_datadist(data)

  fit_default <- rms::orm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    x = TRUE,
    y = TRUE
  )
  fit_scaled <- rms::orm(
    y ~ rms::rcs(time, 3) * tx + yprev + age,
    data = data,
    x = TRUE,
    y = TRUE,
    scale = TRUE
  )
  baseline <- data[!duplicated(data$id), ]

  predict_scaled <- stats::predict(
    fit_scaled,
    newdata = baseline,
    type = "fitted.ind"
  )
  fast_scaled <- markov.misc:::predict_orm_response_markov(
    fit_scaled,
    baseline
  )
  expect_equal(unname(fast_scaled), unname(predict_scaled), tolerance = 1e-12)

  sops_default <- markov.misc::sops(
    fit_default,
    baseline,
    times = 1:4,
    y_levels = fit_default$yunique,
    absorb = "6"
  )
  sops_scaled <- markov.misc::sops(
    fit_scaled,
    baseline,
    times = 1:4,
    y_levels = fit_scaled$yunique,
    absorb = "6"
  )
  expect_equal(sops_scaled$estimate, sops_default$estimate, tolerance = 1e-10)

  robust_default <- rms::robcov(fit_default, cluster = data$id)
  robust_scaled <- rms::robcov(fit_scaled, cluster = data$id)
  avg_default <- markov.misc::avg_sops(
    robust_default,
    refit_data = data,
    variables = "tx",
    times = 1:4,
    y_levels = fit_default$yunique,
    absorb = "6",
    id_var = "id"
  )
  avg_scaled <- markov.misc::avg_sops(
    robust_scaled,
    refit_data = data,
    variables = "tx",
    times = 1:4,
    y_levels = fit_scaled$yunique,
    absorb = "6",
    id_var = "id"
  )
  expect_equal(avg_scaled$estimate, avg_default$estimate, tolerance = 1e-10)

  withr::local_seed(61051)
  mvn_default <- markov.misc::inferences(
    avg_default,
    method = "mvn",
    n_draws = 2,
    return_draws = TRUE
  )
  withr::local_seed(61051)
  mvn_scaled <- markov.misc::inferences(
    avg_scaled,
    method = "mvn",
    n_draws = 2,
    return_draws = TRUE
  )
  expect_equal(mvn_scaled$conf.low, mvn_default$conf.low, tolerance = 1e-8)
  expect_equal(mvn_scaled$conf.high, mvn_default$conf.high, tolerance = 1e-8)

  withr::local_seed(61052)
  score_default <- markov.misc::inferences(
    avg_default,
    method = "score_bootstrap",
    cluster = data$id,
    n_draws = 2,
    return_draws = TRUE
  )
  withr::local_seed(61052)
  score_scaled <- markov.misc::inferences(
    avg_scaled,
    method = "score_bootstrap",
    cluster = data$id,
    n_draws = 2,
    return_draws = TRUE
  )
  expect_equal(
    score_scaled$conf.low,
    score_default$conf.low,
    tolerance = 1e-8
  )
  expect_equal(
    score_scaled$conf.high,
    score_default$conf.high,
    tolerance = 1e-8
  )

  avg_scaled_full <- markov.misc::avg_sops(
    fit_scaled,
    newdata = data,
    variables = "tx",
    times = 1:4,
    y_levels = fit_scaled$yunique,
    absorb = "6",
    id_var = "id"
  )
  withr::local_seed(61053)
  boot_scaled <- markov.misc::inferences(
    avg_scaled_full,
    method = "bootstrap",
    n_draws = 1,
    update_datadist = FALSE,
    return_draws = TRUE
  )
  expect_equal(attr(boot_scaled, "method"), "bootstrap")
  expect_equal(attr(boot_scaled, "n_successful"), 1L)
  expect_false(anyNA(boot_scaled$conf.low))
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
    refit_data = data,
    variables = "tx",
    times = 1:4,
    y_levels = fit$yunique,
    absorb = "6",
    id_var = "id"
  )

  expect_error(
    markov.misc::inferences(
      avg,
      method = "mvn",
      vcov = stats::vcov(fit),
      n_draws = 2
    ),
    "Dimension mismatch"
  )

  withr::local_seed(61031)
  mvn <- markov.misc::inferences(
    avg,
    method = "mvn",
    n_draws = 2,
    return_draws = TRUE
  )
  expect_equal(attr(mvn, "n_successful"), 2L)
  expect_false(anyNA(mvn$conf.low))
  expect_s3_class(markov.misc::get_draws(mvn), "data.frame")

  expect_error(
    markov.misc::inferences(
      avg,
      method = "score_bootstrap",
      n_draws = 2
    ),
    "`cluster` is required"
  )

  withr::local_seed(61032)
  score <- markov.misc::inferences(
    avg,
    method = "score_bootstrap",
    cluster = data$id,
    n_draws = 2,
    return_draws = TRUE
  )
  expect_equal(attr(score, "method"), "score_bootstrap")
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
    y_levels = fit$yunique,
    absorb = "6",
    id_var = "id"
  )

  withr::local_seed(61041)
  boot <- markov.misc::inferences(
    avg,
    method = "bootstrap",
    n_draws = 1,
    update_datadist = FALSE,
    return_draws = TRUE
  )

  expect_equal(attr(boot, "method"), "bootstrap")
  expect_equal(attr(boot, "n_successful"), 1L)
  expect_false(anyNA(boot$conf.low))
})
