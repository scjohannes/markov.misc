test_that("score-bootstrap simulation adds CIs and metadata for avg_sops", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2001)
  data <- make_test_data(n_patients = 60, seed = 2001, follow_up_time = 10)
  baseline_data <- data[data$time == 1, ]
  model <- make_test_model(data, robust = TRUE)

  result <- avg_sops(
    model = model,
    newdata = baseline_data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  ) |>
    inferences(
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 40,
      return_draws = TRUE
    )

  expect_contains(names(result), c("conf.low", "conf.high", "std.error"))
  expect_equal(attr(result, "method"), "simulation")
  expect_equal(attr(result, "engine"), "score_bootstrap")
  expect_equal(attr(result, "score_weight_dist"), "exponential")
  expect_equal(attr(result, "n_sim"), 40)
  expect_true(attr(result, "n_successful") <= 40)
  expect_true(attr(result, "n_successful") > 0)
  expect_all_true(result$std.error >= 0)
  expect_true(any(result$std.error > 0, na.rm = TRUE))
})

test_that("score-bootstrap simulation stores draw-level output with expected structure", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2002)
  data <- make_test_data(n_patients = 60, seed = 2002, follow_up_time = 10)
  baseline_data <- data[data$time == 1, ]
  model <- make_test_model(data, robust = TRUE)

  result <- avg_sops(
    model = model,
    newdata = baseline_data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  ) |>
    inferences(
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 30,
      return_draws = TRUE
    )

  draws <- get_draws(result)

  expect_s3_class(draws, "data.frame")
  expect_contains(names(draws), c("draw_id", "estimate", "time", "state", "tx"))
  expect_equal(length(unique(draws$draw_id)), attr(result, "n_successful"))
  expect_all_true(is.finite(draws$draw_id))
  expect_all_true(is.finite(draws$estimate))
})

test_that("score-bootstrap simulation does not expose draws when return_draws is FALSE", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2003)
  data <- make_test_data(n_patients = 60, seed = 2003, follow_up_time = 10)
  baseline_data <- data[data$time == 1, ]
  model <- make_test_model(data, robust = TRUE)

  result <- avg_sops(
    model = model,
    newdata = baseline_data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  ) |>
    inferences(
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 20,
      return_draws = FALSE
    )

  expect_error(get_draws(result), "No draws found")
})

test_that("score-bootstrap simulation rejects user-supplied vcov", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2004)
  data <- make_test_data(n_patients = 60, seed = 2004, follow_up_time = 10)
  baseline_data <- data[data$time == 1, ]
  model <- make_test_model(data, robust = TRUE)
  custom_vcov <- diag(length(coef(model)))

  avg_result <- avg_sops(
    model = model,
    newdata = baseline_data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  )

  expect_error(
    inferences(
      avg_result,
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 10,
      vcov = custom_vcov
    ),
    "cannot be supplied"
  )
})

test_that("score-bootstrap simulation validates score_weight_dist", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2005)
  data <- make_test_data(n_patients = 60, seed = 2005, follow_up_time = 10)
  baseline_data <- data[data$time == 1, ]
  model <- make_test_model(data, robust = TRUE)

  avg_result <- avg_sops(
    model = model,
    newdata = baseline_data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  )

  expect_error(
    inferences(
      avg_result,
      method = "simulation",
      engine = "score_bootstrap",
      score_weight_dist = "invalid_dist",
      n_sim = 10
    ),
    "exponential"
  )
})

test_that("generate_score_bootstrap_draws returns valid dimensions and normalized weights", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2006)
  data <- make_test_data(n_patients = 35, seed = 2006, follow_up_time = 8)
  model <- make_test_model(data, robust = TRUE)
  baseline_data <- data[data$time == 1, ]

  draws <- generate_score_bootstrap_draws(
    model = model,
    baseline_data = baseline_data,
    id_var = "id",
    n_sim = 25
  )

  expect_true(is.matrix(draws$beta_draws))
  expect_true(is.matrix(draws$baseline_weights))
  expect_equal(nrow(draws$beta_draws), 25)
  expect_equal(ncol(draws$beta_draws), length(model$coefficients))
  expect_equal(nrow(draws$baseline_weights), 25)
  expect_equal(ncol(draws$baseline_weights), nrow(baseline_data))
  expect_equal(
    as.numeric(rowSums(draws$baseline_weights)),
    rep(1, 25),
    tolerance = 1e-10
  )
  expect_true(all(is.finite(draws$beta_draws)))
  expect_true(all(is.finite(draws$baseline_weights)))
  expect_true(all(draws$baseline_weights >= 0))
})

test_that("generate_score_bootstrap_draws is reproducible with a fixed seed", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 30, seed = 2007, follow_up_time = 8)
  model <- make_test_model(data, robust = TRUE)
  baseline_data <- data[data$time == 1, ]

  set.seed(9999)
  draws_a <- generate_score_bootstrap_draws(
    model = model,
    baseline_data = baseline_data,
    id_var = "id",
    n_sim = 15
  )

  set.seed(9999)
  draws_b <- generate_score_bootstrap_draws(
    model = model,
    baseline_data = baseline_data,
    id_var = "id",
    n_sim = 15
  )

  expect_equal(draws_a$beta_draws, draws_b$beta_draws)
  expect_equal(draws_a$baseline_weights, draws_b$baseline_weights)
})

test_that("generate_score_bootstrap_draws errors when baseline IDs do not match clusters", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2008)
  data <- make_test_data(n_patients = 30, seed = 2008, follow_up_time = 8)
  model <- make_test_model(data, robust = TRUE)
  baseline_data <- data[data$time == 1, ]
  baseline_data$id <- as.character(baseline_data$id)
  baseline_data$id[1] <- "missing_cluster_id"

  expect_error(
    generate_score_bootstrap_draws(
      model = model,
      baseline_data = baseline_data,
      id_var = "id",
      n_sim = 10
    ),
    "were not found in the cluster variable"
  )
})

test_that("generate_score_bootstrap_draws errors when required robust components are missing", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(2009)
  data <- make_test_data(n_patients = 30, seed = 2009, follow_up_time = 8)
  model <- make_test_model(data, robust = TRUE)
  baseline_data <- data[data$time == 1, ]

  broken_model <- model
  broken_model$scores <- NULL

  expect_error(
    generate_score_bootstrap_draws(
      model = broken_model,
      baseline_data = baseline_data,
      id_var = "id",
      n_sim = 10
    ),
    "missing required components.*scores"
  )
})
