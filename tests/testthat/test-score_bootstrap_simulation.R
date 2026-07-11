test_that("score-bootstrap simulation adds CIs and metadata for avg_sops", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2001)
  withr::local_seed(2001)
  result <- avg_sops(
    model = case$model,
    refit_data = case$data,
    variables = list(tx = c(0, 1)),
    times = case$times,
    y_levels = case$y_levels,
    absorb = case$absorb,
    id_var = "id"
  ) |>
    inferences(
      method = "score_bootstrap",
      n_draws = 40,
      return_draws = TRUE
    )

  expect_equal(attr(result, "method"), "score_bootstrap")
  expect_equal(attr(result, "n_draws"), 40)
  expect_true(attr(result, "n_successful") <= 40)
  expect_true(attr(result, "n_successful") > 0)
  expect_inference_intervals(result, require_positive_std_error = TRUE)
})

test_that("score-bootstrap simulation stores draw-level output with expected structure", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2002)
  withr::local_seed(2002)
  result <- avg_sops(
    model = case$model,
    refit_data = case$data,
    variables = list(tx = c(0, 1)),
    times = case$times,
    y_levels = case$y_levels,
    absorb = case$absorb,
    id_var = "id"
  ) |>
    inferences(
      method = "score_bootstrap",
      n_draws = 30,
      return_draws = TRUE
    )

  draws <- get_draws(result)

  expect_s3_class(draws, "data.frame")
  expect_contains(names(draws), c("draw_id", "estimate", "time", "state", "tx"))
  expect_equal(length(unique(draws$draw_id)), attr(result, "n_successful"))
  expect_all_true(is.finite(draws$draw_id))
  expect_all_true(is.finite(draws$estimate))
})

test_that("score-bootstrap treats supplied newdata as fixed profiles", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2010)
  profiles <- case$baseline[seq_len(5), , drop = FALSE]
  profiles <- profiles[, setdiff(names(profiles), "id"), drop = FALSE]

  avg <- avg_sops(
    model = case$model,
    newdata = profiles,
    refit_data = case$data,
    variables = list(tx = c(0, 1)),
    times = case$times[1:2],
    y_levels = case$y_levels,
    absorb = case$absorb,
    id_var = "id"
  )

  withr::local_seed(2010)
  expect_warning(
    result <- inferences(
      avg,
      method = "score_bootstrap",
      n_draws = 3,
      return_draws = TRUE
    ),
    "baseline weights have been set to NULL"
  )

  draws <- get_draws(result)

  expect_equal(attr(result, "n_successful"), 3L)
  expect_equal(nrow(attr(avg, "newdata_pred")), nrow(profiles) * 2L)
  expect_contains(names(draws), c("draw_id", "estimate", "time", "state", "tx"))
  expect_false("id" %in% names(draws))
})

test_that("score-bootstrap simulation does not expose draws when return_draws is FALSE", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2003)
  withr::local_seed(2003)
  result <- avg_sops(
    model = case$model,
    refit_data = case$data,
    variables = list(tx = c(0, 1)),
    times = case$times,
    y_levels = case$y_levels,
    absorb = case$absorb,
    id_var = "id"
  ) |>
    inferences(
      method = "score_bootstrap",
      n_draws = 20,
      return_draws = FALSE
    )

  expect_error(get_draws(result), "No draws found")
})

test_that("score-bootstrap simulation rejects user-supplied vcov", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2004)
  custom_vcov <- diag(length(coef(case$model)))

  avg_result <- avg_sops(
    model = case$model,
    refit_data = case$data,
    variables = list(tx = c(0, 1)),
    times = case$times,
    y_levels = case$y_levels,
    absorb = case$absorb,
    id_var = "id"
  )

  expect_error(
    inferences(
      avg_result,
      method = "score_bootstrap",
      n_draws = 10,
      vcov = custom_vcov
    ),
    "cannot be supplied"
  )
})

test_that("score-bootstrap uses the fixed exponential weight distribution", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(seed = 2005)
  avg_result <- avg_sops(
    model = case$model,
    refit_data = case$data,
    variables = list(tx = c(0, 1)),
    times = case$times,
    y_levels = case$y_levels,
    absorb = case$absorb,
    id_var = "id"
  )

  expect_false("score_weight_dist" %in% names(formals(inferences)))
  expect_s3_class(
    inferences(avg_result, method = "score_bootstrap", n_draws = 10),
    "markov_avg_sops"
  )
})

test_that("generate_score_bootstrap_draws returns valid dimensions and normalized weights", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  case <- make_score_bootstrap_case(
    seed = 2006,
    n_patients = 35,
    follow_up_time = 8
  )

  draws <- generate_score_bootstrap_draws(
    model = case$model,
    baseline_data = case$baseline,
    id_var = "id",
    n_draws = 25
  )

  expect_true(is.matrix(draws$beta_draws))
  expect_true(is.matrix(draws$baseline_weights))
  expect_equal(nrow(draws$beta_draws), 25)
  expect_equal(ncol(draws$beta_draws), length(case$model$coefficients))
  expect_equal(nrow(draws$baseline_weights), 25)
  expect_equal(ncol(draws$baseline_weights), nrow(case$baseline))
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

  withr::local_seed(9999)
  draws_a <- generate_score_bootstrap_draws(
    model = model,
    baseline_data = baseline_data,
    id_var = "id",
    n_draws = 15
  )

  withr::local_seed(9999)
  draws_b <- generate_score_bootstrap_draws(
    model = model,
    baseline_data = baseline_data,
    id_var = "id",
    n_draws = 15
  )

  expect_equal(draws_a$beta_draws, draws_b$beta_draws)
  expect_equal(draws_a$baseline_weights, draws_b$baseline_weights)
})

test_that("generate_score_bootstrap_draws errors when baseline IDs do not match clusters", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  withr::local_seed(2008)
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
      n_draws = 10
    ),
    "were not found in the cluster variable"
  )
})

test_that("generate_score_bootstrap_draws errors when required robust components are missing", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  withr::local_seed(2009)
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
      n_draws = 10
    ),
    "missing required components.*scores"
  )
})
