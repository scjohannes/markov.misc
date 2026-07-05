test_that("soprob_markov matches Hmisc::soprobMarkovOrdm for a proportional-odds VGAM model", {
  skip_if_not_installed("Hmisc")
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 18L
  data <- make_test_data(
    n_patients = 60,
    seed = 1234,
    follow_up_time = follow_up
  )
  model <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  expected <- Hmisc::soprobMarkovOrdm(
    model,
    data = data[1, ],
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6"
  )
  actual <- soprob_markov(
    model,
    data = data[1, ],
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6"
  )

  expect_probability_array(actual, expected_dim = c(1L, follow_up, 6L))
  expect_equal(actual[1, , ], expected, tolerance = 1e-10)
})

test_that("soprob_markov gives the same predictions for inline and precomputed spline bases", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 16L
  data <- make_test_data(
    n_patients = 60,
    seed = 234,
    follow_up_time = follow_up
  )
  time_spline <- rms::rcs(data$time, 4)
  data$time_lin <- as.vector(time_spline[, 1])
  data$time_nlin_1 <- as.vector(time_spline[, 2])
  data$time_nlin_2 <- as.vector(time_spline[, 3])
  t_covs <- make_time_covariates(
    data,
    "time",
    "time_lin",
    "time_nlin_1",
    "time_nlin_2"
  )

  inline_model <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  precomputed_model <- VGAM::vglm(
    ordered(y) ~ time_lin + time_nlin_1 + time_nlin_2 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  baseline <- data[data$time == 1, , drop = FALSE][seq_len(6), , drop = FALSE]
  inline_sops <- soprob_markov(
    inline_model,
    data = baseline,
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6"
  )
  precomputed_sops <- soprob_markov(
    precomputed_model,
    data = baseline,
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  expect_probability_array(inline_sops, expected_dim = c(6L, follow_up, 6L))
  expect_probability_array(
    precomputed_sops,
    expected_dim = c(6L, follow_up, 6L)
  )
  expect_equal(inline_sops, precomputed_sops, tolerance = 1e-10)
})

test_that("soprob_markov supports inline splines of numeric previous state", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  raw_data <- sim_trajectories_brownian(
    n_patients = 100,
    follow_up_time = 7,
    n_states = 10,
    thresholds = seq(-3, 3, length.out = 9),
    allowed_start_state = as.integer(2:9),
    absorbing_state = 10,
    x0_sd = 2,
    sigma_rw = 0.25,
    seed = 4242
  )
  data <- prepare_markov_data(
    raw_data,
    absorbing_state = 10,
    factor_previous = FALSE
  )
  baseline <- data[data$time == 1, , drop = FALSE][seq_len(8), , drop = FALSE]
  fit_data <- data[, setdiff(names(data), "id"), drop = FALSE]
  model <- vglm_markov(
    ordered(y) ~ rms::rcs(time, 4) + tx + rms::rcs(yprev, 6),
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = fit_data
  )

  slow <- soprob_markov(
    object = model,
    data = baseline,
    times = 1:7,
    ylevels = 1:10,
    absorb = "10",
    pvarname = "yprev"
  )
  components <- markov.misc:::markov_msm_build(
    model = model,
    data = baseline,
    times = 1:7,
    ylevels = 1:10,
    pvarname = "yprev"
  )
  gamma <- markov.misc:::compute_Gamma(
    stats::coef(model),
    VGAM::constraints(model)
  )
  fast <- markov.misc:::markov_msm_run(
    components = components,
    Gamma = gamma,
    times = 1:7,
    absorb = "10"
  )

  expect_type(data$yprev, "double")
  expect_probability_array(slow, expected_dim = c(8L, 7L, 10L))
  expect_probability_array(fast, expected_dim = c(8L, 7L, 10L))
  expect_equal(unname(fast), unname(slow), tolerance = 1e-10)
})

test_that("soprob_markov handles partial proportional-odds constraints", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 18L
  data <- make_test_data(
    n_patients = 70,
    seed = 456,
    follow_up_time = follow_up
  )
  time_spline <- rms::rcs(data$time, 4)
  data$time_lin <- as.vector(time_spline[, 1])
  data$time_nlin_1 <- as.vector(time_spline[, 2])
  data$time_nlin_2 <- as.vector(time_spline[, 3])
  t_covs <- make_time_covariates(
    data,
    "time",
    "time_lin",
    "time_nlin_1",
    "time_nlin_2"
  )

  n_thresholds <- 5L
  constraints <- list(
    "(Intercept)" = diag(n_thresholds),
    "yprev" = cbind(PO_effect = 1),
    "tx" = cbind(PO_effect = 1),
    "time_lin" = cbind(
      time_lin_1 = 1,
      time_lin_5 = 0:(n_thresholds - 1)
    ),
    "time_nlin_1" = cbind(PO_effect = 1),
    "time_nlin_2" = cbind(PO_effect = 1)
  )
  model <- VGAM::vglm(
    ordered(y) ~ time_lin + time_nlin_1 + time_nlin_2 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ time_lin),
    data = data,
    constraints = constraints
  )
  model@call$constraints <- constraints

  result <- soprob_markov(
    model,
    data = data[1, ],
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  expect_probability_array(result, expected_dim = c(1L, follow_up, 6L))
  expect_absorbing_probability_monotone(result, absorbing_state = 6)
})

test_that("soprob_markov incorporates treatment interactions into predictions", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 18L
  data <- make_test_data(
    n_patients = 70,
    seed = 789,
    follow_up_time = follow_up
  )
  time_spline <- rms::rcs(data$time, 4)
  data$time_lin <- as.vector(time_spline[, 1])
  data$time_nlin_1 <- as.vector(time_spline[, 2])
  data$time_nlin_2 <- as.vector(time_spline[, 3])
  t_covs <- make_time_covariates(
    data,
    "time",
    "time_lin",
    "time_nlin_1",
    "time_nlin_2"
  )

  n_thresholds <- 5L
  constraints <- list(
    "(Intercept)" = diag(n_thresholds),
    "yprev" = cbind(PO_effect = 1),
    "tx" = cbind(PO_effect = 1),
    "time_lin" = cbind(
      time_lin_1 = 1,
      time_lin_5 = 0:(n_thresholds - 1)
    ),
    "time_nlin_1" = cbind(PO_effect = 1),
    "time_nlin_2" = cbind(PO_effect = 1),
    "time_lin:tx" = cbind(PO_effect = 1),
    "time_nlin_1:tx" = cbind(PO_effect = 1),
    "time_nlin_2:tx" = cbind(PO_effect = 1),
    "time_lin:yprev" = cbind(PO_effect = 1)
  )
  model <- VGAM::vglm(
    ordered(y) ~ (time_lin + time_nlin_1 + time_nlin_2) *
      tx +
      yprev +
      yprev:time_lin,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ time_lin),
    data = data,
    constraints = constraints
  )
  model@call$constraints <- constraints

  baseline <- data[data$time == 1, , drop = FALSE][1, , drop = FALSE]
  control <- baseline
  treated <- baseline
  control$tx <- 0
  treated$tx <- 1

  treated_sops <- soprob_markov(
    model,
    data = treated,
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )
  control_sops <- soprob_markov(
    model,
    data = control,
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  expect_probability_array(treated_sops, expected_dim = c(1L, follow_up, 6L))
  expect_probability_array(control_sops, expected_dim = c(1L, follow_up, 6L))
  expect_gt(max(abs(treated_sops - control_sops)), 1e-6)
})

test_that("soprob_markov carries absorbing-state mass forward", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 12L
  data <- make_test_data(
    n_patients = 40,
    seed = 999,
    follow_up_time = follow_up
  )
  model <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  result <- soprob_markov(
    model,
    data = data[seq_len(8), ],
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = "6"
  )

  expect_probability_array(result, expected_dim = c(8L, follow_up, 6L))
  expect_absorbing_probability_monotone(result, absorbing_state = 6)
})

test_that("soprob_markov normalizes transitions when no state is absorbing", {
  skip_if_not_installed("VGAM")

  follow_up <- 10L
  data <- as.data.frame(sim_trajectories_brownian(
    n_patients = 80,
    follow_up_time = follow_up,
    treatment_prob = 0.5,
    absorbing_state = NULL,
    allowed_start_state = as.integer(1:6),
    thresholds = c(-2, -1, 0, 1, 2),
    x0_sd = 2,
    seed = 1000
  ))
  data$yprev <- factor(data$yprev, levels = 1:6)
  baseline <- data[data$time == 1, , drop = FALSE][seq_len(5), , drop = FALSE]
  model <- VGAM::vglm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  result <- soprob_markov(
    model,
    data = baseline,
    time = seq_len(follow_up),
    ylevels = factor(1:6),
    absorb = NULL
  )

  expect_probability_array(result, expected_dim = c(5L, follow_up, 6L))
})

test_that("avg_sops() equals manual G-computation over individual SOPs", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  follow_up <- 6L
  data <- make_test_data(
    n_patients = 40,
    follow_up_time = follow_up,
    seed = 123
  )
  t_covs <- make_time_covariates(data, "time", "time_lin", "time_nlin_1")
  model <- VGAM::vglm(
    ordered(y) ~ (time_lin + time_nlin_1) * tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  robust_model <- robcov_vglm(model, cluster = data$id)
  baseline <- data[data$time == 1, , drop = FALSE]

  baseline_tx1 <- baseline
  baseline_tx1$tx <- 1
  manual_tx1 <- sops(
    robust_model,
    newdata = baseline_tx1,
    times = seq_len(follow_up),
    ylevels = 1:6,
    absorb = "6",
    id_var = "id",
    t_covs = t_covs
  )
  manual_tx1 <- aggregate(
    estimate ~ time + state,
    data = manual_tx1,
    FUN = mean
  )
  manual_tx1$tx <- 1

  baseline_tx0 <- baseline
  baseline_tx0$tx <- 0
  manual_tx0 <- sops(
    robust_model,
    newdata = baseline_tx0,
    times = seq_len(follow_up),
    ylevels = 1:6,
    absorb = "6",
    id_var = "id",
    t_covs = t_covs
  )
  manual_tx0 <- aggregate(
    estimate ~ time + state,
    data = manual_tx0,
    FUN = mean
  )
  manual_tx0$tx <- 0

  manual <- rbind(manual_tx0, manual_tx1)
  manual <- manual[order(manual$tx, manual$state, manual$time), ]
  rownames(manual) <- NULL
  automated <- avg_sops(
    robust_model,
    newdata = baseline,
    variables = "tx",
    times = seq_len(follow_up),
    ylevels = 1:6,
    absorb = "6",
    id_var = "id",
    t_covs = t_covs
  )
  automated <- as.data.frame(automated)
  automated <- automated[order(automated$tx, automated$state, automated$time), ]
  rownames(automated) <- NULL

  expect_equal(manual$estimate, automated$estimate, tolerance = 1e-15)
  expect_equal(manual$tx, automated$tx)
  expect_equal(as.character(manual$state), as.character(automated$state))
  expect_equal(manual$time, automated$time)
})
