test_that("soprobMarkovOrdm and soprob_markov yield identical results - single patient", {
  withr::local_seed(1234)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 100,
    follow_up_time = follow_up,
    seed = 12345,
    treatment_effect = 0
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()

  t_covs <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  baseline_data <- data[data$time == 1, ]

  absorb <- 6
  ylevels <- factor(sort(unique(data$y)), ordered = FALSE)
  times <- 1:follow_up

  fit.b <- VGAM::vglm(
    y ~ rms::rcs(time, 3) *
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  fit.a <- VGAM::vglm(
    y ~
      (time_lin +
        time_nlin_1) *
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Test single patient (first patient)
  a <- Hmisc::soprobMarkovOrdm(
    fit.b,
    baseline_data[1, ],
    times = times,
    ylevels = ylevels,
    absorb = absorb
  )

  b <- soprob_markov(
    object = fit.a,
    data = baseline_data[1, ],
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = "time_lin",
    pvarname = "yprev",
    t_covs = t_covs
  )[1, , ]

  expect_equal(a, b, tolerance = 1e-15)
})

test_that("soprobMarkovOrdm and soprob_markov yield identical results - single patient; partial PO", {
  withr::local_seed(1234)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 100,
    follow_up_time = follow_up,
    seed = 12345,
    treatment_effect = 0
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()

  t_covs <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  baseline_data <- data[data$time == 1, ]

  absorb <- 6
  ylevels <- factor(sort(unique(data$y)), ordered = FALSE)
  times <- 1:follow_up

  fit.b <- VGAM::vglm(
    y ~ rms::rcs(time, 3) +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data
  )

  fit.a <- VGAM::vglm(
    y ~
      time_lin +
      time_nlin_1 +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data
  )

  # Test single patient (first patient)
  a <- Hmisc::soprobMarkovOrdm(
    fit.b,
    baseline_data[1, ],
    times = times,
    ylevels = ylevels,
    absorb = absorb
  )

  b <- soprob_markov(
    object = fit.a,
    data = baseline_data[1, ],
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = "time_lin",
    pvarname = "yprev",
    t_covs = t_covs
  )[1, , ]

  expect_equal(a, b, tolerance = 1e-15)
})

test_that("soprobMarkovOrdm and soprob_markov yield identical results - all patients", {
  withr::local_seed(7890)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 50, # Use smaller sample for full comparison
    follow_up_time = follow_up,
    seed = 12345,
    treatment_effect = 0
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()

  t_covs <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  baseline_data <- data[data$time == 1, ]
  ylevels <- factor(sort(unique(data$y)), ordered = FALSE)
  times <- 1:follow_up
  absorb <- 6

  fit.b <- VGAM::vglm(
    y ~ rms::rcs(time, 3) +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data
  )

  fit.a <- VGAM::vglm(
    y ~
      time_lin +
      time_nlin_1 +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx),
    data = data
  )

  # Test ALL patients
  # Get results from soprobMarkovOrdm (one patient at a time)
  results_hmisc <- array(
    dim = c(nrow(baseline_data), length(times), length(ylevels)),
    dimnames = list(NULL, NULL, levels(ylevels))
  )

  for (i in seq_len(nrow(baseline_data))) {
    results_hmisc[i, , ] <- Hmisc::soprobMarkovOrdm(
      fit.b,
      baseline_data[i, ],
      times = times,
      ylevels = ylevels,
      absorb = absorb
    )
  }

  # Get results from soprob_markov (all at once)
  results_new <- soprob_markov(
    object = fit.a,
    data = baseline_data,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = "time_lin",
    pvarname = "yprev",
    t_covs = t_covs
  )

  # Strip dimnames for fair comparison (soprob_markov adds them, manual loop doesn't)
  dimnames(results_new)[1:2] <- list(NULL, NULL)

  expect_equal(results_hmisc, results_new, tolerance = 1e-15)
})


test_that("Fast path yields same results as slow path for constrained PPO", {
  withr::local_seed(12345)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 25,
    follow_up_time = follow_up,
    seed = 12345
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()

  t_covs <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  absorb <- 6
  ylevels <- factor(sort(unique(data$y)), ordered = FALSE)
  times <- 1:follow_up

  # Define the constraint list
  cons <- list(
    "(Intercept)" = diag(5),

    # linear time: Intercept + linear trend across thresholds
    "time_lin" = cbind(
      time_lin_1 = 1,
      time_lin_5 = 0:4
    ),
    "time_nlin_1" = cbind(PO_effect = 1),
    "tx" = cbind(
      PO_effect = rep(1, 5),
      Deviation_5th = c(rep(0, 4), 1)
    ),
    "yprev" = cbind(PO_effect = 1),
    "age" = cbind(PO_effect = 1),
    "time_lin:tx" = cbind(PO_effect = 1),
    "time_nlin_1:tx" = cbind(PO_effect = 1),
    "time_lin:yprev" = cbind(PO_effect = 1),
    "tx:age" = cbind(PO_effect = 1)
  )

  fit2 <- VGAM::vglm(
    y ~ time_lin +
      time_nlin_1 +
      tx +
      yprev +
      age +
      time_lin:tx +
      time_nlin_1:tx +
      yprev:time_lin +
      tx:age,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ tx + time_lin),
    data = data,
    constraints = cons
  )

  # Setup comparison
  newdata <- data[data$time == 1, ] # Baseline rows

  # --- SLOW PATH (soprob_markov) ---
  res_slow <- soprob_markov(
    object = fit2,
    data = newdata,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = "time_lin",
    pvarname = "yprev",
    t_covs = t_covs
  )

  # --- FAST PATH (markov_msm_build + run) ---
  components <- markov_msm_build(
    model = fit2,
    data = newdata,
    t_covs = t_covs,
    ylevels = ylevels,
    pvarname = "yprev"
  )

  # Get Gamma (Effective coefficients)
  beta <- get_coef(fit2)
  C_list <- VGAM::constraints(fit2)
  Gamma_hat <- compute_Gamma(beta, C_list)

  # Run
  res_fast <- markov_msm_run(
    components = components,
    Gamma = Gamma_hat,
    times = times,
    absorb = absorb
  )

  # Fix dimnames for comparison
  dimnames(res_slow)[1:2] <- list(NULL)

  expect_equal(res_slow, res_fast, tolerance = 1e-15)
})


test_that("Fast path yields same results as slow path for PPO", {
  withr::local_seed(1234567)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 200, # Reduced from 500 for faster tests
    follow_up_time = follow_up,
    seed = 12345,
    treatment_effect = 0
  ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
    dplyr::ungroup()

  t_covs <- data |>
    dplyr::select(time_lin, time_nlin_1) |>
    dplyr::distinct() |>
    as.data.frame()

  fit <- VGAM::vglm(
    y ~ time_lin +
      time_nlin_1 +
      tx +
      yprev +
      age +
      time_lin:tx +
      time_nlin_1:tx +
      yprev:time_lin +
      tx:age,
    family = VGAM::cumulative(
      reverse = TRUE,
      parallel = FALSE ~ tx + time_lin
    ),
    data = data
  )

  # Setup comparison
  newdata <- data[data$time == 1, ] # Baseline rows
  times <- 1:follow_up
  ylevels <- factor(sort(unique(data$y)), ordered = FALSE)
  absorb <- 6

  # --- SLOW PATH (soprob_markov) ---
  res_slow <- soprob_markov(
    object = fit,
    data = newdata,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = "time_lin",
    pvarname = "yprev",
    t_covs = t_covs
  )

  # --- FAST PATH (markov_msm_build + run) ---
  components <- markov_msm_build(
    model = fit,
    data = newdata,
    t_covs = t_covs,
    ylevels = ylevels,
    pvarname = "yprev"
  )

  # Get Gamma (Effective coefficients)
  beta <- get_coef(fit)
  C_list <- VGAM::constraints(fit)
  Gamma_hat <- compute_Gamma(beta, C_list)

  # Run
  res_fast <- markov_msm_run(
    components = components,
    Gamma = Gamma_hat,
    times = times,
    absorb = absorb
  )

  # Fix dimnames for comparison
  dimnames(res_slow)[1:2] <- list(NULL)

  expect_equal(res_slow, res_fast, tolerance = 1e-10)
})
