test_that("soprobMarkovOrdm and soprob_markov yield identical results - single patient", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(1234)

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

  fit.b <- VGAM::vglm(
    y ~ rcs(time, 3) +
      tx +
      yprev +
      age,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  fit.a <- VGAM::vglm(
    y ~
      time_lin +
      time_nlin_1 +
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
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(1234)

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

  fit.b <- VGAM::vglm(
    y ~ rcs(time, 3) +
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
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(7890)

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

  fit.b <- VGAM::vglm(
    y ~ rcs(time, 3) +
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


test_that("Fast path (pre-calculated design matrix) yields same results as slow path (soprob_markov) for constrained PPO", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(12345)

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

  # 2. Fit Complex Model
  # Formula includes:
  # - Main effects: time, t_sq, tx, yprev
  # - Interactions: time:tx, t_sq:tx, yprev:tx, yprev:time

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
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # 2. Define the constraint list
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

  # 3. Setup comparison
  newdata <- data[data$time == 1, ] # Baseline rows
  times <- 1:follow_up
  ylevels <- factor(sort(unique(data$y)), ordered = FALSE)
  absorb <- 6

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
  # Access internal functions
  markov_msm_build <- markov.misc:::markov_msm_build
  markov_msm_run <- markov.misc:::markov_msm_run
  compute_Gamma <- markov.misc:::compute_Gamma

  components <- markov_msm_build(
    model = fit2,
    data = newdata,
    t_covs = t_covs,
    ylevels = ylevels,
    pvarname = "yprev"
  )

  # Get Gamma (Effective coefficients)
  beta <- markov.misc:::get_coef(fit2)
  C_list <- VGAM::constraints(fit2)
  Gamma_hat <- compute_Gamma(beta, C_list)

  # Run
  res_fast <- markov_msm_run(
    components = components,
    Gamma = Gamma_hat,
    times = times,
    absorb = absorb
  )

  # 4. Compare
  # res_slow and res_fast should be identical arrays [n_pat x n_times x n_states]

  # Fix dimnames for comparison
  # soprob_markov adds dimnames for patients and times, markov_msm_run only for states
  dimnames(res_slow)[1:2] <- list(NULL)

  expect_equal(res_slow, res_fast, tolerance = 1e-10)
})


test_that("Fast path (pre-calculated design matrix) yields same results as slow path for PPO", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  set.seed(1234567)

  follow_up <- 30

  data <- make_test_data(
    n_patients = 500,
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

  # 2. Fit Complex Model
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
      parallel = FALSE ~ tx + time_lin + age
    ),
    data = data
  )

  # 3. Setup comparison
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
  # Access internal functions
  markov_msm_build <- markov.misc:::markov_msm_build
  markov_msm_run <- markov.misc:::markov_msm_run
  compute_Gamma <- markov.misc:::compute_Gamma

  components <- markov_msm_build(
    model = fit,
    data = newdata,
    t_covs = t_covs,
    ylevels = ylevels,
    pvarname = "yprev"
  )

  # Get Gamma (Effective coefficients)
  beta <- markov.misc:::get_coef(fit)
  C_list <- VGAM::constraints(fit)
  Gamma_hat <- compute_Gamma(beta, C_list)

  # Run
  res_fast <- markov_msm_run(
    components = components,
    Gamma = Gamma_hat,
    times = times,
    absorb = absorb
  )

  # 4. Compare
  # res_slow and res_fast should be identical arrays [n_pat x n_times x n_states]

  # Fix dimnames for comparison
  # soprob_markov adds dimnames for patients and times, markov_msm_run only for states
  dimnames(res_slow)[1:2] <- list(NULL)

  expect_equal(res_slow, res_fast, tolerance = 1e-15)
})

# test_that("Fast path performance is significantly better than slow path", {
#   skip_on_cran()
#   skip_if_not_installed("VGAM")
#   skip_if_not_installed("rms")

#   # 1. Setup Data & Model (Identical to correctness test)
#   set.seed(1234567)
#   follow_up <- 30

#   # Generate slightly larger n to ensure the speed difference is measurable
#   # but not so large that the test suite hangs.
#   n_patients <- 50000

#   data <- make_test_data(
#     n_patients = n_patients,
#     follow_up_time = follow_up,
#     seed = 12345,
#     treatment_effect = 0
#   ) |>
#     dplyr::group_by(id) |>
#     dplyr::mutate(age = stats::rnorm(1, 50, 10)) |>
#     dplyr::ungroup()

#   t_covs <- data |>
#     dplyr::select(time_lin, time_nlin_1) |>
#     dplyr::distinct() |>
#     as.data.frame()

#   fit <- VGAM::vglm(
#     y ~ time_lin +
#       time_nlin_1 +
#       tx +
#       yprev +
#       age +
#       time_lin:tx +
#       time_nlin_1:tx +
#       yprev:time_lin +
#       tx:age,
#     family = VGAM::cumulative(
#       reverse = TRUE,
#       parallel = TRUE
#     ),
#     data = data
#   )

#   # 2. Setup Comparison Inputs
#   newdata <- data[data$time == 1, ]
#   times <- 1:follow_up
#   ylevels <- factor(sort(unique(data$y)), ordered = FALSE)
#   absorb <- 6

#   # --- MEASURE SLOW PATH ---
#   # We measure the entire execution of the legacy function
#   time_slow <- system.time({
#     res_slow <- soprob_markov(
#       object = fit,
#       data = newdata,
#       times = times,
#       ylevels = ylevels,
#       absorb = absorb,
#       tvarname = "time_lin",
#       pvarname = "yprev",
#       t_covs = t_covs
#     )
#   })[["elapsed"]]

#   # --- MEASURE FAST PATH ---
#   # We measure the TOTAL time (Build + Setup + Run) to be fair.
#   # If the 'Build' phase was excluded, the speedup would be even higher (useful for bootstrapping),
#   # but checking the total time confirms the implementation is efficient end-to-end.

#   markov_msm_build <- markov.misc:::markov_msm_build
#   markov_msm_run <- markov.misc:::markov_msm_run
#   compute_Gamma <- markov.misc:::compute_Gamma
#   get_coef <- markov.misc:::get_coef

#   time_fast <- system.time({
#     # A. Build components
#     components <- markov_msm_build(
#       model = fit,
#       data = newdata,
#       t_covs = t_covs,
#       ylevels = ylevels,
#       pvarname = "yprev"
#     )

#     # B. Calculate Gamma
#     beta <- get_coef(fit)
#     C_list <- VGAM::constraints(fit)
#     Gamma_hat <- compute_Gamma(beta, C_list)

#     # C. Run
#     res_fast <- markov_msm_run(
#       components = components,
#       Gamma = Gamma_hat,
#       times = times,
#       absorb = absorb
#     )
#   })[["elapsed"]]

#   # 3. Assertions

#   # Calculate speedup factor
#   speedup <- time_slow / time_fast

#   # Print results to the test log for visibility
#   message(sprintf(
#     "\nPerformance check:\nSlow path: %.3fs\nFast path: %.3fs\nSpeedup: %.1fx",
#     time_slow,
#     time_fast,
#     speedup
#   ))

#   # Assert that fast path is indeed faster
#   expect_lt(time_fast, time_slow)

#   # Assert that it is *significantly* faster (e.g., at least 2x faster)
#   # Adjust tolerance based on expected performance gains.
#   # If the fast path is 100x faster, a tolerance of 2 is very safe.
#   expect_true(
#     speedup > 10,
#     label = "Fast path should be at least 2x faster than slow path"
#   )
# })
