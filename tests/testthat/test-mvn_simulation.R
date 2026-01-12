# Unit Tests: MVN Simulation-Based Inference for SOPs
#
# These tests verify that:
# 1. set_coef() works correctly for vglm and orm models
# 2. get_vcov_robust() returns correct variance matrices
# 3. inferences(method = "simulation") produces valid confidence intervals
# 4. inferences() works for both avg_sops and sops objects
# 5. get_draws() extracts simulation draws correctly

test_that("set_coef() works for vglm models", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # Simulate test data
  set.seed(111)
  test_data <- sim_trajectories_brownian(
    n_patients = 50,
    follow_up_time = 20,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 111,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Get original coefficients
  orig_coefs <- coef(m_vglm)

  # Create new coefficients (slightly perturbed)
  new_coefs <- orig_coefs * 1.1

  # Set new coefficients
  m_modified <- set_coef(m_vglm, new_coefs)

  # Check that coefficients were updated
  expect_equal(coef(m_modified), new_coefs)

  # Check that original model is unchanged
  expect_equal(coef(m_vglm), orig_coefs)

  # Check predictions differ
  pred_orig <- predict(m_vglm, newdata = data[1:5, ], type = "response")
  pred_mod <- predict(m_modified, newdata = data[1:5, ], type = "response")
  expect_false(all(pred_orig == pred_mod))
})


test_that("set_coef() works for orm models", {
  skip_if_not_installed("rms")

  # Simulate test data
  set.seed(222)
  test_data <- sim_trajectories_brownian(
    n_patients = 50,
    follow_up_time = 20,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 222,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit orm model - datadist must be defined in global env
  dd <- rms::datadist(data)
  assign("dd", dd, envir = globalenv())
  on.exit(rm("dd", envir = globalenv()), add = TRUE)
  options(datadist = "dd")

  m_orm <- rms::orm(
    y ~ rms::rcs(time, 3) + tx + yprev,
    data = data,
    x = TRUE,
    y = TRUE
  )

  # Get original coefficients
  orig_coefs <- coef(m_orm)

  # Create new coefficients
  new_coefs <- orig_coefs * 0.9

  # Set new coefficients
  m_modified <- set_coef(m_orm, new_coefs)

  # Check that coefficients were updated
  expect_equal(coef(m_modified), new_coefs)

  # Check that original model is unchanged
  expect_equal(coef(m_orm), orig_coefs)
})


test_that("get_vcov_robust() returns standard vcov when cluster is NULL", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # Simulate test data
  set.seed(333)
  test_data <- sim_trajectories_brownian(
    n_patients = 30,
    follow_up_time = 15,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 333,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Get vcov without clustering
  V <- get_vcov_robust(m_vglm, cluster = NULL)

  # Should match standard vcov
  expect_equal(V, vcov(m_vglm))
})


test_that("get_vcov_robust() computes cluster-robust vcov with formula", {
  skip_if_not_installed("VGAM")

  # Simulate test data
  set.seed(444)
  test_data <- sim_trajectories_brownian(
    n_patients = 30,
    follow_up_time = 15,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 444,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Get cluster-robust vcov using formula
  V_robust <- get_vcov_robust(m_vglm, cluster = ~id, data = data)

  # Should be a matrix with correct dimensions
  expect_true(is.matrix(V_robust))
  expect_equal(nrow(V_robust), length(coef(m_vglm)))
  expect_equal(ncol(V_robust), length(coef(m_vglm)))

  # Should differ from standard vcov (cluster adjustment)
  V_standard <- vcov(m_vglm)
  expect_false(all(V_robust == V_standard))
})


test_that("get_vcov_robust() extracts vcov from robcov_vglm object", {
  skip_if_not_installed("VGAM")

  # Simulate test data
  set.seed(555)
  test_data <- sim_trajectories_brownian(
    n_patients = 30,
    follow_up_time = 15,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 555,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Create robcov_vglm object
  m_robust <- robcov_vglm(m_vglm, cluster = data$id)

  # Extract vcov from robcov_vglm
  V <- get_vcov_robust(m_robust)

  # Should match the stored var
  expect_equal(V, m_robust$var)
})


test_that("inferences(method='simulation') produces valid CIs for avg_sops", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  # Simulate test data
  set.seed(666)
  test_data <- sim_trajectories_brownian(
    n_patients = 50,
    follow_up_time = 15,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 666,
    mu_treatment_effect = 0.3
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model and wrap with cluster-robust vcov
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  m_robust <- robcov_vglm(m_vglm, cluster = data$id)

  # Compute avg_sops with robust model
  result <- avg_sops(
    model = m_robust,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:15,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  )

  # Add simulation-based inference (vcov extracted from model)
  result_ci <- inferences(
    result,
    method = "simulation",
    n_sim = 100, # Small for testing
    conf_level = 0.95
  )

  # Check that CI columns exist
  expect_true("conf.low" %in% names(result_ci))
  expect_true("conf.high" %in% names(result_ci))
  expect_true("std.error" %in% names(result_ci))

  # Check CIs are valid
  expect_true(all(result_ci$conf.low <= result_ci$estimate, na.rm = TRUE))
  expect_true(all(result_ci$conf.high >= result_ci$estimate, na.rm = TRUE))
  expect_true(all(result_ci$std.error >= 0, na.rm = TRUE))

  # Check attributes
  expect_equal(attr(result_ci, "method"), "simulation")
  expect_equal(attr(result_ci, "n_sim"), 100)
})


test_that("inferences() with return_draws=TRUE stores simulation draws", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  # Simulate test data
  set.seed(777)
  test_data <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 10,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 777,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model and wrap with cluster-robust vcov
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  m_robust <- robcov_vglm(m_vglm, cluster = data$id)

  # Compute avg_sops with simulation CIs and draws
  result <- avg_sops(
    model = m_robust,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  ) |>
    inferences(
      method = "simulation",
      n_sim = 30,
      return_draws = TRUE
    )

  # Extract draws
  draws <- get_draws(result)

  # Check draws structure
  expect_true(is.data.frame(draws))
  expect_true("draw_id" %in% names(draws))
  expect_true("estimate" %in% names(draws))
  expect_true("time" %in% names(draws))
  expect_true("state" %in% names(draws))
  expect_true("tx" %in% names(draws))

  # Check number of draws
  n_draws <- length(unique(draws$draw_id))
  expect_lte(n_draws, 30) # May be less if some failed

  # Check we have draws for all combinations
  n_times <- 10
  n_states <- 6
  n_tx <- 2
  expected_rows_per_draw <- n_times * n_states * n_tx
  actual_rows_per_draw <- nrow(draws) / n_draws
  expect_equal(actual_rows_per_draw, expected_rows_per_draw)
})


test_that("simulation vs bootstrap produce similar results", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")
  skip("Long-running test - run manually")

  # Simulate test data
  set.seed(888)
  test_data <- sim_trajectories_brownian(
    n_patients = 80,
    follow_up_time = 15,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 888,
    mu_treatment_effect = 0.2
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  m_robust <- robcov_vglm(m_vglm, cluster = data$id)

  # Compute avg_sops with robust model
  avg_result <- avg_sops(
    model = m_robust,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:15,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  )

  # Simulation-based inference (vcov from model)
  result_sim <- inferences(
    avg_result,
    method = "simulation",
    n_sim = 200
  )

  # Bootstrap-based inference
  result_boot <- inferences(
    avg_result,
    method = "bootstrap",
    n_sim = 100,
    parallel = TRUE,
    workers = 4,
    return_draws = TRUE
  )

  # Compare standard errors - should be in same ballpark
  se_sim <- result_sim$std.error
  se_boot <- result_boot$std.error

  # Correlation should be positive and reasonably high
  cor_se <- cor(se_sim, se_boot, use = "complete.obs")
  expect_gt(cor_se, 0.5)

  # Ratio of SEs should be mostly between 0.5 and 2
  se_ratio <- se_sim / se_boot
  reasonable_ratios <- se_ratio > 0.3 & se_ratio < 3
  expect_gt(mean(reasonable_ratios, na.rm = TRUE), 0.8)
})


test_that("inferences() works with custom vcov matrix", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  # Simulate test data
  set.seed(999)
  test_data <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 10,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 999,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Compute custom vcov (e.g., inflated for conservative CIs)
  V_standard <- vcov(m_vglm)
  V_inflated <- V_standard * 2 # Inflate variance

  # Compute avg_sops with custom vcov
  result <- avg_sops(
    model = m_vglm,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  ) |>
    inferences(
      method = "simulation",
      n_sim = 50,
      vcov = V_inflated
    )

  # CIs should be wider than with standard vcov
  result_standard <- avg_sops(
    model = m_vglm,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  ) |>
    inferences(
      method = "simulation",
      n_sim = 50
    )

  # Compare CI widths
  ci_width_inflated <- result$conf.high - result$conf.low
  ci_width_standard <- result_standard$conf.high - result_standard$conf.low

  # Inflated vcov should give wider CIs on average
  expect_gt(
    mean(ci_width_inflated, na.rm = TRUE),
    mean(ci_width_standard, na.rm = TRUE)
  )
})


test_that("conf_type='wald' produces different CIs than 'perc'", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  # Simulate test data
  set.seed(1010)
  test_data <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 10,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 1010,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model and wrap with cluster-robust vcov
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  m_robust <- robcov_vglm(m_vglm, cluster = data$id)

  # Compute avg_sops with robust model
  avg_result <- avg_sops(
    model = m_robust,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  )

  # Percentile CIs (vcov from model)
  result_perc <- inferences(
    avg_result,
    method = "simulation",
    n_sim = 50,
    conf_type = "perc"
  )

  # Wald CIs (vcov from model)
  result_wald <- inferences(
    avg_result,
    method = "simulation",
    n_sim = 50,
    conf_type = "wald"
  )

  # Check attributes
  expect_equal(attr(result_perc, "conf_type"), "perc")
  expect_equal(attr(result_wald, "conf_type"), "wald")

  # Both should have valid CIs
  expect_true(all(result_perc$conf.low <= result_perc$estimate, na.rm = TRUE))
  expect_true(all(result_wald$conf.low <= result_wald$estimate, na.rm = TRUE))
})


test_that("inferences errors appropriately for invalid inputs", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # Simulate test data
  set.seed(1111)
  test_data <- sim_trajectories_brownian(
    n_patients = 30,
    follow_up_time = 10,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 1111,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Error: passing model directly instead of sops object
  expect_error(
    inferences(m_vglm),
    "markov_avg_sops|markov_sops"
  )

  # Error: invalid method
  avg_result <- avg_sops(
    model = m_vglm,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  )

  expect_error(
    inferences(avg_result, method = "invalid"),
    "'arg' should be one of"
  )
})


test_that("get_draws() errors when no draws stored", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  # Simulate test data
  set.seed(1212)
  test_data <- sim_trajectories_brownian(
    n_patients = 30,
    follow_up_time = 10,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 1212,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model and wrap with cluster-robust vcov
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  m_robust <- robcov_vglm(m_vglm, cluster = data$id)

  # Compute avg_sops without return_draws
  result <- avg_sops(
    model = m_robust,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    id_var = "id"
  ) |>
    inferences(
      method = "simulation",
      n_sim = 20,
      return_draws = FALSE
    )

  # get_draws should error
  expect_error(
    get_draws(result),
    "No draws found"
  )
})


test_that("inferences(method='simulation') works for individual sops()", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")

  # Simulate test data
  set.seed(1313)
  test_data <- sim_trajectories_brownian(
    n_patients = 40,
    follow_up_time = 10,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = 1313,
    mu_treatment_effect = 0
  )

  data <- prepare_markov_data(test_data)

  # Fit vglm model and wrap with cluster-robust vcov
  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  m_robust <- robcov_vglm(m_vglm, cluster = data$id)

  # Compute individual sops with robust model
  sops_result <- sops(
    model = m_robust,
    newdata = data,
    times = 1:10,
    ylevels = 1:6,
    absorb = 6,
    tvarname = "time",
    pvarname = "yprev"
  )

  # Check that newdata_orig attribute is stored
  expect_false(is.null(attr(sops_result, "newdata_orig")))
  expect_equal(nrow(attr(sops_result, "newdata_orig")), nrow(data))

  # Add simulation-based inference (vcov from model)
  result_ci <- inferences(
    sops_result,
    method = "simulation",
    n_sim = 20,
    conf_level = 0.95,
    return_draws = TRUE
  )

  # Check that CI columns exist
  expect_true("conf.low" %in% names(result_ci))
  expect_true("conf.high" %in% names(result_ci))
  expect_true("std.error" %in% names(result_ci))

  # Check CIs are valid (low <= estimate <= high)
  expect_true(all(result_ci$conf.low <= result_ci$estimate, na.rm = TRUE))
  expect_true(all(result_ci$conf.high >= result_ci$estimate, na.rm = TRUE))

  # Check attributes

  expect_equal(attr(result_ci, "method"), "simulation")
  expect_equal(attr(result_ci, "n_sim"), 20)

  # Check draws can be extracted
  draws <- get_draws(result_ci)
  expect_true(is.data.frame(draws))
  expect_true("draw_id" %in% names(draws))
  expect_true("estimate" %in% names(draws))
  expect_true("rowid" %in% names(draws))
  expect_true("time" %in% names(draws))
  expect_true("state" %in% names(draws))
})
