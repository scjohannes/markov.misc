# Unit Tests: State Occupation Probabilities (SOPs)
#
# These tests verify that:
# 1. soprob_markov() gives same results as Hmisc::soprobMarkovOrdm()
# 2. soprob_markov() gives same results with rcs() vs precomputed basis functions
# 3. Partial proportional odds models work correctly

test_that("soprob_markov matches Hmisc::soprobMarkovOrdm for simple PO model", {
  skip_if_not_installed("Hmisc")
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 100, seed = 1234, follow_up_time = 30)

  # Fit simple proportional odds model
  m1 <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Calculate SOPs using Hmisc
  sops_hmisc <- Hmisc::soprobMarkovOrdm(
    m1,
    data = data[1, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6"
  )

  # Calculate SOPs using markov.misc
  sops_markov <- soprob_markov(
    m1,
    data = data[1, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6"
  )

  # Check dimensions match
  expect_equal(dim(sops_hmisc), dim(sops_markov[1, , ]))

  # Check values are equal
  expect_equal(sops_hmisc, sops_markov[1, , ], tolerance = 1e-10)
})

test_that("soprob_markov works with rcs() and precomputed spline basis", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 100, seed = 234, follow_up_time = 30)

  # Create spline basis with rcs()
  time_spl_m <- rms::rcs(data$time, 4)
  knots <- attr(time_spl_m, "parms")

  # Extract individual spline terms
  data$time_lin <- as.vector(time_spl_m[, 1])
  data$time_nlin_1 <- as.vector(time_spl_m[, 2])
  data$time_nlin_2 <- as.vector(time_spl_m[, 3])

  # Create time covariates data frame
  t_covs <- data.frame(
    time_lin = data$time_lin,
    time_nlin_1 = data$time_nlin_1,
    time_nlin_2 = data$time_nlin_2
  )
  t_covs <- unique(t_covs)
  t_covs <- t_covs[order(t_covs$time_lin), , drop = FALSE]
  rownames(t_covs) <- NULL

  # Fit model with rcs() in formula
  m_rcs <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Fit model with precomputed spline terms
  m_precomp <- VGAM::vglm(
    ordered(y) ~ time_lin + time_nlin_1 + time_nlin_2 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Calculate SOPs with rcs() model
  sops_rcs <- soprob_markov(
    m_rcs,
    data = data[1, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6"
  )

  # Calculate SOPs with precomputed model (requires t_covs)
  sops_precomp <- soprob_markov(
    m_precomp,
    data = data[1, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  # Check dimensions match
  expect_equal(dim(sops_rcs), dim(sops_precomp))

  # Check values are equal
  expect_equal(sops_rcs[1, , ], sops_precomp[1, , ], tolerance = 1e-10)

  # Test with multiple patients
  sops_rcs_multi <- soprob_markov(
    m_rcs,
    data = data[1:10, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6"
  )

  sops_precomp_multi <- soprob_markov(
    m_precomp,
    data = data[1:10, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  # Check all patients match
  for (i in 1:10) {
    expect_equal(
      sops_rcs_multi[i, , ],
      sops_precomp_multi[i, , ],
      tolerance = 1e-10,
      label = paste("Patient", i)
    )
  }
})

test_that("soprob_markov handles partial proportional odds models", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 100, seed = 456, follow_up_time = 30)

  # Create spline basis
  time_spl_m <- rms::rcs(data$time, 4)
  data$time_lin <- as.vector(time_spl_m[, 1])
  data$time_nlin_1 <- as.vector(time_spl_m[, 2])
  data$time_nlin_2 <- as.vector(time_spl_m[, 3])

  t_covs <- data.frame(
    time_lin = data$time_lin,
    time_nlin_1 = data$time_nlin_1,
    time_nlin_2 = data$time_nlin_2
  )
  t_covs <- unique(t_covs)
  t_covs <- t_covs[order(t_covs$time_lin), , drop = FALSE]
  rownames(t_covs) <- NULL

  # Define constraints for partial PO model
  M <- 5
  cons_time <- list(
    "(Intercept)" = diag(M),
    "yprev" = cbind(PO_effect = 1),
    "tx" = cbind(PO_effect = 1),
    "time_lin" = cbind(
      time_lin_1 = 1,
      time_lin_5 = 0:(M - 1)
    ),
    "time_nlin_1" = cbind(PO_effect = 1),
    "time_nlin_2" = cbind(PO_effect = 1)
  )

  # Fit partial PO model
  m_ppo <- VGAM::vglm(
    ordered(y) ~ time_lin + time_nlin_1 + time_nlin_2 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ time_lin),
    data = data,
    constraints = cons_time
  )

  # Bake constraint matrix into model object
  m_ppo@call$constraints <- cons_time

  # Calculate SOPs
  sops_ppo <- soprob_markov(
    m_ppo,
    data = data[1, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  # Check output structure
  expect_equal(dim(sops_ppo), c(1, 30, 6))
  expect_type(sops_ppo, "double")

  # Check probabilities sum to 1 at each time point
  for (t in 1:30) {
    expect_equal(sum(sops_ppo[1, t, ]), 1, tolerance = 1e-10)
  }

  # Check probabilities are between 0 and 1
  expect_true(all(sops_ppo >= 0 & sops_ppo <= 1))

  # Check absorbing state (state 6) is monotonically increasing
  state_6_probs <- sops_ppo[1, , 6]
  expect_true(all(diff(state_6_probs) >= -1e-10)) # Allow tiny numerical errors
})

test_that("soprob_markov handles interaction terms", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 100, seed = 789, follow_up_time = 30)

  # Create spline basis
  time_spl_m <- rms::rcs(data$time, 4)
  data$time_lin <- as.vector(time_spl_m[, 1])
  data$time_nlin_1 <- as.vector(time_spl_m[, 2])
  data$time_nlin_2 <- as.vector(time_spl_m[, 3])

  t_covs <- data.frame(
    time_lin = data$time_lin,
    time_nlin_1 = data$time_nlin_1,
    time_nlin_2 = data$time_nlin_2
  )
  t_covs <- unique(t_covs)
  t_covs <- t_covs[order(t_covs$time_lin), , drop = FALSE]
  rownames(t_covs) <- NULL

  # Define constraints for partial PO model with interactions
  M <- 5
  cons_time_tx <- list(
    "(Intercept)" = diag(M),
    "yprev" = cbind(PO_effect = 1),
    "tx" = cbind(PO_effect = 1),
    "time_lin" = cbind(
      time_lin_1 = 1,
      time_lin_5 = 0:(M - 1)
    ),
    "time_nlin_1" = cbind(PO_effect = 1),
    "time_nlin_2" = cbind(PO_effect = 1),
    "time_lin:tx" = cbind(PO_effect = 1),
    "time_nlin_1:tx" = cbind(PO_effect = 1),
    "time_nlin_2:tx" = cbind(PO_effect = 1),
    "time_lin:yprev" = cbind(PO_effect = 1)
  )

  # Fit model with interactions
  m_interact <- VGAM::vglm(
    ordered(y) ~ (time_lin + time_nlin_1 + time_nlin_2) *
      tx +
      yprev +
      yprev:time_lin,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ time_lin),
    data = data,
    constraints = cons_time_tx
  )

  # Bake constraint matrix into model object
  m_interact@call$constraints <- cons_time_tx

  # Calculate SOPs for treatment group
  data_tx <- data[data$time == 1 & data$tx == 1, , drop = FALSE]
  sops_tx <- soprob_markov(
    m_interact,
    data = data_tx[1, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  # Calculate SOPs for control group
  data_ctrl <- data[data$time == 1 & data$tx == 0, , drop = FALSE]
  sops_ctrl <- soprob_markov(
    m_interact,
    data = data_ctrl[1, ],
    time = 1:30,
    ylevels = factor(1:6),
    absorb = "6",
    t_covs = t_covs
  )

  # Check output structure
  expect_equal(dim(sops_tx), c(1, 30, 6))
  expect_equal(dim(sops_ctrl), c(1, 30, 6))

  # Check probabilities sum to 1 at each time point
  for (t in 1:30) {
    expect_equal(sum(sops_tx[1, t, ]), 1, tolerance = 1e-10)
    expect_equal(sum(sops_ctrl[1, t, ]), 1, tolerance = 1e-10)
  }

  # Check that treatment and control give different results (due to interaction)
  # At least some time points should differ
  differences <- abs(sops_tx[1, , ] - sops_ctrl[1, , ])
  expect_true(any(differences > 1e-6))
})

test_that("soprob_markov respects absorbing state constraint", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 50, seed = 999, follow_up_time = 20)

  # Fit simple model
  m1 <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Calculate SOPs
  sops <- soprob_markov(
    m1,
    data = data[1:10, ],
    time = 1:20,
    ylevels = factor(1:6),
    absorb = "6"
  )

  # Check that state 6 probabilities are monotonically increasing for all patients
  for (i in 1:10) {
    state_6_probs <- sops[i, , 6]
    diffs <- diff(state_6_probs)
    # Allow tiny negative values due to numerical precision
    expect_true(all(diffs >= -1e-10), label = paste("Patient", i))
  }

  # Check that probabilities sum to 1 at each time point for all patients
  for (i in 1:10) {
    for (t in 1:20) {
      prob_sum <- sum(sops[i, t, ])
      expect_equal(
        prob_sum,
        1,
        tolerance = 1e-10,
        label = paste("Patient", i, "Time", t)
      )
    }
  }
})

test_that("soprob_markov works without an absorbing state", {
  skip_if_not_installed("VGAM")

  data <- as.data.frame(sim_trajectories_brownian(
    n_patients = 120,
    follow_up_time = 12,
    treatment_prob = 0.5,
    absorbing_state = NULL,
    allowed_start_state = as.integer(1:6),
    thresholds = c(-2, -1, 0, 1, 2),
    x0_sd = 2,
    seed = 1000
  ))
  data$yprev <- factor(data$yprev, levels = 1:6)

  expect_true(6 %in% as.integer(as.character(data$yprev)))

  baseline_data <- data[data$time == 1, , drop = FALSE]

  model <- VGAM::vglm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  result <- soprob_markov(
    model,
    data = baseline_data[1:5, , drop = FALSE],
    time = 1:12,
    ylevels = factor(1:6),
    absorb = NULL
  )

  expect_equal(dim(result), c(5L, 12L, 6L))
  expect_true(all(is.finite(result)))
  expect_true(all(result >= 0))

  for (i in 1:5) {
    for (t in 1:12) {
      expect_equal(
        sum(result[i, t, ]),
        1,
        tolerance = 1e-10,
        label = paste("Patient", i, "Time", t)
      )
    }
  }
})

# Unit Tests: Equivalence of sops() and avg_sops()
#
# This test verifies that manual marginalization using sops() matches
# G-computation performed by avg_sops().

test_that("marginalization via sops() and avg_sops() result in exactly the same results", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  withr::local_seed(123)

  # 1. Setup Data and Model
  FU <- 10
  data <- make_test_data(n_patients = 50, follow_up_time = FU, seed = 123)

  # Extract t_covs from the data (already contains spline terms)
  t_covs <- data.frame(
    time = data$time,
    time_lin = as.vector(data$time_lin),
    time_nlin_1 = as.vector(data$time_nlin_1)
  )
  t_covs <- unique(t_covs)
  t_covs <- t_covs[order(t_covs$time), , drop = FALSE]
  rownames(t_covs) <- NULL

  # Fit an interaction model
  m_inter <- VGAM::vglm(
    ordered(y) ~ (time_lin + time_nlin_1) * tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  # Use robust model as in the example
  m_robust <- robcov_vglm(m_inter, cluster = data$id)

  # 2. Manual marginalization via sops()

  # Filter to baseline data for G-computation
  baseline_data <- data[data$time == 1, , drop = FALSE]

  baseline_tx1 <- baseline_data
  baseline_tx1$tx <- 1
  sops_1 <- sops(
    m_robust,
    newdata = baseline_tx1,
    times = 1:FU,
    ylevels = 1:6,
    absorb = "6",
    id_var = "id",
    t_covs = t_covs
  )
  sops_1 <- aggregate(estimate ~ time + state, data = sops_1, FUN = mean)
  sops_1$tx <- 1

  baseline_tx0 <- baseline_data
  baseline_tx0$tx <- 0
  sops_0 <- sops(
    m_robust,
    newdata = baseline_tx0,
    times = 1:FU,
    ylevels = 1:6,
    absorb = "6",
    id_var = "id",
    t_covs = t_covs
  )
  sops_0 <- aggregate(estimate ~ time + state, data = sops_0, FUN = mean)
  sops_0$tx <- 0

  sops_manual <- rbind(sops_0, sops_1)
  sops_manual <- sops_manual[order(sops_manual$tx, sops_manual$state, sops_manual$time), ]
  rownames(sops_manual) <- NULL

  # 3. Automated marginalization via avg_sops()
  sops_avg <- avg_sops(
    m_robust,
    newdata = baseline_data,
    variables = "tx",
    times = 1:FU,
    ylevels = 1:6,
    absorb = "6",
    id_var = "id",
    t_covs = t_covs
  )
  sops_avg <- as.data.frame(sops_avg)
  sops_avg <- sops_avg[order(sops_avg$tx, sops_avg$state, sops_avg$time), ]
  rownames(sops_avg) <- NULL

  # 4. Compare results
  expect_equal(sops_manual$estimate, sops_avg$estimate, tolerance = 1e-15)
  expect_equal(sops_manual$tx, sops_avg$tx)
  expect_equal(as.character(sops_manual$state), as.character(sops_avg$state))
  expect_equal(sops_manual$time, sops_avg$time)
})
