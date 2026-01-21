# Unit Tests: Robust Covariance Estimation for vglm
#
# These tests verify that:
# 1. robcov_vglm produces valid sandwich covariance estimates
# 2. Results match rms::robcov for equivalent models
# 3. Clustering works correctly
# 4. Score computations have correct properties
# 5. Model information is preserved
# 6. Z-statistics and p-values are computed correctly

test_that("robcov_vglm produces valid covariance estimates", {
  skip_if_not_installed("VGAM")

  # Generate simple test data
  set.seed(123)
  n <- 200
  x <- rnorm(n)
  lp <- 0.5 * x
  probs <- cbind(
    plogis(-1 - lp),
    plogis(0 - lp) - plogis(-1 - lp),
    plogis(1 - lp) - plogis(0 - lp),
    1 - plogis(1 - lp)
  )
  y <- apply(probs, 1, function(p) sample(1:4, 1, prob = pmax(p, 0.001)))

  test_data <- data.frame(y = ordered(y), x = x)

  # Fit vglm
  m <- VGAM::vglm(
    y ~ x,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = test_data
  )

  # Compute robust covariance
  result <- robcov_vglm(m)

  # Check output structure
  expect_s3_class(result, "robcov_vglm")

  # Check variance and SE components
  expect_true("var" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true("scores" %in% names(result))
  expect_true("bread" %in% names(result))
  expect_true("meat" %in% names(result))

  # Check new components

  expect_true("coefficients" %in% names(result))
  expect_true("z" %in% names(result))
  expect_true("pvalues" %in% names(result))
  expect_true("original_se" %in% names(result))

  # Check model information components
  expect_true("family" %in% names(result))
  expect_true("constraints" %in% names(result))
  expect_true("control" %in% names(result))
  expect_true("call" %in% names(result))
  expect_true("original_call" %in% names(result))

  # Check dimensions
  p <- length(coef(m))
  expect_equal(dim(result$var), c(p, p))
  expect_length(result$se, p)
  expect_length(result$coefficients, p)
  expect_length(result$z, p)
  expect_length(result$pvalues, p)
  expect_equal(dim(result$scores), c(n, p))
  expect_equal(dim(result$bread), c(p, p))
  expect_equal(dim(result$meat), c(p, p))

  # Check variance matrix is symmetric
  expect_equal(result$var, t(result$var), tolerance = 1e-10)

  # Check variance matrix is positive semi-definite
  eigenvalues <- eigen(result$var, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))

  # Check standard errors are positive
  expect_true(all(result$se > 0))

  # Check n matches
  expect_equal(result$n, n)
})


test_that("robcov_vglm computes z-statistics and p-values correctly", {
  skip_if_not_installed("VGAM")

  # Generate simple test data with known effect
  set.seed(456)
  n <- 500
  x <- rnorm(n)
  lp <- -0.5 + 1.5 * x # Strong effect
  y <- rbinom(n, 1, plogis(lp))

  test_data <- data.frame(y = y, x = x)

  # Fit vglm
  m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

  # Compute robust covariance
  result <- robcov_vglm(m)

  # Check z = coef / se
  expected_z <- result$coefficients / result$se
  expect_equal(result$z, expected_z, tolerance = 1e-10)

  # Check p-values = 2 * pnorm(-abs(z))
  expected_p <- 2 * pnorm(-abs(result$z))
  expect_equal(result$pvalues, expected_p, tolerance = 1e-10)

  # Check that strong effect has small p-value
  expect_lt(result$pvalues["x"], 0.05)

  # Check all p-values are in [0, 1]
  expect_true(all(result$pvalues >= 0 & result$pvalues <= 1))
})


test_that("robcov_vglm matches rms::robcov for binary outcome", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # Generate binary outcome data
  set.seed(456)
  n <- 500
  x <- rnorm(n)
  lp <- -0.5 + 0.8 * x
  y <- rbinom(n, 1, plogis(lp))

  test_data <- data.frame(y = y, x = x)

  # Fit orm
  dd <- rms::datadist(test_data)
  old_dd <- options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  on.exit({
    options(old_dd)
    rm("dd", envir = globalenv())
  }, add = TRUE)

  m_orm <- rms::orm(y ~ x, data = test_data, x = TRUE, y = TRUE)

  # Fit vglm with binomialff (standard logistic)
  m_vglm <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

  # Compute robust covariance (unclustered)
  orm_robust <- rms::robcov(m_orm)
  vglm_robust <- robcov_vglm(m_vglm)

  # Compare SEs (should match to numerical precision)
  se_orm <- sqrt(diag(vcov(orm_robust)))
  se_vglm <- vglm_robust$se

  # Use unname() because names differ between packages
  expect_equal(unname(se_vglm), unname(se_orm), tolerance = 1e-5)
})


test_that("robcov_vglm matches rms::robcov for binary outcome with clustering", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # Generate clustered binary outcome data
  set.seed(789)
  n_clusters <- 50
  n_per_cluster <- 10
  n <- n_clusters * n_per_cluster
  cluster <- rep(1:n_clusters, each = n_per_cluster)

  x <- rnorm(n)
  cluster_effect <- rep(rnorm(n_clusters, sd = 0.3), each = n_per_cluster)
  lp <- -0.5 + 0.8 * x + cluster_effect
  y <- rbinom(n, 1, plogis(lp))

  test_data <- data.frame(y = y, x = x, cluster = cluster)

  # Fit orm
  dd <- rms::datadist(test_data)
  old_dd <- options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  on.exit({
    options(old_dd)
    rm("dd", envir = globalenv())
  }, add = TRUE)

  m_orm <- rms::orm(y ~ x, data = test_data, x = TRUE, y = TRUE)

  # Fit vglm
  m_vglm <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

  # Compute cluster-robust covariance (no small-sample adjustment to match rms)
  orm_robust <- rms::robcov(m_orm, cluster = test_data$cluster)
  vglm_robust <- robcov_vglm(
    m_vglm,
    cluster = test_data$cluster,
    adjust = FALSE
  )

  # Compare SEs
  se_orm <- sqrt(diag(vcov(orm_robust)))
  se_vglm <- vglm_robust$se

  # Use unname() because names differ between packages
  expect_equal(unname(se_vglm), unname(se_orm), tolerance = 1e-6)

  # Check number of clusters is correct
  expect_equal(vglm_robust$n_clusters, n_clusters)
})


test_that("robcov_vglm gives similar results to rms::robcov for ordinal outcome", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # Generate ordinal outcome data
  set.seed(101)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  lp <- 0.5 * x1 + 0.3 * x2
  cuts <- c(-1, 0, 1)
  probs <- matrix(NA, nrow = n, ncol = 4)
  for (j in 1:4) {
    if (j == 1) {
      probs[, j] <- plogis(cuts[1] - lp)
    } else if (j == 4) {
      probs[, j] <- 1 - plogis(cuts[3] - lp)
    } else {
      probs[, j] <- plogis(cuts[j] - lp) - plogis(cuts[j - 1] - lp)
    }
  }
  y <- apply(probs, 1, function(p) sample(1:4, 1, prob = pmax(p, 0.001)))

  test_data <- data.frame(y = ordered(y), x1 = x1, x2 = x2)

  # Fit orm
  dd <- rms::datadist(test_data)
  old_dd <- options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  on.exit({
    options(old_dd)
    rm("dd", envir = globalenv())
  }, add = TRUE)

  m_orm <- rms::orm(y ~ x1 + x2, data = test_data, x = TRUE, y = TRUE)

  # Fit vglm
  m_vglm <- VGAM::vglm(
    y ~ x1 + x2,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = test_data
  )

  # Compute robust covariance
  orm_robust <- rms::robcov(m_orm)
  vglm_robust <- robcov_vglm(m_vglm)

  # Compare regression coefficient SEs (allowing ~2% tolerance for ordinal)
  # Extract regression coefficients (not intercepts)
  se_orm_reg <- sqrt(diag(orm_robust$var))[c("x1", "x2")]
  se_vglm_reg <- vglm_robust$se[c("x1", "x2")]

  # Ratios should be close to 1
  ratios <- se_vglm_reg / se_orm_reg
  expect_true(all(ratios > 0.97 & ratios < 1.03))
})



test_that("score contributions sum to approximately zero", {
  skip_if_not_installed("VGAM")

  # Generate test data
  set.seed(202)
  n <- 300
  x <- rnorm(n)
  lp <- 0.5 * x
  probs <- cbind(
    plogis(-1 - lp),
    plogis(0 - lp) - plogis(-1 - lp),
    1 - plogis(0 - lp)
  )
  y <- apply(probs, 1, function(p) sample(1:3, 1, prob = pmax(p, 0.001)))

  test_data <- data.frame(y = ordered(y), x = x)

  # Fit vglm
  m <- VGAM::vglm(
    y ~ x,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = test_data
  )

  # Compute robust covariance
  result <- robcov_vglm(m)

  # Score columns should sum to approximately zero (property of MLE)
  score_sums <- colSums(result$scores)
  expect_true(all(abs(score_sums) < 1e-3))
})


test_that("small-sample adjustment works correctly", {
  skip_if_not_installed("VGAM")

  # Generate clustered data
  set.seed(303)
  n_clusters <- 20
  n_per_cluster <- 15
  n <- n_clusters * n_per_cluster
  cluster <- rep(1:n_clusters, each = n_per_cluster)

  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.5 * x))

  test_data <- data.frame(y = y, x = x, cluster = cluster)

  # Fit vglm
  m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

  # Compute with and without adjustment
  result_noadj <- robcov_vglm(m, cluster = cluster, adjust = FALSE)
  result_adj <- robcov_vglm(m, cluster = cluster, adjust = TRUE)

  # Adjustment factor is G/(G-1) = 20/19 ≈ 1.0526
  expected_ratio <- sqrt(n_clusters / (n_clusters - 1))

  # SEs with adjustment should be larger by this factor
  se_ratio <- result_adj$se / result_noadj$se
  expect_equal(
    unname(se_ratio),
    rep(expected_ratio, length(se_ratio)),
    tolerance = 1e-10
  )
})


test_that("robcov_vglm validates inputs correctly", {
  skip_if_not_installed("VGAM")

  # Generate test data
  set.seed(404)
  n <- 100
  test_data <- data.frame(y = sample(1:3, n, replace = TRUE), x = rnorm(n))

  # Fit vglm
  m <- VGAM::vglm(
    ordered(y) ~ x,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = test_data
  )

  # Should error with non-vglm object
  expect_error(
    robcov_vglm(lm(y ~ x, data = test_data)),
    "must be a vglm object"
  )

  # Should error with wrong cluster length
  expect_error(
    robcov_vglm(m, cluster = 1:10),
    "Length of 'cluster' must equal"
  )
})


test_that("compare_se_orm_vglm produces correct comparison", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  # Generate test data with 4 levels (more typical for ordinal)
  set.seed(505)
  n <- 300
  x <- rnorm(n)
  lp <- 0.5 * x
  probs <- cbind(
    plogis(-1.5 - lp),
    plogis(-0.5 - lp) - plogis(-1.5 - lp),
    plogis(0.5 - lp) - plogis(-0.5 - lp),
    1 - plogis(0.5 - lp)
  )
  y <- apply(probs, 1, function(p) sample(1:4, 1, prob = pmax(p, 0.001)))

  test_data <- data.frame(y = ordered(y), x = x)

  # Fit orm
  dd <- rms::datadist(test_data)
  old_dd <- options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  on.exit({
    options(old_dd)
    rm("dd", envir = globalenv())
  }, add = TRUE)

  m_orm <- rms::orm(y ~ x, data = test_data, x = TRUE, y = TRUE)

  # Fit vglm
  m_vglm <- VGAM::vglm(
    y ~ x,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = test_data
  )

  # Compare
  comparison <- compare_se_orm_vglm(m_orm, m_vglm)

  # Check output structure
  expect_s3_class(comparison, "data.frame")
  expect_true(all(
    c(
      "parameter",
      "se_orm_model",
      "se_vglm_model",
      "se_orm_robust",
      "se_vglm_robust",
      "ratio_model",
      "ratio_robust"
    ) %in%
      names(comparison)
  ))

  # Check that ratios are reasonable (close to 1)
  expect_true(all(comparison$ratio_model > 0.9 & comparison$ratio_model < 1.1))
  expect_true(all(
    comparison$ratio_robust > 0.9 & comparison$ratio_robust < 1.1
  ))

  # Check number of rows matches number of parameters (should have 3 intercepts + 1 coef = 4)
  n_params <- length(coef(m_orm))
  expect_equal(nrow(comparison), n_params)
})



test_that("print and summary methods work", {
  skip_if_not_installed("VGAM")

  # Generate test data
  set.seed(606)
  n <- 100
  test_data <- data.frame(y = sample(0:1, n, replace = TRUE), x = rnorm(n))

  # Fit vglm
  m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

  # Compute robust covariance
  result <- robcov_vglm(m)

  # Test print method
  expect_snapshot(print(result))

  # Test summary method
  expect_snapshot(summary(result))

  # Test with clustering
  cluster <- rep(1:10, each = 10)
  result_cl <- robcov_vglm(m, cluster = cluster)
  expect_snapshot(summary(result_cl))

  # Test with clustering and adjustment
  result_adj <- robcov_vglm(m, cluster = cluster, adjust = TRUE)
  expect_snapshot(summary(result_adj))

  # Test print method returns invisibly
  expect_invisible(print(result))

  # Test summary returns coefficient table
  summary_result <- summary(result)
  expect_true("coefficients" %in% names(summary_result))
  expect_s3_class(summary_result$coefficients, "data.frame")
  expect_true(all(
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in%
      names(summary_result$coefficients)
  ))
})



test_that("coef and vcov methods work correctly", {
  skip_if_not_installed("VGAM")

  # Generate test data
  set.seed(707)
  n <- 100
  test_data <- data.frame(y = sample(0:1, n, replace = TRUE), x = rnorm(n))

  # Fit vglm
  m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

  # Compute robust covariance
  result <- robcov_vglm(m)

  # Test coef method
  expect_equal(coef(result), result$coefficients)
  expect_equal(coef(result), coef(m))
  expect_length(coef(result), 2) # intercept + x

  # Test vcov method
  expect_equal(vcov(result), result$var)
  expect_equal(dim(vcov(result)), c(2, 2))
  expect_equal(vcov(result), t(vcov(result))) # symmetric
})


test_that("robcov_vglm works with multiple predictors", {
  skip_if_not_installed("VGAM")

  # Generate test data with multiple predictors
  set.seed(707)
  n <- 400
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  lp <- 0.3 * x1 + 0.5 * x2 - 0.2 * x3
  probs <- cbind(
    plogis(-1 - lp),
    plogis(0 - lp) - plogis(-1 - lp),
    plogis(1 - lp) - plogis(0 - lp),
    1 - plogis(1 - lp)
  )
  y <- apply(probs, 1, function(p) sample(1:4, 1, prob = pmax(p, 0.001)))

  test_data <- data.frame(y = ordered(y), x1 = x1, x2 = x2, x3 = x3)

  # Fit vglm
  m <- VGAM::vglm(
    y ~ x1 + x2 + x3,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = test_data
  )

  # Compute robust covariance
  result <- robcov_vglm(m)

  # Check we have correct number of parameters (3 intercepts + 3 reg coefs = 6)
  expect_length(result$se, 6)
  expect_length(result$coefficients, 6)
  expect_length(result$z, 6)
  expect_length(result$pvalues, 6)

  # All SEs should be positive
  expect_true(all(result$se > 0))

  # Check parameter names
  expect_true("x1" %in% names(result$se))
  expect_true("x2" %in% names(result$se))
  expect_true("x3" %in% names(result$se))

  # Names should be consistent across all vectors
  expect_equal(names(result$coefficients), names(result$se))
  expect_equal(names(result$coefficients), names(result$z))
  expect_equal(names(result$coefficients), names(result$pvalues))
})


test_that("robcov_vglm preserves model information", {
  skip_if_not_installed("VGAM")

  # Generate test data
  set.seed(808)
  n <- 200
  test_data <- data.frame(
    y = ordered(sample(1:3, n, replace = TRUE)),
    x = rnorm(n)
  )

  # Fit vglm
  m <- VGAM::vglm(
    y ~ x,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = test_data
  )

  # Compute robust covariance
  result <- robcov_vglm(m)

  # Check model information is preserved
  expect_equal(result$coefficients, coef(m))
  expect_equal(result$fitted.values, m@fitted.values)
  expect_equal(result$residuals, m@residuals)
  expect_equal(result$df.residual, m@df.residual)
  expect_equal(result$df.total, m@df.total)

  # Check family is preserved
  expect_equal(class(result$family), class(m@family))

  # Check constraints are preserved
  expect_equal(result$constraints, m@constraints)

  # Check control is preserved
  expect_equal(result$control, m@control)

  # Check original vcov matches
  expect_equal(result$original_vcov, vcov(m))
  expect_equal(result$original_se, sqrt(diag(vcov(m))))
})
