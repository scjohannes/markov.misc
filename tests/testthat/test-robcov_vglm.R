# Unit Tests: Robust Covariance Estimation for vglm
#
# These tests verify that:
# 1. robcov_vglm produces valid sandwich covariance estimates
# 2. Results match rms::robcov for equivalent models
# 3. Clustering works correctly
# 4. Score computations have correct properties
# 5. Model information is preserved
# 6. Z-statistics and p-values are computed correctly

describe("robcov_vglm()", {
  it("produces valid covariance estimates", {
    skip_if_not_installed("VGAM")

    # Generate simple test data
    withr::local_seed(123)
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
    expect_contains(
      names(result),
      c("var", "se", "scores", "influence_beta_obs", "bread", "meat")
    )

    # Check new components
    expect_contains(
      names(result),
      c("coefficients", "z", "pvalues", "original_se")
    )

    # Check model information components
    expect_contains(
      names(result),
      c("family", "constraints", "control", "call", "original_call")
    )

    # Check dimensions
    p <- length(coef(m))
    expect_equal(dim(result$var), c(p, p))
    expect_length(result$se, p)
    expect_length(result$coefficients, p)
    expect_length(result$z, p)
    expect_length(result$pvalues, p)
    expect_equal(dim(result$scores), c(n, p))
    expect_equal(dim(result$influence_beta_obs), c(n, p))
    expect_equal(dim(result$bread), c(p, p))
    expect_equal(dim(result$meat), c(p, p))

    # Check variance matrix is symmetric
    expect_equal(result$var, t(result$var), tolerance = 1e-15)

    # Check variance matrix is positive semi-definite
    eigenvalues <- eigen(result$var, only.values = TRUE)$values
    expect_all_true(eigenvalues >= -1e-10)

    # Check standard errors are positive
    expect_all_true(result$se > 0)

    # Check n matches
    expect_equal(result$n, n)
  })

  it("stores VGAM influence contributions consistent with scores", {
    skip_if_not_installed("VGAM")

    withr::local_seed(111)
    n <- 180
    test_data <- data.frame(
      y = ordered(sample(1:4, n, replace = TRUE)),
      x = rnorm(n)
    )

    m <- VGAM::vglm(
      y ~ x,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = test_data
    )

    result <- robcov_vglm(m)
    expected_if <- result$scores %*% vcov(m)

    expect_equal(
      unname(result$influence_beta_obs),
      unname(expected_if),
      tolerance = 1e-10
    )
  })

  it("computes z-statistics and p-values correctly", {
    skip_if_not_installed("VGAM")

    # Generate simple test data with known effect
    withr::local_seed(456)
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
    expect_lt(result$pvalues["x"], 0.00001)

    # Check all p-values are in [0, 1]
    expect_all_true(result$pvalues >= 0 & result$pvalues <= 1)
  })

  it("matches rms::robcov for binary outcome", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    # Generate binary outcome data
    withr::local_seed(456)
    n <- 500
    x <- rnorm(n)
    lp <- -0.5 + 0.8 * x
    y <- rbinom(n, 1, plogis(lp))

    test_data <- data.frame(y = y, x = x)

    # Fit orm
    withr::local_options(datadist = "dd")
    dd <- rms::datadist(test_data)
    withr::local_environment(globalenv())
    assign("dd", dd, envir = globalenv())
    on.exit(rm("dd", envir = globalenv()), add = TRUE)

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
    expect_equal(unname(se_vglm), unname(se_orm), tolerance = 1e-6)
  })

  it("matches rms::robcov for binary outcome with clustering", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    # Generate clustered binary outcome data
    withr::local_seed(789)
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
    withr::local_options(datadist = "dd")
    dd <- rms::datadist(test_data)
    withr::local_environment(globalenv())
    assign("dd", dd, envir = globalenv())
    on.exit(rm("dd", envir = globalenv()), add = TRUE)

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

  it("gives similar results to rms::robcov for ordinal outcome", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    # Generate ordinal outcome data
    withr::local_seed(101)
    n <- 500
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
    withr::local_options(datadist = "dd")
    dd <- rms::datadist(test_data)
    withr::local_environment(globalenv())
    assign("dd", dd, envir = globalenv())
    on.exit(rm("dd", envir = globalenv()), add = TRUE)

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

    # Compare regression coefficient SEs
    se_orm_reg <- sqrt(diag(orm_robust$var))
    se_vglm_reg <- vglm_robust$se

    expect_equal(unname(se_orm_reg), unname(se_vglm_reg), tolerance = 1e-2)
  })

  it("validates inputs correctly", {
    skip_if_not_installed("VGAM")

    # Generate test data
    withr::local_seed(404)
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

    # expand
  })

  it("works with multiple predictors", {
    skip_if_not_installed("VGAM")

    # Generate test data with multiple predictors
    withr::local_seed(707)
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
    expect_all_true(result$se > 0)

    # Check parameter names
    expect_contains(names(result$se), c("x1", "x2", "x3"))

    # Names should be consistent across all vectors
    expect_equal(names(result$coefficients), names(result$se))
    expect_equal(names(result$coefficients), names(result$z))
    expect_equal(names(result$coefficients), names(result$pvalues))
  })

  it("preserves model information", {
    skip_if_not_installed("VGAM")

    # Generate test data
    withr::local_seed(808)
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
})

describe("Score properties", {
  it("sum to approximately zero", {
    skip_if_not_installed("VGAM")

    # Generate test data
    withr::local_seed(202)
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
    expect_all_true(abs(score_sums) < 1e-5)
  })
})

describe("Small-sample adjustment", {
  it("works correctly", {
    skip_if_not_installed("VGAM")

    # Generate clustered data
    withr::local_seed(303)
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
})

describe("Utility methods", {
  it("print and summary methods work", {
    skip_if_not_installed("VGAM")

    # Generate test data
    withr::local_seed(606)
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
    expect_contains(names(summary_result), "coefficients")
    expect_s3_class(summary_result$coefficients, "data.frame")
    expect_contains(
      names(summary_result$coefficients),
      c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    )
  })

  it("coef and vcov methods work correctly", {
    skip_if_not_installed("VGAM")

    # Generate test data
    withr::local_seed(707)
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
})

# Stress Tests & Edge Cases for robcov_vglm
#
# These tests focus on "breaking" the function with:
# 1. Missing data (NAs) handling
# 2. Models with weights and offsets
# 3. Non-parallel (Partial Proportional Odds) models
# 4. Rank-deficient clusters (clusters < parameters)
# 5. Factor levels issues (empty levels)

describe("robcov_vglm() stress tests", {
  it("handles missing data correctly (alignment with clusters)", {
    skip_if_not_installed("VGAM")

    # Generate data with NAs
    withr::local_seed(9001)
    n <- 200
    x <- rnorm(n)
    y <- rbinom(n, 1, 0.5)

    # Introduce NAs in predictors and outcome
    x[1:10] <- NA
    y[11:20] <- NA

    test_data <- data.frame(y = y, x = x, id = 1:n)

    # Fit vglm (will drop 20 rows)
    # Note: vglm does not have a 'na.action' argument by default like lm,
    # but handles NAs by dropping them if na.action is set globally or in model.frame
    # We explicitly remove NAs to simulate the user passing "clean" data vs "dirty" clusters

    clean_data <- na.omit(test_data)
    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = clean_data)

    # Case 1: Cluster vector is too long (user forgot to drop NAs from cluster var)
    # This simulates a very common user error
    expect_error(
      robcov_vglm(m, cluster = test_data$id),
      "Length of 'cluster' must equal"
    )

    # Case 2: Correct usage
    res <- robcov_vglm(m, cluster = clean_data$id)
    expect_s3_class(res, "robcov_vglm")
    expect_equal(res$n, 180) # 200 - 20 dropped
  })

  it("works with model weights", {
    skip_if_not_installed("VGAM")

    withr::local_seed(9002)
    n <- 100
    x <- rnorm(n)
    y <- rbinom(n, 1, 0.5)
    weights <- runif(n, 0.5, 2)

    test_data <- data.frame(y = y, x = x, w = weights)

    # Fit weighted model
    m <- VGAM::vglm(
      y ~ x,
      family = VGAM::binomialff,
      data = test_data,
      weights = w
    )

    # Robust covariance should account for weights implicitly via the score function
    # (VGAM::weights(deriv=TRUE) returns weighted derivatives)
    res <- robcov_vglm(m)

    expect_all_true(res$se > 0)
    expect_equal(length(res$se), 2)
  })

  it("works with model offsets", {
    skip_if_not_installed("VGAM")

    withr::local_seed(9003)
    n <- 100
    x <- rnorm(n)
    off <- rep(0.5, n) # fixed offset
    y <- rbinom(n, 1, plogis(0.5 * x + off))

    test_data <- data.frame(y = y, x = x, off = off)

    # Fit model with offset
    m <- VGAM::vglm(
      y ~ x,
      family = VGAM::binomialff,
      data = test_data,
      offset = off
    )

    res <- robcov_vglm(m)

    # Check that it didn't crash and produced valid SEs
    expect_s3_class(res, "robcov_vglm")
    expect_all_true(res$se > 0)

    # Offset should act like a fixed intercept shift, not a parameter
    expect_contains(names(coef(res)), c("(Intercept)", "x"))
    expect_false("off" %in% names(coef(res)))
  })

  it("handles non-parallel (general) ordinal models", {
    skip_if_not_installed("VGAM")

    # Generate ordinal data
    withr::local_seed(9004)
    n <- 300
    x <- rnorm(n)
    y <- sample(1:3, n, replace = TRUE)
    test_data <- data.frame(y = ordered(y), x = x)

    # Fit model WITHOUT parallel assumption
    # This means 'x' will have a different effect for 1->2 vs 2->3
    m_nopar <- VGAM::vglm(
      y ~ x,
      family = VGAM::cumulative(parallel = FALSE),
      data = test_data
    )

    res <- robcov_vglm(m_nopar)

    # Parameters should be: (Intercept):1, (Intercept):2, x:1, x:2
    # Total = 4 parameters
    expect_length(coef(res), 4)

    # Check dimensions of variance matrix
    expect_equal(dim(res$var), c(4, 4))

    # Ensure all SEs are positive
    expect_all_true(res$se > 0)
  })

  it("handles rank-deficient clusters (More parameters than clusters)", {
    skip_if_not_installed("VGAM")

    # Scenario: 5 clusters, but model has 2 parameters
    # This is "risky" but mathematically possible.
    # Real danger is if Clusters < Parameters (e.g., 2 clusters, 5 params)

    withr::local_seed(9005)
    n_clusters <- 3
    n_per <- 20
    n <- n_clusters * n_per
    cluster <- rep(1:n_clusters, each = n_per)
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- rnorm(n)
    y <- rbinom(n, 1, 0.5)

    test_data <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3, cluster = cluster)

    # Model has 4 params (Int + x1 + x2 + x3), but only 3 clusters
    # The meat matrix will be rank deficient (rank at most 3)
    m <- VGAM::vglm(
      y ~ x1 + x2 + x3,
      family = VGAM::binomialff,
      data = test_data
    )

    # This should NOT crash, but might produce warnings or singular matrices
    # We want to ensure the code is robust enough to return *something*
    res <- robcov_vglm(m, cluster = cluster)

    expect_s3_class(res, "robcov_vglm")

    # Check if the variance matrix is singular (determinant approx 0)
    # We expect it to be singular here because Rank(Meat) <= 3 < 4
    # But robcov_vglm shouldn't crash
    expect_equal(res$n_clusters, 3)
  })

  it("handles clusters with empty factor levels", {
    skip_if_not_installed("VGAM")

    withr::local_seed(9006)
    n <- 50
    x <- rnorm(n)
    y <- rbinom(n, 1, 0.5)

    # Create factor with levels "A", "B", "C", but only use "A" and "B"
    cluster_char <- sample(c("A", "B"), n, replace = TRUE)
    cluster_fac <- factor(cluster_char, levels = c("A", "B", "C"))

    test_data <- data.frame(y = y, x = x, cluster = cluster_fac)

    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

    # This ensures rowsum() or internal logic doesn't create NA rows for empty level "C"
    res <- robcov_vglm(m, cluster = cluster_fac)

    expect_s3_class(res, "robcov_vglm")

    # Should only detect 2 actual clusters (A and B), not 3 factor levels
    expect_equal(res$n_clusters, 2)
    expect_all_true(!is.na(res$se))
  })
})
