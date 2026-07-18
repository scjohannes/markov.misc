# Unit Tests: Robust Covariance Estimation for vglm
#
# These tests verify that:
# 1. robcov_vglm produces valid sandwich covariance estimates
# 2. Results match rms::robcov for equivalent models
# 3. Clustering works correctly
# 4. Score computations have correct properties
# 5. Model information is preserved
# 6. Z-statistics and p-values are computed correctly

vglm_candidate_at <- function(fit, beta) {
  X_vlm <- stats::model.matrix(fit, type = "vlm")
  n <- stats::nobs(fit, type = "lm")
  M <- VGAM::npred(fit)
  fitted_eta <- matrix(
    drop(X_vlm %*% stats::coef(fit)),
    nrow = n,
    ncol = M,
    byrow = TRUE
  )
  eta <- matrix(
    drop(X_vlm %*% beta),
    nrow = n,
    ncol = M,
    byrow = TRUE
  ) +
    fit@predictors -
    fitted_eta
  candidate <- fit
  candidate@coefficients <- beta
  candidate@predictors <- eta
  candidate@fitted.values <- candidate@family@linkinv(
    eta = eta,
    extra = candidate@extra
  )
  candidate
}

finite_difference_loglik_scores <- function(fit) {
  beta <- stats::coef(fit)
  steps <- 1e-5 * pmax(1, abs(beta))
  vapply(
    seq_along(beta),
    function(j) {
      beta_plus <- beta_minus <- beta
      beta_plus[j] <- beta_plus[j] + steps[j]
      beta_minus[j] <- beta_minus[j] - steps[j]
      ll_plus <- as.numeric(stats::logLik(
        vglm_candidate_at(fit, beta_plus),
        summation = FALSE
      ))
      ll_minus <- as.numeric(stats::logLik(
        vglm_candidate_at(fit, beta_minus),
        summation = FALSE
      ))
      (ll_plus - ll_minus) / (2 * steps[j])
    },
    numeric(stats::nobs(fit, type = "lm"))
  )
}

finite_difference_score_jacobian <- function(fit) {
  evaluator <- markov.misc:::vglm_score_evaluator(fit)
  beta <- stats::coef(fit)
  steps <- 1e-5 * pmax(1, abs(beta))
  vapply(
    seq_along(beta),
    function(j) {
      beta_plus <- beta_minus <- beta
      beta_plus[j] <- beta_plus[j] + steps[j]
      beta_minus[j] <- beta_minus[j] - steps[j]
      (evaluator$evaluate(beta_plus) -
        evaluator$evaluate(beta_minus)) /
        (2 * steps[j])
    },
    numeric(length(beta))
  )
}

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
      c(
        "var",
        "se",
        "scores",
        "influence_beta_obs",
        "bread",
        "bread_type",
        "bread_diagnostics",
        "meat"
      )
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
    expect_equal(result$bread_type, "observed")
    expect_equal(result$type, "HC0")
    expect_identical(result$cadjust, FALSE)

    result_vglm_bread <- robcov_vglm(m, bread = "vglm")
    expect_equal(result_vglm_bread$bread, vcov(m))
    expect_equal(result_vglm_bread$bread_type, "vglm")
    expected_vglm_var <- crossprod(
      result_vglm_bread$scores %*% vcov(m)
    )
    expected_vglm_var <- (expected_vglm_var + t(expected_vglm_var)) / 2
    dimnames(expected_vglm_var) <- dimnames(vcov(m))
    expect_equal(result_vglm_bread$var, expected_vglm_var)

    expected_var <- crossprod(result$scores %*% result$bread)
    expected_var <- (expected_var + t(expected_var)) / 2
    dimnames(expected_var) <- dimnames(result$var)
    expect_equal(result$var, expected_var, tolerance = 1e-15)

    # Check variance matrix is symmetric
    expect_equal(result$var, t(result$var), tolerance = 1e-15)

    # Check variance matrix is positive semi-definite
    eigenvalues <- eigen(result$var, only.values = TRUE)$values
    expect_all_true(eigenvalues >= -1e-10)

    # Check standard errors are positive
    expect_all_true(result$se > 0)

    # Check n matches
    expect_equal(result$n, n)

    result_row_cluster <- robcov_vglm(m, cluster = seq_len(n), adjust = FALSE)
    expect_equal(result_row_cluster$var, result$var, tolerance = 1e-15)
  })

  it("stores influence contributions consistent with the selected bread", {
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
    expected_if <- result$scores %*% result$bread

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

  it("matches rms::robcov for an equivalent ordinal likelihood", {
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

    expect_equal(
      unname(vglm_robust$var),
      unname(orm_robust$var),
      tolerance = 1e-5
    )
    expect_equal(
      unname(vglm_robust$se),
      unname(sqrt(diag(orm_robust$var))),
      tolerance = 1e-5
    )

    cluster <- rep(seq_len(100), each = 5)
    orm_clustered <- rms::robcov(m_orm, cluster = cluster)
    vglm_clustered <- robcov_vglm(
      m_vglm,
      cluster = cluster,
      cadjust = FALSE
    )
    expect_equal(
      unname(vglm_clustered$var),
      unname(orm_clustered$var),
      tolerance = 1e-5
    )
    expect_lt(
      max(abs(vglm_clustered$se / sqrt(diag(orm_clustered$var)) - 1)),
      1e-5
    )
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

  it("rejects non-converged fits and malformed covariance components", {
    skip_if_not_installed("VGAM")

    withr::local_seed(405)
    n <- 120
    data <- data.frame(
      y = stats::rbinom(n, 1, 0.5),
      x = stats::rnorm(n)
    )
    fit <- VGAM::vglm(
      y ~ x,
      family = VGAM::binomialff,
      data = data
    )

    nonconverged <- fit
    nonconverged@iter <- nonconverged@control$maxit
    expect_snapshot(error = TRUE, {
      robcov_vglm(nonconverged)
    })

    testthat::local_mocked_bindings(
      compute_scores_vglm = function(fit) {
        scores <- matrix(
          0,
          nrow = stats::nobs(fit, type = "lm"),
          ncol = length(stats::coef(fit))
        )
        colnames(scores) <- rev(names(stats::coef(fit)))
        scores
      },
      .package = "markov.misc"
    )
    expect_snapshot(error = TRUE, {
      robcov_vglm(fit)
    })
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

  it("matches the loop-equivalent score calculation", {
    skip_if_not_installed("VGAM")

    withr::local_seed(203)
    n <- 160
    test_data <- data.frame(
      y = ordered(sample(1:4, n, replace = TRUE)),
      x1 = rnorm(n),
      x2 = rnorm(n)
    )
    m <- VGAM::vglm(
      y ~ x1 + x2,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = test_data
    )

    scores <- compute_scores_vglm(m)
    deriv_eta <- VGAM::weights(
      m,
      type = "working",
      deriv.arg = TRUE
    )$deriv
    X_vlm <- stats::model.matrix(m, type = "vlm")
    M <- VGAM::npred(m)

    scores_loop <- matrix(0, nrow = n, ncol = length(coef(m)))
    for (i in seq_len(n)) {
      vlm_rows <- ((i - 1) * M + 1):(i * M)
      scores_loop[i, ] <- as.vector(crossprod(
        X_vlm[vlm_rows, , drop = FALSE],
        deriv_eta[i, ]
      ))
    }
    colnames(scores_loop) <- names(coef(m))
    rownames(scores_loop) <- rownames(scores)

    expect_equal(scores, scores_loop)
  })
})

test_that("compare_se_orm_vglm compares model-based and robust SEs", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  withr::local_seed(8128)
  n <- 80
  data <- data.frame(
    y = rbinom(n, 1, 0.5),
    x = rnorm(n),
    cluster = rep(seq_len(20), each = 4)
  )

  dd <- rms::datadist(data)
  withr::local_options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  withr::defer(rm("dd", envir = globalenv()))

  orm_fit <- rms::orm(y ~ x, data = data, x = TRUE, y = TRUE)
  vglm_fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)

  unclustered <- compare_se_orm_vglm(orm_fit, vglm_fit)
  expect_warning(
    clustered <- compare_se_orm_vglm(orm_fit, vglm_fit, cluster = data$cluster),
    "Cluster-robust p-values",
    fixed = TRUE
  )

  expect_named(
    unclustered,
    c(
      "parameter",
      "se_orm_model",
      "se_vglm_model",
      "se_orm_robust",
      "se_vglm_robust",
      "ratio_model",
      "ratio_robust"
    )
  )
  expect_equal(nrow(clustered), nrow(unclustered))
  expect_true(all(is.finite(clustered$ratio_model)))
  expect_true(all(is.finite(clustered$ratio_robust)))
})

test_that("compare_se_orm_vglm warns when model intercept structures differ", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  withr::local_seed(8129)
  n <- 90
  data <- data.frame(
    y_ord = ordered(rep(1:3, length.out = n)),
    y_bin = rep(0:1, length.out = n),
    x = rnorm(n)
  )

  dd <- rms::datadist(data)
  withr::local_options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  withr::defer(rm("dd", envir = globalenv()))

  orm_fit <- rms::orm(y_ord ~ x, data = data, x = TRUE, y = TRUE)
  vglm_fit <- VGAM::vglm(y_bin ~ x, family = VGAM::binomialff, data = data)

  expect_warning(
    comparison <- compare_se_orm_vglm(orm_fit, vglm_fit),
    "Number of intercepts differs",
    fixed = TRUE
  )
  expect_s3_class(comparison, "data.frame")
})

describe("Small-sample adjustment", {
  it("separates HC1 and cluster corrections", {
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
    expect_warning(
      result_noadj <- robcov_vglm(m, cluster = cluster, cadjust = FALSE),
      "fewer than 30 clusters"
    )
    expect_warning(
      result_adj <- robcov_vglm(m, cluster = cluster, cadjust = TRUE),
      "fewer than 30 clusters"
    )
    expect_warning(
      result_hc1 <- robcov_vglm(
        m,
        cluster = cluster,
        type = "HC1",
        cadjust = FALSE
      ),
      "fewer than 30 clusters"
    )
    expect_warning(
      result_both <- robcov_vglm(
        m,
        cluster = cluster,
        type = "HC1",
        cadjust = TRUE
      ),
      "fewer than 30 clusters"
    )

    # Adjustment factor is G/(G-1) = 20/19 ≈ 1.0526
    expected_ratio <- sqrt(n_clusters / (n_clusters - 1))

    # SEs with adjustment should be larger by this factor
    se_ratio <- result_adj$se / result_noadj$se
    expect_equal(
      unname(se_ratio),
      rep(expected_ratio, length(se_ratio)),
      tolerance = 1e-10
    )
    expect_equal(
      result_adj$var,
      suppressWarnings(robcov_vglm(m, cluster = cluster))$var
    )
    hc1_factor <- (n - 1) / (n - length(coef(m)))
    expect_equal(result_hc1$var, result_noadj$var * hc1_factor)
    expect_equal(
      result_both$var,
      result_noadj$var * hc1_factor * n_clusters / (n_clusters - 1)
    )
    expect_equal(result_both$adjustment_factor, hc1_factor * 20 / 19)

    expect_warning(
      result_legacy <- robcov_vglm(m, cluster = cluster, adjust = TRUE),
      "fewer than 30 clusters"
    )
    expect_equal(result_legacy$var, result_adj$var)

    expect_snapshot(error = TRUE, {
      robcov_vglm(
        m,
        cluster = cluster,
        adjust = TRUE,
        cadjust = TRUE
      )
    })
  })
})

test_that("analytic scores and observed bread match independent differences", {
  skip_if_not_installed("VGAM")

  withr::local_seed(204)
  n <- 240
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  eta <- 0.4 * x - 0.2 * z
  cuts <- c(-1, 0, 1)
  cumulative <- vapply(
    cuts,
    function(cut) {
      stats::plogis(cut - eta)
    },
    numeric(n)
  )
  probabilities <- cbind(
    cumulative[, 1],
    cumulative[, 2] - cumulative[, 1],
    cumulative[, 3] - cumulative[, 2],
    1 - cumulative[, 3]
  )
  y <- apply(probabilities, 1, function(probability) {
    sample(seq_len(4), 1, prob = probability)
  })
  data <- data.frame(
    y = ordered(y),
    x = x,
    z = z,
    offset = 0.1 * z,
    weight = sample(1:3, n, replace = TRUE)
  )
  fit <- VGAM::vglm(
    y ~ x + offset(offset),
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    weights = weight,
    x = TRUE,
    y = TRUE
  )

  scores <- compute_scores_vglm(fit)
  numerical_scores <- finite_difference_loglik_scores(fit)
  expect_equal(scores, numerical_scores, tolerance = 1e-6, ignore_attr = TRUE)

  numerical_jacobian <- finite_difference_score_jacobian(fit)
  numerical_information <- -(numerical_jacobian + t(numerical_jacobian)) / 2
  numerical_bread <- solve(numerical_information)
  observed <- robcov_vglm(fit)
  expect_equal(
    unname(observed$bread),
    unname(numerical_bread),
    tolerance = 1e-6
  )
})

test_that("observed bread supports vector-valued inverse links", {
  skip_if_not_installed("VGAM")

  withr::local_seed(206)
  n <- 180
  data <- data.frame(x = stats::rnorm(n))
  data$y <- 0.2 + 0.4 * data$x + stats::rnorm(n)
  fit <- VGAM::vglm(
    y ~ x,
    family = VGAM::uninormal,
    data = data
  )

  inverse_link <- fit@family@linkinv(
    eta = fit@predictors,
    extra = fit@extra
  )
  expect_null(dim(inverse_link))

  evaluator <- markov.misc:::vglm_score_evaluator(fit)
  evaluated_scores <- evaluator$evaluate(stats::coef(fit))
  expected_scores <- colSums(compute_scores_vglm(fit))
  expect_equal(evaluated_scores, expected_scores, tolerance = 1e-12)

  robust <- robcov_vglm(fit)
  expect_s3_class(robust, "robcov_vglm")
  expect_equal(robust$bread_type, "observed")
  expect_equal(all(is.finite(robust$bread)), TRUE)
  expect_equal(all(is.finite(robust$var)), TRUE)
})

test_that("observed bread supports proportional-odds deviation geometries", {
  skip_if_not_installed("VGAM")

  withr::local_seed(205)
  n <- 500
  x <- stats::runif(n, -1, 1)
  z <- stats::rnorm(n)
  group <- factor(sample(c("a", "b"), n, replace = TRUE))
  eta <- 0.5 * x - 0.25 * z + 0.1 * (group == "b")
  cuts <- c(-1.2, -0.3, 0.5, 1.3)
  cumulative <- vapply(
    cuts,
    function(cut) {
      stats::plogis(cut - eta)
    },
    numeric(n)
  )
  probabilities <- cbind(
    cumulative[, 1],
    cumulative[, 2] - cumulative[, 1],
    cumulative[, 3] - cumulative[, 2],
    cumulative[, 4] - cumulative[, 3],
    1 - cumulative[, 4]
  )
  y <- apply(probabilities, 1, function(probability) {
    sample(seq_len(5), 1, prob = probability)
  })
  data <- data.frame(y = ordered(y), x = x, z = z, group = group)

  full_nonparallel <- VGAM::vglm(
    y ~ x + z,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
    data = data,
    x = TRUE,
    y = TRUE
  )
  custom_constraints <- VGAM::constraints(full_nonparallel)
  n_thresholds <- VGAM::npred(full_nonparallel)
  custom_constraints[["x"]] <- cbind(
    PO = rep(1, n_thresholds),
    linear = seq_len(n_thresholds) - mean(seq_len(n_thresholds)),
    terminal = c(rep(0, n_thresholds - 1), 1)
  )
  custom_constraints[["z"]] <- cbind(PO = rep(1, n_thresholds))

  fits <- list(
    proportional = VGAM::vglm(
      y ~ x + z + group,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      x = TRUE,
      y = TRUE
    ),
    selected_nonparallel = VGAM::vglm(
      y ~ x + z + group,
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE ~ x),
      data = data,
      x = TRUE,
      y = TRUE
    ),
    fully_nonparallel = full_nonparallel,
    custom_deviations = VGAM::vglm(
      y ~ x + z,
      family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
      data = data,
      constraints = custom_constraints,
      x = TRUE,
      y = TRUE
    ),
    interactions = VGAM::vglm(
      y ~ x * z + group,
      family = VGAM::cumulative(
        reverse = TRUE,
        parallel = FALSE ~ x + x:z
      ),
      data = data,
      x = TRUE,
      y = TRUE
    )
  )

  for (fit in fits) {
    robust <- robcov_vglm(
      fit,
      cluster = seq_len(n),
      cadjust = FALSE
    )
    expect_equal(robust$bread_type, "observed")
    expect_equal(colnames(robust$scores), names(stats::coef(fit)))
    expect_equal(robust$var, t(robust$var), tolerance = 1e-12)
    expect_gt(min(eigen(robust$bread, symmetric = TRUE)$values), 0)
    expect_gte(min(eigen(robust$var, symmetric = TRUE)$values), -1e-10)
    expect_equal(all(is.finite(robust$se)), TRUE)

    draws <- generate_score_bootstrap_draws(
      robust,
      baseline_data = data.frame(id = seq_len(n)),
      id_var = "id",
      n_draws = 1
    )
    expect_equal(dim(draws$beta_draws), c(1L, length(stats::coef(fit))))
  }
})

describe("Score bootstrap compatibility", {
  it("uses the selected bread directly for one-step score updates", {
    skip_if_not_installed("VGAM")

    withr::local_seed(304)
    n_clusters <- 30
    n_per_cluster <- 4
    n <- n_clusters * n_per_cluster
    cluster <- rep(sprintf("id%02d", seq_len(n_clusters)), each = n_per_cluster)
    x <- rnorm(n)
    y <- rbinom(n, 1, plogis(-0.2 + 0.4 * x))
    test_data <- data.frame(y = y, x = x, id = cluster)

    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)
    robust <- robcov_vglm(m, cluster = test_data$id, adjust = FALSE)
    baseline_data <- data.frame(id = unique(test_data$id))

    cluster_ids <- unique(as.character(robust$cluster))
    scores_clustered <- rowsum(
      robust$scores,
      factor(as.character(robust$cluster), levels = cluster_ids),
      reorder = FALSE
    )

    withr::local_seed(305)
    w_cluster <- stats::rexp(length(cluster_ids), rate = 1)
    U_star <- as.vector(crossprod(w_cluster - 1, scores_clustered))
    expected_beta <- robust$coefficients + as.vector(robust$bread %*% U_star)

    withr::local_seed(305)
    draws <- generate_score_bootstrap_draws(
      robust,
      baseline_data = baseline_data,
      id_var = "id",
      n_draws = 1
    )

    expect_equal(unname(draws$beta_draws[1, ]), unname(expected_beta))
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

    # Test with clustering and no adjustment
    cluster <- rep(1:10, each = 10)
    expect_warning(
      result_cl <- robcov_vglm(m, cluster = cluster, adjust = FALSE),
      "fewer than 30 clusters"
    )
    expect_snapshot(summary(result_cl))

    # Test with clustering and explicit adjustment
    expect_warning(
      result_adj <- robcov_vglm(m, cluster = cluster, adjust = TRUE),
      "fewer than 30 clusters"
    )
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

    withr::local_seed(9001)
    n <- 200
    x <- rnorm(n)
    y <- rbinom(n, 1, 0.5)
    x[1:10] <- NA
    y[11:20] <- NA

    test_data <- data.frame(y = y, x = x, id = 1:n)
    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)
    omitted <- as.integer(unlist(m@na.action, use.names = FALSE))
    aligned_cluster <- test_data$id[-omitted]

    res_original <- robcov_vglm(m, cluster = test_data$id)
    res_aligned <- robcov_vglm(m, cluster = aligned_cluster)

    expect_s3_class(res_original, "robcov_vglm")
    expect_equal(res_original$n, 180) # 200 - 20 dropped
    expect_equal(res_original$cluster, aligned_cluster)
    expect_equal(res_original$var, res_aligned$var)
  })

  it("rejects missing and single-valued clusters", {
    skip_if_not_installed("VGAM")

    withr::local_seed(90011)
    n <- 80
    test_data <- data.frame(
      y = rbinom(n, 1, 0.5),
      x = rnorm(n),
      id = rep(1:8, each = 10)
    )
    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

    cluster_missing <- test_data$id
    cluster_missing[1] <- NA

    expect_error(
      robcov_vglm(m, cluster = cluster_missing),
      "contains missing values"
    )
    expect_error(
      robcov_vglm(m, cluster = rep(1, n)),
      "at least two clusters"
    )
  })

  it("warns for few clusters while retaining z-test p-values", {
    skip_if_not_installed("VGAM")

    withr::local_seed(90012)
    n_clusters <- 5
    n_per <- 20
    n <- n_clusters * n_per
    test_data <- data.frame(
      y = rbinom(n, 1, 0.5),
      x = rnorm(n),
      id = rep(seq_len(n_clusters), each = n_per)
    )
    m <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = test_data)

    expect_warning(
      res <- robcov_vglm(m, cluster = test_data$id),
      "fewer than 30 clusters"
    )
    expect_equal(res$pvalues, 2 * pnorm(-abs(res$z)))
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

    # This should NOT crash, but warns because few clusters give fragile z-tests.
    expect_warning(
      res <- robcov_vglm(m, cluster = cluster),
      "fewer than 30 clusters"
    )

    expect_s3_class(res, "robcov_vglm")

    # Check if the variance matrix is singular (determinant approx 0)
    # We expect it to be singular here because Rank(Meat) <= 3 < 4
    # But robcov_vglm shouldn't crash
    expect_equal(res$n_clusters, 3)
  })

  it("uses a stable PSD construction for ill-conditioned clustered fits", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    trial <- sim_actt2_brownian(
      n_patients = 40,
      follow_up_time = 6,
      seed = 1
    )
    data <- prepare_markov_data(trial, absorbing_state = 8)
    fit <- suppressWarnings(vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    ))
    robust <- robcov_vglm(fit, cluster = data$id)

    cluster_scores <- rowsum(
      robust$scores,
      as.factor(robust$cluster),
      reorder = FALSE
    )
    expected <- crossprod(cluster_scores %*% robust$bread) *
      robust$adjustment_factor
    dimnames(expected) <- dimnames(robust$var)
    expect_equal(robust$var, expected, tolerance = 1e-12)
    expect_gte(
      min(eigen(robust$var, symmetric = TRUE, only.values = TRUE)$values),
      -1e-10
    )
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
    expect_warning(
      res <- robcov_vglm(m, cluster = cluster_fac),
      "fewer than 30 clusters"
    )

    expect_s3_class(res, "robcov_vglm")

    # Should only detect 2 actual clusters (A and B), not 3 factor levels
    expect_equal(res$n_clusters, 2)
    expect_all_true(!is.na(res$se))
  })
})
