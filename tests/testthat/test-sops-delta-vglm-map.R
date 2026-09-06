test_that("VGLM coefficient maps match native expansion on every basis vector", {
  skip_if_not_installed("VGAM")
  set.seed(7307)
  data <- data.frame(
    x = stats::rnorm(240),
    group = factor(rep(letters[1:3], 80))
  )
  data$y <- ordered(cut(
    0.3 * data$x + stats::rlogis(240),
    c(-Inf, -1, 0, 1, Inf)
  ))
  models <- list(
    VGAM::vglm(
      y ~ rms::rcs(x, 3) * group,
      VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    ),
    VGAM::vglm(
      y ~ group + I(x^2) + x,
      VGAM::cumulative(reverse = TRUE, parallel = TRUE, Thresh = "equid"),
      data = data
    )
  )
  # A non-unit common-slope basis is still proportional odds. Reparameterize
  # the fitted model without changing its effective coefficients.
  scaled <- models[[2L]]
  scaled@constraints[["x"]] <- -2 * scaled@constraints[["x"]]
  scaled@coefficients["x"] <- scaled@coefficients["x"] / -2
  models[[3L]] <- scaled

  for (model in models) {
    map <- get_effective_coef_map(model)
    beta <- stats::coef(model)
    expect_identical(colnames(map), names(beta))
    expect_identical(
      attr(map, "gamma_dimnames")[[2L]],
      rownames(stats::coef(model, matrix = TRUE))
    )
    for (j in seq_along(beta)) {
      basis <- model
      basis@coefficients[] <- 0
      basis@coefficients[j] <- 1
      expect_equal(
        unname(map[, j]),
        as.vector(t(stats::coef(basis, matrix = TRUE))),
        tolerance = 1e-12
      )
    }
  }

  zero_columns <- scaled
  zero_columns@constraints[["x"]] <- matrix(0, nrow = 3, ncol = 0)
  expect_error(get_effective_coef_map(zero_columns), "no coefficient columns")
})
