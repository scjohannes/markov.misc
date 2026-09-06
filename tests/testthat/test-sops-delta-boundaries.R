test_that("zero or tiny uncertainty does not certify a structural boundary", {
  estimate <- rep(c(0, 1), 3)
  standard_error <- rep(c(0, 1e-20, 1e-7), each = 2)
  warning <- NULL
  bounds <- withCallingHandlers(
    delta_interval_bounds(estimate, standard_error, 0.95, "logit"),
    warning = function(condition) {
      warning <<- conditionMessage(condition)
      invokeRestart("muffleWarning")
    }
  )
  expect_match(warning, "nonstructural boundary", fixed = TRUE)
  expect_identical(bounds$conf.low, rep(NA_real_, 6))
  expect_identical(bounds$conf.high, rep(NA_real_, 6))
  wald <- delta_interval_bounds(estimate, standard_error, 0.95, "wald")
  expect_equal(wald$conf.low, estimate - stats::qnorm(0.975) * standard_error)
  expect_equal(wald$conf.high, estimate + stats::qnorm(0.975) * standard_error)
})

test_that("absorption during prediction does not certify numerical saturation", {
  data <- make_test_data(n_patients = 60, follow_up_time = 8, seed = 8118)
  model <- suppressWarnings(vglm_markov(
    ordered(y) ~ tx,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  ))
  states <- seq_len(nrow(get_effective_coefs(model)) + 1L)
  baseline <- data[!is.na(data$y) & data$yprev != max(states), , drop = FALSE]
  baseline <- baseline[!duplicated(baseline$id), , drop = FALSE][1:4, ]
  individual <- sops(
    model,
    newdata = baseline,
    times = 1:3,
    y_levels = states,
    absorb = max(states)
  )
  expect_gt(min(individual$estimate), 0)
  expect_lt(max(individual$estimate), 1)
  absorbing <- individual[individual$state == max(states), ]
  expect_equal(
    all(vapply(
      split(absorbing, absorbing$rowid),
      function(profile) {
        all(diff(profile$estimate[order(profile$time)]) > 0)
      },
      logical(1)
    )),
    TRUE
  )

  coefficients <- get_coef(model)
  coefficients[grepl("Intercept", names(coefficients))] <-
    coefficients[grepl("Intercept", names(coefficients))] + 1000
  saturated_model <- set_coef(model, coefficients)
  points <- list(
    sops(
      saturated_model,
      newdata = baseline,
      times = 1:3,
      y_levels = states,
      absorb = max(states)
    ),
    avg_sops(
      saturated_model,
      newdata = baseline,
      times = 1:3,
      variables = list(tx = c(0, 1)),
      y_levels = states,
      absorb = max(states)
    )
  )
  for (point in points) {
    expect_equal(all(point$estimate %in% c(0, 1)), TRUE)
    warning <- NULL
    result <- withCallingHandlers(
      inferences(point, method = "delta", vcov = stats::vcov(model)),
      warning = function(condition) {
        warning <<- conditionMessage(condition)
        invokeRestart("muffleWarning")
      }
    )
    expect_match(warning, "nonstructural boundary", fixed = TRUE)
    expect_equal(result$std.error, rep(0, nrow(result)))
    expect_identical(result$conf.low, rep(NA_real_, nrow(result)))
    expect_identical(result$conf.high, rep(NA_real_, nrow(result)))
  }
})
