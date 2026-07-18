delta_vglm_case <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      data <- make_test_data(
        n_patients = 100,
        follow_up_time = 8,
        seed = 7301
      )
      model <- make_test_model(data)
      y_levels <- seq_len(nrow(get_effective_coefs(model)) + 1L)
      baseline <- data[!duplicated(data$id), , drop = FALSE][
        1:4,
        ,
        drop = FALSE
      ]
      plan <- compile_sop_execution_plan(
        model = model,
        newdata = baseline,
        times = 1:5,
        y_levels = y_levels,
        absorb = max(y_levels)
      )
      value <<- list(
        data = data,
        model = model,
        baseline = baseline,
        y_levels = y_levels,
        plan = plan
      )
    }
    value
  }
})

delta_orm_case <- local({
  value <- NULL
  function() {
    if (is.null(value)) {
      data <- make_test_data(
        n_patients = 100,
        follow_up_time = 8,
        seed = 7311
      )
      model <- rms::orm(
        y ~ time + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      )
      baseline <- data[!duplicated(data$id), , drop = FALSE][
        1:3,
        ,
        drop = FALSE
      ]
      plan <- compile_sop_execution_plan(
        model = model,
        newdata = baseline,
        times = 1:4,
        y_levels = model$yunique,
        absorb = max(model$yunique)
      )
      value <<- list(
        data = data,
        model = model,
        baseline = baseline,
        y_levels = model$yunique,
        plan = plan
      )
    }
    value
  }
})

central_sop_jacobian <- function(plan, model) {
  coef <- get_coef(model)
  point <- run_sop_execution_plan(plan, get_effective_coefs(model))
  out <- array(0, dim = c(dim(point), length(coef)))
  for (coefficient in seq_along(coef)) {
    step <- .Machine$double.eps^(1 / 3) * max(1, abs(coef[coefficient]))
    upper <- coef
    lower <- coef
    upper[coefficient] <- upper[coefficient] + step
    lower[coefficient] <- lower[coefficient] - step
    upper_model <- set_coef(model, upper)
    lower_model <- set_coef(model, lower)
    upper_probability <- run_sop_execution_plan(
      plan,
      get_effective_coefs(upper_model)
    )
    lower_probability <- run_sop_execution_plan(
      plan,
      get_effective_coefs(lower_model)
    )
    out[,,, coefficient] <-
      (upper_probability - lower_probability) / (2 * step)
  }
  dimnames(out) <- list(NULL, NULL, dimnames(point)[[3L]], names(coef))
  out
}

make_delta_factor_visit_case <- function(
  n_patients = 80,
  n_visits = 4,
  seed = 7321
) {
  set.seed(seed)
  visits <- as.character(seq_len(n_visits))
  data <- expand.grid(
    id = seq_len(n_patients),
    time = visits,
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  data$time <- factor(data$time, levels = visits)
  treatment <- stats::rbinom(n_patients, 1, 0.5)
  data$tx <- treatment[data$id]
  data$yprev <- factor(
    sample(seq_len(4), nrow(data), replace = TRUE),
    levels = as.character(seq_len(4))
  )
  eta <-
    -0.35 *
    data$tx +
    0.25 * as.integer(data$time) +
    0.2 * as.integer(as.character(data$yprev))
  probability <- stats::plogis(eta - mean(eta))
  y <- pmin(
    4L,
    pmax(1L, 1L + stats::rbinom(nrow(data), 3L, probability))
  )
  data$y <- ordered(y, levels = as.character(seq_len(4)))
  data
}

test_that("VGLM analytical SOP Jacobians match raw-coefficient differences", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  case <- delta_vglm_case()

  analytical <- run_sop_delta_plan(case$plan, case$model)
  point <- run_sop_execution_plan(
    case$plan,
    get_effective_coefs(case$model)
  )
  numeric <- central_sop_jacobian(case$plan, case$model)

  expect_named(
    analytical,
    c("probabilities", "jacobian", "coef", "coef_names")
  )
  expect_identical(
    dim(analytical$jacobian),
    c(dim(analytical$probabilities), length(get_coef(case$model)))
  )
  expect_lt(
    max(abs(unname(analytical$probabilities) - unname(point))),
    1e-12
  )
  expect_lt(max(abs(unname(analytical$jacobian) - unname(numeric))), 1e-7)
  expect_lt(
    max(abs(apply(analytical$probabilities, c(1L, 2L), sum) - 1)),
    1e-10
  )
  expect_lt(
    max(abs(apply(analytical$jacobian, c(1L, 2L, 4L), sum))),
    1e-10
  )
  expect_identical(analytical$coef_names, names(get_coef(case$model)))
  expect_identical(names(analytical$coef), analytical$coef_names)
})

test_that("ORM analytical SOP Jacobians match raw-coefficient differences", {
  skip_if_not_installed("rms")
  case <- delta_orm_case()

  analytical <- run_sop_delta_plan(case$plan, case$model)
  point <- run_sop_execution_plan(
    case$plan,
    get_effective_coefs(case$model)
  )
  numeric <- central_sop_jacobian(case$plan, case$model)

  expect_lt(
    max(abs(unname(analytical$probabilities) - unname(point))),
    1e-12
  )
  expect_lt(max(abs(unname(analytical$jacobian) - unname(numeric))), 1e-7)
  expect_lt(
    max(abs(apply(analytical$probabilities, c(1L, 2L), sum) - 1)),
    1e-10
  )
  expect_lt(
    max(abs(apply(analytical$jacobian, c(1L, 2L, 4L), sum))),
    1e-10
  )
})

test_that("effective coefficient maps exactly preserve named raw coefficients", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  for (case in list(delta_vglm_case(), delta_orm_case())) {
    map <- get_effective_coef_map(case$model)
    coef <- get_coef(case$model)
    gamma <- get_effective_coefs(case$model)

    expect_identical(colnames(map), names(coef))
    expect_identical(anyDuplicated(colnames(map)), 0L)
    expect_identical(attr(map, "gamma_dim"), dim(gamma))
    expect_identical(attr(map, "gamma_dimnames"), dimnames(gamma))
    expect_equal(
      unname(drop(map %*% coef)),
      unname(as.vector(gamma)),
      tolerance = 1e-12
    )
  }
})

test_that("absorbing states propagate probability and derivative mass", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  case <- delta_vglm_case()
  analytical <- run_sop_delta_plan(case$plan, case$model)
  absorbing_state <- length(case$y_levels)

  absorbing_probability <- analytical$probabilities[,, absorbing_state]
  expect_true(all(apply(absorbing_probability, 1L, diff) >= -1e-12))
  expect_gt(max(abs(analytical$jacobian[, -1L, absorbing_state, ])), 0)
  expect_lt(
    max(abs(apply(analytical$jacobian, c(1L, 2L, 4L), sum))),
    1e-10
  )

  wrapper <- structure(
    list(
      coefficients = get_coef(case$model),
      vglm_fit = case$model
    ),
    class = "robcov_vglm"
  )
  robust_result <- run_sop_delta_plan(case$plan, wrapper)
  expect_equal(
    unname(robust_result$jacobian),
    unname(analytical$jacobian),
    tolerance = 0
  )
})

test_that("factor visit designs use the same analytical recursion", {
  skip_if_not_installed("VGAM")
  data <- make_delta_factor_visit_case()
  model <- VGAM::vglm(
    y ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  baseline <- data[!duplicated(data$id), , drop = FALSE][1:5, , drop = FALSE]
  plan <- compile_sop_execution_plan(
    model = model,
    newdata = baseline,
    times = levels(data$time),
    y_levels = seq_len(4),
    absorb = 4
  )

  analytical <- run_sop_delta_plan(plan, model)
  oracle <- run_sop_delta_r_oracle(plan, model)
  point <- run_sop_execution_plan(plan, get_effective_coefs(model))

  expect_equal(
    unname(analytical$probabilities),
    unname(oracle$probabilities),
    tolerance = 1e-12
  )
  expect_equal(
    unname(analytical$jacobian),
    unname(oracle$jacobian),
    tolerance = 1e-12
  )
  expect_lt(
    max(abs(unname(analytical$probabilities) - unname(point))),
    1e-12
  )
  expect_lt(
    max(abs(apply(analytical$jacobian, c(1L, 2L, 4L), sum))),
    1e-10
  )
})

test_that("native recursion reproduces the R analytic oracle", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  for (case in list(delta_vglm_case(), delta_orm_case())) {
    native <- run_sop_delta_plan(case$plan, case$model)
    oracle <- run_sop_delta_r_oracle(case$plan, case$model)

    expect_equal(
      unname(native$probabilities),
      unname(oracle$probabilities),
      tolerance = 1e-12
    )
    expect_equal(
      unname(native$jacobian),
      unname(oracle$jacobian),
      tolerance = 1e-12
    )
  }
})

test_that("grouped recursion retains only the required analytical state", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  case <- delta_vglm_case()
  plan <- case$plan
  plan$components$X_init <- plan$components$X_init[
    rep(seq_len(4L), 2L),
    ,
    drop = FALSE
  ]
  plan$components$n_pat <- 8L
  plan$components$X_transition <- lapply(
    plan$components$X_transition,
    function(design) {
      if (is.null(design)) {
        return(NULL)
      }
      origin_count <- nrow(design) / 4L
      rows <- unlist(lapply(seq_len(origin_count), function(origin) {
        source <- (origin - 1L) * 4L + seq_len(4L)
        c(source, source)
      }))
      design[rows, , drop = FALSE]
    }
  )

  full <- run_sop_delta_plan(plan, case$model)
  grouped <- run_sop_delta_plan(
    plan,
    case$model,
    average_group_size = 4L
  )
  population <- run_sop_delta_plan(
    plan,
    case$model,
    average_group_size = 4L,
    retain_individual_probabilities = TRUE
  )

  expected_probability <- apply(full$probabilities, c(2L, 3L), mean)
  expected_jacobian <- apply(full$jacobian, c(2L, 3L, 4L), mean)
  expect_identical(
    dim(grouped$probabilities),
    c(2L, 5L, length(case$y_levels))
  )
  expect_identical(
    dim(grouped$jacobian),
    c(2L, 5L, length(case$y_levels), length(get_coef(case$model)))
  )
  expect_null(grouped$individual_probabilities)
  expect_equal(
    unname(grouped$probabilities[1L, , ]),
    unname(expected_probability),
    tolerance = 1e-15
  )
  expect_equal(
    unname(grouped$jacobian[1L, , , ]),
    unname(expected_jacobian),
    tolerance = 1e-15
  )
  expect_equal(
    unname(grouped$probabilities[2L, , ]),
    unname(expected_probability),
    tolerance = 1e-15
  )
  expect_identical(
    dim(population$individual_probabilities),
    dim(full$probabilities)
  )
  expect_equal(
    unname(population$individual_probabilities),
    unname(full$probabilities),
    tolerance = 0
  )

  individual_bytes <- sop_delta_preflight(plan, length(get_coef(case$model)))
  grouped_bytes <- sop_delta_preflight(
    plan,
    length(get_coef(case$model)),
    average_group_size = 4L
  )
  population_bytes <- sop_delta_preflight(
    plan,
    length(get_coef(case$model)),
    average_group_size = 4L,
    retain_individual_probabilities = TRUE
  )
  expect_lt(grouped_bytes, population_bytes)
  expect_lt(population_bytes, individual_bytes)
})

test_that("analytical SOP validation rejects unsupported model structures", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  local_reproducible_output(width = 80)
  case <- delta_vglm_case()

  second_order <- case$plan
  second_order$recursion_order <- 2L
  expect_snapshot(
    run_sop_delta_plan(second_order, case$model),
    error = TRUE
  )

  unsupported_link <- case$plan
  unsupported_link$link <- "probit"
  expect_snapshot(
    run_sop_delta_plan(unsupported_link, case$model),
    error = TRUE
  )

  unsupported_model <- stats::lm(mpg ~ wt, data = mtcars)
  expect_snapshot(
    run_sop_delta_plan(case$plan, unsupported_model),
    error = TRUE
  )

  partial_model <- VGAM::vglm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
    data = case$data
  )
  partial_plan <- compile_sop_execution_plan(
    model = partial_model,
    newdata = case$baseline,
    times = 1:3,
    y_levels = seq_len(nrow(get_effective_coefs(partial_model)) + 1L),
    absorb = max(case$y_levels)
  )
  expect_snapshot(
    run_sop_delta_plan(partial_plan, partial_model),
    error = TRUE
  )
})

test_that("crossed raw ordinal probabilities are rejected without clipping", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  local_reproducible_output(width = 80)
  case <- delta_vglm_case()
  constraints <- VGAM::constraints(case$model)
  intercept <- grep("Intercept", names(constraints), fixed = TRUE)
  preceding <- if (intercept == 1L) {
    0L
  } else {
    sum(vapply(constraints[seq_len(intercept - 1L)], ncol, integer(1)))
  }
  intercept_coef <- preceding + seq_len(ncol(constraints[[intercept]]))
  coef <- get_coef(case$model)
  coef[intercept_coef] <- seq(-4, 4, length.out = length(intercept_coef))
  crossed_model <- set_coef(case$model, coef)

  condition <- rlang::catch_cnd(
    run_sop_delta_plan(case$plan, crossed_model)
  )
  expect_s3_class(condition, "markov_misc_delta_crossed_probabilities")
  expect_lt(condition$minimum_probability, 0)
  expect_snapshot(
    run_sop_delta_plan(case$plan, crossed_model),
    error = TRUE
  )
})

test_that("analytical allocation guard is typed and configurable", {
  local_reproducible_output(width = 80)
  withr::local_options(markov.misc.delta_max_bytes = 1)

  condition <- rlang::catch_cnd(
    delta_assert_bytes(16, "Test allocation")
  )
  expect_s3_class(condition, "markov_misc_delta_too_large")
  expect_identical(condition$required_bytes, 16)
  expect_identical(condition$max_bytes, 1)
  expect_snapshot(
    delta_assert_bytes(16, "Test allocation"),
    error = TRUE
  )
})
