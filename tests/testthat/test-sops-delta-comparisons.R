make_delta_comparison_factor_data <- function(
  n_patients = 45,
  n_visits = 4,
  seed = 8421
) {
  set.seed(seed)
  visits <- paste0("v", seq_len(n_visits))
  data <- expand.grid(
    id = seq_len(n_patients),
    time = visits,
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  data$time <- factor(data$time, levels = visits)
  treatment <- stats::rbinom(n_patients, 1, 0.5)
  sex <- sample(c("F", "M"), n_patients, replace = TRUE)
  data$tx <- treatment[data$id]
  data$sex <- factor(sex[data$id])
  data$yprev <- factor(
    sample(seq_len(4), nrow(data), replace = TRUE),
    levels = as.character(seq_len(4))
  )
  data$yprev2 <- factor(
    sample(seq_len(4), nrow(data), replace = TRUE),
    levels = as.character(seq_len(4))
  )
  eta <-
    -0.4 *
    data$tx +
    0.22 * as.integer(data$time) +
    0.18 * as.integer(as.character(data$yprev))
  probability <- stats::plogis(eta - mean(eta))
  outcome <- pmin(
    4L,
    pmax(1L, 1L + stats::rbinom(nrow(data), 3L, probability))
  )
  data$y <- ordered(outcome, levels = as.character(seq_len(4)))
  data
}

delta_comparison_factor_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(value)) {
      data <- make_delta_comparison_factor_data()
      old_options <- options(datadist = NULL)
      on.exit(options(old_options), add = TRUE)
      fit <- suppressWarnings(orm_markov(
        y ~ time + tx + yprev,
        data = data,
        id_var = "id",
        opt_method = "LM",
        scale = TRUE,
        maxit = 100
      ))
      visits <- levels(data$time)
      value <<- list(
        data = data,
        fit = fit,
        visits = visits,
        time_map = stats::setNames(c(3, 7, 14, 21), visits),
        target_times = 0:21,
        y_levels = fit$yunique,
        covariance = fit$var
      )
    }
    value
  }
})

delta_factor_time_point <- function(model, case) {
  avg_comparisons(
    model,
    variables = list(tx = c(0, 1)),
    estimand = "time_in_state",
    state_sets = list(lower = c("1", "2")),
    comparison = "difference",
    times = case$visits,
    y_levels = case$y_levels,
    time_map = case$time_map,
    origin_time = 0,
    target_times = case$target_times,
    origin = "empirical_baseline",
    time_unit = "days"
  )
}

central_delta_comparison_jacobian <- function(model, point_function) {
  coefficient <- get_coef(model)
  point <- point_function(model)$estimate
  jacobian <- matrix(
    NA_real_,
    nrow = length(point),
    ncol = length(coefficient),
    dimnames = list(NULL, names(coefficient))
  )
  for (j in seq_along(coefficient)) {
    step <- .Machine$double.eps^(1 / 3) * max(1, abs(coefficient[j]))
    upper <- lower <- coefficient
    upper[j] <- upper[j] + step
    lower[j] <- lower[j] - step
    jacobian[, j] <- (point_function(set_coef(model, upper))$estimate -
      point_function(set_coef(model, lower))$estimate) /
      (2 * step)
  }
  jacobian
}

manual_factor_time_auc <- function(avg, case, states, level) {
  baseline <- attr(avg, "newdata_orig")
  baseline <- baseline[!duplicated(baseline$id), , drop = FALSE]
  anchor <- mean(as.character(baseline$yprev) %in% as.character(states))
  keep <-
    as.character(avg$state) %in%
    as.character(states) &
    as.character(avg$tx) == as.character(level)
  visits <- stats::aggregate(
    estimate ~ time,
    data = as.data.frame(avg)[keep, , drop = FALSE],
    FUN = sum
  )
  visits <- visits[
    match(case$visits, as.character(visits$time)),
    ,
    drop = FALSE
  ]
  interpolated <- stats::approx(
    x = c(0, unname(case$time_map)),
    y = c(anchor, visits$estimate),
    xout = case$target_times,
    rule = 1,
    ties = "ordered"
  )$y
  list(
    anchor = anchor,
    auc = sum(
      diff(case$target_times) *
        (utils::head(interpolated, -1L) + utils::tail(interpolated, -1L)) /
        2
    )
  )
}

test_that("joint average SOP Jacobians retain cross-arm covariance", {
  case <- delta_comparison_factor_case()
  variables <- list(tx = c(0, 1, 2))
  state_sets <- list(lower = c("1", "2"))
  times <- case$visits[1:2]

  avg <- avg_sops(
    case$fit,
    variables = variables,
    times = times,
    y_levels = case$y_levels
  )
  avg <- inferences(
    avg,
    method = "delta",
    target = "empirical",
    vcov = case$covariance,
    conf_type = "wald"
  )
  comparison <- avg_comparisons(
    case$fit,
    variables = variables,
    estimand = "sop",
    state_sets = state_sets,
    comparison = "difference",
    times = times,
    y_levels = case$y_levels
  )
  comparison <- inferences(
    comparison,
    method = "delta",
    target = "empirical",
    vcov = case$covariance
  )

  average_jacobian <- get_jacobian(avg)
  comparison_jacobian <- get_jacobian(comparison)
  expected <- matrix(
    0,
    nrow = nrow(comparison),
    ncol = ncol(average_jacobian),
    dimnames = list(NULL, colnames(average_jacobian))
  )
  arm_jacobians <- vector("list", nrow(comparison))
  for (i in seq_len(nrow(comparison))) {
    common <-
      as.character(avg$time) == as.character(comparison$time[i]) &
      as.character(avg$state) %in% state_sets[[comparison$state_set[i]]]
    high <- common &
      as.character(avg$tx) == as.character(comparison$comparison_level[i])
    reference <- common &
      as.character(avg$tx) == as.character(comparison$reference_level[i])
    high_jacobian <- colSums(average_jacobian[high, , drop = FALSE])
    reference_jacobian <- colSums(
      average_jacobian[reference, , drop = FALSE]
    )
    expected[i, ] <- high_jacobian - reference_jacobian
    arm_jacobians[[i]] <- list(
      high = high_jacobian,
      reference = reference_jacobian
    )
  }

  discrepancy <- max(abs(comparison_jacobian - expected))
  expect_lt(discrepancy, 1e-12)

  covariance <- case$covariance[
    colnames(comparison_jacobian),
    colnames(comparison_jacobian),
    drop = FALSE
  ]
  correct <- naive <- numeric(nrow(comparison))
  cross_covariance <- numeric(nrow(comparison))
  for (i in seq_len(nrow(comparison))) {
    high <- arm_jacobians[[i]]$high
    reference <- arm_jacobians[[i]]$reference
    correct[i] <- drop((high - reference) %*% covariance %*% (high - reference))
    naive[i] <- drop(high %*% covariance %*% high) +
      drop(reference %*% covariance %*% reference)
    cross_covariance[i] <- drop(high %*% covariance %*% reference)
  }
  expect_equal(correct, comparison$std.error^2, tolerance = 1e-12)
  expect_gt(max(abs(cross_covariance)), 1e-8)
  expect_gt(max(abs(correct - naive)), 1e-8)
})

test_that("factor-time ORM differences match public finite differences", {
  case <- delta_comparison_factor_case()
  point <- delta_factor_time_point(case$fit, case)
  inferred <- inferences(
    point,
    method = "delta",
    target = "empirical",
    vcov = case$covariance
  )

  numeric_jacobian <- central_delta_comparison_jacobian(
    case$fit,
    function(model) delta_factor_time_point(model, case)
  )
  analytical_jacobian <- get_jacobian(inferred)
  discrepancy <- max(abs(analytical_jacobian - numeric_jacobian))
  expect_lt(discrepancy, 1e-7)

  avg <- avg_sops(
    case$fit,
    variables = list(tx = c(0, 1)),
    times = case$visits,
    y_levels = case$y_levels
  )
  operator <- delta_comparison_operator(
    point,
    avg,
    attr(point, "comparison_args"),
    attr(point, "avg_args")
  )
  expect_lt(
    max(abs(drop(operator %*% avg$estimate) - point$estimate)),
    1e-12
  )

  stale <- point
  stale$estimate[1] <- stale$estimate[1] + 5e-12
  condition <- tryCatch(
    delta_comparison_operator(
      stale,
      avg,
      attr(stale, "comparison_args"),
      attr(stale, "avg_args")
    ),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(
    conditionMessage(condition),
    "did not reproduce the stored point estimates",
    fixed = TRUE
  )

  reference <- manual_factor_time_auc(avg, case, c("1", "2"), 0)
  treatment <- manual_factor_time_auc(avg, case, c("1", "2"), 1)
  expect_equal(treatment$anchor - reference$anchor, 0)
  expect_equal(treatment$auc - reference$auc, point$estimate, tolerance = 1e-12)

  covariance <- case$covariance[
    colnames(analytical_jacobian),
    colnames(analytical_jacobian),
    drop = FALSE
  ]
  expected_se <- sqrt(drop(
    analytical_jacobian %*% covariance %*% t(analytical_jacobian)
  ))
  expect_equal(inferred$std.error, expected_se, tolerance = 1e-12)
})

test_that("superpopulation comparisons retain transformed stacked influence", {
  case <- delta_comparison_factor_case()
  point <- avg_comparisons(
    case$fit,
    variables = list(tx = c(0, 1)),
    estimand = "sop",
    state_sets = list(lower = c("1", "2")),
    comparison = "difference",
    times = case$visits[1:3],
    y_levels = case$y_levels
  )
  inferred <- inferences(
    point,
    method = "delta",
    target = "superpopulation"
  )
  analytical <- attr(inferred, "analytical")
  rows <- c(1L, nrow(inferred))
  expected <- stats::cov(analytical$influence[, rows, drop = FALSE]) /
    nrow(analytical$influence)

  expect_identical(analytical$representation, "influence")
  expect_equal(
    unname(stats::vcov(inferred, rows = rows)),
    unname(expected),
    tolerance = 1e-12
  )
  expect_equal(
    inferred$std.error[rows]^2,
    unname(diag(expected)),
    tolerance = 1e-12
  )

  expect_snapshot(
    inferences(
      point,
      method = "delta",
      target = "superpopulation",
      vcov = case$covariance
    ),
    error = TRUE
  )
})

test_that("analytical comparison scope rejects deferred estimands", {
  case <- delta_comparison_factor_case()
  local_reproducible_output(width = 80)
  base <- avg_comparisons(
    case$fit,
    variables = list(tx = c(0, 1)),
    estimand = "sop",
    state_sets = "1",
    comparison = "difference",
    times = case$visits[1:2],
    y_levels = case$y_levels
  )

  ratio <- base
  ratio$comparison <- "ratio"
  ratio_args <- attr(ratio, "comparison_args")
  ratio_args$comparison <- "ratio"
  attr(ratio, "comparison_args") <- ratio_args
  expect_snapshot(inferences(ratio, method = "delta"), error = TRUE)

  time_benefit <- base
  time_benefit$estimand <- "time_benefit"
  time_benefit_args <- attr(time_benefit, "comparison_args")
  time_benefit_args$estimand <- "time_benefit"
  attr(time_benefit, "comparison_args") <- time_benefit_args
  expect_snapshot(inferences(time_benefit, method = "delta"), error = TRUE)

  stratified <- base
  stratified_args <- attr(stratified, "avg_args")
  stratified_args$by <- "sex"
  attr(stratified, "avg_args") <- stratified_args
  expect_snapshot(inferences(stratified, method = "delta"), error = TRUE)

  second_order <- base
  attr(second_order, "p2_var") <- "yprev2"
  expect_snapshot(
    inferences(
      second_order,
      method = "delta",
      vcov = case$covariance
    ),
    error = TRUE
  )
})

test_that("analytical comparisons reject partial proportional odds", {
  skip_if_not_installed("VGAM")
  case <- delta_comparison_factor_case()
  local_reproducible_output(width = 80)
  partial_model <- VGAM::vglm(
    y ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
    data = case$data
  )
  baseline <- case$data[!duplicated(case$data$id), , drop = FALSE]
  partial <- avg_comparisons(
    partial_model,
    newdata = baseline,
    variables = list(tx = c(0, 1)),
    estimand = "sop",
    state_sets = "1",
    comparison = "difference",
    times = case$visits[1:2],
    y_levels = case$y_levels
  )
  expect_snapshot(
    inferences(
      partial,
      method = "delta",
      vcov = stats::vcov(partial_model)
    ),
    error = TRUE
  )
})
