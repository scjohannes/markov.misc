make_factor_visit_case <- function(n_patients = 80, n_visits = 4, seed = 2026) {
  set.seed(seed)
  visits <- as.character(seq_len(n_visits))
  ids <- seq_len(n_patients)
  tx <- stats::rbinom(n_patients, 1, 0.5)

  data <- expand.grid(
    id = ids,
    time = visits,
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  data$time <- factor(data$time, levels = visits)
  data$tx <- tx[data$id]
  data$yprev <- factor(
    sample(seq_len(4), nrow(data), replace = TRUE),
    levels = as.character(seq_len(4))
  )
  eta <- -0.35 * data$tx + 0.25 * as.integer(data$time) +
    0.2 * as.integer(as.character(data$yprev))
  p_worse <- stats::plogis(eta - mean(eta))
  y_num <- pmin(
    4L,
    pmax(1L, 1L + stats::rbinom(nrow(data), 3L, p_worse))
  )
  data$y <- ordered(y_num, levels = as.character(seq_len(4)))
  data
}

make_interpolation_case <- function() {
  x <- expand.grid(
    time = factor(c("v1", "v2"), levels = c("v1", "v2")),
    state = factor(c("1", "2"), levels = c("1", "2")),
    tx = 0:1,
    KEEP.OUT.ATTRS = FALSE
  )
  x$estimate <- mapply(
    function(time, state, tx) {
      if (tx == 0 && time == "v1" && state == "1") return(0.8)
      if (tx == 0 && time == "v1" && state == "2") return(0.2)
      if (tx == 0 && time == "v2" && state == "1") return(0.4)
      if (tx == 0 && time == "v2" && state == "2") return(0.6)
      if (tx == 1 && time == "v1" && state == "1") return(0.7)
      if (tx == 1 && time == "v1" && state == "2") return(0.3)
      if (tx == 1 && time == "v2" && state == "1") return(0.5)
      0.5
    },
    as.character(x$time),
    as.character(x$state),
    x$tx
  )

  class(x) <- c("markov_avg_sops", class(x))
  attr(x, "pvarname") <- "yprev"
  attr(x, "ylevels") <- factor(1:2)
  attr(x, "avg_args") <- list(
    variables = list(tx = c(0, 1)),
    by = NULL,
    times = factor(c("v1", "v2"), levels = c("v1", "v2")),
    id_var = "id"
  )
  attr(x, "newdata_orig") <- data.frame(
    id = seq_len(4),
    yprev = factor(c("1", "2", "2", "2"), levels = c("1", "2")),
    tx = c(0, 1, 0, 1)
  )
  x
}

add_interpolation_draws <- function(x) {
  draws <- x[rep(seq_len(nrow(x)), each = 3L), , drop = FALSE]
  draws$draw_id <- rep(seq_len(3), times = nrow(x))
  offset <- c(-0.1, 0, 0.1)[draws$draw_id]
  draws$estimate <- draws$estimate + ifelse(draws$state == "1", offset, -offset)
  draws$estimate <- pmin(pmax(draws$estimate, 0), 1)
  rownames(draws) <- NULL

  attr(x, "simulation_draws") <- draws
  attr(x, "conf_level") <- 0.8
  attr(x, "conf_type") <- "perc"
  x
}

add_interpolation_intervals <- function(x) {
  x$conf.low <- pmax(x$estimate - 0.05, 0)
  x$conf.high <- pmin(x$estimate + 0.05, 1)
  x$std.error <- 0.025
  x
}

test_that("factor visit time works through SOP APIs and the fast path", {
  skip_if_not_installed("VGAM")

  data <- make_factor_visit_case()
  fit <- VGAM::vglm(
    y ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  baseline <- data[!duplicated(data$id), , drop = FALSE]

  direct <- soprob_markov(
    fit,
    data = baseline,
    times = 1:4,
    ylevels = factor(1:4),
    pvarname = "yprev"
  )
  expect_equal(dim(direct), c(nrow(baseline), 4L, 4L))
  expect_equal(dimnames(direct)[[2]], as.character(1:4))

  individual <- sops(
    fit,
    newdata = baseline,
    times = NULL,
    ylevels = factor(1:4),
    pvarname = "yprev"
  )
  expect_equal(sort(unique(as.character(individual$time))), as.character(1:4))

  averaged <- avg_sops(
    fit,
    newdata = baseline,
    variables = list(tx = c(0, 1)),
    times = NULL,
    ylevels = factor(1:4),
    pvarname = "yprev"
  )
  expect_equal(as.character(attr(averaged, "avg_args")$times), as.character(1:4))

  components <- markov.misc:::markov_msm_build(
    model = fit,
    data = baseline,
    times = 1:4,
    ylevels = factor(1:4),
    pvarname = "yprev"
  )
  fast <- markov.misc:::markov_msm_run(
    components,
    get_effective_coefs(fit),
    times = attr(averaged, "avg_args")$times
  )
  expect_equal(unname(fast), unname(direct), tolerance = 1e-10)
})

test_that("factor visit gap requires explicit numeric t_covs", {
  skip_if_not_installed("VGAM")

  data <- make_factor_visit_case(n_patients = 40)
  fit <- VGAM::vglm(
    y ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  baseline <- data[!duplicated(data$id), , drop = FALSE]

  expect_error(
    soprob_markov(
      fit,
      data = baseline,
      times = 1:4,
      ylevels = factor(1:4),
      pvarname = "yprev",
      gap = "gap"
    ),
    "Factor visit time with `gap` requires numeric gap values"
  )

  out <- soprob_markov(
    fit,
    data = baseline,
    times = 1:4,
    ylevels = factor(1:4),
    pvarname = "yprev",
    gap = "gap",
    t_covs = data.frame(gap = c(3, 4, 7, 14))
  )
  expect_equal(dim(out), c(nrow(baseline), 4L, 4L))
})

test_that("factor visit simulation and bootstrap inference smoke-test", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("mvtnorm")

  data <- make_factor_visit_case(n_patients = 60, seed = 2027)
  fit <- VGAM::vglm(
    y ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  baseline <- data[!duplicated(data$id), , drop = FALSE]

  avg_baseline <- avg_sops(
    fit,
    newdata = baseline,
    variables = list(tx = c(0, 1)),
    times = NULL,
    ylevels = factor(1:4),
    pvarname = "yprev"
  )
  sim <- inferences(avg_baseline, n_sim = 2, return_draws = TRUE)
  expect_s3_class(sim, "markov_avg_sops")
  expect_false(is.null(attr(sim, "simulation_draws")))

  avg_full <- avg_sops(
    fit,
    newdata = data,
    variables = list(tx = c(0, 1)),
    times = NULL,
    ylevels = factor(1:4),
    pvarname = "yprev"
  )
  boot <- inferences(avg_full, method = "bootstrap", n_sim = 1, return_draws = TRUE)
  expect_s3_class(boot, "markov_avg_sops")
  expect_false(is.null(attr(boot, "bootstrap_draws")))
})

test_that("interpolate_sops maps visits, anchors baseline, and normalizes states", {
  x <- make_interpolation_case()
  out <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    xout = c(0, 3, 5, 7),
    origin_time = 0
  )

  mid <- out[out$tx == 0 & out$time == 5 & out$state == "1", ]
  expect_equal(mid$estimate, 0.6, tolerance = 1e-10)

  origin <- out[out$time == 0, c("tx", "state", "estimate")]
  expect_equal(
    origin[origin$tx == 0, c("state", "estimate")],
    origin[origin$tx == 1, c("state", "estimate")],
    ignore_attr = TRUE
  )
  expect_equal(origin$estimate[origin$state == "1" & origin$tx == 0], 0.25)
  expect_equal(origin$estimate[origin$state == "2" & origin$tx == 0], 0.75)

  sums <- stats::aggregate(estimate ~ tx + time, out, sum)
  expect_equal(sums$estimate, rep(1, nrow(sums)), tolerance = 1e-10)
})

test_that("interpolate_sops interpolates stored draws from deterministic origin anchor", {
  x <- add_interpolation_draws(make_interpolation_case())
  out <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    xout = 0:7,
    origin_time = 0
  )
  draws <- attr(out, "simulation_draws")

  expect_s3_class(draws, "data.frame")
  expect_equal(sort(unique(draws$time)), 0:7)

  origin <- out[out$time == 0 & out$tx == 0 & out$state == "1", ]
  expect_equal(origin$estimate, 0.25, tolerance = 1e-10)
  expect_equal(origin$conf.low, 0.25, tolerance = 1e-10)
  expect_equal(origin$conf.high, 0.25, tolerance = 1e-10)
  expect_equal(origin$std.error, 0, tolerance = 1e-10)

  day_1 <- out[out$time == 1 & out$tx == 0 & out$state == "1", ]
  expect_false(is.na(day_1$conf.low))
  expect_false(is.na(day_1$conf.high))
  expect_gt(day_1$conf.high, day_1$conf.low)

  day_1_draws <- draws[draws$time == 1 & draws$tx == 0 & draws$state == "1", ]
  expect_equal(
    day_1$conf.low,
    unname(stats::quantile(day_1_draws$estimate, probs = 0.1)),
    tolerance = 1e-10
  )
  expect_equal(
    day_1$conf.high,
    unname(stats::quantile(day_1_draws$estimate, probs = 0.9)),
    tolerance = 1e-10
  )
})

test_that("interpolate_sops warns when interval-bearing objects do not store draws", {
  x <- add_interpolation_intervals(make_interpolation_case())
  expect_warning(
    out <- interpolate_sops(
      x,
      time_map = c(v1 = 3, v2 = 7),
      xout = 1:7,
      origin_time = 0
    ),
    "no stored draws"
  )

  day_1 <- out[out$time == 1 & out$tx == 0 & out$state == "1", ]
  expect_false(is.na(day_1$conf.low))
  expect_false(is.na(day_1$conf.high))

  attr(x, "method") <- "posterior"
  expect_warning(
    interpolate_sops(
      x,
      time_map = c(v1 = 3, v2 = 7),
      xout = 1:7,
      origin_time = 0
    ),
    "posterior draws"
  )
})

test_that("interpolate_sops handles blrm-style stored posterior draws", {
  x <- add_interpolation_intervals(make_interpolation_case())
  draws <- attr(add_interpolation_draws(make_interpolation_case()), "simulation_draws")
  attr(x, "draws") <- draws
  attr(x, "method") <- "posterior"
  attr(x, "conf_level") <- 0.8

  expect_warning(
    out <- interpolate_sops(
      x,
      time_map = c(v1 = 3, v2 = 7),
      xout = 0:7,
      origin_time = 0
    ),
    NA
  )

  posterior_draws <- attr(out, "draws")
  expect_s3_class(posterior_draws, "data.frame")
  expect_equal(sort(unique(posterior_draws$time)), 0:7)

  origin <- out[out$time == 0 & out$tx == 0 & out$state == "1", ]
  expect_equal(origin$conf.low, 0.25, tolerance = 1e-10)
  expect_equal(origin$conf.high, 0.25, tolerance = 1e-10)

  day_1 <- out[out$time == 1 & out$tx == 0 & out$state == "1", ]
  day_1_draws <- posterior_draws[
    posterior_draws$time == 1 &
      posterior_draws$tx == 0 &
      posterior_draws$state == "1",
  ]
  expect_equal(
    day_1$conf.low,
    unname(stats::quantile(day_1_draws$estimate, probs = 0.1)),
    tolerance = 1e-10
  )
  expect_equal(
    day_1$conf.high,
    unname(stats::quantile(day_1_draws$estimate, probs = 0.9)),
    tolerance = 1e-10
  )
})

test_that("plot_sops bars work with derived facet labels on interpolated draws", {
  x <- add_interpolation_draws(make_interpolation_case())
  out <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    xout = 1:7,
    origin_time = 0
  )
  label_tx <- function(x) {
    x <- as.character(x)
    factor(
      ifelse(x %in% c("1", "Treatment"), "Treatment", "Control"),
      levels = c("Control", "Treatment")
    )
  }
  out$tx_label <- label_tx(out$tx)
  draws <- attr(out, "simulation_draws")
  draws$tx_label <- label_tx(draws$tx)
  attr(out, "simulation_draws") <- draws

  expect_equal(anyDuplicated(out[c("time", "state", "tx_label")]), 0L)
  expect_s3_class(plot_sops(out, geom = "bar", facet_var = "tx_label"), "ggplot")
})

test_that("interpolate_sops guardrails are clear", {
  x <- make_interpolation_case()
  expect_error(
    interpolate_sops(x, time_map = c(v1 = 3)),
    "`time_map` is missing entries"
  )
  expect_error(
    interpolate_sops(x, time_map = c(v1 = 3, v2 = 7), xout = 8),
    "`xout` must stay within the supported time range"
  )
})

test_that("time_in_state uses trapezoidal AUC on mapped real time", {
  x <- make_interpolation_case()

  auc <- time_in_state(
    x,
    target_states = "1",
    time_map = c(v1 = 3, v2 = 7),
    origin_time = 0
  )
  auc <- auc[order(auc$tx), , drop = FALSE]

  expect_equal(auc$total_time[auc$tx == 0], 3.975, tolerance = 1e-10)
  expect_equal(auc$total_time[auc$tx == 1], 3.825, tolerance = 1e-10)

  day_1_auc <- time_in_state(
    x,
    target_states = "1",
    time_map = c(v1 = 3, v2 = 7),
    origin_time = 0,
    xout = 1:7
  )
  day_1_auc <- day_1_auc[order(day_1_auc$tx), , drop = FALSE]

  expect_equal(day_1_auc$total_time[day_1_auc$tx == 0], 3.633333, tolerance = 1e-6)
  expect_equal(day_1_auc$total_time[day_1_auc$tx == 1], 3.5, tolerance = 1e-10)

  interpolated <- interpolate_sops(
    x,
    time_map = c(v1 = 3, v2 = 7),
    xout = 1:7,
    origin_time = 0
  )
  interpolated_auc <- time_in_state(interpolated, target_states = "1")
  interpolated_auc <- interpolated_auc[order(interpolated_auc$tx), , drop = FALSE]

  expect_equal(interpolated_auc, day_1_auc)
})
