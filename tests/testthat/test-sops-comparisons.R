make_mock_avg_sops <- function(by = NULL, tx_levels = c(0, 1, 2)) {
  if (is.null(by)) {
    out <- expand.grid(
      time = 1:2,
      state = c("1", "2"),
      tx = tx_levels,
      KEEP.OUT.ATTRS = FALSE
    )
    values <- c(
      "0" = 0.2,
      "1" = 0.5,
      "2" = 0.1
    )
    out$estimate <- values[as.character(out$tx)]
    out$estimate[out$time == 2 & out$tx == 0] <- 0.4
    out$estimate[out$time == 2 & out$tx == 1] <- 0.7
    out$estimate[out$time == 2 & out$tx == 2] <- 0.2
    out$estimate[out$state == "2"] <- 1 - out$estimate[out$state == "2"]
  } else {
    out <- expand.grid(
      time = 1,
      state = c("1", "2"),
      tx = c(0, 1),
      sex = c("f", "m"),
      KEEP.OUT.ATTRS = FALSE
    )
    out$estimate <- c(0.2, 0.8, 0.6, 0.4, 0.3, 0.7, 0.9, 0.1)
  }
  class(out) <- c("markov_avg_sops", class(out))
  attr(out, "avg_args") <- list(
    variables = list(tx = tx_levels),
    by = by,
    times = sort(unique(out$time)),
    id_var = "id"
  )
  attr(out, "y_levels") <- c("1", "2")
  attr(out, "call_args") <- list(times = sort(unique(out$time)))
  attr(out, "time_var") <- "time"
  attr(out, "p_var") <- "yprev"
  attr(out, "newdata_orig") <- data.frame(
    id = 1:2,
    yprev = c("1", "2")
  )
  attr(out, "newdata_supplied") <- TRUE
  out
}

test_that("avg_comparisons() computes SOP differences and multiple levels", {
  avg <- make_mock_avg_sops()

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1, 2)),
        estimand = "sop",
        state_sets = "1",
        times = 1:2,
        comparison = "difference"
      )
    }
  )

  out <- out[order(out$comparison_level, out$time), ]
  rownames(out) <- NULL

  expect_s3_class(out, "markov_avg_comparisons")
  expect_equal(out$estimate, c(0.3, 0.3, -0.1, -0.2))
  expect_equal(out$reference_level, c(0, 0, 0, 0))
  expect_equal(out$comparison_level, c(1, 1, 2, 2))
})

test_that("avg_comparisons() computes time-in-state ratios for lumped states", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        comparison = "ratio"
      )
    }
  )

  expect_equal(out$estimate, (0.5 + 0.7) / (0.2 + 0.4))
  expect_equal(out$state_set, "1")
  expect_equal(out$contrast, "1 / 0")
})

test_that("avg_comparisons() preserves by strata", {
  avg <- make_mock_avg_sops(by = "sex", tx_levels = c(0, 1))

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        estimand = "sop",
        state_sets = "1",
        times = 1,
        by = "sex"
      )
    }
  )

  out <- out[order(out$sex), ]
  rownames(out) <- NULL

  expect_equal(as.character(out$sex), c("f", "m"))
  expect_equal(out$estimate, c(0.4, 0.6))
})

test_that("avg_comparisons() uses real-time AUC for time-in-state", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        time_map = c("1" = 1, "2" = 3)
      )
    }
  )

  expect_equal(out$estimate, 0.6)
})

test_that("real-time comparisons use baseline support without integrating it", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      explicit <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        time_map = c("1" = 7, "2" = 14),
        baseline_time = 0,
        target_times = 1:14
      )
      disabled <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        time_map = c("1" = 7, "2" = 14),
        baseline_time = NULL
      )
    }
  )

  difference <- stats::approx(
    x = c(0, 7, 14),
    y = c(0, 0.3, 0.3),
    xout = 1:14
  )$y
  expected <- sum(
    utils::head(difference, -1L) + utils::tail(difference, -1L)
  ) /
    2
  expect_equal(explicit$estimate, expected)
  expect_equal(disabled$estimate, 2.1)
})

test_that("avg_comparisons() forwards model-structure args for marginal metrics", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))
  time_covariates <- data.frame(visit = 1:2, spline = c(0, 1))
  captured <- new.env(parent = emptyenv())

  with_mocked_bindings(
    avg_sops = function(...) {
      captured$args <- list(...)
      avg
    },
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        estimand = "time_in_state",
        state_sets = "1",
        times = 1:2,
        y_levels = c("1", "2"),
        absorb = "2",
        time_var = "visit",
        p_var = "prev",
        p2_var = "preprev",
        gap_var = "gap",
        time_covariates = time_covariates
      )
    }
  )

  expect_s3_class(out, "markov_avg_comparisons")
  expect_equal(captured$args$y_levels, c("1", "2"))
  expect_equal(captured$args$absorb, "2")
  expect_equal(captured$args$time_var, "visit")
  expect_equal(captured$args$p_var, "prev")
  expect_equal(captured$args$p2_var, "preprev")
  expect_equal(captured$args$gap, "gap")
  expect_equal(captured$args$time_covariates, time_covariates)
})

test_that("avg_comparisons() computes time benefit from paired profiles", {
  newdata <- data.frame(
    id = 1:2,
    tx = c(0, 1),
    yprev = c(2, 2),
    prev = c(2, 2),
    preprev = c(1, 1)
  )
  time_covariates <- data.frame(visit = 1, spline = 0)
  captured <- new.env(parent = emptyenv())

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    sops = function(model, newdata, times, ...) {
      captured$args <- list(...)
      out <- expand.grid(
        rowid = seq_len(nrow(newdata)),
        time = 1,
        state = c("1", "2"),
        KEEP.OUT.ATTRS = FALSE
      )
      out$tx <- rep(newdata$tx, times = 2)
      out$yprev <- rep(newdata$yprev, times = 2)
      out$prev <- rep(newdata$prev, times = 2)
      out$estimate <- c(0, 0, 1, 0, 1, 1, 0, 1)
      class(out) <- c("markov_sops", class(out))
      attr(out, "call_args") <- list(times = 1)
      attr(out, "y_levels") <- c("1", "2")
      attr(out, "time_var") <- captured$args$time_var
      attr(out, "p_var") <- captured$args$p_var
      attr(out, "p2_var") <- captured$args$p2_var
      out
    },
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model"),
        newdata = newdata,
        variables = list(tx = c(0, 1)),
        estimand = "time_benefit",
        times = 1,
        y_levels = c("1", "2"),
        time_var = "visit",
        p_var = "prev",
        p2_var = "preprev",
        gap_var = "gap",
        time_covariates = time_covariates
      )
    }
  )

  expect_equal(out$estimate, 0.5)
  expect_equal(out$estimand, "time_benefit")
  expect_equal(captured$args$time_var, "visit")
  expect_equal(captured$args$p_var, "prev")
  expect_equal(captured$args$p2_var, "preprev")
  expect_equal(captured$args$gap, "gap")
  expect_equal(captured$args$time_covariates, time_covariates)
})

test_that("time-benefit array reducer preserves scenario blocks", {
  setup <- list(
    n_cf = 2,
    n_each = 2,
    variables = list(tx = c(0, 1)),
    y_levels = c("1", "2"),
    times = 1,
    by = NULL,
    baseline_data = data.frame(id = 1:2)
  )
  sops_array <- array(NA_real_, dim = c(4, 1, 2))
  sops_array[1, 1, ] <- c(0, 1)
  sops_array[2, 1, ] <- c(0, 1)
  sops_array[3, 1, ] <- c(1, 0)
  sops_array[4, 1, ] <- c(0, 1)

  out <- markov.misc:::reduce_time_benefit_array_for_setup(
    sops_array = sops_array,
    setup = setup,
    comparison = "difference",
    weights = NULL,
    time_unit = NULL
  )

  expect_equal(out$estimate, 0.5)
})

test_that("linear ordinal-benefit reducer matches the dense score definition", {
  reference <- rbind(
    c(0.2, 0.3, 0.5),
    c(0.7, 0.2, 0.1)
  )
  comparison <- rbind(
    c(0.4, 0.4, 0.2),
    c(0.1, 0.3, 0.6)
  )
  states <- seq_len(ncol(reference))
  score <- outer(states, states, function(ref, cmp) sign(ref - cmp))
  expected <- rowSums((reference %*% score) * comparison)

  expect_equal(
    markov.misc:::ordinal_pairwise_benefit(reference, comparison),
    expected,
    tolerance = 1e-15
  )
})

test_that("time-benefit array reducer uses real-time AUC when mapped", {
  baseline <- data.frame(id = 1, tx = 0, yprev = 2, rowid = 1)
  grid <- data.frame(tx = c(0, 1))
  setup <- list(
    n_cf = 2,
    n_each = 1,
    variables = list(tx = c(0, 1)),
    y_levels = c("1", "2"),
    times = 1:2,
    by = NULL,
    baseline_data = baseline,
    newdata_orig = baseline,
    newdata_pred = markov.misc:::create_counterfactual_data(
      baseline,
      grid,
      list(tx = c(0, 1))
    ),
    grid = grid,
    id_var = "id",
    time_var = "time",
    p_var = "yprev",
    newdata_supplied = TRUE
  )
  sops_array <- array(NA_real_, dim = c(2, 2, 2))
  sops_array[1, 1, ] <- c(0, 1)
  sops_array[1, 2, ] <- c(0, 1)
  sops_array[2, 1, ] <- c(0, 1)
  sops_array[2, 2, ] <- c(1, 0)

  out <- markov.misc:::reduce_time_benefit_array_for_setup(
    sops_array = sops_array,
    setup = setup,
    comparison = "difference",
    weights = NULL,
    time_unit = NULL,
    time_map = c("1" = 0, "2" = 4),
    baseline_time = NULL
  )

  expect_equal(out$estimate, 2)
})

test_that("time-benefit simulation inference uses real-time comparison scale", {
  baseline <- data.frame(id = 1, tx = 0, yprev = 2, rowid = 1)
  grid <- data.frame(tx = c(0, 1))
  newdata_pred <- markov.misc:::create_counterfactual_data(
    baseline,
    grid,
    list(tx = c(0, 1))
  )
  object <- data.frame(
    estimand = "time_benefit",
    term = "tx",
    reference_level = 0,
    comparison_level = 1,
    contrast = "1 - 0",
    comparison = "difference",
    estimate = 2
  )
  class(object) <- c("markov_avg_comparisons", class(object))
  attr(object, "model") <- structure(
    list(coefficients = c(a = 0, b = 0)),
    class = "mock_model"
  )
  attr(object, "avg_args") <- list(
    variables = list(tx = c(0, 1)),
    by = NULL,
    times = 1:2,
    id_var = "id"
  )
  attr(object, "comparison_args") <- list(
    estimand = "time_benefit",
    comparison = "difference",
    time_map = c("1" = 0, "2" = 4),
    baseline_time = NULL,
    target_times = NULL,
    time_unit = NULL
  )
  attr(object, "newdata_orig") <- baseline
  attr(object, "newdata_pred") <- newdata_pred
  attr(object, "newdata_supplied") <- TRUE
  attr(object, "comparison_baseline_data") <- baseline
  attr(object, "comparison_grid") <- grid
  attr(object, "comparison_n_each") <- 1L
  attr(object, "call_args") <- list(times = 1:2, y_levels = c("1", "2"))
  attr(object, "time_var") <- "time"
  attr(object, "p_var") <- "yprev"
  attr(object, "y_levels") <- c("1", "2")

  sops_array <- array(NA_real_, dim = c(2, 2, 2))
  sops_array[1, 1, ] <- c(0, 1)
  sops_array[1, 2, ] <- c(0, 1)
  sops_array[2, 1, ] <- c(0, 1)
  sops_array[2, 2, ] <- c(1, 0)

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2) * 1e-8,
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) sops_array,
    {
      withr::local_seed(1)
      out <- markov.misc:::inferences_avg_comparisons_time_benefit_simulation(
        object = object,
        engine = "mvn",
        score_weight_dist = "exponential",
        n_draws = 2,
        vcov = NULL,
        cluster = NULL,
        workers = NULL,
        conf_level = 0.95,
        conf_type = "perc",
        return_draws = TRUE
      )
    }
  )

  expect_equal(out$conf.low, 2)
  expect_equal(out$conf.high, 2)
  expect_equal(attr(out, "draws")$estimate, c(2, 2))
})

test_that("avg_comparisons() validates v1 interface boundaries", {
  avg <- make_mock_avg_sops()

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model"),
          variables = list(tx = c(0, 1))
        ),
        "`times` must be supplied to `avg_comparisons()`.",
        fixed = TRUE
      )
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model"),
          variables = list(tx = c(0, 1), z = c(0, 1)),
          times = 1
        ),
        "exactly one named variable",
        fixed = TRUE
      )
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model"),
          variables = list(tx = c(0, 1)),
          estimand = "time_benefit",
          times = 1,
          state_sets = "1"
        ),
        "`state_sets` is not used",
        fixed = TRUE
      )
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model"),
          variables = list(tx = c(0, 1)),
          estimand = "time_benefit",
          times = 1,
          comparison = "ratio"
        ),
        "only supports",
        fixed = TRUE
      )
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model"),
          variables = list(tx = c(0, 1)),
          times = 1,
          origin_time = 0
        ),
        "`origin_time` is no longer supported; use `baseline_time` instead.",
        fixed = TRUE
      )
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model"),
          variables = list(tx = c(0, 1)),
          times = 1,
          origin = "none"
        ),
        "`origin` is no longer supported; use `baseline_time` instead.",
        fixed = TRUE
      )
    }
  )
})

test_that("comparison inference reduces paired draw-level SOPs", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))
  draws <- avg[rep(seq_len(nrow(avg)), times = 2), ]
  draws$draw_id <- rep(1:2, each = nrow(avg))
  draws$estimate <- avg$estimate
  draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "1"] <-
    draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "1"] +
    0.1
  draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "2"] <-
    draws$estimate[draws$draw_id == 1 & draws$tx == 1 & draws$state == "2"] -
    0.1
  draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "1"] <-
    draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "1"] -
    0.1
  draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "2"] <-
    draws$estimate[draws$draw_id == 2 & draws$tx == 1 & draws$state == "2"] +
    0.1

  inferred_avg <- avg
  attr(inferred_avg, "draws") <- draws
  attr(inferred_avg, "conf_level") <- 0.95
  attr(inferred_avg, "conf_type") <- "perc"
  attr(inferred_avg, "method") <- "mvn"

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      cmp <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        estimand = "sop",
        times = 1:2,
        state_sets = "1"
      )
    }
  )

  with_mocked_bindings(
    avg_sops = function(...) avg,
    inferences = function(...) inferred_avg,
    {
      out <- markov.misc:::inferences_avg_comparisons_linear(
        object = cmp,
        method = "mvn",
        score_weight_dist = "exponential",
        n_draws = 2,
        vcov = NULL,
        cluster = NULL,
        workers = NULL,
        conf_level = 0.95,
        conf_type = "perc",
        return_draws = TRUE,
        update_datadist = TRUE,
        use_coefstart = FALSE
      )
    }
  )

  reduced <- attr(out, "draws")
  reduced <- reduced[order(reduced$draw_id, reduced$time), ]
  rownames(reduced) <- NULL

  expect_contains(names(out), c("conf.low", "conf.high", "std.error"))
  expect_equal(reduced$estimate, c(0.4, 0.4, 0.2, 0.2))
  expect_equal(reduced$draw_id, c(1, 1, 2, 2))
})
