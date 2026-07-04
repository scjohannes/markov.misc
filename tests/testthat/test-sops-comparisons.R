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
  attr(out, "ylevels") <- c("1", "2")
  attr(out, "call_args") <- list(times = sort(unique(out$time)))
  attr(out, "tvarname") <- "time"
  attr(out, "pvarname") <- "yprev"
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
        metric = "sop",
        states = "1",
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
        metric = "time_in_state",
        states = "1",
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
        metric = "sop",
        states = "1",
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
        metric = "time_in_state",
        states = "1",
        times = 1:2,
        time_map = c("1" = 0, "2" = 2)
      )
    }
  )

  expect_equal(out$estimate, 0.6)
})

test_that("avg_comparisons() forwards model-structure args for marginal metrics", {
  avg <- make_mock_avg_sops(tx_levels = c(0, 1))
  t_covs <- data.frame(visit = 1:2, spline = c(0, 1))
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
        metric = "time_in_state",
        states = "1",
        times = 1:2,
        ylevels = c("1", "2"),
        absorb = "2",
        tvarname = "visit",
        pvarname = "prev",
        p2varname = "preprev",
        gap = "gap",
        t_covs = t_covs
      )
    }
  )

  expect_s3_class(out, "markov_avg_comparisons")
  expect_equal(captured$args$ylevels, c("1", "2"))
  expect_equal(captured$args$absorb, "2")
  expect_equal(captured$args$tvarname, "visit")
  expect_equal(captured$args$pvarname, "prev")
  expect_equal(captured$args$p2varname, "preprev")
  expect_equal(captured$args$gap, "gap")
  expect_equal(captured$args$t_covs, t_covs)
})

test_that("avg_comparisons() computes time benefit from paired profiles", {
  newdata <- data.frame(
    id = 1:2,
    tx = c(0, 1),
    yprev = c(2, 2),
    prev = c(2, 2),
    preprev = c(1, 1)
  )
  t_covs <- data.frame(visit = 1, spline = 0)
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
      attr(out, "ylevels") <- c("1", "2")
      attr(out, "tvarname") <- captured$args$tvarname
      attr(out, "pvarname") <- captured$args$pvarname
      attr(out, "p2varname") <- captured$args$p2varname
      out
    },
    {
      out <- avg_comparisons(
        structure(list(), class = "mock_model"),
        newdata = newdata,
        variables = list(tx = c(0, 1)),
        metric = "time_benefit",
        times = 1,
        ylevels = c("1", "2"),
        tvarname = "visit",
        pvarname = "prev",
        p2varname = "preprev",
        gap = "gap",
        t_covs = t_covs
      )
    }
  )

  expect_equal(out$estimate, 0.5)
  expect_equal(out$metric, "time_benefit")
  expect_equal(captured$args$tvarname, "visit")
  expect_equal(captured$args$pvarname, "prev")
  expect_equal(captured$args$p2varname, "preprev")
  expect_equal(captured$args$gap, "gap")
  expect_equal(captured$args$t_covs, t_covs)
})

test_that("time-benefit array reducer preserves scenario blocks", {
  setup <- list(
    n_cf = 2,
    n_each = 2,
    variables = list(tx = c(0, 1)),
    ylevels = c("1", "2"),
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
          metric = "time_benefit",
          times = 1,
          states = "1"
        ),
        "`states` is not used",
        fixed = TRUE
      )
      expect_error(
        avg_comparisons(
          structure(list(), class = "mock_model"),
          variables = list(tx = c(0, 1)),
          metric = "time_benefit",
          times = 1,
          comparison = "ratio"
        ),
        "only supports",
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
  attr(inferred_avg, "simulation_draws") <- draws
  attr(inferred_avg, "conf_level") <- 0.95
  attr(inferred_avg, "conf_type") <- "perc"
  attr(inferred_avg, "method") <- "simulation"
  attr(inferred_avg, "engine") <- "mvn"

  with_mocked_bindings(
    avg_sops = function(...) avg,
    {
      cmp <- avg_comparisons(
        structure(list(), class = "mock_model"),
        variables = list(tx = c(0, 1)),
        metric = "sop",
        times = 1:2,
        states = "1"
      )
    }
  )

  with_mocked_bindings(
    avg_sops = function(...) avg,
    inferences = function(...) inferred_avg,
    {
      out <- markov.misc:::inferences_avg_comparisons_linear(
        object = cmp,
        method = "simulation",
        engine = "mvn",
        score_weight_dist = "exponential",
        n_sim = 2,
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

  reduced <- attr(out, "simulation_draws")
  reduced <- reduced[order(reduced$draw_id, reduced$time), ]
  rownames(reduced) <- NULL

  expect_contains(names(out), c("conf.low", "conf.high", "std.error"))
  expect_equal(reduced$estimate, c(0.4, 0.4, 0.2, 0.2))
  expect_equal(reduced$draw_id, c(1, 1, 2, 2))
})
