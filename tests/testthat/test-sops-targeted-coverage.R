test_that("sops() covers validation and stratified aggregation branches", {
  model <- structure(list(), class = "mock_model")
  newdata <- data.frame(
    id = 1:3,
    time = c(1, 2, 1),
    yprev = c(1, 2, 1),
    tx = c(0, 1, 1)
  )
  sops_array <- array(
    c(
      0.1, 0.2, 0.3,
      0.4, 0.5, 0.6,
      0.9, 0.8, 0.7,
      0.6, 0.5, 0.4
    ),
    dim = c(3, 2, 2)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    soprob_markov = function(object, data, times, ylevels, absorb, tvarname, pvarname, gap, t_covs) {
      expect_equal(times, c(1, 2))
      expect_equal(ylevels, c(1, 2))
      sops_array
    },
    {
      expect_error(
        sops(model, newdata = NULL, times = 1, ylevels = 1:2),
        "Please provide `newdata`",
        fixed = TRUE
      )
      expect_error(
        sops(model, newdata = data.frame(id = 1), ylevels = 1:2),
        "`times` must be specified",
        fixed = TRUE
      )
      expect_error(
        sops(model, newdata = newdata, times = 1:2, ylevels = 1:2, by = 1),
        "'by' must be a character vector",
        fixed = TRUE
      )
      expect_error(
        sops(model, newdata = newdata, times = 1:2, ylevels = 1:2, by = "missing"),
        "not found in newdata",
        fixed = TRUE
      )

      result <- sops(
        model,
        newdata = newdata,
        times = 1:2,
        ylevels = 1:2,
        by = "tx"
      )

      inferred_times <- sops(
        model,
        newdata = newdata,
        times = NULL,
        ylevels = 1:2
      )
    }
  )

  expect_equal(sort(unique(inferred_times$time)), c(1, 2))
  expected <- aggregate(
    estimate ~ time + state + tx,
    data = data.frame(
      tx = rep(newdata$tx, times = 4),
      time = rep(rep(1:2, each = 3), times = 2),
      state = rep(1:2, each = 6),
      estimate = as.vector(sops_array)
    ),
    FUN = mean
  )
  result <- result[order(result$time, result$state, result$tx), ]
  expected <- expected[order(expected$time, expected$state, expected$tx), ]
  rownames(result) <- NULL
  rownames(expected) <- NULL

  expect_equal(
    as.data.frame(result[, c("time", "state", "tx", "estimate")]),
    expected[, c("time", "state", "tx", "estimate")]
  )
  expect_equal(attr(result, "by"), "tx")
})

test_that("sops() infers ylevels from supported model containers", {
  newdata <- data.frame(id = 1, time = 1, yprev = 1)
  sops_array <- array(c(0.25, 0.75), dim = c(1, 1, 2))

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    soprob_markov = function(object, data, times, ylevels, absorb, tvarname, pvarname, gap, t_covs) {
      expect_equal(ylevels, c("1", "2"))
      sops_array
    },
    {
      vglm_model <- methods::new("vglm")
      vglm_model@extra <- list(colnames.y = c("1", "2"))
      vglm_model@family <- methods::new("vglmff")
      vglm_model@call <- quote(VGAM::vglm(y ~ x, family = VGAM::cumulative(reverse = TRUE)))

      robust_model <- list(extra = list(colnames.y = c("1", "2")))
      class(robust_model) <- "robcov_vglm"

      orm_model <- list(yunique = c("1", "2"))
      class(orm_model) <- "orm"

      expect_equal(nrow(sops(vglm_model, newdata = newdata, times = 1)), 2)
      expect_equal(nrow(sops(robust_model, newdata = newdata, times = 1)), 2)
      expect_equal(nrow(sops(orm_model, newdata = newdata, times = 1)), 2)

      vglm_model@extra <- list(colnames.y = NULL)
      expect_error(
        sops(vglm_model, newdata = newdata, times = 1),
        "`ylevels` cannot be NULL",
        fixed = TRUE
      )

      robust_model$extra <- list(colnames.y = NULL)
      expect_error(
        sops(robust_model, newdata = newdata, times = 1),
        "`ylevels` cannot be NULL",
        fixed = TRUE
      )

      orm_model$yunique <- NULL
      expect_error(
        sops(orm_model, newdata = newdata, times = 1),
        "`ylevels` cannot be NULL",
        fixed = TRUE
      )

      expect_error(
        sops(structure(list(), class = "other"), newdata = newdata, times = 1),
        "`ylevels` cannot be NULL",
        fixed = TRUE
      )
    }
  )
})

test_that("avg_sops() covers validation and by branches", {
  model <- structure(list(), class = "mock_model")
  newdata <- data.frame(
    id = c(1, 2),
    tx = c(0, 1),
    subgroup = c("a", "b"),
    yprev = c(1, 2)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    sops = function(model, newdata, times, ...) {
      out <- expand.grid(
        time = times,
        state = c(1, 2),
        tx = c(0, 1),
        subgroup = c("a", "b"),
        KEEP.OUT.ATTRS = FALSE
      )
      out$estimate <- seq_len(nrow(out)) / 100
      attr(out, "model") <- model
      attr(out, "call_args") <- list(times = times)
      attr(out, "tvarname") <- "time"
      attr(out, "pvarname") <- "yprev"
      attr(out, "ylevels") <- c(1, 2)
      attr(out, "absorb") <- NULL
      attr(out, "gap") <- NULL
      attr(out, "t_covs") <- NULL
      out
    },
    {
      expect_error(
        avg_sops(model, newdata = newdata, variables = NULL, times = 1),
        "`variables` is required",
        fixed = TRUE
      )
      expect_error(
        avg_sops(model, newdata = NULL, variables = "tx", times = 1),
        "Provide newdata",
        fixed = TRUE
      )
      expect_error(
        avg_sops(model, newdata = data.frame(tx = 0), variables = "tx", times = 1),
        "ID variable 'id' not found",
        fixed = TRUE
      )
      expect_error(
        avg_sops(model, newdata = newdata, variables = list(missing = 1), times = 1),
        "Variables not in data",
        fixed = TRUE
      )

      result <- avg_sops(
        model,
        newdata = newdata,
        variables = "tx",
        by = "subgroup",
        times = 1:2
      )
    }
  )

  expect_true("subgroup" %in% names(result))
  expect_equal(attr(result, "avg_args")$by, "subgroup")
})

test_that("avg_sops() reports missing grouping variables from individual SOPs", {
  model <- structure(list(), class = "mock_model")
  newdata <- data.frame(id = 1, tx = 0, yprev = 1)

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    sops = function(model, newdata, times, ...) {
      out <- data.frame(time = 1, state = 1, tx = 0, estimate = 0.5)
      attr(out, "model") <- model
      attr(out, "call_args") <- list(times = times)
      out
    },
    {
      expect_error(
        avg_sops(
          model,
          newdata = newdata,
          variables = "tx",
          by = "subgroup",
          times = 1
        ),
        "Grouping variables missing",
        fixed = TRUE
      )
    }
  )
})

test_that("get_draws() covers invalid input and no-key fallback", {
  expect_error(
    get_draws(data.frame(x = 1)),
    "get_draws() requires an object from inferences",
    fixed = TRUE
  )

  object <- data.frame(time = 1, state = 1, estimate = 0.5)
  class(object) <- c("markov_avg_sops", class(object))
  attr(object, "draws") <- data.frame(draw_id = 1:2, estimate = c(0.4, 0.6))

  draws <- get_draws(object)

  expect_equal(names(draws), c("draw_id", "draw"))
  expect_equal(draws$draw, c(0.4, 0.6))
})

test_that("get_draws() joins metadata with base merge", {
  object <- data.frame(time = 1:2, state = 1, estimate = c(0.5, 0.6))
  class(object) <- c("markov_avg_sops", class(object))
  attr(object, "simulation_draws") <- data.frame(
    draw_id = c(1, 1),
    time = 1:2,
    state = 1,
    estimate = c(0.4, 0.7)
  )

  draws <- get_draws(object)

  draws <- draws[order(draws$time), ]
  rownames(draws) <- NULL

  expect_equal(draws$draw, c(0.4, 0.7))
  expect_equal(draws$estimate, c(0.5, 0.6))
})

test_that("markov_msm_build() covers default columns and missing design columns", {
  model <- stats::lm(y ~ x, data = data.frame(y = 1:3, x = 1:3))
  data <- data.frame(id = 1:2, x = c(0, 1))

  with_mocked_bindings(
    get_effective_coefs = function(model) {
      out <- matrix(0, nrow = 1, ncol = 1)
      colnames(out) <- "missing_column"
      out
    },
    {
      expect_error(
        components <- markov.misc:::markov_msm_build(
          model = model,
          data = data,
          ylevels = 1:2,
          pvarname = "yprev"
        ),
        "subscript out of bounds"
      )
    }
  )
})

test_that("markov_msm_run() covers null design matrix and no absorbing states", {
  components <- list(
    X_0 = NULL,
    X_slopes = list(),
    X_prev = list(NULL, NULL),
    X_interactions = list(list(), list()),
    y_base_idx = c(1L, 2L),
    t_covs = NULL,
    n_pat = 2L,
    n_states = 2L,
    M = 1L,
    ylevels = 1:2,
    col_names = "x"
  )
  gamma <- matrix(0, nrow = 1, ncol = 1)
  colnames(gamma) <- "x"

  out <- markov.misc:::markov_msm_run(
    components = components,
    Gamma = gamma,
    times = 1:2,
    absorb = NULL
  )

  expect_equal(dim(out), c(2L, 2L, 2L))
  expect_equal(out[, 1, 1], c(0.5, 0.5))
  expect_equal(out[, 2, 1], c(0.5, 0.5))
})
