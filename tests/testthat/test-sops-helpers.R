test_that("sops() validates inputs and supports stratified aggregation", {
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
    soprob_markov = function(object, data, times, ylevels, absorb, tvarname, pvarname, gap, t_covs, ...) {
      expect_equal(times, c(1, 2))
      expect_equal(ylevels, c(1, 2))
      sops_array
    },
    {
      expect_error(
        sops(model, newdata = NULL, times = 1, ylevels = 1:2),
        "Provide newdata",
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
    soprob_markov = function(object, data, times, ylevels, absorb, tvarname, pvarname, gap, t_covs, ...) {
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

test_that("avg_sops() validates inputs and preserves grouping metadata", {
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

test_that("get_draws() rejects plain data frames and falls back without join keys", {
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

test_that("get_draws() preserves covariates for individual SOP draws", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("dplyr")

  withr::local_seed(123)
  data <- make_test_data(n_patients = 20, seed = 123, follow_up_time = 30)
  data$age <- stats::rnorm(nrow(data), 60, 10)

  m_vglm <- VGAM::vglm(
    ordered(y) ~ rms::rcs(time, 3) + tx + yprev + age,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  expect_warning(
    m_vglm_rob <- robcov_vglm(m_vglm, cluster = data$id),
    "fewer than 30 clusters"
  )

  sops_result <- sops(
    model = m_vglm_rob,
    newdata = dplyr::filter(data, time == 1),
    times = 1:5,
    ylevels = 1:6,
    absorb = 6,
    tvarname = "time",
    pvarname = "yprev"
  )

  result_ci <- inferences(
    sops_result,
    method = "simulation",
    n_sim = 10,
    return_draws = TRUE
  )

  draws <- get_draws(result_ci)
  sample_row <- draws[1, ]
  orig_row <- data[sample_row$rowid, ]

  expect_s3_class(draws, "data.frame")
  expect_contains(names(draws), c("draw_id", "rowid", "time", "state", "estimate"))
  expect_contains(names(draws), c("tx", "age"))
  expect_equal(sample_row$tx, orig_row$tx)
  expect_equal(sample_row$age, orig_row$age)
})

test_that("markov_msm_build() reports missing design columns", {
  model <- VGAM::vglm(
    y ~ x,
    family = VGAM::binomialff,
    data = data.frame(
      y = c(0, 1, 1, 0, 0, 1),
      x = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    )
  )
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
        "Could not construct a design matrix"
      )
    }
  )
})

test_that("markov_msm_run() propagates probabilities without absorbing states", {
  components <- list(
    X_init = matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x")),
    X_transition = list(
      list(
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x")),
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x"))
      ),
      list(
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x")),
        matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "x"))
      )
    ),
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

test_that("standardize_sops() validates inputs and marginalizes arrays", {
  model <- structure(list(), class = "vglm")

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    {
      expect_error(
        standardize_sops(model, data = NULL),
        "No data provided",
        fixed = TRUE
      )
      expect_error(
        standardize_sops(
          model,
          data = data.frame(id = 1, time = 1, tx = 0),
          ylevels = 1:2
        ),
        "Previous-state variable",
        fixed = TRUE
      )
    }
  )

  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    yprev = c(1, 1, 2, 2),
    tx = c(0, 0, 1, 1)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    soprob_markov = function(object, data, times, ylevels, absorb, tvarname, pvarname, gap, t_covs, ...) {
      if (all(data$tx == 1)) {
        array(
          c(0.8, 0.6, 0.2, 0.4, 0.7, 0.5, 0.3, 0.5),
          dim = c(2, 2, 2)
        )
      } else {
        array(
          c(0.4, 0.2, 0.6, 0.8, 0.3, 0.1, 0.7, 0.9),
          dim = c(2, 2, 2)
        )
      }
    },
    {
      out <- standardize_sops(
        model,
        data = data,
        times = 1:2,
        ylevels = 1:2,
        absorb = NULL
      )
    }
  )

  expect_equal(dim(out$sop_tx), c(2L, 2L))
  expect_equal(colnames(out$sop_tx), c("1", "2"))
  expect_equal(out$sop_tx[, 1], c(mean(c(0.8, 0.6)), mean(c(0.2, 0.4))))
  expect_equal(out$sop_ctrl[, 2], c(mean(c(0.3, 0.1)), mean(c(0.7, 0.9))))
})

test_that("standardize_sops() supports model data and draw arrays", {
  model <- structure(
    list(
      x = data.frame(
        id = c(1, 2),
        time = c(1, 1),
        yprev = c(1, 2),
        tx = c(0, 1)
      ),
      y = c(1, 2)
    ),
    class = "vglm"
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    soprob_markov = function(...) {
      arr <- array(
        seq_len(2 * 2 * 2 * 2) / 100,
        dim = c(2, 2, 2, 2),
        dimnames = list(NULL, NULL, NULL, c("1", "2"))
      )
      arr
    },
    {
      out <- standardize_sops(
        model,
        data = NULL,
        times = 1:2,
        ylevels = 1:2,
        absorb = NULL
      )
    }
  )

  expect_equal(dim(out$sop_tx), c(2L, 2L, 2L))
  expect_equal(dimnames(out$sop_tx)[[3]], c("1", "2"))
})

test_that("standardize_sops() rejects unexpected SOP array dimensions", {
  model <- structure(list(), class = "vglm")
  data <- data.frame(id = 1, time = 1, yprev = 1, tx = 0)

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    soprob_markov = function(...) matrix(1, nrow = 1, ncol = 1),
    {
      expect_error(
        standardize_sops(model, data = data, times = 1, ylevels = 1),
        "Unexpected array dimensions",
        fixed = TRUE
      )
    }
  )
})

test_that("time_in_state() handles bootstrap data frames", {
  boot <- data.frame(
    boot_id = rep(1:2, each = 4),
    time = rep(1:2, times = 4),
    tx = rep(c(1, 0), each = 2, times = 2),
    state_1 = c(0.8, 0.7, 0.5, 0.4, 0.6, 0.5, 0.3, 0.2),
    state_2 = c(0.2, 0.3, 0.5, 0.6, 0.4, 0.5, 0.7, 0.8)
  )

  one_state <- time_in_state(boot, target_states = 1)
  two_states <- time_in_state(boot, target_states = c(1, 2))

  expect_equal(one_state$SOP_tx, c(1.5, 1.1))
  expect_equal(one_state$SOP_ctrl, c(0.9, 0.5))
  expect_equal(two_states$delta, c(0, 0))

  expect_error(
    time_in_state(data.frame(boot_id = 1, time = 1)),
    "must contain 'boot_id', 'time', and 'tx'",
    fixed = TRUE
  )
  expect_error(
    time_in_state(data.frame(boot_id = 1, time = 1, tx = 0)),
    "does not contain any columns starting with 'state_'",
    fixed = TRUE
  )
  expect_error(
    time_in_state(boot, target_states = 3),
    "target state columns were not found",
    fixed = TRUE
  )
})

test_that("time_in_state() handles arrays and validation branches", {
  arr3 <- array(
    c(
      0.8, 0.6,
      0.2, 0.4,
      0.7, 0.5,
      0.3, 0.5
    ),
    dim = c(2, 2, 2),
    dimnames = list(NULL, NULL, c("1", "2"))
  )
  arr4 <- array(
    seq_len(2 * 2 * 2 * 2) / 100,
    dim = c(2, 2, 2, 2),
    dimnames = list(NULL, NULL, NULL, c("1", "2"))
  )
  unnamed <- array(rep(c(0.2, 0.8), 4), dim = c(2, 2, 2))

  expect_equal(time_in_state(arr3, target_states = "1"), c(1.0, 1.0))
  expect_equal(dim(time_in_state(arr4, target_states = c("1", "2"))), c(2L, 2L))
  expect_equal(time_in_state(unnamed, target_states = 1), c(0.4, 1.6))
  expect_error(time_in_state(1:3), "must be an array or data frame", fixed = TRUE)
  expect_error(time_in_state(matrix(1, 2, 2)), "must be 3D or 4D", fixed = TRUE)
  expect_error(
    time_in_state(arr3, target_states = "3"),
    "Target states not found",
    fixed = TRUE
  )
})

test_that("low-level SOP helpers validate model families and state metadata", {
  skip_if_not_installed("VGAM")
  loadNamespace("VGAM")

  expect_error(
    markov.misc:::validate_markov_model(lm(mpg ~ wt, data = mtcars)),
    NA
  )

  model <- methods::new("vglm")
  model@family <- methods::new("vglmff", vfamily = "binomialff")
  expect_error(
    markov.misc:::validate_markov_model(model),
    "only the cumulative family",
    fixed = TRUE
  )

  model@family <- methods::new("vglmff", vfamily = "cumulative")
  model@call <- quote(VGAM::vglm(y ~ x, family = VGAM::cumulative()))
  expect_error(
    markov.misc:::validate_markov_model(model),
    "argument is missing|Cannot determine"
  )

  model@call <- quote(VGAM::vglm(y ~ x, family = VGAM::cumulative(reverse = FALSE)))
  expect_error(
    markov.misc:::validate_markov_model(model),
    "must use reverse = TRUE",
    fixed = TRUE
  )

  offset_data <- data.frame(
    y = ordered(rep(1:3, length.out = 30)),
    x = rep(c(0, 1), length.out = 30),
    off = seq(-0.2, 0.2, length.out = 30)
  )
  offset_model <- suppressWarnings(
    VGAM::vglm(
      y ~ x + offset(off),
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = offset_data
    )
  )
  expect_error(
    markov.misc:::validate_markov_model(offset_model),
    "Model offsets are not supported",
    fixed = TRUE
  )

  null_offset_model <- suppressWarnings(
    VGAM::vglm(
      y ~ x,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = offset_data,
      offset = NULL
    )
  )
  expect_error(
    markov.misc:::validate_markov_model(null_offset_model),
    NA
  )

  fit_with_optional_offset <- function(offset = NULL) {
    suppressWarnings(
      VGAM::vglm(
        y ~ x,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = offset_data,
        offset = offset
      )
    )
  }
  expect_error(
    markov.misc:::validate_markov_model(fit_with_optional_offset()),
    NA
  )

  expect_error(
    markov.misc:::score_bootstrap_components(list()),
    "supports robcov_vglm and orm",
    fixed = TRUE
  )

  expect_error(
    markov.misc:::marginalize_sops_array(
      array(1, dim = c(1, 1, 1)),
      grid = data.frame(tx = 0),
      times = 1,
      ylevels = 1,
      variables = list(tx = 0),
      n_cf = 1,
      n_each = 1,
      weights = 0
    ),
    "positive finite sum",
    fixed = TRUE
  )

  expect_warning(
    no_rowid <- markov.misc:::array_to_df_individual(
      array(1, dim = c(1, 1, 1)),
      times = 1,
      ylevels = 1,
      newdata = data.frame(id = 10),
      by = "missing"
    ),
    "Variables in 'by' not found",
    fixed = TRUE
  )
  expect_equal(no_rowid$rowid, 1)
})

test_that("low-level SOP helpers reject orm offsets", {
  skip_if_not_installed("rms")

  offset_data <- data.frame(
    y = ordered(rep(1:3, length.out = 45)),
    x = rep(c(0, 1, 2), length.out = 45),
    off = seq(-0.2, 0.2, length.out = 45)
  )
  dd <- rms::datadist(offset_data)
  assign("dd_offset_test", dd, envir = globalenv())
  withr::defer(rm("dd_offset_test", envir = globalenv()))
  withr::local_options(datadist = "dd_offset_test")

  offset_model <- rms::orm(
    y ~ x + offset(off),
    data = offset_data,
    x = TRUE,
    y = TRUE
  )

  expect_error(
    markov.misc:::validate_markov_model(offset_model),
    "Model offsets are not supported",
    fixed = TRUE
  )
})
