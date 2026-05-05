test_that("taooh validates model class and data", {
  expect_error(
    taooh(lm(mpg ~ wt, data = mtcars), data = mtcars),
    "model must be an orm object",
    fixed = TRUE
  )

  model <- structure(list(), class = "vglm")
  expect_error(
    taooh(model, data = NULL),
    "Provide data used for model fitting",
    fixed = TRUE
  )
})

test_that("taooh computes treatment-control time in target states", {
  model <- structure(list(), class = "vglm")
  data <- data.frame(
    id = c("ctrl", "ctrl", "tx", "tx"),
    time = c(1, 2, 1, 2),
    y = c(1, 2, 1, 2),
    yprev = c(1, 1, 1, 1),
    tx = c(0, 0, 1, 1)
  )

  with_mocked_bindings(
    soprobMarkovOrdm = function(object, data, times, ylevels, absorb, tvarname, pvarname, gap) {
      if (data$tx == 1) {
        matrix(c(0.6, 0.4, 0.0, 0.3, 0.7, 0.0), nrow = 2, byrow = TRUE)
      } else {
        matrix(c(0.2, 0.8, 0.0, 0.1, 0.9, 0.0), nrow = 2, byrow = TRUE)
      }
    },
    {
      expect_warning(
        result <- taooh(
          model,
          data = data,
          times = 1:2,
          ylevels = 1:3,
          absorb = 3,
          target_states = c(1, 2)
        ),
        "Variable 'yprev' should be a factor"
      )
    }
  )

  expect_equal(result[["SOP_tx"]], 2)
  expect_equal(result[["SOP_ctrl"]], 2)
  expect_equal(result[["delta_taooh"]], 0)
})

test_that("taooh can compute a single target state", {
  model <- structure(list(), class = "vglm")
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(1, 2, 1, 2),
    yprev = factor(c(1, 1, 1, 1)),
    tx = c(0, 0, 1, 1)
  )

  with_mocked_bindings(
    soprobMarkovOrdm = function(object, data, times, ylevels, absorb, tvarname, pvarname, gap) {
      if (data$tx == 1) {
        matrix(c(0.8, 0.2, 0.6, 0.4), nrow = 2, byrow = TRUE)
      } else {
        matrix(c(0.3, 0.7, 0.2, 0.8), nrow = 2, byrow = TRUE)
      }
    },
    {
      result <- taooh(
        model,
        data = data,
        ylevels = 1:2,
        absorb = 2,
        target_states = 1
      )
    }
  )

  expect_equal(result[["SOP_tx"]], 1.4)
  expect_equal(result[["SOP_ctrl"]], 0.5)
  expect_equal(result[["delta_taooh"]], 0.9)
})

test_that("taooh_bootstrap validates models and returns bootstrap records", {
  expect_error(
    taooh_bootstrap(lm(mpg ~ wt, data = mtcars), data = mtcars, n_boot = 1),
    "model must be an orm object",
    fixed = TRUE
  )

  model <- structure(list(), class = "vglm")
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = factor(c(1, 2, 1, 2)),
    yprev = factor(c(1, 1, 1, 1)),
    tx = c(0, 0, 1, 1)
  )

  with_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      list(
        data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1),
        data.frame(original_id = c(2, 1), new_id = c("2_1", "1_1"), boot_id = 2)
      )
    },
    apply_to_bootstrap = function(boot_samples, analysis_fn, data, id_var, parallel, workers, packages, globals) {
      expect_equal(parallel, FALSE)
      expect_null(workers)
      list(0.2, 0.4)
    },
    {
      result <- taooh_bootstrap(model, data = data, n_boot = 2, parallel = FALSE)
    }
  )

  expect_equal(result$id, c("Bootstrap01", "Bootstrap02"))
  expect_equal(result$models, list(0.2, 0.4))
})

test_that("taooh_bootstrap2 validates and summarizes bootstrap draws", {
  expect_error(
    taooh_bootstrap2(lm(mpg ~ wt, data = mtcars), data = mtcars, n_boot = 1),
    "model must be an orm object",
    fixed = TRUE
  )

  model <- structure(list(), class = "vglm")
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = factor(c(1, 2, 1, 2)),
    yprev = factor(c(1, 1, 1, 1)),
    tx = c(0, 0, 1, 1)
  )

  with_mocked_bindings(
    plan = function(...) NULL,
    future_map_dbl = function(.x, .f, .options = NULL) c(0.1, NA_real_, 0.5),
    furrr_options = function(...) list(...),
    {
      withr::local_seed(123)
      result <- taooh_bootstrap2(
        model,
        data = data,
        n_boot = 3,
        workers = 1,
        parallel = FALSE,
        quantiles = 0.8
      )
    }
  )

  expect_equal(result$t_est, NA_real_)
  expect_equal(
    unname(result$ci),
    unname(stats::quantile(c(0.1, NA, 0.5), c(0.2, 0.8), na.rm = TRUE))
  )
  expect_equal(result$errors, 1)
  expect_equal(result$all_results, c(0.1, NA_real_, 0.5))
})
