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

test_that("taooh_bootstrap stores NA for failed bootstrap fits", {
  model <- structure(list(coefficients = c(tx = 1)), class = "vglm")
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = factor(c(1, 2, 1, 2)),
    yprev = factor(c(1, 1, 1, 1)),
    tx = c(0, 0, 1, 1)
  )

  call_id <- 0
  with_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      replicate(
        n_boot,
        data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1")),
        simplify = FALSE
      )
    },
    apply_to_bootstrap = function(boot_samples, analysis_fn, data, id_var, parallel, workers, packages, globals) {
      lapply(boot_samples, function(sample) analysis_fn(data))
    },
    bootstrap_analysis_wrapper = function(...) {
      call_id <<- call_id + 1
      if (call_id == 2) {
        list(model = NULL, data = data, ylevels = 1:2, absorb = 2)
      } else {
        list(model = model, data = data, ylevels = 1:2, absorb = 2)
      }
    },
    taooh = function(...) c(delta_taooh = 0.25),
    {
      result <- taooh_bootstrap(model, data = data, n_boot = 2, parallel = FALSE)
    }
  )

  expect_equal(as.numeric(result$models[[1]]), 0.25)
  expect_true(is.na(result$models[[2]]))
})

test_that("taooh_bootstrap records taooh calculation failures", {
  model <- structure(list(coefficients = c(tx = 1)), class = "vglm")
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = factor(c(1, 2, 1, 2)),
    yprev = factor(c(1, 1, 1, 1)),
    tx = c(0, 0, 1, 1)
  )

  with_mocked_bindings(
    fast_group_bootstrap = function(...) {
      list(data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1")))
    },
    apply_to_bootstrap = function(boot_samples, analysis_fn, data, id_var, parallel, workers, packages, globals) {
      list(analysis_fn(data))
    },
    bootstrap_analysis_wrapper = function(...) {
      list(model = model, data = data, ylevels = 1:2, absorb = 2)
    },
    taooh = function(...) stop("taooh failed"),
    {
      expect_warning(
        result <- taooh_bootstrap(model, data = data, n_boot = 1, parallel = FALSE),
        "taooh calculation failed",
        fixed = TRUE
      )
    }
  )

  expect_true(is.na(result$models[[1]]))
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

test_that("taooh_bootstrap2 records missing-state, fit, and taooh failures", {
  withr::local_options(datadist = NULL)
  model <- structure(list(coefficients = c(tx = 1)), class = "orm")
  data <- data.frame(
    id = rep(1, 6),
    time = rep(1:2, times = 3),
    y = factor(c(1, 2, 1, 2, 1, 2), levels = 1:3),
    yprev = factor(c(1, 1, 1, 1, 2, 2), levels = 1:3),
    tx = c(0, 0, 1, 1, 0, 0)
  )
  withr::defer(
    if (exists("dd", envir = globalenv())) rm("dd", envir = globalenv())
  )

  update_call <- 0
  taooh_call <- 0
  with_mocked_bindings(
    plan = function(...) NULL,
    future_map_dbl = function(.x, .f, .options = NULL) vapply(.x, .f, numeric(1)),
    furrr_options = function(...) list(...),
    update = function(object, data, ...) {
      update_call <<- update_call + 1
      if (update_call == 2) stop("fit failed")
      object
    },
    taooh = function(...) {
      taooh_call <<- taooh_call + 1
      if (taooh_call == 2) stop("taooh failed")
      c(delta_taooh = 0.2)
    },
    {
      warnings <- character()
      result <- withCallingHandlers(
        taooh_bootstrap2(
          model,
          data = data,
          n_boot = 3,
          workers = 1,
          parallel = FALSE,
          ylevels = 1:3,
          absorb = 3,
          quantiles = 0.75
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    }
  )

  expect_true(any(grepl("not present", warnings, fixed = TRUE)))
  expect_true(any(grepl("failed", warnings, fixed = TRUE)))
  expect_equal(unname(result$all_results), c(0.2, NA_real_, NA_real_))
  expect_equal(result$errors, 2)
})

test_that("taooh_bootstrap2 preserves results when all states are present", {
  withr::local_options(datadist = NULL)
  model <- structure(list(coefficients = c(tx = 1)), class = "vglm")
  data <- data.frame(
    id = rep(1, 6),
    time = rep(1:2, times = 3),
    y = factor(c(1, 2, 2, 3, 1, 3), levels = 1:3),
    yprev = factor(c(1, 1, 2, 2, 3, 3), levels = 1:3),
    tx = c(0, 0, 1, 1, 0, 0)
  )

  with_mocked_bindings(
    plan = function(...) NULL,
    future_map_dbl = function(.x, .f, .options = NULL) vapply(.x, .f, numeric(1)),
    furrr_options = function(...) list(...),
    update = function(object, data, ...) object,
    taooh = function(...) c(delta_taooh = 0.4),
    {
      result <- taooh_bootstrap2(
        model,
        data = data,
        n_boot = 1,
        workers = 1,
        parallel = FALSE,
        ylevels = 1:3,
        absorb = 3
      )
    }
  )

  expect_equal(unname(result$all_results), 0.4)
  expect_equal(result$errors, 0)
})

test_that("taooh_bootstrap2 handles default workers, parallel plans, and numeric states", {
  model <- structure(list(coefficients = c(tx = 1)), class = "vglm")
  data <- data.frame(
    id = rep(1, 4),
    time = rep(1:2, times = 2),
    y = c(1, 3, 1, 3),
    yprev = c(1, 1, 3, 3),
    tx = c(0, 0, 1, 1)
  )

  plan_calls <- list()
  with_mocked_bindings(
    plan = function(...) {
      plan_calls[[length(plan_calls) + 1L]] <<- list(...)
    },
    future_map_dbl = function(.x, .f, .options = NULL) vapply(.x, .f, numeric(1)),
    furrr_options = function(...) list(...),
    update = function(object, data, ...) object,
    taooh = function(...) c(delta_taooh = 0.6),
    {
      withr::local_seed(124)
      result <- taooh_bootstrap2(
        model,
        data = data,
        n_boot = 1,
        workers = NULL,
        parallel = TRUE,
        ylevels = 1:3,
        absorb = 3
      )
    }
  )

  expect_equal(unname(result$all_results), 0.6)
  expect_true(length(plan_calls) >= 1)
})
