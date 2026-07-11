test_that("fast_group_bootstrap samples group ids and creates unique ids", {
  data <- data.frame(
    id = c("a", "a", "b", "c", "c"),
    value = 1:5
  )

  withr::local_seed(42)
  boot <- fast_group_bootstrap(data, id_var = "id", n_boot = 3)

  expect_length(boot, 3)
  expect_equal(vapply(boot, nrow, integer(1)), rep(3L, 3))
  expect_named(boot[[1]], c("original_id", "new_id", "boot_id"))
  expect_all_true(boot[[1]]$original_id %in% unique(data$id))
  expect_equal(boot[[2]]$boot_id, rep(2, 3))
  expect_equal(length(unique(boot[[1]]$new_id)), nrow(boot[[1]]))
})

test_that("fast_group_bootstrap validates the grouping column", {
  data <- data.frame(id = 1:3)

  expect_error(
    fast_group_bootstrap(data, id_var = "patient", n_boot = 1),
    "id_var 'patient' not found in data",
    fixed = TRUE
  )
})

test_that("materialize_bootstrap_sample preserves sampled copies and row data", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(0, 1, 1, 0)
  )
  boot_ids <- data.frame(
    original_id = c(2, 1, 2),
    new_id = c("2_1", "1_1", "2_2"),
    boot_id = 1
  )

  result <- materialize_bootstrap_sample(boot_ids, data, "id")

  expect_equal(nrow(result), 6)
  expect_equal(result$new_id, rep(c("2_1", "1_1", "2_2"), each = 2))
  expect_equal(result$time, c(1, 2, 1, 2, 1, 2))
  expect_equal(result$y, c(1, 0, 0, 1, 1, 0))
})

test_that("fractional bootstrap weights are mean-one cluster weights", {
  data <- data.frame(
    id = c("a", "a", "b", "c", "c"),
    value = 1:5
  )

  withr::local_seed(123)
  boot_a <- markov.misc:::generate_fwb_bootstrap_weights(
    data,
    id_var = "id",
    n_boot = 2
  )
  withr::local_seed(123)
  boot_b <- markov.misc:::generate_fwb_bootstrap_weights(
    data,
    id_var = "id",
    n_boot = 2
  )

  expect_equal(boot_a, boot_b)
  expect_length(boot_a, 2)
  expect_named(boot_a[[1]], c("original_id", "fwb_weight", "boot_id"))
  expect_equal(boot_a[[1]]$original_id, c("a", "b", "c"))
  expect_true(all(boot_a[[1]]$fwb_weight > 0))
  expect_equal(mean(boot_a[[1]]$fwb_weight), 1, tolerance = 1e-12)
  expect_equal(boot_a[[2]]$boot_id, rep(2, 3))
})

test_that("fractional bootstrap materialization expands cluster weights to rows", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(1, 2, 1, 2),
    y = c(0, 1, 1, 0)
  )
  fwb_weights <- data.frame(
    original_id = c(1, 2),
    fwb_weight = c(0.25, 1.75),
    boot_id = 1
  )

  result <- markov.misc:::materialize_fwb_bootstrap_sample(
    fwb_weights,
    data,
    "id"
  )

  expect_equal(nrow(result), nrow(data))
  expect_equal(result$fwb_weight, c(0.25, 0.25, 1.75, 1.75))
  expect_equal(result$boot_id, rep(1, 4))
  expect_equal(
    markov.misc:::fwb_baseline_weights(
      fwb_weights,
      data[!duplicated(data$id), ],
      "id"
    ),
    c(0.25, 1.75)
  )
})

test_that("apply_to_bootstrap materializes samples before sequential analysis", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    y = c(1, 2, 10, 20)
  )
  boot_samples <- list(
    data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1),
    data.frame(original_id = c(2, 2), new_id = c("2_1", "2_2"), boot_id = 2)
  )

  result <- apply_to_bootstrap(
    boot_samples = boot_samples,
    analysis_fn = function(boot_data) {
      c(rows = nrow(boot_data), total = sum(boot_data$y))
    },
    data = data,
    id_var = "id",
    workers = 1
  )

  expect_equal(result[[1]], c(rows = 4, total = 33))
  expect_equal(result[[2]], c(rows = 4, total = 60))
})

test_that("apply_to_fwb_bootstrap materializes weights before analysis", {
  data <- data.frame(
    id = c(1, 1, 2, 2),
    y = c(1, 2, 10, 20)
  )
  fwb_samples <- list(
    data.frame(original_id = c(1, 2), fwb_weight = c(0.5, 1.5), boot_id = 1),
    data.frame(original_id = c(1, 2), fwb_weight = c(1.2, 0.8), boot_id = 2)
  )

  result <- markov.misc:::apply_to_fwb_bootstrap(
    fwb_samples = fwb_samples,
    analysis_fn = function(boot_data, fwb_weights) {
      c(
        rows = nrow(boot_data),
        total = sum(boot_data$y * boot_data$fwb_weight),
        n_weights = nrow(fwb_weights)
      )
    },
    data = data,
    id_var = "id",
    workers = 1
  )

  expect_equal(result[[1]], c(rows = 4, total = 46.5, n_weights = 2))
  expect_equal(result[[2]], c(rows = 4, total = 27.6, n_weights = 2))
})

test_that("apply_to_fwb_bootstrap works with parallel workers", {
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))

  data <- data.frame(
    id = c(1, 1, 2, 2),
    y = c(1, 2, 10, 20)
  )
  fwb_samples <- list(
    data.frame(original_id = c(1, 2), fwb_weight = c(0.5, 1.5), boot_id = 1)
  )

  result <- suppressWarnings(
    markov.misc:::apply_to_fwb_bootstrap(
      fwb_samples = fwb_samples,
      analysis_fn = function(boot_data, fwb_weights) {
        c(
          total = sum(boot_data$y * boot_data$fwb_weight),
          n_weights = nrow(fwb_weights)
        )
      },
      data = data,
      id_var = "id",
      workers = 2,
      packages = character(0)
    )
  )

  expect_equal(result[[1]], c(total = 46.5, n_weights = 2))
  expect_equal(future::nbrOfWorkers(), 1L)
})

test_that("apply_to_bootstrap maps deprecated parallel argument to workers", {
  data <- data.frame(id = c(1, 2), y = c(3, 4))
  boot_samples <- list(
    data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1)
  )

  expect_warning(
    result <- apply_to_bootstrap(
      boot_samples = boot_samples,
      analysis_fn = function(boot_data) sum(boot_data$y),
      data = data,
      id_var = "id",
      parallel = FALSE
    ),
    "parallel.*deprecated"
  )
  expect_equal(result[[1]], 7)
})

test_that("apply_to_bootstrap restores future plan after parallel success", {
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))

  data <- data.frame(id = c(1, 2), y = c(3, 4))
  boot_samples <- list(
    data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1)
  )

  result <- suppressWarnings(
    apply_to_bootstrap(
      boot_samples = boot_samples,
      analysis_fn = function(boot_data) sum(boot_data$y),
      data = data,
      id_var = "id",
      workers = 2,
      packages = character(0)
    )
  )

  expect_equal(result[[1]], 7)
  expect_equal(future::nbrOfWorkers(), 1L)
})

test_that("apply_to_bootstrap restores future plan after parallel errors", {
  future::plan(future::sequential)
  withr::defer(future::plan(future::sequential))

  data <- data.frame(id = c(1, 2), y = c(3, 4))
  boot_samples <- list(
    data.frame(original_id = c(1, 2), new_id = c("1_1", "2_1"), boot_id = 1)
  )

  expect_error(
    suppressWarnings(
      apply_to_bootstrap(
        boot_samples = boot_samples,
        analysis_fn = function(boot_data) stop("boom"),
        data = data,
        id_var = "id",
        workers = 2,
        packages = character(0)
      )
    ),
    "boom"
  )
  expect_equal(future::nbrOfWorkers(), 1L)
})

test_that("bootstrap_analysis_wrapper relevels data and refits models", {
  original_data <- data.frame(
    y = factor(c(1, 2, 3)),
    x = c(0, 1, 2)
  )
  boot_data <- data.frame(
    y = factor(c(1, 3, 3), levels = c(1, 2, 3)),
    x = c(0, 2, 3)
  )
  model <- lm(as.numeric(y) ~ x, data = original_data)

  result <- bootstrap_analysis_wrapper(
    boot_data = boot_data,
    model = model,
    factor_cols = "y",
    original_data = original_data,
    y_levels = 1:3,
    absorb = 3,
    update_datadist = FALSE
  )

  expect_s3_class(result$model, "lm")
  expect_equal(levels(result$data$y), c("1", "2"))
  expect_equal(result$ylevels, 1:2)
  expect_equal(unname(result$absorb), 2)
  expect_equal(result$missing_states, "2")
})

test_that("relevel_factors_consecutive preserves numeric previous state", {
  original_data <- data.frame(
    y = ordered(c(1, 2, 3, 4)),
    yprev = c(1, 2, 3, 4)
  )
  boot_data <- data.frame(
    y = ordered(c(1, 3, 4), levels = 1:4),
    yprev = c(1, 3, 4)
  )

  result <- relevel_factors_consecutive(
    data = boot_data,
    factor_cols = c("y", "yprev"),
    original_data = original_data,
    y_levels = 1:4,
    absorb = 4
  )

  expect_s3_class(result$data$y, "ordered")
  expect_type(result$data$yprev, "double")
  expect_equal(result$data$yprev, c(1, 2, 3))
  expect_equal(result$ylevels, 1:3)
  expect_equal(unname(result$absorb), 3)
  expect_equal(result$missing_levels, "2")
})

test_that("bootstrap_analysis_wrapper reports model fitting failures", {
  boot_data <- data.frame(y = factor(c(1, 2)), x = c(1, 2))
  model <- structure(list(), class = "not_updateable")

  expect_warning(
    result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = "y",
      original_data = boot_data,
      update_datadist = FALSE
    ),
    "Bootstrap model fitting failed"
  )

  expect_null(result$model)
  expect_equal(result$missing_states, character(0))
})

test_that("bootstrap_analysis_wrapper scopes orm datadist updates", {
  skip_if_not_installed("rms")

  update.orm <- function(object, data, ...) {
    structure(list(data = data), class = "orm")
  }
  assign("update.orm", update.orm, envir = globalenv())
  withr::defer(rm("update.orm", envir = globalenv()))

  old_dd <- structure(list(source = "existing"), class = "datadist")
  assign("dd", old_dd, envir = globalenv())
  options(datadist = "existing_dd")
  withr::defer(options(datadist = NULL))
  withr::defer(rm("dd", envir = globalenv()))

  boot_data <- data.frame(
    y = factor(c(1, 2, 3)),
    x = c(0, 1, 2),
    boot_id = 1,
    fwb_weight = c(0.5, 1, 1.5),
    .markov_misc_fit_weight = c(0.5, 1, 1.5)
  )
  result <- expect_no_warning(
    bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = structure(list(), class = "orm"),
      factor_cols = "y",
      original_data = boot_data,
      update_datadist = TRUE
    )
  )

  expect_s3_class(result$model, "orm")
  expect_equal(getOption("datadist"), "existing_dd")
  expect_identical(get("dd", envir = globalenv()), old_dd)
})

test_that("update_bootstrap_model suppresses only known orm weight warning", {
  update.orm <- function(object, data, weights = NULL, ...) {
    warning(
      "currently weights are ignored in model validation and bootstrapping orm fits",
      call. = FALSE
    )
    structure(list(data = data), class = "orm")
  }
  assign("update.orm", update.orm, envir = globalenv())
  withr::defer(rm("update.orm", envir = globalenv()))

  boot_data <- data.frame(y = 1:3, .markov_misc_fit_weight = c(0.5, 1, 1.5))

  result <- expect_no_warning(
    markov.misc:::update_bootstrap_model(
      structure(list(), class = "orm"),
      boot_data,
      fit_weights = boot_data$.markov_misc_fit_weight
    )
  )
  expect_s3_class(result, "orm")

  update.orm <- function(object, data, weights = NULL, ...) {
    warning("important warning", call. = FALSE)
    structure(list(data = data), class = "orm")
  }
  assign("update.orm", update.orm, envir = globalenv())

  expect_warning(
    markov.misc:::update_bootstrap_model(
      structure(list(), class = "orm"),
      boot_data,
      fit_weights = boot_data$.markov_misc_fit_weight
    ),
    "important warning"
  )
})

test_that("bootstrap_analysis_wrapper falls back when coefstart update fails", {
  coef.vglm <- function(object, ...) object$coefficients
  update.vglm <- function(object, data, coefstart = NULL, ...) {
    if (!is.null(coefstart)) {
      stop("coefstart failed")
    }
    structure(
      list(data = data, coefficients = object$coefficients),
      class = "vglm"
    )
  }
  assign("coef.vglm", coef.vglm, envir = globalenv())
  assign("update.vglm", update.vglm, envir = globalenv())
  withr::defer(rm("coef.vglm", envir = globalenv()))
  withr::defer(rm("update.vglm", envir = globalenv()))

  boot_data <- data.frame(y = factor(c(1, 2)), x = c(0, 1))
  model <- structure(list(coefficients = c(a = 1)), class = "vglm")

  result <- bootstrap_analysis_wrapper(
    boot_data = boot_data,
    model = model,
    factor_cols = "y",
    original_data = boot_data,
    update_datadist = FALSE,
    use_coefstart = TRUE
  )

  expect_s3_class(result$model, "vglm")
  expect_equal(result$model$data, boot_data)
})

test_that("bootstrap_analysis_wrapper multiplies original and fractional weights", {
  update.mock_weighted <- function(object, data, weights = NULL, ...) {
    weight_expr <- substitute(weights)
    structure(
      list(
        data = data,
        weights = eval(weight_expr, data, parent.frame())
      ),
      class = "mock_weighted"
    )
  }
  assign("update.mock_weighted", update.mock_weighted, envir = globalenv())
  withr::defer(rm("update.mock_weighted", envir = globalenv()))

  original_data <- data.frame(
    y = factor(c(1, 2, 3)),
    x = c(0, 1, 2),
    w = c(2, 3, 4)
  )
  boot_data <- original_data
  model <- structure(
    list(call = quote(mock_fit(y ~ x, data = original_data, weights = w))),
    class = "mock_weighted"
  )

  result <- bootstrap_analysis_wrapper(
    boot_data = boot_data,
    model = model,
    factor_cols = "y",
    original_data = original_data,
    update_datadist = FALSE,
    fit_weights = c(0.5, 1, 1.5)
  )

  expect_s3_class(result$model, "mock_weighted")
  expect_equal(result$model$weights, c(1, 3, 6))
  expect_equal(result$model$data$.markov_misc_fit_weight, c(1, 3, 6))
})

test_that("bootstrap_analysis_wrapper recovers stored fit weights", {
  make_weighted_lm <- function() {
    original_data <- data.frame(
      y = c(1, 2, 3, 4),
      x = c(0, 1, 2, 3)
    )
    original_weights <- c(2, 3, 4, 5)

    list(
      model = stats::lm(
        y ~ x,
        data = original_data,
        weights = original_weights
      ),
      data = original_data
    )
  }

  weighted_fit <- make_weighted_lm()

  result <- expect_no_warning(
    bootstrap_analysis_wrapper(
      boot_data = weighted_fit$data,
      model = weighted_fit$model,
      factor_cols = character(0),
      original_data = weighted_fit$data,
      update_datadist = FALSE,
      fit_weights = c(0.5, 1, 1.5, 2)
    )
  )

  expect_s3_class(result$model, "lm")
  expect_equal(result$model$weights, c(1, 3, 6, 10))
  expect_equal(
    result$model$model[["(weights)"]],
    c(1, 3, 6, 10)
  )
})
