test_that("bootstrap_model_coefs validates inputs", {
  expect_error(
    bootstrap_model_coefs(
      lm(mpg ~ wt, data = mtcars),
      data = mtcars,
      n_boot = 1
    ),
    "model must be an orm object",
    fixed = TRUE
  )

  model <- structure(list(), class = "vglm")
  expect_error(
    bootstrap_model_coefs(model, data = NULL, n_boot = 1),
    "No data provided",
    fixed = TRUE
  )

  expect_error(
    bootstrap_model_coefs(model, data = data.frame(patient = 1:2), n_boot = 1),
    "id_var 'id' not found in data",
    fixed = TRUE
  )

  expect_warning(
    expect_error(
      bootstrap_model_coefs(
        lm(mpg ~ wt, data = mtcars),
        data = mtcars,
        n_boot = 1,
        parallel = TRUE
      ),
      "model must be an orm object",
      fixed = TRUE
    ),
    "parallel.*deprecated"
  )
})

test_that("bootstrap_model_coefs orchestrates bootstrap coefficient extraction", {
  model <- structure(
    list(coefficients = c("(Intercept)" = 1, tx = -0.2)),
    class = "vglm"
  )
  data <- data.frame(
    id = c(1, 1, 2, 2),
    y = c(0, 1, 0, 1),
    tx = c(0, 0, 1, 1),
    yprev = factor(c(0, 0, 0, 0))
  )

  with_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      expect_equal(id_var, "id")
      replicate(
        n_boot,
        data.frame(
          original_id = c(1, 2),
          new_id = c("1_1", "2_1"),
          boot_id = 1
        ),
        simplify = FALSE
      )
    },
    apply_to_bootstrap = function(
      boot_samples,
      analysis_fn,
      data,
      id_var,
      workers,
      packages,
      globals
    ) {
      expect_length(boot_samples, 2)
      expect_equal(workers, 1)
      list(
        list("(Intercept)" = 1, tx = -0.1),
        list("(Intercept)" = 2, tx = -0.3)
      )
    },
    {
      expect_warning(
        result <- bootstrap_model_coefs(
          model,
          data = data,
          n_boot = 2,
          workers = 1,
          parallel = FALSE
        ),
        "parallel.*deprecated"
      )
    }
  )

  expect_equal(result$boot_id, 1:2)
  expect_equal(result$tx, c(-0.1, -0.3))
})

test_that("bootstrap_model_coefs preserves numeric yprev before bootstrapping", {
  model <- structure(
    list(coefficients = c("(Intercept)" = 1, tx = -0.2)),
    class = "vglm"
  )
  data <- data.frame(
    id = c(1, 2),
    y = c(0, 1),
    tx = c(0, 1),
    yprev = c(0, 1)
  )

  observed_is_factor <- NULL
  with_mocked_bindings(
    fast_group_bootstrap = function(data, id_var, n_boot) {
      observed_is_factor <<- is.factor(data$yprev)
      list(data.frame(
        original_id = c(1, 2),
        new_id = c("1_1", "2_1"),
        boot_id = 1
      ))
    },
    apply_to_bootstrap = function(...) list(list(tx = 1)),
    {
      result <- bootstrap_model_coefs(model, data = data, n_boot = 1)
    }
  )

  expect_false(observed_is_factor)
  expect_equal(result$tx, 1)
})

test_that("bootstrap_model_coefs records coefficients and failed refits", {
  model <- structure(
    list(coefficients = c("(Intercept)" = 1, tx = 2)),
    class = "orm"
  )
  data <- data.frame(
    id = c(1, 1, 2, 2),
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
    apply_to_bootstrap = function(
      boot_samples,
      analysis_fn,
      data,
      id_var,
      workers,
      packages,
      globals
    ) {
      lapply(boot_samples, function(sample) analysis_fn(data))
    },
    bootstrap_analysis_wrapper = function(...) {
      call_id <<- call_id + 1
      if (call_id == 2) {
        list(model = NULL)
      } else {
        list(model = model)
      }
    },
    {
      result <- bootstrap_model_coefs(
        model,
        data = data,
        n_boot = 2,
        workers = 1,
        use_coefstart = TRUE
      )
    }
  )

  expect_equal(result$tx, c(2, NA))
})

test_that("tidy_bootstrap_coefs summarizes coefficient draws", {
  boot_coefs <- data.frame(
    boot_id = 1:4,
    tx = c(0.1, 0.2, 0.3, 0.4),
    age = c(1, 2, NA, 4)
  )

  result <- tidy_bootstrap_coefs(
    boot_coefs,
    probs = c(0.75, 0.25),
    estimate = "mean",
    na.rm = TRUE
  )

  expect_equal(result$term, c("tx", "age"))
  expect_equal(result$estimate, c(0.25, mean(c(1, 2, 4))))
  expect_equal(result$n_iter, c(4, 4))

  no_estimate <- tidy_bootstrap_coefs(boot_coefs, estimate = NULL)
  expect_false("estimate" %in% names(no_estimate))
})

test_that("tidy_bootstrap_coefs validates inputs", {
  expect_error(
    tidy_bootstrap_coefs(list()),
    "boot_coefs must be a data frame",
    fixed = TRUE
  )
  expect_error(
    tidy_bootstrap_coefs(data.frame(iter = 1, tx = 1)),
    "id_col 'boot_id' not found",
    fixed = TRUE
  )
  expect_error(
    tidy_bootstrap_coefs(data.frame(boot_id = 1, tx = 1), probs = 0.1),
    "probs must be a numeric vector of length 2",
    fixed = TRUE
  )
  expect_error(
    tidy_bootstrap_coefs(data.frame(boot_id = 1, tx = 1), estimate = "mode"),
    "estimate must be",
    fixed = TRUE
  )
  expect_error(
    tidy_bootstrap_coefs(
      data.frame(boot_id = 1, tx = 1),
      na.rm = c(TRUE, FALSE)
    ),
    "na.rm must be a single logical value",
    fixed = TRUE
  )
})
