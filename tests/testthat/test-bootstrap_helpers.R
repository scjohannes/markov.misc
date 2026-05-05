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
    ylevels = 1:3,
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
    ylevels = 1:4,
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
