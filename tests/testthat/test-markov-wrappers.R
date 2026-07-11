test_that("vglm_markov stores fitting data and returns robust wrapper with id_var", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1001)

  fit <- vglm_markov(
    ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id"
  )

  expect_s3_class(fit, "robcov_vglm")
  expect_s4_class(fit$vglm_fit, "vglm")
  expect_equal(attr(fit, "markov_id_var"), "id")
  expect_equal(attr(fit$vglm_fit, "markov_id_var"), "id")
  expect_equal(nrow(attr(fit, "markov_data")), nrow(data))
})

test_that("vglm_markov aligns stored fitting data after subset and NA drops", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1007)
  data$time_lin[2] <- NA_real_

  fit <- vglm_markov(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id",
    subset = time <= 5
  )

  stored_data <- attr(fit, "markov_data")
  expected_data <- data[data$time <= 5 & !is.na(data$time_lin), , drop = FALSE]

  expect_s3_class(fit, "robcov_vglm")
  expect_equal(
    stored_data[c("id", "time", "time_lin")],
    expected_data[c("id", "time", "time_lin")]
  )
  expect_equal(
    attr(fit$vglm_fit, "markov_data")[c("id", "time")],
    expected_data[c("id", "time")]
  )
  expect_equal(nrow(fit$vglm_fit@x), nrow(stored_data))
})

test_that("orm_markov stores fitting data and applies rms robust covariance", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1002)
  dd <- rms::datadist(data)
  old_dd_exists <- exists("dd", envir = globalenv(), inherits = FALSE)
  old_dd <- if (old_dd_exists) get("dd", envir = globalenv()) else NULL
  assign("dd", dd, envir = globalenv())
  old_options <- options(datadist = "dd")
  on.exit(
    {
      options(old_options)
      if (old_dd_exists) {
        assign("dd", old_dd, envir = globalenv())
      } else if (exists("dd", envir = globalenv(), inherits = FALSE)) {
        remove(list = "dd", envir = globalenv())
      }
    },
    add = TRUE
  )

  fit <- orm_markov(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id"
  )

  expect_s3_class(fit, "orm")
  expect_equal(attr(fit, "markov_id_var"), "id")
  expect_equal(nrow(attr(fit, "markov_data")), nrow(data))
  expect_false(is.null(fit$orig.var))

  subset_fit <- expect_no_error(
    orm_markov(
      ordered(y) ~ time + tx + yprev,
      data = data,
      id_var = "id",
      subset = time <= 4
    )
  )
  expect_equal(nrow(attr(subset_fit, "markov_data")), sum(data$time <= 4))
  expect_true(all(attr(subset_fit, "markov_data")$time <= 4))
})

test_that("blrm_markov stores fitting data without sampling when requested", {
  skip_if_not_installed("rms")
  skip_if_not_installed("rmsb")

  data <- make_test_data(n_patients = 12, follow_up_time = 6, seed = 1005)
  dd <- rms::datadist(data)
  old_dd_exists <- exists("dd", envir = globalenv(), inherits = FALSE)
  old_dd <- if (old_dd_exists) get("dd", envir = globalenv()) else NULL
  assign("dd", dd, envir = globalenv())
  old_options <- options(datadist = "dd")
  on.exit(
    {
      options(old_options)
      if (old_dd_exists) {
        assign("dd", old_dd, envir = globalenv())
      } else if (exists("dd", envir = globalenv(), inherits = FALSE)) {
        remove(list = "dd", envir = globalenv())
      }
    },
    add = TRUE
  )

  fit <- blrm_markov(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id",
    subset = time <= 3,
    standata = TRUE
  )

  expect_type(fit, "list")
  expect_equal(attr(fit, "markov_id_var"), "id")
  expect_equal(nrow(attr(fit, "markov_data")), sum(data$time <= 3))
  expect_true(all(attr(fit, "markov_data")$time <= 3))

  expect_error(
    blrm_markov(
      ordered(y) ~ offset(time) + tx + yprev,
      data = data,
      id_var = "id",
      standata = TRUE
    ),
    "offset"
  )
})

test_that("sops uses stored full data and extracts earliest prediction rows", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1003)
  data <- data[!(data$id == 1 & data$time == min(data$time)), , drop = FALSE]
  fit <- vglm_markov(
    ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id"
  )

  baseline_idx <- vapply(
    split(seq_len(nrow(data)), data$id),
    function(idx) idx[which.min(data$time[idx])],
    integer(1)
  )
  manual_baseline <- data[baseline_idx, , drop = FALSE]
  manual <- sops(
    fit,
    newdata = manual_baseline,
    times = 1:3,
    y_levels = 1:6,
    absorb = 6
  )
  automatic <- sops(
    fit,
    times = 1:3,
    y_levels = 1:6,
    absorb = 6
  )

  compare_cols <- c("rowid", "time", "state", "estimate", "id", "tx", "yprev")
  expect_equal(
    automatic[, compare_cols],
    manual[, compare_cols],
    tolerance = 1e-8
  )
  prediction_data <- attr(automatic, "newdata_pred")
  expect_equal(nrow(prediction_data), length(unique(data$id)))
  expect_equal(
    prediction_data[c("id", "time", "yprev")],
    manual_baseline[c("id", "time", "yprev")]
  )
  expect_equal(prediction_data$time[prediction_data$id == 1], 2)
  expect_equal(nrow(attr(automatic, "refit_data")), nrow(data))
  expect_equal(attr(automatic, "id_var"), "id")
})

test_that("sops treats supplied newdata rows as fixed prediction profiles", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1008)
  fit <- vglm_markov(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id"
  )

  profiles <- data[data$time == 1, , drop = FALSE][seq_len(3), ]
  profiles$id <- c(1, 1, 1)
  profiles$time <- c(1, 2, 3)
  profiles$rowid <- c(99, 99, 99)
  profiles_no_id <- profiles[, setdiff(names(profiles), "id"), drop = FALSE]

  out <- sops(
    fit,
    newdata = profiles_no_id,
    times = 1:2,
    y_levels = 1:6,
    absorb = 6
  )
  prediction_data <- attr(out, "newdata_pred")

  expect_equal(nrow(prediction_data), nrow(profiles_no_id))
  expect_equal(prediction_data$rowid, seq_len(nrow(profiles_no_id)))
  expect_equal(sort(unique(out$rowid)), seq_len(nrow(profiles_no_id)))
  expect_equal(attr(out, "id_var"), "id")
  expect_equal(nrow(attr(out, "refit_data")), nrow(data))

  withr::local_seed(10081)
  expect_warning(
    inferred <- inferences(
      out,
      method = "fwb",
      n_draws = 1,
      return_draws = TRUE,
      use_coefstart = FALSE
    ),
    "fixed prediction profiles",
    fixed = TRUE
  )
  draws <- get_draws(inferred)

  expect_equal(attr(inferred, "n_successful"), 1L)
  expect_equal(sort(unique(draws$rowid)), seq_len(nrow(profiles_no_id)))
})

test_that("avg_sops treats supplied newdata rows as standardization profiles", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1009)
  fit <- vglm_markov(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id"
  )

  profiles <- data[data$time == 1, , drop = FALSE][seq_len(4), ]
  profiles$id <- c(1, 1, 1, 1)
  profiles$time <- c(1, 2, 3, 4)
  profiles$rowid <- 100 + seq_len(nrow(profiles))
  profiles_no_id <- profiles[, setdiff(names(profiles), "id"), drop = FALSE]

  out <- avg_sops(
    fit,
    newdata = profiles_no_id,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    y_levels = 1:6,
    absorb = 6
  )
  prediction_data <- attr(out, "newdata_pred")

  expect_equal(nrow(prediction_data), nrow(profiles_no_id) * 2L)
  expect_equal(
    prediction_data$rowid,
    rep(seq_len(nrow(profiles_no_id)), times = 2L)
  )
  expect_equal(attr(out, "id_var"), "id")
  expect_equal(nrow(attr(out, "refit_data")), nrow(data))

  withr::local_seed(10091)
  expect_warning(
    inferred <- inferences(
      out,
      method = "fwb",
      n_draws = 1,
      return_draws = TRUE,
      use_coefstart = FALSE
    ),
    "baseline weights have been set to NULL"
  )
  expect_equal(attr(inferred, "n_successful"), 1L)
})

test_that("Markov wrappers warn about duplicate id-time rows", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1010)
  data <- rbind(data, data[1, , drop = FALSE])

  expect_warning(
    vglm_markov(
      ordered(y) ~ time_lin + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data,
      id_var = "id"
    ),
    "duplicate `id` x `time` combinations"
  )
})

test_that("individual sops supports FWB but rejects standard bootstrap", {
  skip_if_not_installed("VGAM")

  data <- make_test_data(n_patients = 35, follow_up_time = 6, seed = 1004)
  fit <- vglm_markov(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data,
    id_var = "id"
  )

  object <- sops(fit, times = 1:2, y_levels = 1:6, absorb = 6)

  expect_error(
    inferences(
      object,
      method = "bootstrap",
      n_draws = 1
    ),
    "Standard refit bootstrap is not supported"
  )

  withr::local_seed(10041)
  out <- inferences(
    object,
    method = "fwb",
    n_draws = 2,
    return_draws = TRUE,
    use_coefstart = FALSE
  )

  expect_s3_class(out, "markov_sops")
  expect_equal(attr(out, "method"), "fwb")
  expect_equal(attr(out, "n_successful"), 2)
  expect_true(all(c("conf.low", "conf.high", "std.error") %in% names(out)))
  expect_false(is.null(attr(out, "draws")))
})

test_that("orm_markov supports FWB for grouped and ungrouped sops", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 45, follow_up_time = 6, seed = 1006)
  ids <- unique(data$id)
  sex_by_id <- stats::setNames(
    rep(c("female", "male"), length.out = length(ids)),
    ids
  )
  data$sex <- factor(sex_by_id[as.character(data$id)])

  dd <- rms::datadist(data)
  old_dd_exists <- exists("dd", envir = globalenv(), inherits = FALSE)
  old_dd <- if (old_dd_exists) get("dd", envir = globalenv()) else NULL
  assign("dd", dd, envir = globalenv())
  old_options <- options(datadist = "dd")
  on.exit(
    {
      options(old_options)
      if (old_dd_exists) {
        assign("dd", old_dd, envir = globalenv())
      } else if (exists("dd", envir = globalenv(), inherits = FALSE)) {
        remove(list = "dd", envir = globalenv())
      }
    },
    add = TRUE
  )

  fit <- orm_markov(
    ordered(y) ~ time + tx + yprev,
    data = data,
    id_var = "id"
  )

  grouped <- sops(
    fit,
    by = "sex",
    times = 1:2,
    y_levels = fit$yunique,
    absorb = "6"
  )

  withr::local_seed(10061)
  grouped_ci <- inferences(
    grouped,
    method = "fwb",
    n_draws = 2,
    return_draws = TRUE,
    update_datadist = FALSE
  )

  expect_equal(attr(grouped_ci, "n_successful"), 2L)
  expect_setequal(grouped_ci$sex, levels(data$sex))
  expect_equal(anyNA(grouped_ci$conf.low), FALSE)
  expect_contains(
    names(attr(grouped_ci, "draws")),
    c("sex", "draw_id")
  )

  individual <- sops(
    fit,
    times = 1:2,
    y_levels = fit$yunique,
    absorb = "6"
  )

  individual_ci <- inferences(
    individual,
    method = "fwb",
    n_draws = 2,
    return_draws = TRUE,
    update_datadist = FALSE
  )
  draws <- get_draws(individual_ci)

  expect_equal(attr(individual_ci, "n_successful"), 2L)
  expect_equal(length(unique(draws$rowid)), length(unique(data$id)))
  expect_contains(names(draws), c("draw_id", "rowid", "sex", "draw"))
  expect_equal(anyNA(draws$draw), FALSE)
})
