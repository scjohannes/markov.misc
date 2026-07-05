test_that("tidy_po handles glm alternatives and excludes intercepts", {
  data <- data.frame(
    y = c(0, 0, 1, 1, 0, 1),
    tx = c(0, 0, 1, 1, 1, 0),
    x = c(1, 2, 1, 2, 3, 3)
  )
  fit <- glm(y ~ tx + x, data = data, family = binomial())

  two_sided <- tidy_po(fit)
  greater <- tidy_po(fit, alternative = "greater")
  less <- tidy_po(fit, alternative = "less")

  expect_equal(two_sided$term, c("tx", "x"))
  expect_equal(unname(two_sided$estimate), unname(coef(fit)[c("tx", "x")]))
  expect_equal(two_sided$statistic, two_sided$estimate / two_sided$std_error)
  expect_true(all(is.finite(two_sided$conf_low)))
  expect_true(all(is.finite(two_sided$conf_high)))
  expect_equal(unname(greater$conf_high), rep(Inf, 2))
  expect_equal(unname(less$conf_low), rep(-Inf, 2))
})

test_that("sample_from_arrow samples balanced treatment groups", {
  skip_if_not_installed("arrow")

  path <- withr::local_tempdir()
  data <- data.frame(
    id = rep(1:8, each = 2),
    tx = rep(rep(c(0, 1), each = 4), each = 2),
    y = seq_len(16)
  )
  arrow::write_dataset(data, path)

  sampled <- sample_from_arrow(
    data_path = path,
    sample_size = 4,
    allocation_ratio = 0.5,
    seed = 1
  )

  expect_equal(length(unique(sampled$id)), 4)
  expect_equal(as.integer(table(unique(sampled[c("id", "tx")])$tx)), c(2L, 2L))
})

test_that("sample_from_arrow samples with replacement using synthetic IDs", {
  skip_if_not_installed("arrow")

  path <- withr::local_tempdir()
  data <- data.frame(
    id = rep(1:4, each = 2),
    tx = rep(rep(c(0, 1), each = 2), each = 2),
    y = seq_len(8)
  )
  arrow::write_dataset(data, path)

  sampled <- sample_from_arrow(
    data_path = path,
    sample_size = 4,
    allocation_ratio = 0.5,
    seed = 2,
    replace = TRUE
  )

  expect_equal(nrow(sampled), 8)
  expect_equal(unique(sampled$id), c("3_1", "3_2", "2_1", "2_2"))
  expect_equal(as.integer(table(sampled$id)), rep(2L, 4))
  arm_by_id <- tapply(sampled$tx, sampled$id, unique)
  expect_equal(
    as.vector(arm_by_id[c("3_1", "3_2", "2_1", "2_2")]),
    c(1L, 1L, 0L, 0L)
  )
  expect_false("new_id" %in% names(sampled))
  expect_false("boot_id" %in% names(sampled))
})

test_that("sample_from_arrow samples with replacement beyond observed IDs", {
  skip_if_not_installed("arrow")

  path <- withr::local_tempdir()
  data <- data.frame(
    id = rep(1:4, each = 2),
    tx = rep(rep(c(0, 1), each = 2), each = 2),
    y = seq_len(8)
  )
  arrow::write_dataset(data, path)

  sampled <- sample_from_arrow(
    data_path = path,
    sample_size = 6,
    allocation_ratio = 0.5,
    seed = 1,
    replace = TRUE
  )

  expect_equal(length(unique(sampled$id)), 6)
  expect_equal(sum(sampled$tx == 1), 6)
  expect_equal(sum(sampled$tx == 0), 6)
  expect_false("new_id" %in% names(sampled))
  expect_false("boot_id" %in% names(sampled))
})

test_that("sample_from_arrow errors without replacement beyond observed IDs", {
  skip_if_not_installed("arrow")

  path <- withr::local_tempdir()
  data <- data.frame(
    id = rep(1:4, each = 2),
    tx = rep(rep(c(0, 1), each = 2), each = 2),
    y = seq_len(8)
  )
  arrow::write_dataset(data, path)

  expect_error(
    sample_from_arrow(
      data_path = path,
      sample_size = 6,
      allocation_ratio = 0.5,
      replace = FALSE
    ),
    "Requested 3 treatment patients, but only 2 are available.",
    fixed = TRUE
  )
})

test_that("sample_from_arrow validates sampling controls", {
  skip_if_not_installed("arrow")

  path <- withr::local_tempdir()
  data <- data.frame(id = 1:2, tx = c(0, 1), y = 1:2)
  arrow::write_dataset(data, path)

  expect_error(
    sample_from_arrow(path, sample_size = 1.5),
    "`sample_size` must be a single positive integer.",
    fixed = TRUE
  )
  expect_error(
    sample_from_arrow(path, sample_size = 2, allocation_ratio = 1.5),
    "`allocation_ratio` must be a single number between 0 and 1.",
    fixed = TRUE
  )
  expect_error(
    sample_from_arrow(path, sample_size = 2, replace = NA),
    "`replace` must be `TRUE` or `FALSE`.",
    fixed = TRUE
  )
})

test_that("sample_from_arrow supports custom ID and treatment variables", {
  skip_if_not_installed("arrow")

  path <- withr::local_tempdir()
  data <- data.frame(
    patient_id = rep(1:8, each = 2),
    arm = factor(rep(rep(c("usual_care", "active"), each = 4), each = 2)),
    y = seq_len(16)
  )
  arrow::write_dataset(data, path)

  expect_error(
    sample_from_arrow(
      data_path = path,
      sample_size = 4,
      allocation_ratio = 0.5,
      id_var = "patient_id",
      tx_var = "arm"
    ),
    "pass `control_value` and `treatment_value`",
    fixed = TRUE
  )

  sampled <- sample_from_arrow(
    data_path = path,
    sample_size = 4,
    allocation_ratio = 0.5,
    seed = 1,
    id_var = "patient_id",
    tx_var = "arm",
    control_value = "usual_care",
    treatment_value = "active"
  )

  expect_equal(length(unique(sampled$patient_id)), 4)
  expect_equal(
    as.integer(table(unique(sampled[c("patient_id", "arm")])$arm)),
    c(2L, 2L)
  )
})

test_that("tidy_po handles vglm, orm, clm, and coxph models", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  skip_if_not_installed("ordinal")
  skip_if_not_installed("survival")

  data <- data.frame(
    y_ord = ordered(rep(1:3, times = 8)),
    y_bin = rep(c(0, 1), times = 12),
    tx = rep(c(0, 1), each = 12),
    x = rep(c(-1, 0, 1), times = 8),
    time = seq_len(24) + 1,
    status = rep(c(1, 0, 1, 1), length.out = 24)
  )

  vglm_fit <- suppressWarnings(
    VGAM::vglm(
      y_ord ~ tx + x,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  )
  vglm_tidy <- tidy_po(vglm_fit)
  expect_setequal(vglm_tidy$term, c("tx", "x"))

  dd <- rms::datadist(data)
  withr::local_options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  withr::defer(rm("dd", envir = globalenv()))

  orm_fit <- rms::orm(y_ord ~ tx + x, data = data, x = TRUE, y = TRUE)
  orm_tidy <- tidy_po(orm_fit)
  orm_robust_tidy <- tidy_po(rms::robcov(orm_fit))
  expect_setequal(orm_tidy$term, c("tx", "x"))
  expect_setequal(orm_robust_tidy$term, c("tx", "x"))

  clm_fit <- suppressWarnings(ordinal::clm(y_ord ~ tx + x, data = data))
  clm_tidy <- tidy_po(clm_fit)
  expect_setequal(clm_tidy$term, c("tx", "x"))

  cox_fit <- suppressWarnings(
    survival::coxph(survival::Surv(time, status) ~ tx + x, data = data)
  )
  cox_tidy <- tidy_po(cox_fit)
  expect_setequal(cox_tidy$term, c("tx", "x"))
})

test_that("tidy_po warns for rank-deficient coxph coefficients", {
  skip_if_not_installed("survival")

  data <- data.frame(
    time = 1:8,
    status = c(1, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 1, 0, 1, 0, 1, 0, 1)
  )
  fit <- survival::coxph(
    survival::Surv(time, status) ~ x + duplicate_x,
    data = transform(data, duplicate_x = x),
    singular.ok = TRUE
  )

  expect_warning(
    out <- tidy_po(fit),
    "Some coefficients in coxph model are NA",
    fixed = TRUE
  )
  expect_true(any(is.na(out$estimate)))
})

test_that("tidy_po rejects unsupported model objects", {
  expect_error(
    tidy_po(list()),
    "fit must be a vglm, orm, clm, coxph, or glm object",
    fixed = TRUE
  )
})

test_that("assess_operating_characteristics combines summaries and saves details", {
  output_path <- withr::local_tempdir()
  sampled <- data.frame(
    id = rep(1:4, each = 2),
    tx = rep(c(0, 1), each = 4),
    y = seq_len(8)
  )

  with_mocked_bindings(
    sample_from_arrow = function(
      data_path,
      sample_size,
      allocation_ratio,
      seed = NULL,
      replace = FALSE,
      id_var = "id",
      tx_var = "tx",
      control_value = 0,
      treatment_value = 1
    ) {
      expect_equal(data_path, "mock-path")
      expect_equal(sample_size, 4)
      expect_equal(allocation_ratio, 0.5)
      expect_equal(replace, FALSE)
      expect_equal(id_var, "id")
      expect_equal(tx_var, "tx")
      expect_equal(control_value, 0)
      expect_equal(treatment_value, 1)
      sampled
    },
    {
      result <- assess_operating_characteristics(
        iter_num = 2,
        data_paths = list(mock = "mock-path"),
        output_path = output_path,
        fit_functions = list(
          mock = function(data, iter) {
            list(
              data.frame(
                iter = iter,
                analysis = "mock",
                se_type = "naive",
                term = "tx",
                estimate = mean(data$tx),
                std_error = 0.1,
                statistic = 5,
                p_value = 0.01,
                conf_low = 0.2,
                conf_high = 0.8
              ),
              data.frame(saved = seq_len(2))
            )
          }
        ),
        sample_size = 4
      )
    }
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$iter, 2)
  expect_equal(result$analysis, "mock")
  detail_file <- file.path(output_path, "details", "mock", "mock_iter_2.rds")
  expect_true(file.exists(detail_file))
  expect_equal(readRDS(detail_file)[[1]]$saved, 1:2)
})

test_that("assess_operating_characteristics passes custom variable names to sampler", {
  output_path <- withr::local_tempdir()
  sampled <- data.frame(
    patient_id = rep(1:4, each = 2),
    arm = rep(c("usual_care", "active"), each = 4),
    y = seq_len(8)
  )

  with_mocked_bindings(
    sample_from_arrow = function(
      data_path,
      sample_size,
      allocation_ratio,
      seed = NULL,
      replace = FALSE,
      id_var = "id",
      tx_var = "tx",
      control_value = 0,
      treatment_value = 1
    ) {
      expect_equal(replace, FALSE)
      expect_equal(id_var, "patient_id")
      expect_equal(tx_var, "arm")
      expect_equal(control_value, "usual_care")
      expect_equal(treatment_value, "active")
      sampled
    },
    {
      result <- assess_operating_characteristics(
        iter_num = 2,
        data_paths = list(mock = "mock-path"),
        output_path = output_path,
        fit_functions = list(
          mock = function(data, iter) {
            list(data.frame(
              iter = iter,
              analysis = "mock",
              se_type = "naive",
              term = "arm",
              estimate = mean(data$arm == "active"),
              std_error = 0.1,
              statistic = 5,
              p_value = 0.01,
              conf_low = 0.2,
              conf_high = 0.8
            ))
          }
        ),
        sample_size = 4,
        id_var = "patient_id",
        tx_var = "arm",
        control_value = "usual_care",
        treatment_value = "active"
      )
    }
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$term, "arm")
})

test_that("assess_operating_characteristics passes replacement sampling to sampler", {
  output_path <- withr::local_tempdir()
  sampled <- data.frame(
    id = rep(c("1_1", "1_2", "2_1", "2_2"), each = 2),
    tx = rep(c(0, 1), each = 4),
    y = seq_len(8)
  )

  with_mocked_bindings(
    sample_from_arrow = function(
      data_path,
      sample_size,
      allocation_ratio,
      seed = NULL,
      replace = FALSE,
      id_var = "id",
      tx_var = "tx",
      control_value = 0,
      treatment_value = 1
    ) {
      expect_equal(replace, TRUE)
      sampled
    },
    {
      result <- assess_operating_characteristics(
        iter_num = 1,
        data_paths = list(mock = "mock-path"),
        output_path = output_path,
        fit_functions = list(
          mock = function(data, iter) {
            list(data.frame(
              iter = iter,
              analysis = "mock",
              se_type = "naive",
              term = "tx",
              estimate = mean(data$tx),
              std_error = 0.1,
              statistic = 5,
              p_value = 0.01,
              conf_low = 0.2,
              conf_high = 0.8
            ))
          }
        ),
        sample_size = 4,
        replace = TRUE
      )
    }
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$analysis, "mock")
})

test_that("assess_operating_characteristics can rerandomize sampled treatment", {
  output_path <- withr::local_tempdir()
  sampled <- data.frame(
    id = rep(1:6, each = 2),
    tx = "usual_care",
    y = seq_len(12)
  )

  observed_tx <- NULL
  with_mocked_bindings(
    sample_from_arrow = function(...) sampled,
    {
      result <- assess_operating_characteristics(
        iter_num = 1,
        data_paths = list(mock = "mock-path"),
        output_path = output_path,
        fit_functions = list(
          mock = function(data, iter) {
            observed_tx <<- tapply(data$tx, data$id, unique)
            list(data.frame(
              iter = iter,
              analysis = "mock",
              se_type = "naive",
              term = "tx",
              estimate = mean(data$tx == "active"),
              std_error = 1,
              statistic = 0,
              p_value = 1,
              conf_low = -1,
              conf_high = 1
            ))
          }
        ),
        sample_size = 6,
        allocation_ratio = 0.5,
        rerandomize = TRUE,
        control_value = "usual_care",
        treatment_value = "active"
      )
    }
  )

  expect_equal(nrow(result), 1)
  expect_equal(sum(observed_tx == "active"), 3)
  expect_setequal(unique(observed_tx), c("usual_care", "active"))
})

test_that("assess_operating_characteristics validates output and inputs", {
  expect_error(
    assess_operating_characteristics(
      iter_num = 1,
      data_paths = list(),
      output_path = NULL,
      fit_functions = list()
    ),
    "output_path must be provided",
    fixed = TRUE
  )

  with_mocked_bindings(
    sample_from_arrow = function(...) data.frame(id = 1, tx = 0),
    {
      expect_error(
        assess_operating_characteristics(
          iter_num = 1,
          data_paths = list(mock = "mock-path"),
          output_path = withr::local_tempdir(),
          fit_functions = list(mock = function(data, iter) data.frame(x = 1))
        ),
        "Fitting function must return a list",
        fixed = TRUE
      )
    }
  )
})

test_that("assess_operating_characteristics warns and skips missing data paths", {
  expect_warning(
    result <- assess_operating_characteristics(
      iter_num = 1,
      data_paths = list(),
      output_path = withr::local_tempdir(),
      fit_functions = list(missing = function(data, iter) {
        list(data.frame(x = 1))
      })
    ),
    "No data path provided for analysis: missing"
  )

  expect_equal(nrow(result), 0)
})

test_that("summarize_oc_results computes operating characteristics with truth", {
  results <- data.frame(
    analysis = rep(c("a", "b"), each = 3),
    term = "tx",
    estimate = c(0.9, 1.1, 1.2, 1.8, 2.1, 2.2),
    std_error = c(0.2, 0.3, 0.4, 0.4, 0.5, 0.6),
    p_value = c(0.01, 0.20, 0.03, 0.20, 0.01, 0.04),
    conf_low = c(0.5, 0.8, 0.9, 1.2, 1.7, 1.9),
    conf_high = c(1.3, 1.4, 1.5, 2.4, 2.5, 2.6)
  )

  summary <- summarize_oc_results(results, true_value = 1, alpha = 0.05)

  expect_equal(summary$analysis, c("a", "b"))
  expect_equal(
    summary$mean_estimate,
    c(mean(c(0.9, 1.1, 1.2)), mean(c(1.8, 2.1, 2.2)))
  )
  expect_equal(summary$bias, summary$mean_estimate - 1)
  expect_equal(summary$power, c(2 / 3, 2 / 3))
  expect_equal(summary$coverage, c(1, 0))
  expect_equal(summary$n_iters, c(3L, 3L))
})

test_that("summarize_oc_results can omit truth-dependent metrics", {
  results <- data.frame(
    analysis = "a",
    term = "tx",
    estimate = c(1, 2),
    std_error = c(0.2, 0.4),
    p_value = c(0.01, 0.2),
    conf_low = c(0, 1),
    conf_high = c(2, 3)
  )

  summary <- summarize_power_results(results)

  expect_equal(summary$coverage, NA_real_)
  expect_false("bias" %in% names(summary))
})
