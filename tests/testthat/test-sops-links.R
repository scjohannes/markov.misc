test_that("Markov workflows reject non-logit vglm links", {
  skip_if_not_installed("VGAM")

  data <- suppressWarnings(
    make_test_data(n_patients = 80, seed = 841, follow_up_time = 5)
  )
  fit <- suppressWarnings(
    vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(
        reverse = TRUE,
        parallel = TRUE,
        link = "probitlink"
      ),
      data = data
    )
  )

  error <- tryCatch(
    markov.misc:::validate_markov_model(fit),
    error = identity
  )
  expect_s3_class(error, "markov_misc_unsupported_link")
  expect_equal(
    conditionMessage(error),
    paste0(
      "Only cumulative-logit models are supported; ",
      "refit the model with a logit link."
    )
  )
  expect_equal(error$link, "probitlink")
})

test_that("Markov workflows reject non-logistic orm families", {
  skip_if_not_installed("rms")

  data <- suppressWarnings(
    make_test_data(n_patients = 80, seed = 842, follow_up_time = 5)
  )
  dd_name <- ".markov_misc_link_test_dd"
  assign(dd_name, rms::datadist(data), envir = globalenv())
  old <- options(datadist = dd_name)
  on.exit(
    {
      options(old)
      rm(list = dd_name, envir = globalenv())
    },
    add = TRUE
  )
  fit <- suppressWarnings(
    orm_markov(
      ordered(y) ~ time + tx + yprev,
      family = "probit",
      data = data
    )
  )

  error <- tryCatch(
    markov.misc:::validate_markov_model(fit),
    error = identity
  )
  expect_s3_class(error, "markov_misc_unsupported_link")
  expect_equal(error$link, "probit")
})

test_that("Markov workflows accept logit vglm and orm models", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- suppressWarnings(
    make_test_data(n_patients = 80, seed = 843, follow_up_time = 5)
  )
  vglm_fit <- suppressWarnings(
    vglm_markov(
      ordered(y) ~ time + tx + yprev,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  )
  dd_name <- ".markov_misc_link_accept_dd"
  assign(dd_name, rms::datadist(data), envir = globalenv())
  old <- options(datadist = dd_name)
  on.exit(
    {
      options(old)
      rm(list = dd_name, envir = globalenv())
    },
    add = TRUE
  )
  orm_fit <- suppressWarnings(
    orm_markov(
      ordered(y) ~ time + tx + yprev,
      data = data
    )
  )

  expect_null(markov.misc:::validate_markov_model(vglm_fit))
  expect_null(markov.misc:::validate_markov_model(orm_fit))
})

test_that("Markov workflows reject raw and manually robust backend fits", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- suppressWarnings(
    make_test_data(n_patients = 40, seed = 844, follow_up_time = 5)
  )
  raw_fit <- VGAM::vglm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )
  robust_fit <- robcov_vglm(raw_fit, cluster = data$id)
  raw_orm_fit <- rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    x = TRUE,
    y = TRUE
  )

  expect_error(
    markov.misc:::validate_markov_model(raw_fit),
    "require fits created by `vglm_markov()`",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::validate_markov_model(robust_fit),
    "require fits created by `vglm_markov()`",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::validate_markov_model(raw_orm_fit),
    "require fits created by `orm_markov()`",
    fixed = TRUE
  )
})
