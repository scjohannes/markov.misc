test_that("Markov workflows reject non-logit vglm links", {
  skip_if_not_installed("VGAM")

  data <- suppressWarnings(
    make_test_data(n_patients = 80, seed = 841, follow_up_time = 5)
  )
  fit <- VGAM::vglm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(
      reverse = TRUE,
      parallel = TRUE,
      link = "probitlink"
    ),
    data = data
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
  fit <- rms::orm(
    ordered(y) ~ time + tx + yprev,
    family = "probit",
    data = data,
    x = TRUE,
    y = TRUE
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
  vglm_fit <- VGAM::vglm(
    ordered(y) ~ time + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
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
  orm_fit <- rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = data,
    x = TRUE,
    y = TRUE
  )

  expect_null(markov.misc:::validate_markov_model(vglm_fit))
  expect_null(markov.misc:::validate_markov_model(orm_fit))
})
