test_that("vglm.markov() supplies rms formula helpers for %ia%", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2468, follow_up_time = 12)
  baseline <- data[data$time == 1, , drop = FALSE]
  stripped_env <- baseenv()

  form <- stats::as.formula(
    "ordered(y) ~ rcs(time, 4) + tx + yprev + time %ia% yprev",
    env = stripped_env
  )
  expect_equal(exists("%ia%", envir = environment(form), inherits = TRUE), FALSE)
  expect_equal(exists("rcs", envir = environment(form), inherits = TRUE), FALSE)

  fit <- suppressWarnings(
    markov.misc::vglm.markov(
      form,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  )

  ia_terms <- grep("%ia%", names(fit@constraints), value = TRUE, fixed = TRUE)
  expect_length(ia_terms, nlevels(data$yprev) - 1)
  expect_equal(grep("time'", ia_terms, fixed = TRUE), integer(0))
  expect_equal(grep("time * yprev=", ia_terms, fixed = TRUE), seq_along(ia_terms))

  out <- markov.misc::soprob_markov(
    fit,
    data = baseline[1:3, , drop = FALSE],
    times = 1:4,
    ylevels = 1:6,
    absorb = "6"
  )

  expect_equal(dim(out), c(3L, 4L, 6L))
  expect_equal(unname(rowSums(out[, 4, ])), rep(1, 3), tolerance = 1e-10)
})

test_that("vglm.markov() keeps explicit rms::%ia% calls working", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2469, follow_up_time = 12)
  form <- stats::as.formula(
    paste(
      "ordered(y) ~ rms::rcs(time, 4) + tx + yprev +",
      "rms::`%ia%`(time, yprev)"
    ),
    env = baseenv()
  )

  fit <- suppressWarnings(
    markov.misc::vglm.markov(
      form,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  )

  ia_terms <- grep("rms::`%ia%`", names(fit@constraints), value = TRUE, fixed = TRUE)
  expect_length(ia_terms, nlevels(data$yprev) - 1)
  expect_equal(grep("time'", ia_terms, fixed = TRUE), integer(0))
})

test_that("vglm.markov() preserves rms %ia% semantics for spline bases", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2470, follow_up_time = 12)
  form <- stats::as.formula(
    "ordered(y) ~ rcs(time, 4) + tx + yprev + rcs(time, 4) %ia% yprev",
    env = baseenv()
  )

  fit <- suppressWarnings(
    markov.misc::vglm.markov(
      form,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  )

  ia_terms <- grep("%ia%", names(fit@constraints), value = TRUE, fixed = TRUE)
  expect_gt(length(ia_terms), nlevels(data$yprev) - 1)
  expect_gt(length(grep("time'", ia_terms, fixed = TRUE)), 0)
})

test_that("vglm.markov() accepts family names and formula-environment data", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3),
    w = c(1, 2, 1, 2, 1, 2, 1, 2)
  )
  env <- list2env(data, parent = baseenv())
  form <- stats::as.formula("y ~ x", env = env)

  fit_from_env <- suppressWarnings(
    vglm.markov(form, family = "binomialff")
  )
  weighted <- suppressWarnings(
    vglm.markov(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      weights = w,
      smart = FALSE
    )
  )

  expect_s4_class(fit_from_env, "vglm")
  expect_s4_class(weighted, "vglm")
  expect_false(isTRUE(weighted@smart.prediction$smart.arg))
  expect_equal(nrow(weighted@prior.weights), nrow(data))
  expect_equal(as.vector(weighted@prior.weights), data$w)
})

test_that("vglm.markov() rejects invalid families, weights, and form2 subsets", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1),
    x = c(0, 0, 1, 1),
    off = c(0.1, 0.2, 0.3, 0.4),
    bad_w = c(1, 1, -1, 1)
  )

  expect_error(
    vglm.markov(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      weights = bad_w
    ),
    "negative weights",
    fixed = TRUE
  )
  expect_error(
    vglm.markov(y ~ x, family = list(), data = data),
    "is not a VGAM family function",
    fixed = TRUE
  )
  expect_error(
    vglm.markov(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      form2 = ~x,
      subset = data$x == 1
    ),
    "subset",
    fixed = TRUE
  )
  expect_error(
    vglm.markov(
      y ~ x + offset(off),
      family = VGAM::binomialff,
      data = data
    ),
    "Model offsets are not supported",
    fixed = TRUE
  )
  expect_error(
    vglm.markov(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      offset = off
    ),
    "Model offsets are not supported",
    fixed = TRUE
  )
})

test_that("vglm.markov() records NA actions after dropping incomplete rows", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1, NA, 1),
    x = c(0, 0, 1, 1, 2, 2)
  )

  expect_error(
    suppressWarnings(
      vglm.markov(y ~ 0, family = VGAM::binomialff, data = data)
    ),
    "model matrix"
  )
  with_na <- suppressWarnings(
    vglm.markov(y ~ x, family = VGAM::binomialff, data = data)
  )

  expect_s4_class(with_na, "vglm")
  expect_length(with_na@na.action, 1)
})

test_that("vglm.markov() stores form2 design slots", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3),
    z = c(1, 2, 1, 2, 1, 2, 1, 2)
  )

  fit <- suppressWarnings(
    vglm.markov(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      form2 = y ~ z,
      x.arg = TRUE,
      y.arg = TRUE,
      qr.arg = TRUE
    )
  )

  expect_s4_class(fit, "vglm")
  expect_true(length(fit@Xm2) > 0)
  expect_true(length(fit@Ym2) > 0)
  expect_equal(deparse(fit@misc$form2), deparse(y ~ z))
  expect_true(length(fit@callXm2) > 0)
})

test_that("vglm.markov() validates form2 row alignment", {
  skip_if_not_installed("VGAM")

  data <- data.frame(
    y = c(0, 1, 0, 1, 0, 1, 0, 1),
    x = c(0, 0, 1, 1, 2, 2, 3, 3)
  )
  env <- list2env(
    list(z = 1:3, z_y = c(0, 1, 0)),
    parent = baseenv()
  )

  expect_error(
    vglm.markov(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      form2 = stats::as.formula("z_y ~ z", env = env)
    ),
    "number of rows of 'y' and 'Ym2' are unequal",
    fixed = TRUE
  )
  expect_error(
    vglm.markov(
      y ~ x,
      family = VGAM::binomialff,
      data = data,
      form2 = stats::as.formula("~ z", env = env)
    ),
    "number of rows of 'x' and 'Xm2' are unequal",
    fixed = TRUE
  )
})

test_that("add_rms_formula_helpers() returns non-formulas unchanged", {
  expect_equal(markov.misc:::add_rms_formula_helpers(1), 1)

  form <- y ~ x
  environment(form) <- NULL
  out <- markov.misc:::add_rms_formula_helpers(form)

  expect_s3_class(out, "formula")
  expect_true(exists("rcs", envir = environment(out), inherits = TRUE))
})
