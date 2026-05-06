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
