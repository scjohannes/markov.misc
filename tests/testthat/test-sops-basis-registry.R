test_that("RMS basis registry compiles rcs and lsp bases", {
  skip_if_not_installed("rms")
  skip_if_not_installed("Hmisc")

  handlers <- markov.misc:::rms_basis_registry()
  expect_equal(names(handlers), c("rcs", "lsp"))

  values <- c(-1, 0, 1, 3, 6, NA_real_)

  lsp <- markov.misc:::new_compiled_rms_basis(
    handler = "lsp",
    variable = "time",
    parameters = c(1, 4),
    columns = c("time", "time'", "time''")
  )
  actual_lsp <- markov.misc:::evaluate_compiled_rms_basis(lsp, values)
  expected_lsp <- cbind(values, pmax(values - 1, 0), pmax(values - 4, 0))
  expect_equal(as.vector(actual_lsp), as.vector(expected_lsp))
  expect_equal(attr(actual_lsp, "nonlinear"), c(FALSE, TRUE, TRUE))

  knots <- c(0, 1, 4, 5)
  rcs <- markov.misc:::new_compiled_rms_basis(
    handler = "rcs",
    variable = "time",
    parameters = knots,
    columns = c("time", "time'", "time''")
  )
  actual_rcs <- markov.misc:::evaluate_compiled_rms_basis(rcs, values)
  expected_rcs <- Hmisc::rcspline.eval(values, knots = knots, inclx = TRUE)
  expect_equal(
    as.vector(actual_rcs),
    as.vector(expected_rcs),
    tolerance = 1e-13
  )
})

test_that("RMS basis registry recognizes fitted Design metadata", {
  expect_equal(
    markov.misc:::rms_basis_handler_for_design("rcspline", 4L)$id,
    "rcs"
  )
  expect_equal(
    markov.misc:::rms_basis_handler_for_design("lspline", 3L)$id,
    "lsp"
  )
  expect_null(markov.misc:::rms_basis_handler_for_design("polynomial", 2L))
})

test_that("formula helper injection includes every registered RMS basis", {
  skip_if_not_installed("rms")

  form <- stats::as.formula(
    "y ~ rcs(time, 4) + lsp(age, c(30, 60))",
    env = baseenv()
  )
  out <- markov.misc:::add_rms_formula_helpers(form)
  out_env <- environment(out)

  expect_equal(exists("rcs", envir = out_env, inherits = TRUE), TRUE)
  expect_equal(exists("lsp", envir = out_env, inherits = TRUE), TRUE)
  expect_equal(exists("%ia%", envir = out_env, inherits = TRUE), TRUE)
})

test_that("vglm_markov splits lsp basis columns and records metadata", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 140, seed = 2471, follow_up_time = 12)
  form <- stats::as.formula(
    "ordered(y) ~ lsp(time, c(3, 7)) + tx + yprev",
    env = baseenv()
  )
  fit <- suppressWarnings(
    markov.misc::vglm_markov(
      form,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  )

  basis_terms <- attr(fit, "markov_basis_terms", exact = TRUE)
  expect_length(basis_terms, 1L)
  expect_equal(basis_terms[[1L]]$handlers, "lsp")
  expect_length(basis_terms[[1L]]$columns, 3L)
  expect_equal(attr(fit, "markov_split_assign", exact = TRUE), TRUE)

  lsp_constraints <- grep("lsp\\(", names(fit@constraints), value = TRUE)
  expect_length(lsp_constraints, 3L)

  baseline <- data[data$time == 1, , drop = FALSE][1:8, , drop = FALSE]
  fast <- markov.misc::soprob_markov(
    fit,
    newdata = baseline,
    times = 1:8,
    y_levels = 1:6,
    absorb = 6
  )
  reference <- markov.misc:::soprob_markov_reference(
    fit,
    newdata = baseline,
    times = 1:8,
    y_levels = 1:6,
    absorb = 6
  )
  expect_equal(fast, reference, tolerance = 1e-11)
})

test_that("orm lsp terms use the optimized execution plan", {
  skip_if_not_installed("rms")

  data <- make_test_data(n_patients = 100, seed = 2472, follow_up_time = 12)
  dd_name <- ".markov_misc_lsp_dd"
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
    ordered(y) ~ rms::lsp(time, c(3, 7)) + tx + yprev,
    data = data,
    x = TRUE,
    y = TRUE
  )
  baseline <- data[data$time == 1, , drop = FALSE][1:8, , drop = FALSE]

  plan <- markov.misc:::compile_sop_execution_plan(
    fit,
    newdata = baseline,
    times = 1:8,
    y_levels = 1:6,
    absorb = 6,
    builder = "batched"
  )
  fast <- markov.misc:::run_sop_execution_plan(
    plan,
    markov.misc:::get_effective_coefs(fit)
  )
  reference <- markov.misc:::soprob_markov_reference(
    fit,
    newdata = baseline,
    times = 1:8,
    y_levels = 1:6,
    absorb = 6
  )

  expect_s3_class(plan, "markov_sop_exec_plan")
  expect_equal(unname(fast), unname(reference), tolerance = 1e-11)
})
