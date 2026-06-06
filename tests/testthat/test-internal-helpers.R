test_that("relevel_factors_consecutive preserves types and updates absorbing states", {
  int_data <- data.frame(y = c(1L, 2L), yprev = c(1L, 2L))
  expect_warning(
    int_out <- relevel_factors_consecutive(
      int_data,
      original_data = data.frame(y = factor(1:3), yprev = factor(1:3)),
      ylevels = 1:3,
      absorb = c(2, 3)
    ),
    "were dropped",
    fixed = TRUE
  )
  expect_equal(int_out$absorb, 2)
  expect_equal(int_out$ylevels, 1:2)

  expect_warning(
    no_absorb <- relevel_factors_consecutive(
      data.frame(y = c(1, 2), yprev = c(1, 2)),
      original_data = data.frame(y = 1:3, yprev = 1:3),
      ylevels = 1:3,
      absorb = 3
    ),
    "Set to NULL",
    fixed = TRUE
  )
  expect_null(no_absorb$absorb)

  char_out <- relevel_factors_consecutive(
    data.frame(y = c("a", "b"), yprev = c("a", "b")),
    original_data = data.frame(y = c("a", "b", "c"), yprev = c("a", "b", "c")),
    ylevels = NULL,
    absorb = NULL
  )
  expect_equal(char_out$ylevels, NULL)
  expect_type(char_out$data$y, "character")

  inferred <- relevel_factors_consecutive(
    data.frame(y = c(1, 2), yprev = c(1, 2)),
    original_data = NULL,
    ylevels = NULL,
    absorb = NULL
  )
  expect_equal(inferred$missing_levels, character(0))
})

test_that("calc_time_in_state_diff handles factor inputs", {
  trajectories <- data.frame(
    id = rep(1:4, each = 2),
    y = factor(c("home", "home", "home", "icu", "icu", "home", "icu", "icu")),
    tx = factor(rep(c("control", "treated"), each = 4)),
    time = rep(1:2, times = 4)
  )

  out <- calc_time_in_state_diff(
    trajectories,
    target_state = factor("home"),
    reference_level = factor("control")
  )

  expect_equal(out$state, "home")
  expect_equal(out$reference_level, "control")
})

test_that("tidy_bootstrap_coefs can summarize with median estimates", {
  out <- tidy_bootstrap_coefs(
    data.frame(boot_id = 1:3, tx = c(3, 1, 2)),
    estimate = "median"
  )

  expect_equal(out$estimate, 2)
})

test_that("validate_coef_vcov reports structural mismatches", {
  beta <- c(a = 1, b = 2)

  expect_error(
    markov.misc:::validate_coef_vcov(beta, c(1, 2)),
    "must be a matrix",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::validate_coef_vcov(beta, diag(3)),
    "Dimension mismatch",
    fixed = TRUE
  )
  bad_rows <- diag(2)
  rownames(bad_rows) <- c("b", "a")
  colnames(bad_rows) <- names(beta)
  expect_error(
    markov.misc:::validate_coef_vcov(beta, bad_rows),
    "Coefficient names do not match row names",
    fixed = TRUE
  )
  bad_cols <- diag(2)
  rownames(bad_cols) <- names(beta)
  colnames(bad_cols) <- c("b", "a")
  expect_error(
    markov.misc:::validate_coef_vcov(beta, bad_cols),
    "Coefficient names do not match column names",
    fixed = TRUE
  )
})

test_that("set_coef() updates nested robcov_vglm fits and vgam aliases", {
  skip_if_not_installed("VGAM")

  data <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = c(0, 0, 1, 1, 2, 2))
  fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)

  expect_error(
    set_coef(fit, 1),
    "Length of new_coefs",
    fixed = TRUE
  )

  new_coefs <- stats::coef(fit) + 0.1
  robust <- structure(
    list(coefficients = stats::coef(fit), vglm_fit = fit),
    class = "robcov_vglm"
  )
  robust2 <- set_coef(robust, new_coefs)
  alias <- markov.misc:::set_coef.vgam(fit, new_coefs)

  expect_equal(robust2$coefficients, new_coefs)
  expect_equal(stats::coef(robust2$vglm_fit), new_coefs)
  expect_equal(stats::coef(alias), new_coefs)
})

test_that("get_effective_coefs() reports unsupported and mismatched coefficients", {
  skip_if_not_installed("VGAM")

  expect_error(
    get_effective_coefs(list()),
    "only supports vglm and orm",
    fixed = TRUE
  )

  data <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = c(0, 0, 1, 1, 2, 2))
  fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)
  expect_warning(
    gamma <- markov.misc:::get_effective_coefs_vglm(
      fit,
      beta = c(stats::coef(fit), extra = 1)
    ),
    "Mismatch in coefficient consumption",
    fixed = TRUE
  )
  expect_equal(ncol(gamma), length(VGAM::constraints(fit)))

  orm <- structure(list(coefficients = c("y>=2" = 0, x = 1), non.slopes = NULL), class = "orm")
  expect_error(
    markov.misc:::get_effective_coefs_orm(orm),
    "Cannot determine",
    fixed = TRUE
  )
  orm$non.slopes <- 2
  expect_error(
    markov.misc:::get_effective_coefs_orm(orm, beta = 1),
    "shorter than",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::get_effective_coefs_orm(orm, beta = c(1, 2, 3)),
    "Length of beta",
    fixed = TRUE
  )
})

test_that("plot validation branches are covered", {
  expect_error(
    plot_results(data.frame(x = 1, g = "a"), y = x, x = x, group = missing_group),
    "not found in data",
    fixed = TRUE
  )

  plot_data <- data.frame(time = 1:2, y = 1:2, tx = 0)
  p <- plot_sops(plot_data, geom = "line", facet_var = "tx")
  expect_s3_class(p, "ggplot")
})

test_that("relevel_factors_consecutive can relevel from ylevels alone", {
  out <- relevel_factors_consecutive(
    data = data.frame(y = c(1L, 3L), yprev = c(1L, 3L)),
    factor_cols = c("y", "yprev"),
    original_data = NULL,
    ylevels = 1:3,
    absorb = 3
  )

  expect_equal(out$missing_levels, "2")
  expect_type(out$data$y, "integer")
  expect_equal(out$data$y, c(1L, 2L))
  expect_equal(unname(out$absorb), 2)
})

test_that("previous-state coercion respects scalar storage modes", {
  expect_equal(
    markov.misc:::coerce_previous_state_values(c("1", "2"), 1L, "yprev"),
    c(1L, 2L)
  )
  expect_error(
    markov.misc:::coerce_previous_state_values("home", 1L, "yprev"),
    "integer-compatible",
    fixed = TRUE
  )

  expect_equal(
    markov.misc:::coerce_previous_state_values(c("1.5", "2"), 1, "yprev"),
    c(1.5, 2)
  )
  expect_error(
    markov.misc:::coerce_previous_state_values("home", 1, "yprev"),
    "numeric-compatible",
    fixed = TRUE
  )

  expect_equal(
    markov.misc:::coerce_previous_state_values(c(1, 2), "1", "yprev"),
    c("1", "2")
  )
  expect_equal(
    markov.misc:::coerce_previous_state_values(c(1, 2), list(), "yprev"),
    c("1", "2")
  )
})

test_that("orm prediction helpers validate metadata and categorical levels", {
  skip_if_not_installed("rms")

  expect_error(
    markov.misc:::orm_model_matrix(list(), data.frame(x = 1)),
    "requires an orm model",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::orm_model_matrix(
      structure(list(Design = NULL, sformula = NULL), class = "orm"),
      data.frame(x = 1)
    ),
    "stored Design metadata",
    fixed = TRUE
  )

  data <- data.frame(
    y = ordered(c(1, 2, 3, 1, 2, 3)),
    z = factor(c("a", "b", "a", "b", "a", "b"))
  )
  dd <- rms::datadist(data)
  withr::local_options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  withr::defer(rm("dd", envir = globalenv()))

  fit <- rms::orm(y ~ z, data = data, x = TRUE, y = TRUE)
  X <- markov.misc:::orm_model_matrix(
    fit,
    data.frame(y = ordered(c(1, 1), levels = 1:3), z = c("a", "b")),
    include_intercept = TRUE
  )

  expect_equal(nrow(X), 2)
  expect_true("(Intercept)" %in% colnames(X))
  expect_error(
    markov.misc:::orm_model_matrix(
      fit,
      data.frame(y = ordered(1, levels = 1:3), z = "c")
    ),
    "not among fitted orm levels",
    fixed = TRUE
  )

  bad_fit <- fit
  bad_fit$Design$colnames <- c(bad_fit$Design$colnames, "extra")
  expect_error(
    markov.misc:::orm_model_matrix(
      bad_fit,
      data.frame(y = ordered(1, levels = 1:3), z = "a")
    ),
    "Could not reconstruct",
    fixed = TRUE
  )

  expect_equal(
    markov.misc:::normalize_orm_prediction_data(
      fit,
      data.frame(y = ordered(1, levels = 1:3))
    ),
    data.frame(y = ordered(1, levels = 1:3))
  )

  strat_env <- new.env(parent = baseenv())
  strat_env$strat <- function(x) x
  strat_fit <- structure(
    list(
      Design = list(name = character(0), assume.code = integer(0), parms = list(), colnames = "x"),
      sformula = stats::as.formula("y ~ x + strat(z)", env = strat_env)
    ),
    class = "orm"
  )
  expect_error(
    markov.misc:::orm_model_matrix(
      strat_fit,
      data.frame(
        y = ordered(c(1, 1), levels = 1:3),
        x = 1:2,
        z = factor(c("a", "b"))
      )
    ),
    "Could not reconstruct",
    fixed = TRUE
  )
})

test_that("orm score-bootstrap helpers cover validation and alignment branches", {
  skip_if_not_installed("rms")

  expect_error(
    markov.misc:::compute_scores_orm(list()),
    "must be an orm object",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::compute_scores_orm(structure(list(x = NULL, y = NULL), class = "orm")),
    "x = TRUE, y = TRUE",
    fixed = TRUE
  )

  data <- data.frame(
    y = ordered(c(1, 2, 3, 1, 2, 3, 1, 2, 3)),
    x = c(-1, 0, 1, -0.5, 0.5, 1.5, -1.5, 0.2, 1.2)
  )
  dd <- rms::datadist(data)
  withr::local_options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  withr::defer(rm("dd", envir = globalenv()))

  fit <- rms::orm(y ~ x, data = data, x = TRUE, y = TRUE)
  scores <- markov.misc:::compute_scores_orm(fit)
  expect_equal(nrow(scores), nrow(data))
  expect_equal(colnames(scores), names(stats::coef(fit)))

  bad_fit <- fit
  bad_fit$y <- as.character(bad_fit$y)
  bad_fit$y[1] <- 99
  expect_error(
    markov.misc:::compute_scores_orm(bad_fit),
    "Could not match orm response",
    fixed = TRUE
  )

  expect_error(
    markov.misc:::align_cluster_orm(fit, ~id, nrow(scores)),
    "row-level vector",
    fixed = TRUE
  )
  expect_equal(markov.misc:::align_cluster_orm(fit, seq_len(nrow(scores)), nrow(scores)), seq_len(nrow(scores)))

  na_fit <- fit
  na_fit$na.action <- c(2L, 4L)
  expect_equal(
    markov.misc:::align_cluster_orm(na_fit, letters[1:(nrow(scores) + 2)], nrow(scores)),
    letters[-c(2, 4)][1:nrow(scores)]
  )
  expect_error(
    markov.misc:::align_cluster_orm(fit, 1:2, nrow(scores)),
    "does not match",
    fixed = TRUE
  )

  expect_true(is.matrix(markov.misc:::get_orm_model_vcov(fit)))
})

test_that("score bootstrap draw helper reports malformed robust components", {
  baseline <- data.frame(id = c("a", "b"))
  robust <- structure(
    list(
      coefficients = c(a = 0, b = 0),
      bread = diag(2),
      scores = matrix(0, nrow = 2, ncol = 2),
      cluster = c("a", "b")
    ),
    class = "robcov_vglm"
  )

  expect_error(
    markov.misc:::generate_score_bootstrap_draws(robust, baseline, "missing", 1),
    "ID variable",
    fixed = TRUE
  )

  bad_scores <- robust
  bad_scores$scores <- data.frame(a = 1)
  expect_error(
    markov.misc:::generate_score_bootstrap_draws(bad_scores, baseline, "id", 1),
    "must be a matrix",
    fixed = TRUE
  )

  bad_bread <- robust
  bad_bread$bread <- c(1, 2)
  expect_error(
    markov.misc:::generate_score_bootstrap_draws(bad_bread, baseline, "id", 1),
    "`model$bread` must be a matrix",
    fixed = TRUE
  )

  bad_cluster <- robust
  bad_cluster$cluster <- "a"
  expect_error(
    markov.misc:::generate_score_bootstrap_draws(bad_cluster, baseline, "id", 1),
    "does not match",
    fixed = TRUE
  )

  bad_coef <- robust
  bad_coef$coefficients <- c(a = 0)
  expect_error(
    markov.misc:::generate_score_bootstrap_draws(bad_coef, baseline, "id", 1),
    "Coefficient length",
    fixed = TRUE
  )

  expect_error(
    markov.misc:::generate_score_bootstrap_draws(robust, baseline[0, , drop = FALSE], "id", 1),
    "Invalid baseline weights",
    fixed = TRUE
  )
})

test_that("robcov_vglm and score extraction cover malformed internals", {
  skip_if_not_installed("VGAM")

  data <- data.frame(y = c(0, 1, 0, 1, 0, 1), x = c(0, 0, 1, 1, 2, 2))
  fit <- VGAM::vglm(y ~ x, family = VGAM::binomialff, data = data)
  n <- stats::nobs(fit, type = "lm")
  p <- length(stats::coef(fit))

  expect_error(
    compute_scores_vglm(lm(y ~ x, data = data)),
    "must be a vglm object",
    fixed = TRUE
  )

  scores <- compute_scores_vglm(fit)
  expect_equal(dim(scores), c(n, p))

  with_mocked_bindings(
    compute_scores_vglm = function(fit) matrix(0, nrow = n + 1, ncol = p),
    {
      expect_error(
        robcov_vglm(fit),
        "score rows",
        fixed = TRUE
      )
    }
  )
  with_mocked_bindings(
    compute_scores_vglm = function(fit) matrix(0, nrow = n, ncol = p + 1),
    {
      expect_error(
        robcov_vglm(fit),
        "score columns",
        fixed = TRUE
      )
    }
  )
})

test_that("print.robcov_vglm reports clustered fits", {
  obj <- structure(
    list(
      call = NULL,
      n = 4,
      n_clusters = 2,
      coefficients = c(a = 1),
      se = c(a = 0.5)
    ),
    class = "robcov_vglm"
  )

  expect_output(print(obj), "Number of clusters")
})

test_that("get_vcov_robust validates clusters before dispatching by model class", {
  fit <- lm(mpg ~ wt, data = mtcars)
  expect_error(
    markov.misc:::get_vcov_robust(fit, cluster = ~cyl, data = mtcars),
    "Unsupported model class",
    fixed = TRUE
  )

  data_na <- mtcars
  data_na$wt[1] <- NA
  fit_na <- lm(mpg ~ wt, data = data_na)
  expect_error(
    markov.misc:::get_vcov_robust(fit_na, cluster = ~cyl, data = data_na),
    "Unsupported model class",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::get_vcov_robust(
      fit,
      cluster = ~id,
      data = data.frame(id = seq_len(nrow(mtcars) + 1L))
    ),
    "Unsupported model class",
    fixed = TRUE
  )
  expanded <- rbind(mtcars, mtcars[1, , drop = FALSE])
  expanded$id <- seq_len(nrow(expanded))
  expect_error(
    markov.misc:::get_vcov_robust(fit, cluster = ~id, data = expanded),
    "Unsupported model class",
    fixed = TRUE
  )

  formula.weird_model <- function(x, ...) y ~ x
  nobs.weird_model <- function(object, ...) 2L
  assign("formula.weird_model", formula.weird_model, envir = globalenv())
  assign("nobs.weird_model", nobs.weird_model, envir = globalenv())
  withr::defer(rm("formula.weird_model", envir = globalenv()))
  withr::defer(rm("nobs.weird_model", envir = globalenv()))

  weird <- structure(list(), class = "weird_model")
  expect_error(
    markov.misc:::get_vcov_robust(
      weird,
      cluster = ~id,
      data = data.frame(id = 1:2)
    ),
    "Unsupported model class",
    fixed = TRUE
  )
  expect_error(
    markov.misc:::get_vcov_robust(
      weird,
      cluster = ~id,
      data = data.frame(id = 1:3, x = c(1, NA, 3))
    ),
    "Unsupported model class",
    fixed = TRUE
  )

  skip_if_not_installed("rms")
  data <- data.frame(
    y = ordered(rep(1:3, each = 4)),
    x = seq_len(12),
    id = rep(1:4, each = 3)
  )
  dd <- rms::datadist(data)
  withr::local_options(datadist = "dd")
  assign("dd", dd, envir = globalenv())
  withr::defer(rm("dd", envir = globalenv()))
  orm_fit <- rms::orm(y ~ x, data = data, x = TRUE, y = TRUE)
  expect_true(is.matrix(markov.misc:::get_vcov_robust(orm_fit, cluster = data$id)))
})

test_that("standardize_sops and inferences validate dispatch-only branches", {
  with_mocked_bindings(
    validate_markov_model = function(model) NULL,
    {
      expect_error(
        standardize_sops(structure(list(), class = "not_supported"), data = data.frame()),
        "model must be an orm",
        fixed = TRUE
      )
    }
  )

  object <- data.frame(estimate = 1)
  class(object) <- c("markov_avg_sops", class(object))
  expect_error(
    inferences(object, method = "bootstrap", engine = "score_bootstrap"),
    "only used when `method = \"simulation\"`",
    fixed = TRUE
  )
})

test_that("markov_msm_build validates model classes and design alignment", {
  with_mocked_bindings(
    get_effective_coefs = function(model, beta = NULL) {
      matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
    },
    {
      expect_error(
        markov.misc:::markov_msm_build(
          structure(list(), class = "not_supported"),
          data.frame(yprev = 1)
        ),
        "supports only vglm and orm",
        fixed = TRUE
      )
    }
  )
  expect_error(
    markov.misc:::markov_msm_build(
      structure(list(), class = "orm"),
      data.frame(yprev = 1),
      t_covs = data.frame(z = 1:2),
      times = 1
    ),
    "one row per requested time point",
    fixed = TRUE
  )

  calls <- 0
  with_mocked_bindings(
    get_effective_coefs = function(model, beta = NULL) {
      matrix(0, nrow = 1, ncol = 2, dimnames = list(NULL, c("(Intercept)", "x")))
    },
    orm_model_matrix = function(model, newdata, include_intercept = TRUE) {
      calls <<- calls + 1
      if (calls == 1) {
        matrix(1, nrow = nrow(newdata), ncol = 2, dimnames = list(NULL, c("(Intercept)", "x")))
      } else {
        matrix(1, nrow = nrow(newdata), ncol = 1, dimnames = list(NULL, "(Intercept)"))
      }
    },
    {
      expect_error(
        markov.misc:::markov_msm_build(
          structure(list(), class = "orm"),
          data.frame(id = 1, time = 2, yprev = 1),
          times = NULL,
          ylevels = 1:2
        ),
        "Design matrix columns missing",
        fixed = TRUE
      )
    }
  )
})

test_that("plot_results reports unrecognizable grouped expressions", {
  expect_error(
    plot_results(
      data.frame(x = 1, y = 1, g = factor("a")),
      y,
      x = x,
      group = c(g, x)
    ),
    "not recognizable",
    fixed = TRUE
  )
})
