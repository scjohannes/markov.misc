make_avg_sops_object <- function(model, newdata, by = NULL) {
  out <- expand.grid(
    time = 1:2,
    state = 1:2,
    tx = c(0, 1),
    KEEP.OUT.ATTRS = FALSE
  )
  if (!is.null(by)) {
    out[[by]] <- "a"
  }
  out$estimate <- 0.5
  class(out) <- c("markov_avg_sops", class(out))
  attr(out, "model") <- model
  attr(out, "newdata_orig") <- newdata
  attr(out, "call_args") <- list(times = 1:2)
  attr(out, "avg_args") <- list(
    variables = list(tx = c(0, 1)),
    by = by,
    times = 1:2,
    id_var = "id"
  )
  attr(out, "time_var") <- "time"
  attr(out, "p_var") <- "yprev"
  attr(out, "y_levels") <- 1:2
  attr(out, "absorb") <- NULL
  attr(out, "time_covariates") <- NULL
  out
}

make_individual_sops_object <- function(model, newdata, by = NULL) {
  out <- expand.grid(
    rowid = seq_len(nrow(newdata)),
    time = 1:2,
    state = 1:2,
    KEEP.OUT.ATTRS = FALSE
  )
  if (!is.null(by)) {
    out[[by]] <- newdata[[by]][out$rowid]
  }
  out$estimate <- 0.5
  class(out) <- c("markov_sops", class(out))
  attr(out, "model") <- model
  attr(out, "newdata_orig") <- newdata
  attr(out, "call_args") <- list(times = 1:2, by = by)
  attr(out, "time_var") <- "time"
  attr(out, "p_var") <- "yprev"
  attr(out, "y_levels") <- 1:2
  attr(out, "absorb") <- NULL
  attr(out, "time_covariates") <- NULL
  attr(out, "by") <- by
  out
}

test_that("inferences_simulation() validates stored models and simulation engines", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "mock_model")
  newdata <- data.frame(id = 1, tx = 0, yprev = 1, time = 1)
  object <- make_avg_sops_object(model, newdata)

  attr(object, "model") <- NULL
  expect_error(
    markov.misc:::inferences_simulation(
      object,
      engine = "mvn",
      score_weight_dist = "exponential",
      n_draws = 1,
      vcov = NULL,
      cluster = NULL,
      conf_level = 0.95,
      conf_type = "perc",
      workers = NULL,
      return_draws = FALSE
    ),
    "Model not stored",
    fixed = TRUE
  )

  attr(object, "model") <- model
  expect_error(
    markov.misc:::inferences_simulation(
      object,
      engine = "other",
      score_weight_dist = "exponential",
      n_draws = 1,
      vcov = NULL,
      cluster = NULL,
      conf_level = 0.95,
      conf_type = "perc",
      workers = NULL,
      return_draws = FALSE
    ),
    "Unknown simulation engine",
    fixed = TRUE
  )
})

test_that("inferences() validates the homogenized method interface", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  newdata <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    yprev = 1,
    tx = 0
  )
  object <- make_avg_sops_object(model, newdata)

  expect_error(inferences(object, method = "other", n_draws = 1), "arg")
  expect_false("engine" %in% names(formals(inferences)))
  expect_false("score_weight_dist" %in% names(formals(inferences)))
  expect_false("..." %in% names(formals(inferences)))

  seen_engine <- NULL
  with_mocked_bindings(
    inferences_bootstrap = function(object, engine, ...) {
      seen_engine <<- engine
      attr(object, "engine") <- engine
      object
    },
    {
      out <- inferences(
        object,
        method = "fwb",
        n_draws = 1
      )
    }
  )

  expect_equal(seen_engine, "fwb")
  expect_equal(attr(out, "method"), "fwb")
})

test_that("inferences() rejects obsolete and misspelled arguments", {
  object <- structure(
    data.frame(estimate = 0.5),
    class = c("markov_sops", "data.frame")
  )

  expect_error(
    inferences(object, n_boot = 2),
    "unused argument (n_boot = 2)",
    fixed = TRUE
  )
})

test_that("inferences() reports repeated R fallbacks only once", {
  object <- structure(
    data.frame(estimate = 0.5),
    class = c("markov_sops", "data.frame")
  )

  with_mocked_bindings(
    inferences_impl = function(...) {
      for (i in seq_len(1000L)) {
        notify_sop_reference_fallback(paste("draw", i))
      }
      object
    },
    {
      messages <- testthat::capture_messages(inferences(object))
    }
  )

  expect_length(messages, 1L)
  expect_match(messages, "Compiled C\\+\\+ SOP calculations were not used")
})

test_that("inferences_simulation() falls back when fast-path setup fails", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "vglm")
  newdata <- data.frame(id = 1:2, tx = c(0, 1), yprev = c(1, 2), time = 1)
  object <- make_avg_sops_object(model, newdata)

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2),
    validate_coef_vcov = function(beta, Sigma, arg = "vcov") Sigma,
    compile_sop_execution_plan = function(...) stop("build failed"),
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) {
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5,
          0.7,
          0.6,
          0.3,
          0.4,
          0.5,
          0.4,
          0.5,
          0.6
        ),
        dim = c(4, 2, 2)
      )
    },
    {
      withr::local_seed(1)
      expect_message(
        out <- markov.misc:::inferences_simulation(
          object,
          engine = "mvn",
          score_weight_dist = "exponential",
          n_draws = 2,
          vcov = NULL,
          cluster = NULL,
          conf_level = 0.8,
          conf_type = "wald",
          workers = NULL,
          return_draws = TRUE
        ),
        "Compiled C\\+\\+ SOP calculations were not used"
      )
    }
  )

  expect_s3_class(out, "markov_avg_sops")
  expect_equal(attr(out, "n_successful"), 2L)
  expect_equal(attr(out, "conf_type"), "wald")
  expect_s3_class(attr(out, "draws"), "data.frame")
})

test_that("inferences_simulation() errors when all fast-path draws fail", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "vglm")
  newdata <- data.frame(id = 1:2, tx = c(0, 1), yprev = c(1, 2), time = 1)
  object <- make_avg_sops_object(model, newdata)

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2),
    validate_coef_vcov = function(beta, Sigma, arg = "vcov") Sigma,
    compile_sop_execution_plan = function(...) {
      structure(list(ok = TRUE), class = "markov_sop_exec_plan")
    },
    get_effective_coefs = function(model, beta = NULL) {
      matrix(c(0, 0), nrow = 1)
    },
    run_sop_execution_plan = function(...) stop("run failed"),
    {
      withr::local_seed(2)
      expect_warning(
        expect_error(
          markov.misc:::inferences_simulation(
            object,
            engine = "mvn",
            score_weight_dist = "exponential",
            n_draws = 1,
            vcov = NULL,
            cluster = NULL,
            conf_level = 0.95,
            conf_type = "perc",
            workers = NULL,
            return_draws = FALSE
          ),
          "All simulation draws failed",
          fixed = TRUE
        ),
        "markov_msm_run failed"
      )
    }
  )
})

test_that("inferences_simulation() preserves individual slow-path grouping", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "mock_model")
  newdata <- data.frame(
    id = 1:2,
    tx = c(0, 1),
    yprev = c(1, 2),
    time = 1,
    group = c("a", "b")
  )
  object <- make_individual_sops_object(model, newdata, by = "group")

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2),
    validate_coef_vcov = function(beta, Sigma, arg = "vcov") Sigma,
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) {
      array(
        c(0.8, 0.2, 0.6, 0.4, 0.7, 0.3, 0.5, 0.5),
        dim = c(2, 2, 2)
      )
    },
    {
      withr::local_seed(3)
      out <- markov.misc:::inferences_simulation(
        object,
        engine = "mvn",
        score_weight_dist = "exponential",
        n_draws = 2,
        vcov = diag(2),
        cluster = NULL,
        conf_level = 0.95,
        conf_type = "perc",
        workers = NULL,
        return_draws = TRUE
      )
    }
  )

  expect_s3_class(out, "markov_sops")
  expect_equal(attr(out, "n_successful"), 2L)
  expect_true("group" %in% names(out))
  expect_s3_class(attr(out, "draws"), "data.frame")
})

test_that("inferences_simulation() rejects invalid custom vcov objects", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "mock_model")
  newdata <- data.frame(
    id = 1:2,
    tx = c(0, 1),
    yprev = c(1, 2),
    time = 1
  )
  object <- make_avg_sops_object(model, newdata)

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    {
      expect_error(
        markov.misc:::inferences_simulation(
          object,
          engine = "mvn",
          score_weight_dist = "exponential",
          n_draws = 2,
          vcov = data.frame(a = c(1, 0), b = c(0, 1)),
          cluster = NULL,
          conf_level = 0.95,
          conf_type = "perc",
          workers = NULL,
          return_draws = FALSE
        ),
        "`vcov` must be a matrix.",
        fixed = TRUE
      )
    }
  )
})

test_that("inferences_simulation() applies score-bootstrap slow-path weights", {
  model <- structure(
    list(
      coefficients = c(a = 0, b = 0),
      bread = diag(2),
      scores = matrix(0, nrow = 2, ncol = 2),
      cluster = c(1, 2)
    ),
    class = "robcov_vglm"
  )
  newdata <- data.frame(id = 1:2, tx = c(0, 1), yprev = c(1, 2), time = 1)
  object <- make_avg_sops_object(model, newdata)

  with_mocked_bindings(
    generate_score_bootstrap_draws = function(...) {
      list(
        beta_draws = matrix(c(0, 0, 0, 0), nrow = 2, byrow = TRUE),
        baseline_weights = matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
      )
    },
    compile_sop_execution_plan = function(...) stop("build failed"),
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) {
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5,
          0.7,
          0.6,
          0.3,
          0.4,
          0.5,
          0.4,
          0.5,
          0.6
        ),
        dim = c(4, 2, 2)
      )
    },
    {
      out <- markov.misc:::inferences_simulation(
        object,
        engine = "score_bootstrap",
        score_weight_dist = "exponential",
        n_draws = 2,
        vcov = NULL,
        cluster = NULL,
        conf_level = 0.8,
        conf_type = "perc",
        workers = NULL,
        return_draws = FALSE
      )
    }
  )

  expect_s3_class(out, "markov_avg_sops")
  expect_equal(attr(out, "n_successful"), 2L)
})

test_that("inferences_simulation() applies score-bootstrap by weights", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  newdata <- data.frame(
    id = 1:4,
    time = 1,
    yprev = 1,
    tx = 0,
    grp = c("a", "a", "b", "b")
  )
  object <- expand.grid(
    time = 1,
    state = 1,
    tx = 0,
    grp = c("a", "b"),
    KEEP.OUT.ATTRS = FALSE
  )
  object$estimate <- 0.5
  class(object) <- c("markov_avg_sops", class(object))
  attr(object, "model") <- model
  attr(object, "newdata_orig") <- newdata
  attr(object, "call_args") <- list(times = 1)
  attr(object, "avg_args") <- list(
    variables = list(tx = 0),
    by = "grp",
    times = 1,
    id_var = "id"
  )
  attr(object, "time_var") <- "time"
  attr(object, "p_var") <- "yprev"
  attr(object, "y_levels") <- 1
  attr(object, "absorb") <- NULL
  attr(object, "time_covariates") <- NULL

  with_mocked_bindings(
    generate_score_bootstrap_draws = function(...) {
      list(
        beta_draws = matrix(0, nrow = 1, ncol = 1),
        baseline_weights = matrix(c(1, 3, 2, 6), nrow = 1)
      )
    },
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) {
      array(c(0.2, 0.8, 0.1, 0.5), dim = c(4, 1, 1))
    },
    {
      out <- markov.misc:::inferences_simulation(
        object,
        engine = "score_bootstrap",
        score_weight_dist = "exponential",
        n_draws = 1,
        vcov = NULL,
        cluster = NULL,
        conf_level = 0.8,
        conf_type = "perc",
        workers = NULL,
        return_draws = TRUE
      )
    }
  )

  draws <- attr(out, "draws")
  draws <- draws[order(draws$grp), , drop = FALSE]

  expect_equal(draws$estimate, c(0.65, 0.4))
  expect_false("score_weight" %in% names(draws))
  expect_equal(attr(out, "draw_weight_omission_reason"), "averaged_sops")
})

test_that("inferences_simulation() supports score-bootstrap individual sops", {
  model <- structure(
    list(coefficients = c(a = 0, b = 0)),
    class = "robcov_vglm"
  )
  newdata <- data.frame(id = 1:2, yprev = c(1, 2), time = 1)
  object <- make_individual_sops_object(model, newdata)
  attr(object, "id_var") <- "id"
  attr(object, "newdata_supplied") <- FALSE

  with_mocked_bindings(
    generate_score_bootstrap_draws = function(...) {
      list(
        beta_draws = matrix(c(0, 0, 0, 0), nrow = 2, byrow = TRUE),
        baseline_weights = matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
      )
    },
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) {
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5
        ),
        dim = c(2, 2, 2)
      )
    },
    {
      out <- markov.misc:::inferences_simulation(
        object,
        engine = "score_bootstrap",
        score_weight_dist = "exponential",
        n_draws = 2,
        vcov = NULL,
        cluster = NULL,
        conf_level = 0.8,
        conf_type = "perc",
        workers = NULL,
        return_draws = TRUE
      )
    }
  )

  draws <- attr(out, "draws")
  public_draws <- markov.misc::get_draws(out)

  expect_s3_class(out, "markov_sops")
  expect_equal(attr(out, "method"), "score_bootstrap")
  expect_equal(attr(out, "draw_weights_attached"), TRUE)
  expect_equal(attr(out, "draw_weight_col"), "score_weight")
  expect_equal(sort(unique(draws$score_weight)), c(0.3, 0.4, 0.6, 0.7))
  expect_equal(
    sort(unique(public_draws$score_weight)),
    c(0.3, 0.4, 0.6, 0.7)
  )
  expect_equal(public_draws$estimate, rep(0.5, nrow(public_draws)))
})

test_that("inferences_simulation() rejects reserved score weight columns", {
  model <- structure(
    list(coefficients = c(a = 0, b = 0)),
    class = "robcov_vglm"
  )
  newdata <- data.frame(
    id = 1:2,
    yprev = c(1, 2),
    time = 1,
    score_weight = c(10, 20)
  )
  object <- make_individual_sops_object(model, newdata)
  attr(object, "id_var") <- "id"
  attr(object, "newdata_supplied") <- FALSE

  with_mocked_bindings(
    generate_score_bootstrap_draws = function(...) {
      list(
        beta_draws = matrix(c(0, 0), nrow = 1),
        baseline_weights = matrix(c(0.7, 0.3), nrow = 1)
      )
    },
    {
      expect_error(
        markov.misc:::inferences_simulation(
          object,
          engine = "score_bootstrap",
          score_weight_dist = "exponential",
          n_draws = 1,
          vcov = NULL,
          cluster = NULL,
          conf_level = 0.8,
          conf_type = "perc",
          workers = NULL,
          return_draws = TRUE
        ),
        "reserved for bootstrap draw weights",
        fixed = TRUE
      )
    }
  )
})

test_that("inferences_simulation() reports slow-path SOP failures", {
  model <- structure(list(coefficients = c(a = 0, b = 0)), class = "mock_model")
  newdata <- data.frame(id = 1:2, tx = c(0, 1), yprev = c(1, 2), time = 1)
  object <- make_avg_sops_object(model, newdata)

  with_mocked_bindings(
    get_coef = function(model) c(a = 0, b = 0),
    get_vcov_robust = function(model) diag(2),
    validate_coef_vcov = function(beta, Sigma, arg = "vcov") Sigma,
    set_coef = function(model, new_coefs) model,
    soprob_markov = function(...) stop("sop failed"),
    {
      withr::local_seed(4)
      expect_warning(
        expect_error(
          markov.misc:::inferences_simulation(
            object,
            engine = "mvn",
            score_weight_dist = "exponential",
            n_draws = 1,
            vcov = diag(2),
            cluster = NULL,
            conf_level = 0.95,
            conf_type = "perc",
            workers = NULL,
            return_draws = FALSE
          ),
          "All simulation draws failed",
          fixed = TRUE
        ),
        "soprob_markov failed"
      )
    }
  )
})

test_that("soprob_markov() validates model classes and propagates Bayesian draws", {
  expect_error(
    soprob_markov(
      model = structure(list(), class = "other"),
      newdata = data.frame(yprev = 1),
      times = 1,
      y_levels = 1
    ),
    "Object class not supported",
    fixed = TRUE
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    {
      expect_error(
        soprob_markov(
          model = structure(list(vglm_fit = NULL), class = "robcov_vglm"),
          newdata = data.frame(yprev = 1),
          times = 1,
          y_levels = 1:2
        ),
        "does not contain the original vglm fit",
        fixed = TRUE
      )
      expect_error(
        soprob_markov(
          model = structure(list(), class = "vglm"),
          newdata = data.frame(x = 1),
          times = 1,
          y_levels = 1:2
        ),
        "Previous-state variable",
        fixed = TRUE
      )
    }
  )

  draws <- matrix(stats::qlogis(0.4), nrow = 2, ncol = 1)
  colnames(draws) <- "y>=2"
  blrm <- structure(
    list(
      draws = draws,
      non.slopes = 1L,
      pppo = 0L,
      tauInfo = data.frame(name = character()),
      ylevels = 1:2,
      yname = "y"
    ),
    class = "blrm"
  )

  with_mocked_bindings(
    blrm_design_matrix = function(model, newdata, second = FALSE) {
      matrix(numeric(0), nrow = nrow(newdata), ncol = 0)
    },
    {
      out <- soprob_markov(
        model = blrm,
        newdata = data.frame(
          id = 1:2,
          yprev = factor(c(1, 1)),
          gap_var = 0,
          z = 0
        ),
        times = 1:2,
        y_levels = 1:2,
        absorb = 2,
        gap_var = "gap_var",
        time_covariates = data.frame(z = c(1, 2))
      )
    }
  )

  expect_equal(dim(out), c(2L, 2L, 2L, 2L))
  expect_equal(unname(out[,, 2, 2]), array(0.64, dim = c(2, 2)))
})

test_that("inferences_bootstrap() validates inputs and records callback outcomes", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  newdata <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    y = factor(c(1, 2, 1, 2), levels = 1:3),
    yprev = factor(c(1, 1, 1, 1), levels = 1:3),
    tx = c(0, 0, 1, 1),
    grp = "a"
  )
  object <- make_avg_sops_object(model, newdata, by = "grp")
  attr(object, "avg_args")$variables <- list(tx = c(0, 1), grp = "a")
  object <- expand.grid(
    time = 1:2,
    state = 1:3,
    tx = c(0, 1),
    grp = "a",
    KEEP.OUT.ATTRS = FALSE
  )
  object$estimate <- 0.5
  class(object) <- c("markov_avg_sops", class(object))
  attr(object, "model") <- model
  attr(object, "newdata_orig") <- newdata
  attr(object, "avg_args") <- list(
    variables = list(tx = c(0, 1), grp = "a"),
    by = "grp",
    times = 1:2,
    id_var = "id"
  )
  attr(object, "call_args") <- list(times = 1:2)
  attr(object, "time_var") <- "time"
  attr(object, "p_var") <- "yprev"
  attr(object, "y_levels") <- 1:3
  attr(object, "absorb") <- 3
  attr(object, "time_covariates") <- NULL

  expect_error(
    markov.misc:::inferences_bootstrap(
      structure(data.frame(), class = "markov_sops"),
      engine = "standard",
      n_boot = 1,
      workers = NULL,
      conf_level = 0.95,
      return_draws = FALSE,
      update_datadist = TRUE,
      use_coefstart = FALSE
    ),
    "Standard refit bootstrap is not supported",
    fixed = TRUE
  )

  no_data <- object
  attr(no_data, "newdata_orig") <- NULL
  expect_error(
    markov.misc:::inferences_bootstrap(
      no_data,
      engine = "standard",
      n_boot = 1,
      workers = NULL,
      conf_level = 0.95,
      return_draws = FALSE,
      update_datadist = TRUE,
      use_coefstart = FALSE
    ),
    "Full refit data not stored",
    fixed = TRUE
  )

  missing_id <- object
  attr(missing_id, "newdata_orig") <- data.frame(time = 1:2)
  expect_error(
    markov.misc:::inferences_bootstrap(
      missing_id,
      engine = "standard",
      n_boot = 1,
      workers = NULL,
      conf_level = 0.95,
      return_draws = FALSE,
      update_datadist = TRUE,
      use_coefstart = FALSE
    ),
    "ID variable 'id' not found",
    fixed = TRUE
  )

  baseline_only <- object
  attr(baseline_only, "newdata_orig") <- newdata[newdata$time == 1, ]
  expect_error(
    markov.misc:::inferences_bootstrap(
      baseline_only,
      engine = "standard",
      n_boot = 1,
      workers = NULL,
      conf_level = 0.95,
      return_draws = FALSE,
      update_datadist = TRUE,
      use_coefstart = FALSE
    ),
    "requires full longitudinal data",
    fixed = TRUE
  )

  wrapper_call <- 0
  sop_call <- 0
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
      wrapper_call <<- wrapper_call + 1
      if (wrapper_call == 2) {
        return(list(
          model = NULL,
          data = transform(newdata, new_id = id),
          y_levels = 1:3,
          absorb = 3,
          missing_states = character(0)
        ))
      }
      list(
        model = model,
        data = transform(newdata, new_id = id),
        y_levels = 1:2,
        absorb = 2,
        missing_states = 3
      )
    },
    soprob_markov = function(...) {
      sop_call <<- sop_call + 1
      if (sop_call == 2) {
        stop("sop failed")
      }
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5,
          0.7,
          0.6,
          0.3,
          0.4,
          0.5,
          0.4,
          0.5,
          0.6
        ),
        dim = c(4, 2, 2)
      )
    },
    {
      expect_warning(
        out <- markov.misc:::inferences_bootstrap(
          object,
          engine = "standard",
          n_boot = 3,
          workers = NULL,
          conf_level = 0.8,
          return_draws = TRUE,
          update_datadist = FALSE,
          use_coefstart = TRUE
        ),
        "soprob_markov failed"
      )
    }
  )

  expect_s3_class(out, "markov_avg_sops")
  expect_equal(attr(out, "n_successful"), 1L)
  expect_s3_class(attr(out, "draws"), "data.frame")
  expect_equal(as.character(out$grp), rep("a", nrow(out)))

  with_mocked_bindings(
    fast_group_bootstrap = function(...) {
      list(data.frame(original_id = 1, new_id = "1_1"))
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
      list(NULL)
    },
    {
      expect_error(
        markov.misc:::inferences_bootstrap(
          object,
          engine = "standard",
          n_boot = 1,
          workers = NULL,
          conf_level = 0.95,
          return_draws = FALSE,
          update_datadist = FALSE,
          use_coefstart = FALSE
        ),
        "All bootstrap iterations failed",
        fixed = TRUE
      )
    }
  )
})

test_that("inferences_bootstrap() supports fractional weighted refits", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  newdata <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    y = factor(c(1, 2, 1, 2), levels = 1:2),
    yprev = factor(c(1, 1, 1, 1), levels = 1:2),
    tx = c(0, 0, 1, 1)
  )
  object <- make_avg_sops_object(model, newdata)

  fit_weights_seen <- NULL
  with_mocked_bindings(
    generate_fwb_bootstrap_weights = function(data, id_var, n_boot) {
      list(data.frame(
        original_id = c(1, 2),
        fwb_weight = c(0.5, 1.5),
        boot_id = 1
      ))
    },
    bootstrap_analysis_wrapper = function(
      boot_data,
      model,
      factor_cols,
      original_data,
      y_levels,
      absorb,
      update_datadist,
      use_coefstart,
      fit_weights
    ) {
      fit_weights_seen <<- fit_weights
      list(
        model = model,
        data = boot_data,
        y_levels = 1:2,
        absorb = NULL,
        missing_states = character(0)
      )
    },
    soprob_markov = function(...) {
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5,
          0.7,
          0.6,
          0.3,
          0.4,
          0.5,
          0.4,
          0.5,
          0.6
        ),
        dim = c(4, 2, 2)
      )
    },
    {
      out <- markov.misc:::inferences_bootstrap(
        object,
        engine = "fwb",
        n_boot = 1,
        workers = NULL,
        conf_level = 0.95,
        return_draws = TRUE,
        update_datadist = FALSE,
        use_coefstart = FALSE
      )
    }
  )

  expect_equal(fit_weights_seen, c(0.5, 0.5, 1.5, 1.5))
  expect_equal(attr(out, "method"), "fwb")
  expect_equal(attr(out, "fwb_weight_type"), "exponential")
  expect_equal(attr(out, "fwb_weight_scale"), "cluster_mean_1")
  expect_s3_class(attr(out, "draws"), "data.frame")
})

test_that("inferences_bootstrap() applies FWB by weights", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  refit_data <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, times = 4),
    y = factor(rep(c(1, 2), times = 4), levels = 1:2),
    yprev = factor(1, levels = 1:2),
    tx = 0,
    grp = rep(c("a", "a", "b", "b"), each = 2)
  )
  object <- expand.grid(
    time = 1,
    state = 1,
    tx = 0,
    grp = c("a", "b"),
    KEEP.OUT.ATTRS = FALSE
  )
  object$estimate <- 0.5
  class(object) <- c("markov_avg_sops", class(object))
  attr(object, "model") <- model
  attr(object, "newdata_orig") <- refit_data
  attr(object, "avg_args") <- list(
    variables = list(tx = 0),
    by = "grp",
    times = 1,
    id_var = "id"
  )
  attr(object, "call_args") <- list(times = 1)
  attr(object, "time_var") <- "time"
  attr(object, "p_var") <- "yprev"
  attr(object, "y_levels") <- 1
  attr(object, "absorb") <- NULL
  attr(object, "time_covariates") <- NULL

  with_mocked_bindings(
    generate_fwb_bootstrap_weights = function(data, id_var, n_boot) {
      list(data.frame(
        original_id = 1:4,
        fwb_weight = c(1, 3, 2, 6),
        boot_id = 1
      ))
    },
    apply_to_fwb_bootstrap = function(
      fwb_samples,
      analysis_fn,
      data,
      id_var,
      workers,
      packages,
      globals
    ) {
      list(analysis_fn(data, fwb_samples[[1]]))
    },
    bootstrap_analysis_wrapper = function(
      boot_data,
      model,
      factor_cols,
      original_data,
      y_levels,
      absorb,
      update_datadist,
      use_coefstart,
      fit_weights
    ) {
      list(
        model = model,
        data = boot_data,
        y_levels = 1,
        absorb = NULL,
        missing_states = character(0)
      )
    },
    soprob_markov = function(...) {
      array(c(0.2, 0.8, 0.1, 0.5), dim = c(4, 1, 1))
    },
    {
      out <- markov.misc:::inferences_bootstrap(
        object,
        engine = "fwb",
        n_boot = 1,
        workers = NULL,
        conf_level = 0.95,
        return_draws = TRUE,
        update_datadist = FALSE,
        use_coefstart = FALSE
      )
    }
  )

  draws <- attr(out, "draws")
  draws <- draws[order(draws$grp), , drop = FALSE]

  expect_equal(draws$estimate, c(0.65, 0.4))
  expect_false("fwb_weight" %in% names(draws))
})

test_that("inferences_bootstrap_sops_fwb() attaches stored-data draw weights", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  prediction_data <- data.frame(id = 1:2, time = 1, yprev = 1)
  refit_data <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    y = factor(c(1, 2, 1, 2), levels = 1:2),
    yprev = factor(c(1, 1, 1, 1), levels = 1:2)
  )
  object <- make_individual_sops_object(model, prediction_data)
  attr(object, "newdata_pred") <- prediction_data
  attr(object, "refit_data") <- refit_data
  attr(object, "id_var") <- "id"
  attr(object, "newdata_supplied") <- FALSE

  with_mocked_bindings(
    generate_fwb_bootstrap_weights = function(data, id_var, n_boot) {
      list(data.frame(
        original_id = c(1, 2),
        fwb_weight = c(0.5, 1.5),
        boot_id = 1
      ))
    },
    bootstrap_analysis_wrapper = function(
      boot_data,
      model,
      factor_cols,
      original_data,
      y_levels,
      absorb,
      update_datadist,
      use_coefstart,
      fit_weights
    ) {
      list(
        model = model,
        data = boot_data,
        y_levels = 1:2,
        absorb = NULL,
        missing_states = character(0)
      )
    },
    soprob_markov = function(...) {
      array(
        c(
          0.8,
          0.7,
          0.2,
          0.3,
          0.6,
          0.5,
          0.4,
          0.5
        ),
        dim = c(2, 2, 2)
      )
    },
    {
      out <- markov.misc:::inferences_bootstrap_sops_fwb(
        object,
        n_boot = 1,
        workers = NULL,
        conf_level = 0.95,
        return_draws = TRUE,
        update_datadist = FALSE,
        use_coefstart = FALSE
      )
    }
  )

  draws <- attr(out, "draws")
  public_draws <- markov.misc::get_draws(out)

  expect_equal(attr(out, "draw_weights_attached"), TRUE)
  expect_equal(attr(out, "draw_weight_col"), "fwb_weight")
  expect_equal(draws$fwb_weight[1:2], c(0.5, 1.5))
  expect_equal(sort(unique(public_draws$fwb_weight)), c(0.5, 1.5))
  expect_equal(public_draws$estimate, rep(0.5, nrow(public_draws)))
})

test_that("inferences_bootstrap_sops_fwb() rejects reserved draw weight columns", {
  model <- structure(list(coefficients = c(a = 0)), class = "mock_model")
  prediction_data <- data.frame(
    id = 1:2,
    time = 1,
    yprev = 1,
    fwb_weight = c(10, 20)
  )
  refit_data <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(1:2, times = 2),
    y = factor(c(1, 2, 1, 2), levels = 1:2),
    yprev = factor(c(1, 1, 1, 1), levels = 1:2)
  )
  object <- make_individual_sops_object(model, prediction_data)
  attr(object, "newdata_pred") <- prediction_data
  attr(object, "refit_data") <- refit_data
  attr(object, "id_var") <- "id"
  attr(object, "newdata_supplied") <- FALSE

  expect_error(
    markov.misc:::inferences_bootstrap_sops_fwb(
      object,
      n_boot = 1,
      workers = NULL,
      conf_level = 0.95,
      return_draws = TRUE,
      update_datadist = FALSE,
      use_coefstart = FALSE
    ),
    "reserved for bootstrap draw weights",
    fixed = TRUE
  )
})
