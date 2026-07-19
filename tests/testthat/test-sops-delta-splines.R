local_spline_delta_datadist <- function(data, env = parent.frame()) {
  name <- "spline_delta_test_dd"
  old_exists <- exists(name, envir = globalenv(), inherits = FALSE)
  old_value <- if (old_exists) {
    get(name, envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
  assign(name, rms::datadist(data), envir = globalenv())
  old_options <- options(datadist = name)
  withr::defer(
    {
      options(old_options)
      if (old_exists) {
        assign(name, old_value, envir = globalenv())
      } else if (exists(name, envir = globalenv(), inherits = FALSE)) {
        rm(list = name, envir = globalenv())
      }
    },
    envir = env
  )
  invisible(NULL)
}

local_spline_delta_orm_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(value)) {
      state_levels <- as.character(seq_len(8L))
      raw_data <- sim_trajectories_brownian(
        n_patients = 48,
        follow_up_time = 6,
        n_states = length(state_levels),
        thresholds = stats::qnorm(seq(0.05, 0.95, length.out = 7L)),
        allowed_start_state = as.integer(2:7),
        absorbing_state = NULL,
        drift_change_times = NULL,
        latent_dist = "normal",
        seed = 2701
      )
      data <- prepare_markov_data(
        raw_data,
        absorbing_state = NULL,
        factor_previous = FALSE
      )
      data$y <- ordered(data$y, levels = state_levels)

      local_spline_delta_datadist(data)
      model <- orm_markov(
        y ~ rms::rcs(time, 3) + tx + rms::rcs(yprev, 4),
        data = data,
        id_var = "id",
        opt_method = "LM",
        scale = TRUE,
        penalty = 0.1
      )
      value <<- list(data = data, model = model)
    }
    value
  }
})

replay_spline_avg_sops <- function(object, model) {
  avg_args <- attr(object, "avg_args")
  avg_sops(
    model = model,
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    y_levels = attr(object, "y_levels"),
    absorb = attr(object, "absorb"),
    id_var = avg_args$id_var,
    time_var = attr(object, "time_var") %||% "time",
    p_var = attr(object, "p_var") %||% "yprev",
    p2_var = attr(object, "p2_var"),
    gap_var = attr(object, "gap_var"),
    time_covariates = attr(object, "time_covariates")
  )
}

central_spline_avg_jacobian <- function(object, model) {
  coefficient <- get_coef(model)
  avg_args <- attr(object, "avg_args")
  key_columns <- c("time", "state", names(avg_args$variables), avg_args$by)
  object_key <- sop_draw_cell_key(object, key_columns)
  jacobian <- matrix(
    NA_real_,
    nrow = nrow(object),
    ncol = length(coefficient),
    dimnames = list(NULL, names(coefficient))
  )

  for (column in seq_along(coefficient)) {
    step <- .Machine$double.eps^(1 / 3) *
      max(1, abs(coefficient[column]))
    upper <- coefficient
    lower <- coefficient
    upper[column] <- upper[column] + step
    lower[column] <- lower[column] - step
    upper_result <- replay_spline_avg_sops(object, set_coef(model, upper))
    lower_result <- replay_spline_avg_sops(object, set_coef(model, lower))
    upper_index <- match(
      object_key,
      sop_draw_cell_key(upper_result, key_columns)
    )
    lower_index <- match(
      object_key,
      sop_draw_cell_key(lower_result, key_columns)
    )
    if (anyNA(upper_index) || anyNA(lower_index)) {
      stop("Finite-difference spline SOP cells did not align.")
    }
    jacobian[, column] <-
      (upper_result$estimate[upper_index] -
        lower_result$estimate[lower_index]) /
      (2 * step)
  }

  jacobian
}

test_that("backend orm spline covariance labels preserve coefficient order", {
  fit <- local_spline_delta_orm_case()$model
  coefficient_names <- names(stats::coef(fit))
  q <- length(coefficient_names)
  expected <- diag(seq_len(q))

  unnamed_fit <- fit
  unnamed_fit$var <- expected
  unnamed_fit$orig.var <- expected
  unnamed <- get_delta_cluster_vcov(unnamed_fit)
  expect_identical(rownames(unnamed$vcov), coefficient_names)
  expect_identical(colnames(unnamed$vcov), coefficient_names)
  expect_identical(unname(diag(unnamed$vcov)), as.numeric(seq_len(q)))

  rms_names <- c(
    coefficient_names[seq_len(fit$non.slopes)],
    fit$Design$mmcolnames
  )
  aliased_fit <- fit
  aliased_fit$var <- expected
  aliased_fit$orig.var <- expected
  dimnames(aliased_fit$var) <- list(rms_names, rms_names)
  dimnames(aliased_fit$orig.var) <- list(rms_names, rms_names)
  aliased <- get_delta_cluster_vcov(aliased_fit)
  expect_identical(rownames(aliased$vcov), coefficient_names)
  expect_identical(colnames(aliased$vcov), coefficient_names)
  expect_identical(unname(diag(aliased$vcov)), as.numeric(seq_len(q)))

  condition <- tryCatch(
    get_delta_cluster_vcov(fit, vcov = aliased_fit$var),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(conditionMessage(condition), "must be uniquely named")
})

test_that("penalized orm previous-state splines support delta inference", {
  case <- local_spline_delta_orm_case()
  fit <- case$model
  avg <- avg_sops(
    fit,
    variables = list(tx = c(0, 1)),
    times = 1:3,
    y_levels = fit$yunique,
    p_var = "yprev"
  )

  plan <- delta_compile_plan(
    object = avg,
    model_plan = fit,
    newdata = attr(avg, "newdata_pred"),
    times = attr(avg, "avg_args")$times,
    output = "average"
  )
  native <- run_sop_delta_plan(plan, fit)
  oracle <- run_sop_delta_r_oracle(plan, fit)
  expect_equal(native$probabilities, oracle$probabilities, tolerance = 1e-12)
  expect_equal(native$jacobian, oracle$jacobian, tolerance = 1e-12)

  inferred <- inferences(avg, method = "delta", target = "empirical")
  numerical <- central_spline_avg_jacobian(avg, fit)
  analytical <- get_jacobian(inferred)

  expect_identical(colnames(analytical), names(stats::coef(fit)))
  expect_lt(max(abs(analytical - numerical)), 1e-7)
  expect_all_true(is.finite(inferred$std.error))
  expect_all_true(inferred$std.error >= 0)

  condition <- tryCatch(
    inferences(avg, method = "delta", target = "superpopulation"),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(
    conditionMessage(condition),
    "does not currently support penalized orm likelihoods",
    fixed = TRUE
  )
})

test_that("superpopulation delta rejects weighted orm score construction", {
  case <- local_spline_delta_orm_case()
  data <- case$data
  local_spline_delta_datadist(data)
  weighted <- suppressWarnings(orm_markov(
    y ~ rms::rcs(time, 3) + tx + rms::rcs(yprev, 4),
    data = data,
    weights = rep(c(1, 2), length.out = nrow(data)),
    id_var = "id",
    opt_method = "LM",
    scale = TRUE
  ))
  avg <- avg_sops(
    weighted,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    y_levels = weighted$yunique,
    p_var = "yprev"
  )

  condition <- tryCatch(
    inferences(avg, method = "delta", target = "superpopulation"),
    error = identity
  )
  expect_s3_class(condition, "error")
  expect_match(
    conditionMessage(condition),
    "does not currently support orm fits with non-unit case weights",
    fixed = TRUE
  )
})
