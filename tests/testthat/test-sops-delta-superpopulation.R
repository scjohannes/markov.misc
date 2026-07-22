local_delta_vglm_fit <- local({
  fit <- NULL
  function() {
    skip_if_not_installed("VGAM")
    if (is.null(fit)) {
      data <- make_test_data(
        n_patients = 35,
        follow_up_time = 6,
        seed = 2401
      )
      fit <<- vglm_markov(
        ordered(y) ~ time_lin + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data,
        id_var = "id"
      )
    }
    fit
  }
})

local_delta_orm_fit <- local({
  fit <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(fit)) {
      data <- make_test_data(
        n_patients = 35,
        follow_up_time = 6,
        seed = 2402
      )
      fit <<- orm_markov(
        ordered(y) ~ time + tx + yprev,
        data = data,
        id_var = "id"
      )
    }
    fit
  }
})

delta_weighted_sop_estimate <- function(
  model,
  coefficient,
  prediction_data,
  profile_weights,
  variables,
  times,
  y_levels,
  absorb
) {
  weighted_model <- set_coef(model, coefficient)
  probabilities <- soprob_markov(
    model = weighted_model,
    newdata = prediction_data,
    times = times,
    y_levels = y_levels,
    absorb = absorb
  )
  grid <- do.call(expand.grid, variables)
  result <- marginalize_sops_array(
    sops_array = probabilities,
    grid = grid,
    times = times,
    y_levels = y_levels,
    variables = variables,
    n_cf = nrow(grid),
    n_each = length(profile_weights),
    weights = profile_weights,
    newdata = prediction_data
  )
  result$estimate
}

delta_weighted_refit_oracle <- function(
  model,
  fit_data,
  profile_data,
  refit,
  patient_ids,
  variables = list(tx = c(0, 1)),
  times = 1:2,
  y_levels = 1:6,
  absorb = 6,
  step = 1e-2
) {
  point <- avg_sops(
    model,
    variables = variables,
    times = times,
    y_levels = y_levels,
    absorb = absorb
  )
  inferred <- inferences(
    point,
    method = "delta",
    target = "superpopulation"
  )
  analytical <- attr(inferred, "analytical")
  profile_ids <- as.character(profile_data$id)
  n <- length(profile_ids)
  prediction_data <- attr(point, "newdata_pred")
  base_weights <- rep(1, n)
  components <- get_delta_score_components(model)

  numerical_coefficient <- matrix(
    NA_real_,
    nrow = length(patient_ids),
    ncol = length(stats::coef(model)),
    dimnames = list(patient_ids, names(stats::coef(model)))
  )
  numerical_sop <- matrix(
    NA_real_,
    nrow = length(patient_ids),
    ncol = nrow(point),
    dimnames = list(patient_ids, NULL)
  )
  profile_only_sop <- numerical_sop

  for (patient_id in patient_ids) {
    fit_patient <- as.character(fit_data$id) == patient_id
    profile_patient <- profile_ids == patient_id
    upper_fit_weights <- lower_fit_weights <- rep(1, nrow(fit_data))
    upper_fit_weights[fit_patient] <- 1 + step
    lower_fit_weights[fit_patient] <- 1 - step
    upper_profile_weights <- lower_profile_weights <- base_weights
    upper_profile_weights[profile_patient] <- 1 + step
    lower_profile_weights[profile_patient] <- 1 - step

    upper_coefficient <- refit(upper_fit_weights)
    lower_coefficient <- refit(lower_fit_weights)
    numerical_coefficient[patient_id, ] <- n *
      (upper_coefficient - lower_coefficient) /
      (2 * step)

    upper_sop <- delta_weighted_sop_estimate(
      model = model,
      coefficient = upper_coefficient,
      prediction_data = prediction_data,
      profile_weights = upper_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    lower_sop <- delta_weighted_sop_estimate(
      model = model,
      coefficient = lower_coefficient,
      prediction_data = prediction_data,
      profile_weights = lower_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    numerical_sop[patient_id, ] <- n *
      (upper_sop - lower_sop) /
      (2 * step)

    upper_profile_only <- delta_weighted_sop_estimate(
      model = model,
      coefficient = stats::coef(model),
      prediction_data = prediction_data,
      profile_weights = upper_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    lower_profile_only <- delta_weighted_sop_estimate(
      model = model,
      coefficient = stats::coef(model),
      prediction_data = prediction_data,
      profile_weights = lower_profile_weights,
      variables = variables,
      times = times,
      y_levels = y_levels,
      absorb = absorb
    )
    profile_only_sop[patient_id, ] <- n *
      (upper_profile_only - lower_profile_only) /
      (2 * step)
  }

  score_rows <- match(patient_ids, components$ids)
  expected_coefficient <- components$scores[score_rows, , drop = FALSE] %*%
    t(components$Ainv)
  rownames(expected_coefficient) <- patient_ids
  influence_rows <- match(patient_ids, analytical$profile_ids)

  list(
    numerical_coefficient = numerical_coefficient,
    expected_coefficient = expected_coefficient,
    numerical_sop = numerical_sop,
    profile_only_sop = profile_only_sop,
    expected_sop = analytical$influence[influence_rows, , drop = FALSE]
  )
}

test_that("patient clusters resolve explicitly or from stored fitting metadata", {
  fit <- local_delta_vglm_fit()
  fitted_data <- attr(fit, "markov_data")

  stored <- resolve_delta_cluster(fit)
  expect_equal(stored$cluster, fitted_data$id)
  expect_equal(stored$ids, unique(as.character(fitted_data$id)))
  expect_equal(stored$metadata$source, "stored_id_var")
  expect_equal(stored$metadata$cluster_level, "patient")

  explicit <- paste0("patient-", fitted_data$id)
  resolved_explicit <- resolve_delta_cluster(fit, cluster = explicit)
  expect_equal(resolved_explicit$cluster, explicit)
  expect_equal(resolved_explicit$metadata$source, "explicit_vector")

  resolved_formula <- resolve_delta_cluster(fit, cluster = ~id)
  expect_equal(resolved_formula$cluster, fitted_data$id)
  expect_equal(resolved_formula$metadata$source, "explicit_formula")
})

test_that("patient clusters never fall back to fitting rows", {
  fit <- local_delta_vglm_fit()
  attr(fit, "markov_data") <- NULL
  attr(fit, "markov_id_var") <- NULL
  attr(fit$vglm_fit, "markov_data") <- NULL
  attr(fit$vglm_fit, "markov_id_var") <- NULL

  expect_error(
    resolve_delta_cluster(fit),
    "requires patient clustering",
    fixed = TRUE
  )
})

test_that("explicit coefficient covariance has precedence and strict validation", {
  fit <- local_delta_vglm_fit()
  attr(fit, "markov_data") <- NULL
  attr(fit, "markov_id_var") <- NULL
  attr(fit$vglm_fit, "markov_data") <- NULL
  attr(fit$vglm_fit, "markov_id_var") <- NULL
  coefficient_names <- names(fit$coefficients)
  reversed_names <- rev(coefficient_names)
  covariance <- diag(seq_along(reversed_names))
  dimnames(covariance) <- list(reversed_names, reversed_names)

  resolved <- get_delta_cluster_vcov(fit, vcov = covariance)
  expect_identical(rownames(resolved$vcov), coefficient_names)
  expect_identical(colnames(resolved$vcov), coefficient_names)
  expect_equal(resolved$metadata$source, "explicit")
  expect_null(resolved$cluster_info)

  unnamed <- covariance
  dimnames(unnamed) <- NULL
  expect_error(
    get_delta_cluster_vcov(fit, vcov = unnamed),
    "must be uniquely named",
    fixed = TRUE
  )

  nonfinite <- covariance
  nonfinite[1L, 1L] <- Inf
  expect_error(
    get_delta_cluster_vcov(fit, vcov = nonfinite),
    "contains non-finite values",
    fixed = TRUE
  )

  nonsymmetric <- covariance
  nonsymmetric[1L, 2L] <- 1
  expect_error(
    get_delta_cluster_vcov(fit, vcov = nonsymmetric),
    "not numerically symmetric",
    fixed = TRUE
  )

  non_psd <- covariance
  non_psd[1L, 1L] <- -1
  expect_error(
    get_delta_cluster_vcov(fit, vcov = non_psd),
    "not numerically positive semidefinite",
    fixed = TRUE
  )
})

test_that("vglm covariance and score components preserve backend conventions", {
  fit <- local_delta_vglm_fit()
  fitted_data <- attr(fit, "markov_data")
  covariance <- get_delta_cluster_vcov(fit)

  expect_equal(covariance$vcov, fit$var)
  expect_equal(covariance$metadata$type, fit$type)
  expect_identical(covariance$metadata$cadjust, fit$cadjust)
  expect_equal(
    covariance$metadata$adjustment_factor,
    fit$adjustment_factor
  )

  components <- get_delta_score_components(fit)
  cluster_ids <- unique(as.character(fitted_data$id))
  expected_scores <- rowsum(
    fit$scores,
    factor(as.character(fitted_data$id), levels = cluster_ids),
    reorder = FALSE
  )
  rownames(expected_scores) <- cluster_ids

  expect_equal(components$ids, cluster_ids)
  expect_equal(components$scores, expected_scores)
  expect_equal(components$n, length(cluster_ids))
  expect_equal(components$Ainv, length(cluster_ids) * fit$bread)
  expect_false(isTRUE(all.equal(components$Ainv, fit$var)))
  expect_equal(components$metadata$bread_source, "robcov_vglm_observed")
  expect_equal(components$metadata$backend_hc_type_ignored, fit$type)
  expect_identical(
    components$metadata$backend_cadjust_ignored,
    fit$cadjust
  )
  expect_equal(components$metadata$n_patients, length(cluster_ids))
})

test_that("raw vglm models compute patient covariance and score components", {
  fit <- local_delta_vglm_fit()
  raw_fit <- fit$vglm_fit
  fitted_data <- attr(fit, "markov_data")
  cluster_ids <- unique(as.character(fitted_data$id))

  covariance <- get_delta_cluster_vcov(raw_fit, cluster = fitted_data$id)
  components <- get_delta_score_components(raw_fit, cluster = fitted_data$id)
  expected_scores <- rowsum(
    fit$scores,
    factor(as.character(fitted_data$id), levels = cluster_ids),
    reorder = FALSE
  )
  rownames(expected_scores) <- cluster_ids

  expect_equal(covariance$vcov, fit$var, tolerance = 1e-10)
  expect_equal(covariance$metadata$source, "computed_robcov_vglm")
  expect_equal(components$scores, expected_scores)
  expect_equal(components$Ainv, length(cluster_ids) * fit$bread)
})

test_that("orm covariance and score components use inverse total information", {
  fit <- local_delta_orm_fit()
  fitted_data <- attr(fit, "markov_data")
  covariance <- get_delta_cluster_vcov(fit)

  expect_equal(covariance$vcov, fit$var)
  expect_equal(covariance$metadata$type, "HC0")
  expect_false(covariance$metadata$cadjust)

  recomputed <- get_delta_cluster_vcov(fit, cluster = fitted_data$id)
  expect_equal(recomputed$vcov, fit$var, tolerance = 1e-10)
  expect_equal(recomputed$metadata$source, "computed_rms_robcov")

  components <- get_delta_score_components(fit)
  row_scores <- compute_scores_orm(fit)
  cluster_ids <- unique(as.character(fitted_data$id))
  expected_scores <- rowsum(
    row_scores,
    factor(as.character(fitted_data$id), levels = cluster_ids),
    reorder = FALSE
  )
  rownames(expected_scores) <- cluster_ids
  expected_bread <- get_orm_model_vcov(fit)

  expect_equal(components$ids, cluster_ids)
  expect_equal(components$scores, expected_scores)
  expect_equal(components$Ainv, length(cluster_ids) * expected_bread)
  expect_equal(
    components$metadata$bread_source,
    "rms_inverse_total_information"
  )
})

test_that("vglm patient-weight refits validate stacked score orientation", {
  fit <- local_delta_vglm_fit()
  fit_data <- attr(fit, "markov_data")
  profile_data <- attr(fit, "markov_starting_profile_data")
  components <- get_delta_score_components(fit)
  score_size <- rowSums(abs(components$scores))
  patient_ids <- components$ids[order(score_size, decreasing = TRUE)[1:3]]

  oracle <- delta_weighted_refit_oracle(
    model = fit,
    fit_data = fit_data,
    profile_data = profile_data,
    refit = function(weights) {
      weighted <- VGAM::vglm(
        ordered(y) ~ time_lin + tx + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = fit_data,
        weights = weights,
        coefstart = stats::coef(fit),
        control = VGAM::vglm.control(epsilon = 1e-10, maxit = 100)
      )
      stats::coef(weighted)
    },
    patient_ids = patient_ids
  )

  expect_equal(
    oracle$numerical_coefficient,
    oracle$expected_coefficient,
    tolerance = 2e-3
  )
  expect_equal(
    oracle$numerical_sop,
    oracle$expected_sop,
    tolerance = 2e-3
  )
  expect_gt(max(abs(oracle$numerical_sop - oracle$profile_only_sop)), 1e-3)
})

test_that("orm patient-weight refits validate stacked score orientation", {
  fit <- local_delta_orm_fit()
  fit_data <- attr(fit, "markov_data")
  profile_data <- attr(fit, "markov_starting_profile_data")
  components <- get_delta_score_components(fit)
  score_size <- rowSums(abs(components$scores))
  patient_ids <- components$ids[order(score_size, decreasing = TRUE)[1:3]]

  oracle <- delta_weighted_refit_oracle(
    model = fit,
    fit_data = fit_data,
    profile_data = profile_data,
    refit = function(weights) {
      weighted <- suppressWarnings(rms::orm(
        ordered(y) ~ time + tx + yprev,
        data = fit_data,
        weights = weights,
        x = TRUE,
        y = TRUE
      ))
      stats::coef(weighted)
    },
    patient_ids = patient_ids
  )

  expect_equal(
    oracle$numerical_coefficient,
    oracle$expected_coefficient,
    tolerance = 2e-3
  )
  expect_equal(
    oracle$numerical_sop,
    oracle$expected_sop,
    tolerance = 2e-3
  )
  expect_gt(max(abs(oracle$numerical_sop - oracle$profile_only_sop)), 1e-3)
})

test_that("superpopulation score machinery rejects partial proportional odds", {
  skip_if_not_installed("VGAM")
  data <- make_test_data(
    n_patients = 35,
    follow_up_time = 6,
    seed = 2403
  )
  fit <- VGAM::vglm(
    ordered(y) ~ time_lin + tx + yprev,
    family = VGAM::cumulative(
      reverse = TRUE,
      parallel = FALSE ~ tx
    ),
    data = data
  )

  expect_error(
    resolve_delta_cluster(fit, cluster = data$id),
    "requires a full proportional-odds vglm model",
    fixed = TRUE
  )
})

test_that("stacked influence matches fitted patients and retains the cross term", {
  values <- rbind(
    p1 = c(cell_a = 0.2, cell_b = 0.8),
    p2 = c(cell_a = 0.5, cell_b = 0.5)
  )
  jacobian <- rbind(
    cell_a = c(beta_1 = 1, beta_2 = 0.5),
    cell_b = c(beta_1 = -1, beta_2 = -0.5)
  )
  scores <- rbind(
    p2 = c(beta_1 = -0.3, beta_2 = 0.4),
    p1 = c(beta_1 = 0.2, beta_2 = -0.1)
  )
  Ainv <- diag(c(2, 1))
  dimnames(Ainv) <- list(c("beta_1", "beta_2"), c("beta_1", "beta_2"))
  components <- list(
    ids = rownames(scores),
    scores = scores,
    Ainv = Ainv,
    n = nrow(scores)
  )

  result <- delta_stacked_influence(
    individual_values = values,
    average_jacobian = jacobian,
    profile_ids = rownames(values),
    score_components = components
  )

  aligned_scores <- rbind(
    p1 = scores["p1", ],
    p2 = scores["p2", ]
  )
  coefficient_influence <- aligned_scores %*% t(Ainv) %*% t(jacobian)
  profile_component <- sweep(values, 2L, colMeans(values), "-")
  expected_influence <- profile_component + coefficient_influence
  expected_covariance <- stats::cov(expected_influence) / nrow(values)
  centered <- sweep(
    expected_influence,
    2L,
    colMeans(expected_influence),
    "-"
  )
  centered_hc0 <- crossprod(centered) / nrow(values)^2
  raw_uncentered <- crossprod(expected_influence) / nrow(values)^2
  influence_mean <- colMeans(expected_influence)

  expect_equal(result$scores, aligned_scores)
  expect_equal(result$Ainv, Ainv)
  expect_equal(result$coefficient_influence, coefficient_influence)
  expect_equal(result$influence, expected_influence)
  expect_equal(result$covariance, expected_covariance)
  expect_equal(result$std.error, sqrt(diag(expected_covariance)))
  expect_equal(
    result$covariance,
    crossprod(centered) / (nrow(values) * (nrow(values) - 1L))
  )
  expect_equal(
    result$covariance,
    centered_hc0 * nrow(values) / (nrow(values) - 1L)
  )
  expect_equal(
    result$covariance,
    nrow(values) /
      (nrow(values) - 1L) *
      raw_uncentered -
      tcrossprod(influence_mean) / (nrow(values) - 1L)
  )
  expect_equal(
    result$metadata$finite_sample_method,
    "patient_level_sample_covariance"
  )
  expect_equal(result$metadata$finite_sample_factor, 2)
  expect_equal(result$metadata$n_patients, 2)
  expect_equal(result$metadata$n_score_clusters, 2)
  expect_null(result$metadata$missing_score_profiles)
  expect_null(result$metadata$sensitivity_scale)

  additive_covariance <- stats::cov(profile_component) /
    nrow(values) +
    stats::cov(coefficient_influence) / nrow(values)
  expect_gt(max(abs(result$covariance - additive_covariance)), 1e-8)
})

test_that("stacked influence rejects profile-only patients", {
  values <- rbind(
    p1 = c(cell = 0.2),
    p2 = c(cell = 0.4),
    p3 = c(cell = 0.6)
  )
  jacobian <- matrix(1, nrow = 1L, dimnames = list("cell", "beta"))
  scores <- matrix(
    c(0.1, -0.1),
    ncol = 1L,
    dimnames = list(c("p1", "p2"), "beta")
  )

  expect_error(
    delta_stacked_influence(
      values,
      jacobian,
      rownames(values),
      list(
        ids = rownames(scores),
        scores = scores,
        Ainv = matrix(1, dimnames = list("beta", "beta")),
        n = 2L
      )
    ),
    "Target-profile patients have no usable likelihood transition: p3",
    fixed = TRUE
  )
})

test_that("stacked influence rejects score patients absent from profiles", {
  values <- rbind(
    p1 = c(cell = 0.2),
    p2 = c(cell = 0.4)
  )
  jacobian <- matrix(
    1,
    nrow = 1L,
    dimnames = list("cell", "beta")
  )
  scores <- matrix(
    c(0.1, -0.1),
    ncol = 1L,
    dimnames = list(c("p1", "p3"), "beta")
  )
  Ainv <- matrix(1, dimnames = list("beta", "beta"))

  expect_error(
    delta_stacked_influence(
      values,
      jacobian,
      rownames(values),
      list(ids = rownames(scores), scores = scores, Ainv = Ainv, n = 2L)
    ),
    "Likelihood-score patients are absent from the target profiles: p3",
    fixed = TRUE
  )
})

test_that("empirical covariance propagates named raw coefficients", {
  jacobian <- rbind(
    first = c(beta_2 = 2, beta_1 = 1),
    second = c(beta_2 = -1, beta_1 = 0.5)
  )
  covariance <- matrix(
    c(4, 1, 1, 2),
    nrow = 2L,
    dimnames = list(c("beta_1", "beta_2"), c("beta_1", "beta_2"))
  )
  expected <- jacobian %*%
    covariance[colnames(jacobian), colnames(jacobian)] %*%
    t(jacobian)

  expect_equal(delta_empirical_covariance(jacobian, covariance), expected)
})
