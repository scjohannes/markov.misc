local_robcov_orm_case <- local({
  value <- NULL
  function() {
    skip_if_not_installed("rms")
    if (is.null(value)) {
      data <- suppressWarnings(make_test_data(
        n_patients = 28,
        follow_up_time = 5,
        seed = 2811
      ))
      fit <- suppressWarnings(rms::orm(
        ordered(y) ~ time + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      ))
      value <<- list(data = data, fit = fit)
    }
    value
  }
})

orm_test_row_loglik <- function(model, row, coefficient) {
  gamma <- get_effective_coefs(model, beta = coefficient)
  design <- cbind("(Intercept)" = 1, as.matrix(model$x))
  design <- design[, colnames(gamma), drop = FALSE]
  eta <- design[row, , drop = FALSE] %*% t(gamma)
  probabilities <- lp_to_probs(eta, nrow(gamma))
  response_index <- match(
    as_state_labels(model$y[row]),
    as_state_labels(model$yunique)
  )
  weight <- orm_case_weights(model, nrow(model$x))[row]
  weight * log(probabilities[1L, response_index])
}

test_that("orm HC0 covariance equals a directly assembled cluster sandwich", {
  case <- local_robcov_orm_case()
  fit <- case$fit
  cluster <- case$data$id
  robust <- suppressWarnings(robcov_orm(
    fit,
    cluster = cluster,
    type = "HC0",
    cadjust = FALSE
  ))

  scores <- compute_scores_orm(fit)
  cluster_ids <- unique(as.character(cluster))
  clustered_scores <- rowsum(
    scores,
    factor(as.character(cluster), levels = cluster_ids),
    reorder = FALSE
  )
  bread <- orm_model_bread(fit)$bread
  expected <- crossprod(clustered_scores %*% bread)
  dimnames(expected) <- dimnames(robust$var)

  expect_equal(robust$var, expected, tolerance = 1e-12)
  expect_equal(
    attr(robust, "markov_robust_covariance")$adjustment_factor,
    1
  )
})

test_that("orm HC1 and cluster corrections are exact scalar adjustments", {
  case <- local_robcov_orm_case()
  fit <- case$fit
  cluster <- case$data$id
  n <- nrow(fit$x)
  p <- length(stats::coef(fit))
  n_clusters <- length(unique(cluster))

  hc0 <- robcov_orm(fit, cluster, type = "HC0", cadjust = FALSE)
  hc1 <- robcov_orm(fit, cluster, type = "HC1", cadjust = FALSE)
  cluster_adjusted <- robcov_orm(
    fit,
    cluster,
    type = "HC0",
    cadjust = TRUE
  )
  combined <- robcov_orm(fit, cluster, type = "HC1", cadjust = TRUE)

  hc1_factor <- (n - 1) / (n - p)
  cluster_factor <- n_clusters / (n_clusters - 1)
  expect_equal(hc1$var, hc0$var * hc1_factor, tolerance = 1e-12)
  expect_equal(
    cluster_adjusted$var,
    hc0$var * cluster_factor,
    tolerance = 1e-12
  )
  expect_equal(
    combined$var,
    hc0$var * hc1_factor * cluster_factor,
    tolerance = 1e-12
  )
})

test_that("orm robust covariance aggregates repeated rows by represented cluster", {
  case <- local_robcov_orm_case()
  fit <- case$fit
  cluster <- factor(case$data$id, levels = c(unique(case$data$id), 9999))
  robust <- robcov_orm(fit, cluster, type = "HC0", cadjust = FALSE)
  metadata <- attr(robust, "markov_robust_covariance")

  scores <- compute_scores_orm(fit)
  cluster_ids <- unique(as.character(cluster))
  clustered_scores <- rowsum(
    scores,
    factor(as.character(cluster), levels = cluster_ids),
    reorder = FALSE
  )

  expect_lt(nrow(clustered_scores), nrow(scores))
  expect_equal(metadata$n_clusters, length(unique(case$data$id)))
  expect_false("9999" %in% rownames(clustered_scores))
})

test_that("orm_markov stores package robust covariance metadata by default", {
  case <- local_robcov_orm_case()
  fit <- orm_markov(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    id_var = "id"
  )
  metadata <- attr(fit, "markov_robust_covariance")
  n_clusters <- length(unique(case$data$id))

  expect_s3_class(fit, "orm")
  expect_equal(metadata$implementation, "markov.misc::robcov_orm")
  expect_equal(metadata$backend, "orm")
  expect_equal(metadata$type, "HC0")
  expect_true(metadata$cadjust)
  expect_equal(metadata$n_clusters, n_clusters)
  expect_equal(metadata$adjustment_factor, n_clusters / (n_clusters - 1))
  expect_identical(metadata$covariance_identity, fit$var)
})

test_that("weighted orm scores equal numerical weighted row likelihood derivatives", {
  case <- local_robcov_orm_case()
  weights <- rep(c(0.5, 1, 2), length.out = nrow(case$data))
  weighted <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    weights = weights,
    x = TRUE,
    y = TRUE
  ))
  scores <- compute_scores_orm(weighted)
  coefficient <- stats::coef(weighted)
  rows <- unique(c(1L, nrow(weighted$x) %/% 2L, nrow(weighted$x)))
  numerical <- matrix(
    NA_real_,
    nrow = length(rows),
    ncol = length(coefficient),
    dimnames = list(NULL, names(coefficient))
  )
  rownames(numerical) <- rownames(scores)[rows]

  for (i in seq_along(rows)) {
    for (j in seq_along(coefficient)) {
      step <- .Machine$double.eps^(1 / 3) * max(1, abs(coefficient[j]))
      upper <- lower <- coefficient
      upper[j] <- upper[j] + step
      lower[j] <- lower[j] - step
      numerical[i, j] <- (orm_test_row_loglik(weighted, rows[i], upper) -
        orm_test_row_loglik(weighted, rows[i], lower)) /
        (2 * step)
    }
  }

  expect_equal(scores[rows, ], numerical, tolerance = 1e-7)

  robust <- robcov_orm(
    weighted,
    cluster = case$data$id,
    type = "HC0",
    cadjust = TRUE
  )
  expect_true(attr(robust, "markov_robust_covariance")$weighted)
  expect_true(all(is.finite(robust$var)))
  expect_true(all(diag(robust$var) >= 0))

  weighted_data <- case$data
  weighted_data$case_weight <- weights
  wrapped <- suppressWarnings(orm_markov(
    ordered(y) ~ time + tx + yprev,
    data = weighted_data,
    weights = case_weight,
    id_var = "id"
  ))
  average <- avg_sops(
    wrapped,
    variables = list(tx = c(0, 1)),
    times = 1:2,
    y_levels = wrapped$yunique,
    absorb = 6,
    p_var = "yprev"
  )
  inferred <- inferences(average, method = "delta", vcov = "conditional")
  expect_true(all(is.finite(inferred$std.error)))
  expect_true(all(inferred$std.error >= 0))
})

test_that("orm robust covariance rejects invalid fits and inputs", {
  case <- local_robcov_orm_case()
  fit <- case$fit
  cluster <- case$data$id

  missing_cluster <- cluster
  missing_cluster[1L] <- NA
  expect_error(robcov_orm(fit, integer()), "does not match")
  expect_error(robcov_orm(fit, missing_cluster), "contains missing values")
  expect_error(
    robcov_orm(fit, rep(1, length(cluster))),
    "at least two clusters"
  )

  invalid_weights <- fit
  invalid_weights$weights <- rep(1, nrow(fit$x))
  invalid_weights$weights[1L] <- -1
  expect_error(robcov_orm(invalid_weights, cluster), "finite and nonnegative")

  missing_x <- fit
  missing_x$x <- NULL
  expect_error(
    robcov_orm(missing_x, cluster),
    "x = TRUE, y = TRUE",
    fixed = TRUE
  )
  missing_y <- fit
  missing_y$y <- NULL
  expect_error(
    robcov_orm(missing_y, cluster),
    "x = TRUE, y = TRUE",
    fixed = TRUE
  )

  non_logit <- fit
  non_logit$family <- "probit"
  expect_error(robcov_orm(non_logit, cluster), "only ordinary ordinal logistic")

  too_small <- fit
  p <- length(stats::coef(too_small))
  too_small$x <- too_small$x[seq_len(p), , drop = FALSE]
  too_small$y <- too_small$y[seq_len(p)]
  too_small$weights <- rep(1, p)
  expect_error(
    robcov_orm(too_small, cluster = seq_len(p), type = "HC1"),
    "requires more observations than parameters"
  )

  malformed_bread <- fit
  malformed_bread$orig.var <- diag(2)
  expect_error(robcov_orm(malformed_bread, cluster), "orm model bread")
})

test_that("stored orm covariance integrity detects mutation", {
  case <- local_robcov_orm_case()
  robust <- robcov_orm(case$fit, case$data$id)
  expect_true(orm_stored_covariance_valid(robust, case$data$id))

  mutated <- robust
  mutated$var[1L, 1L] <- mutated$var[1L, 1L] + 1
  expect_false(orm_stored_covariance_valid(mutated, case$data$id))
  expect_false(orm_stored_covariance_valid(robust, rev(case$data$id)))
})

test_that("penalized sandwich orm uses its dedicated inverse-sensitivity bread", {
  case <- local_robcov_orm_case()
  penalized <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = case$data,
    penalty = 0.1,
    var.penalty = "sandwich",
    x = TRUE,
    y = TRUE
  ))
  expect_false(is.null(penalized$var.from.info.matrix))

  robust <- robcov_orm(
    penalized,
    cluster = case$data$id,
    type = "HC0",
    cadjust = FALSE
  )
  expected_bread <- delta_normalize_orm_backend_matrix(
    penalized$var.from.info.matrix,
    penalized,
    names(stats::coef(penalized))
  )
  metadata <- attr(robust, "markov_robust_covariance")

  expect_equal(robust$orig.var, expected_bread, tolerance = 1e-12)
  expect_equal(metadata$bread_convention, "var.from.info.matrix")

  missing_bread <- penalized
  missing_bread$var.from.info.matrix <- NULL
  expect_error(
    robcov_orm(missing_bread, cluster = case$data$id),
    "does not retain the required inverse-sensitivity matrix",
    fixed = TRUE
  )

  malformed_bread <- penalized
  malformed_bread$var.from.info.matrix <- diag(2)
  expect_error(
    robcov_orm(malformed_bread, cluster = case$data$id),
    "inverse-sensitivity matrix stored in `var.from.info.matrix` is invalid",
    fixed = TRUE
  )
})

test_that("get_vcov_robust inherits and explicitly overrides stored orm corrections", {
  case <- local_robcov_orm_case()
  stored <- robcov_orm(
    case$fit,
    cluster = case$data$id,
    type = "HC1",
    cadjust = FALSE
  )

  expect_identical(get_vcov_robust(stored), stored$var)
  expect_equal(
    get_vcov_robust(stored, cluster = case$data$id),
    stored$var,
    tolerance = 1e-12
  )

  expected_type_override <- robcov_orm(
    case$fit,
    cluster = case$data$id,
    type = "HC0",
    cadjust = FALSE
  )$var
  expect_equal(
    get_vcov_robust(
      stored,
      cluster = case$data$id,
      type = "HC0"
    ),
    expected_type_override,
    tolerance = 1e-12
  )

  expected_cadjust_override <- robcov_orm(
    case$fit,
    cluster = case$data$id,
    type = "HC1",
    cadjust = TRUE
  )$var
  expect_equal(
    get_vcov_robust(
      stored,
      cluster = case$data$id,
      cadjust = TRUE
    ),
    expected_cadjust_override,
    tolerance = 1e-12
  )
})

test_that("zero-weight rows and clusters do not enter orm corrections", {
  case <- local_robcov_orm_case()
  weights <- rep(1, nrow(case$data))
  zero_cluster_ids <- unique(case$data$id)[1:2]
  weights[case$data$id %in% zero_cluster_ids] <- 0
  weights[which(!case$data$id %in% zero_cluster_ids)[1L]] <- 0
  weighted_data <- case$data
  weighted_data$case_weight <- weights
  weighted <- suppressWarnings(rms::orm(
    ordered(y) ~ time + tx + yprev,
    data = weighted_data,
    weights = case_weight,
    x = TRUE,
    y = TRUE
  ))

  fitted_weights <- orm_case_weights(weighted, nrow(weighted$x))
  represented <- fitted_weights > 0
  represented_clusters <- unique(case$data$id[represented])
  n_represented <- sum(represented)
  n_clusters <- length(represented_clusters)
  p <- length(stats::coef(weighted))

  hc0 <- robcov_orm(
    weighted,
    cluster = case$data$id,
    type = "HC0",
    cadjust = FALSE
  )
  corrected <- robcov_orm(
    weighted,
    cluster = case$data$id,
    type = "HC1",
    cadjust = TRUE
  )
  hc0_metadata <- attr(hc0, "markov_robust_covariance")
  corrected_metadata <- attr(corrected, "markov_robust_covariance")
  expected_factor <-
    (n_represented - 1) / (n_represented - p) * n_clusters / (n_clusters - 1)

  expect_equal(hc0_metadata$n_rows, nrow(weighted$x))
  expect_equal(hc0_metadata$n_represented_rows, n_represented)
  expect_equal(hc0_metadata$n_clusters, n_clusters)
  expect_identical(hc0_metadata$cluster, case$data$id)
  expect_equal(corrected_metadata$adjustment_factor, expected_factor)
  expect_equal(corrected$var, hc0$var * expected_factor, tolerance = 1e-12)
  expect_false(any(zero_cluster_ids %in% represented_clusters))
})
