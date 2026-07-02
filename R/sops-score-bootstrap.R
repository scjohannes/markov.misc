# SOP score-bootstrap draws.

#' Generate One-Step Score-Bootstrap Draws
#'
#' Computes simulation draws using cluster-level score perturbations and a
#' one-step Newton approximation. Also returns patient-level averaging weights
#' (shared with score perturbations) for `avg_sops` marginalization.
#'
#' @param model A `robcov_vglm` object with stored `scores`, `bread`, and
#'   `cluster` components, or an `orm` object fitted with `x = TRUE, y = TRUE`.
#'   For `orm`, `cluster` must be supplied.
#' @param baseline_data Optional baseline data (one row per patient) used only
#'   when patient-level averaging weights should be returned.
#' @param id_var Name of patient ID variable in `baseline_data`.
#' @param n_sim Number of simulation draws.
#' @param score_weight_dist Cluster weight distribution. Currently only
#'   `"exponential"` is supported.
#' @param cluster Optional row-level cluster vector for `orm` models.
#' @param return_baseline_weights Logical. If `TRUE`, return normalized
#'   patient-level averaging weights aligned to `baseline_data`. If `FALSE`,
#'   only coefficient perturbations are generated.
#'
#' @return A list with:
#'   \itemize{
#'     \item `beta_draws`: Matrix `[n_sim x p]` of perturbed coefficients.
#'     \item `baseline_weights`: Matrix `[n_sim x n_patients]` of normalized
#'       patient averaging weights.
#'   }
#'
#' @keywords internal
generate_score_bootstrap_draws <- function(
  model,
  baseline_data,
  id_var,
  n_sim,
  score_weight_dist = "exponential",
  cluster = NULL,
  return_baseline_weights = TRUE
) {
  if (!inherits(model, "robcov_vglm") && !inherits(model, "orm")) {
    stop(
      "`engine = \"score_bootstrap\"` requires a 'robcov_vglm' model ",
      "or an 'orm' model with `cluster` supplied. Wrap vglm fits with ",
      "robcov_vglm(fit, cluster = <id>), or call inferences(..., ",
      "engine = \"score_bootstrap\", cluster = <id>) for orm fits."
    )
  }

  score_weight_dist <- match.arg(score_weight_dist, choices = "exponential")

  if (isTRUE(return_baseline_weights)) {
    if (is.null(baseline_data) || !is.data.frame(baseline_data)) {
      stop("`baseline_data` must be supplied as a data frame.")
    }
    if (!id_var %in% names(baseline_data)) {
      stop("ID variable '", id_var, "' not found in baseline data.")
    }
  }

  components <- score_bootstrap_components(model, cluster = cluster)
  beta_hat <- components$coefficients
  bread <- components$bread
  scores <- components$scores
  cluster <- components$cluster

  if (!is.matrix(scores)) {
    stop("`model$scores` must be a matrix for score bootstrap.")
  }
  if (!is.matrix(bread)) {
    stop("`model$bread` must be a matrix for score bootstrap.")
  }
  if (length(cluster) != nrow(scores)) {
    stop(
      "Length of `model$cluster` (",
      length(cluster),
      ") does not match ",
      "the number of score rows (",
      nrow(scores),
      ")."
    )
  }
  if (length(beta_hat) != ncol(scores)) {
    stop(
      "Coefficient length (",
      length(beta_hat),
      ") does not match ",
      "score columns (",
      ncol(scores),
      ")."
    )
  }

  # Fix cluster ordering to first occurrence and aggregate score contributions.
  cluster_chr <- as.character(cluster)
  cluster_ids <- unique(cluster_chr)
  cluster_fac <- factor(cluster_chr, levels = cluster_ids)
  scores_clustered <- rowsum(scores, cluster_fac, reorder = FALSE)

  cluster_idx <- NULL
  n_pat <- 0L
  if (isTRUE(return_baseline_weights)) {
    baseline_ids <- as.character(baseline_data[[id_var]])
    cluster_idx <- match(baseline_ids, cluster_ids)
    if (anyNA(cluster_idx)) {
      missing_ids <- unique(baseline_ids[is.na(cluster_idx)])
      stop(
        "Some baseline IDs used in avg_sops() were not found in the cluster ",
        "variable stored in robcov_vglm: ",
        paste(utils::head(missing_ids, 5), collapse = ", "),
        if (length(missing_ids) > 5) " ..." else "",
        ". Ensure avg_sops(id_var = ...) matches the clustering used in ",
        "robcov_vglm(cluster = ...)."
      )
    }
    n_pat <- length(baseline_ids)
  }

  n_clusters <- length(cluster_ids)
  p <- length(beta_hat)

  beta_draws <- matrix(NA_real_, nrow = n_sim, ncol = p)
  colnames(beta_draws) <- names(beta_hat)
  baseline_weights <- if (isTRUE(return_baseline_weights)) {
    matrix(NA_real_, nrow = n_sim, ncol = n_pat)
  } else {
    NULL
  }

  for (i in seq_len(n_sim)) {
    # Exponential multipliers are centered at 1; centered form drives score perturbation.
    w_cluster <- stats::rexp(n_clusters, rate = 1)
    u_cluster <- w_cluster - 1

    # multiple each centered patient weight by the patient-level score
    # contribution, then sum across patients to get the overall score
    # perturbation for the draw.

    # Before applying the weights, colSums of the scores would be (by
    # definition) zero (because we're looking at the first derivative at the MLE
    #).

    # U_star is the perturbed total score vector for this draw. It
    # represents the total score perturbation under resampled patient weights.
    U_star <- as.vector(crossprod(u_cluster, scores_clustered))

    # One-Step Newton-Raphson Update:
    # If we ignore covariance: if var(beta_x) very low, then even if U_star is
    # large, then the update will be small.
    #
    # Reminder to self bread %*% U_star is equivalent to sum(bread[i, ] *
    # U_star) for each row of bread.

    # We can use Newton-Raphson to find root of the score equation, which gives
    # us the new MLE estimates for this draw.
    # We want to find where $f(x) = 0$ (our score)
    # $x_{new} = x_{old} - \frac{f(x_{old})}{f'(x_{old})}$
    # $x$ is is beta
    # $f(x)$ is the score equation, which is zero at the MLE
    # $f'(x)$ is the Hessian
    # $$\beta_{new} = \beta_{old} - \frac{U(\beta_{old})}{H(\beta_{old})}$$

    # Dividing by a matrix is the same as multiplying by the inverse, so we can
    # write:
    # $\beta_{new} = \beta_{old} - H^{-1} \times U(\beta_{old})$
    # vcov (V) is defined as the inverse of the negative Hessian ($V = (-H)^{-1}$), so we can rewrite:
    # $\beta_{new} = \beta_{old} + V \times U(\beta_{old})$
    # The new beta (x1) is where the tangent line at the old beta (x0) intersects the x-axis (where score = 0).

    # The bread component is vcov(fit), so this directly applies the
    # one-step update from the total score perturbation.
    beta_draws[i, ] <- beta_hat + as.vector(bread %*% U_star)

    if (isTRUE(return_baseline_weights)) {
      w_patient <- w_cluster[cluster_idx]
      w_sum <- sum(w_patient)
      if (!is.finite(w_sum) || w_sum <= 0) {
        stop("Invalid baseline weights generated for score bootstrap.")
      }
      baseline_weights[i, ] <- w_patient / w_sum
    }
  }

  list(beta_draws = beta_draws, baseline_weights = baseline_weights)
}

score_bootstrap_components <- function(model, cluster = NULL) {
  if (inherits(model, "robcov_vglm")) {
    required <- c("coefficients", "bread", "scores", "cluster")
    missing_required <- required[vapply(
      required,
      function(x) is.null(model[[x]]),
      logical(1)
    )]
    if (length(missing_required) > 0) {
      stop(
        "The supplied robcov_vglm model is missing required components for ",
        "score bootstrap: ",
        paste(missing_required, collapse = ", "),
        ". Refit with robcov_vglm()."
      )
    }

    return(list(
      coefficients = model$coefficients,
      bread = model$bread,
      scores = model$scores,
      cluster = model$cluster
    ))
  }

  if (!inherits(model, "orm")) {
    stop("score_bootstrap_components() supports robcov_vglm and orm models.")
  }

  if (is.null(cluster)) {
    stop(
      "`cluster` is required for orm score-bootstrap inference. ",
      "Supply the row-level patient ID vector used to fit the orm model, ",
      "for example `cluster = data$id`."
    )
  }

  scores <- compute_scores_orm(model)
  cluster <- align_cluster_orm(model, cluster, nrow(scores))
  bread <- get_orm_model_vcov(model)

  list(
    coefficients = stats::coef(model),
    bread = bread,
    scores = scores,
    cluster = cluster
  )
}

compute_scores_orm <- function(model) {
  if (!inherits(model, "orm")) {
    stop("'model' must be an orm object.")
  }
  if (is.null(model$x) || is.null(model$y)) {
    stop(
      "orm score bootstrap requires the model to be fitted with ",
      "`x = TRUE, y = TRUE`."
    )
  }

  beta <- stats::coef(model)
  Gamma <- get_effective_coefs(model, beta = beta)
  X <- cbind("(Intercept)" = 1, as.matrix(model$x))
  X <- X[, colnames(Gamma), drop = FALSE]

  eta <- X %*% t(Gamma)
  cum_probs <- stats::plogis(eta)
  probs <- lp_to_probs(eta, nrow(Gamma))
  probs <- pmax(probs, .Machine$double.eps)

  y_levels <- as_state_labels(model$yunique)
  y_idx <- match(as_state_labels(model$y), y_levels)
  if (anyNA(y_idx)) {
    stop("Could not match orm response values to fitted outcome levels.")
  }

  n <- nrow(X)
  M <- nrow(Gamma)
  d_eta <- matrix(0, nrow = n, ncol = M)

  for (i in seq_len(n)) {
    k <- y_idx[i]
    if (k == 1L) {
      d_eta[i, 1L] <- -cum_probs[i, 1L] * (1 - cum_probs[i, 1L]) / probs[i, 1L]
    } else if (k == M + 1L) {
      d_eta[i, M] <- cum_probs[i, M] * (1 - cum_probs[i, M]) / probs[i, M + 1L]
    } else {
      d_eta[i, k - 1L] <- cum_probs[i, k - 1L] *
        (1 - cum_probs[i, k - 1L]) /
        probs[i, k]
      d_eta[i, k] <- -cum_probs[i, k] *
        (1 - cum_probs[i, k]) /
        probs[i, k]
    }
  }

  intercept_scores <- d_eta
  slope_scores <- as.matrix(model$x) * rowSums(d_eta)
  scores <- cbind(intercept_scores, slope_scores)
  colnames(scores) <- names(beta)
  scores
}

align_cluster_orm <- function(model, cluster, n_scores) {
  if (inherits(cluster, "formula")) {
    stop("`cluster` for orm score bootstrap must be a row-level vector.")
  }

  if (length(cluster) == n_scores) {
    return(cluster)
  }

  na_action <- model$na.action
  if (!is.null(na_action) && length(cluster) > n_scores) {
    keep <- seq_along(cluster)
    keep <- keep[-as.integer(na_action)]
    if (length(keep) == n_scores) {
      return(cluster[keep])
    }
  }

  stop(
    "Length of `cluster` (",
    length(cluster),
    ") does not match the number of orm score rows (",
    n_scores,
    "). Supply the row-level patient ID vector used to fit the model."
  )
}

get_orm_model_vcov <- function(model) {
  beta <- stats::coef(model)

  if (!is.null(model$orig.var)) {
    return(validate_coef_vcov(beta, model$orig.var, arg = "orm model vcov"))
  }

  if (requireNamespace("rms", quietly = TRUE)) {
    robust <- rms::robcov(model)
    if (!is.null(robust$orig.var)) {
      return(validate_coef_vcov(beta, robust$orig.var, arg = "orm model vcov"))
    }
  }

  stop(
    "Could not obtain a full model-based covariance matrix for the orm fit. ",
    "Refit with rms::orm(..., x = TRUE, y = TRUE)."
  )
}
