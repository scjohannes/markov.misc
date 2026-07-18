#' Robust (Sandwich) Covariance Matrix Estimation for vglm Objects
#'
#' Computes the Huber-White sandwich covariance matrix estimator for vglm
#' objects, optionally with cluster-robust standard errors. This provides
#' functionality similar to `rms::robcov()` for vglm models.
#'
#' @param fit A fitted vglm object from the VGAM package.
#' @param cluster Optional vector of cluster identifiers. If provided, computes
#'   cluster-robust standard errors. Must have the same length as the fitted
#'   observations, or the same length as the original pre-NA data when the fitted
#'   model records omitted rows in `fit@na.action`.
#' @param adjust Deprecated compatibility alias for `cadjust`. A non-missing
#'   logical value is accepted so existing calls continue to work. Do not
#'   supply both `adjust` and `cadjust`.
#' @param bread Character string selecting the sandwich bread. `"observed"`
#'   (the default) numerically differentiates the summed analytic scores.
#'   `"vglm"` uses the inverse final VGAM IRLS working information returned by
#'   `stats::vcov(fit)`.
#' @param type Character string selecting an HC correction. `"HC0"` applies no
#'   observation degrees-of-freedom correction. `"HC1"` multiplies the meat by
#'   `(n - 1) / (n - p)`.
#' @param cadjust Logical. Apply the cluster correction `G / (G - 1)`? The
#'   default is `TRUE` when `cluster` is supplied and `FALSE` otherwise. When
#'   explicitly `TRUE` without `cluster`, each fitted row is treated as its own
#'   cluster for this scalar correction.
#'
#' @return A list with class "robcov_vglm" containing all original model
#'   information plus robust covariance estimates:
#'
#'   **Model information (copied from original vglm fit):**
#'   \item{coefficients}{Model coefficients}
#'   \item{fitted.values}{Fitted values from the model}
#'   \item{residuals}{Model residuals}
#'   \item{family}{The VGAM family object}
#'   \item{predictors}{Linear predictors (eta)}
#'   \item{prior.weights}{Prior weights}
#'   \item{df.residual}{Residual degrees of freedom}
#'   \item{df.total}{Total degrees of freedom}
#'   \item{rank}{Model rank}
#'   \item{criterion}{Convergence criterion values}
#'   \item{constraints}{Constraint matrices}
#'   \item{control}{Control parameters}
#'   \item{extra}{Extra information from family}
#'   \item{misc}{Miscellaneous model info}
#'   \item{model}{Model frame (if available)}
#'   \item{x}{Design matrix (if available)}
#'   \item{y}{Response (if available)}
#'   \item{terms}{Model terms}
#'   \item{assign}{Assignment vector}
#'   \item{xlevels}{Factor levels}
#'   \item{na.action}{NA handling}
#'   \item{call}{Original model call}
#'   \item{original_call}{Call to robcov_vglm}
#'
#'   **Robust covariance estimates:**
#'   \item{var}{The robust variance-covariance matrix}
#'   \item{se}{Robust standard errors (sqrt of diagonal)}
#'   \item{z}{Z-statistics (coefficients / robust SE)}
#'   \item{pvalues}{Two-sided p-values based on z-statistics}
#'   \item{original_vcov}{The original model-based variance-covariance matrix}
#'   \item{original_se}{Original model-based standard errors}
#'
#'   **Sandwich components:**
#'   \item{scores}{Matrix of observation-level score contributions (n x p)}
#'   \item{influence_beta_obs}{Matrix of observation-level influence
#'     contributions for \eqn{\beta} (n x p), computed as the scores times the
#'     selected bread}
#'   \item{bread}{The selected covariance-scale bread}
#'   \item{bread_type}{The requested bread type}
#'   \item{meat}{The crossproduct of raw observation scores, or cluster-summed
#'     scores when `cluster` is supplied}
#'
#'   **Clustering information:**
#'   \item{cluster}{The cluster variable (if provided)}
#'   \item{n}{Number of observations}
#'   \item{n_clusters}{Number of clusters (if clustered)}
#'   \item{type}{The HC correction type}
#'   \item{cadjust}{Whether the cluster correction was applied}
#'   \item{adjustment_factor}{The combined HC and cluster scalar correction}
#'   \item{adjust}{Compatibility alias for `cadjust`}
#'
#' @details
#' The sandwich estimator has the form:
#' \deqn{V = B \cdot M \cdot B}
#' where B is the selected covariance-scale bread and M is the score
#' crossproduct. By default, B is the inverse observed information obtained from
#' the numerical Jacobian of the summed analytic scores. This keeps the meat and
#' bread tied to the same estimating equations and provides misspecification-
#' robust inference. Set `bread = "vglm"` to use VGAM's final fitted
#' working/Fisher information instead.
#'
#' Numerical differentiation operates on VGAM's effective VLM coefficients, so
#' proportional-odds, fully nonparallel, and custom partial-proportional-odds
#' constraint matrices use the same implementation. Perturbation steps are
#' reduced adaptively when a candidate produces invalid fitted values.
#'
#' For unclustered data, the meat matrix is:
#' \deqn{M = \sum_{i=1}^{n} \psi_i \psi_i'}
#' where \eqn{\psi_i} is the score (gradient of log-likelihood) for observation i.
#'
#' For clustered data, scores are summed within clusters before computing the
#' meat matrix:
#' \deqn{M = \sum_{g=1}^{G} \left(\sum_{i \in g} \psi_i\right)
#'       \left(\sum_{i \in g} \psi_i\right)'}
#'
#' The `type` and `cadjust` corrections are separate, matching the controls used
#' by `sandwich::meatCL()`. For clustered data, `cadjust = TRUE` multiplies the
#' meat by G/(G-1). `type = "HC1"` independently multiplies it by
#' (n-1)/(n-p). Set `cadjust = FALSE` and `type = "HC0"` when matching
#' unadjusted `rms::robcov()` results.
#'
#' **Z-statistics and p-values**: The returned object includes z-statistics
#' computed as coefficients divided by robust standard errors, and two-sided
#' p-values from the standard normal distribution. This approximation is
#' asymptotic in the number of independent observations or clusters. With few
#' clusters, p-values may be anti-conservative.
#'
#' **Weights and aggregated responses**: Scores are computed at the fitted-row
#' level. With prior or frequency weights, or with aggregated ordinal count
#' responses, this treats each fitted row as the independent score contribution.
#' That can differ from a sandwich estimator computed after expanding data to
#' one row per independent individual.
#'
#' **Comparison with rms::robcov()**: For equivalent models, results should be
#' compared after accounting for parameterization, cutpoint direction, response
#' coding, weights, and missing-data handling.
#'
#' @examples
#' \dontrun{
#' library(VGAM)
#'
#' # Fit a proportional odds model
#' fit <- vglm(y ~ x1 + x2,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = mydata)
#'
#' # Compute robust (unclustered) standard errors
#' robust_vcov <- robcov_vglm(fit)
#' robust_vcov$se
#'
#' # Compute cluster-robust standard errors, with G/(G-1) correction by default
#' robust_vcov_cl <- robcov_vglm(fit, cluster = mydata$cluster_id)
#' robust_vcov_cl$se
#'
#' # Without small-sample correction
#' robust_vcov_noadj <- robcov_vglm(
#'   fit,
#'   cluster = mydata$cluster_id,
#'   cadjust = FALSE
#' )
#' }
#'
#' @seealso \code{\link{compare_se_orm_vglm}} for comparing standard errors
#'   between orm and vglm fits
#'
#' @importFrom stats vcov nobs coef pnorm fitted residuals symnum
#'
#' @export
robcov_vglm <- function(
  fit,
  cluster = NULL,
  adjust = NULL,
  bread = c("observed", "vglm"),
  type = c("HC0", "HC1"),
  cadjust = NULL
) {
  if (!inherits(fit, "vglm")) {
    stop("'fit' must be a vglm object")
  }

  bread_type <- match.arg(bread)
  type <- match.arg(type)
  correction <- resolve_vglm_corrections(
    cluster = cluster,
    adjust = adjust,
    cadjust = cadjust,
    adjust_missing = missing(adjust),
    cadjust_missing = missing(cadjust)
  )
  cadjust <- correction$cadjust

  # --- 1. Extract key quantities from the vglm object ---
  n <- nobs(fit, type = "lm")
  p <- length(coef(fit))
  coefficients <- coef(fit)

  if (
    is.null(names(coefficients)) ||
      anyNA(names(coefficients)) ||
      any(names(coefficients) == "") ||
      anyDuplicated(names(coefficients))
  ) {
    stop("vglm coefficients must have unique, non-missing names.")
  }
  if (any(!is.finite(coefficients))) {
    stop("vglm coefficients must all be finite.")
  }

  # Use VGAM's covariance-scale bread. For cumulative proportional-odds
  # models, VGAM's vcov() is based on final IRLS working/Fisher information,
  # which can differ slightly from observed-Hessian bread used by sandwich
  # methods for ordinal::clm() and MASS::polr().
  model_vcov <- vcov(fit)
  validate_vglm_matrix(
    model_vcov,
    names(coefficients),
    "stats::vcov(fit)",
    positive_definite = TRUE
  )
  original_se <- sqrt(diag(model_vcov))
  names(original_se) <- names(coefficients)

  # --- 2. Compute observation-level score contributions ---
  scores <- compute_scores_vglm(fit)

  # Validate scores dimensions

  if (nrow(scores) != n) {
    stop("Number of score rows does not match number of observations")
  }
  if (ncol(scores) != p) {
    stop("Number of score columns does not match number of parameters")
  }
  if (!identical(colnames(scores), names(coefficients))) {
    stop("Score columns must exactly match vglm coefficient names and order.")
  }
  if (any(!is.finite(scores))) {
    stop("Observation-level vglm scores must all be finite.")
  }

  convergence <- validate_vglm_convergence(fit, scores, model_vcov)

  bread_result <- if (bread_type == "observed") {
    compute_observed_bread_vglm(fit)
  } else {
    list(
      bread = model_vcov,
      information = solve(model_vcov),
      jacobian_asymmetry = NA_real_,
      min_information_eigenvalue = min(
        eigen(
          solve(model_vcov),
          symmetric = TRUE,
          only.values = TRUE
        )$values
      )
    )
  }
  bread <- bread_result$bread
  validate_vglm_matrix(
    bread,
    names(coefficients),
    paste0("selected `", bread_type, "` bread"),
    positive_definite = TRUE
  )
  influence_beta_obs <- scores %*% bread
  dimnames(influence_beta_obs) <- list(
    rownames(scores),
    names(coefficients)
  )

  # --- 3. Compute the meat matrix ---
  if (!is.null(cluster)) {
    cluster <- align_cluster_vglm(fit, cluster, n)

    if (anyNA(cluster)) {
      stop(
        "'cluster' contains missing values after alignment with the fitted ",
        "data."
      )
    }

    # Convert to factor and aggregate scores by cluster
    cluster_group <- as.factor(cluster)

    # Sum scores within clusters using efficient rowsum
    scores_clustered <- rowsum(scores, cluster_group, reorder = FALSE)

    # Count actual clusters present in data (not just factor levels)
    n_clusters <- nrow(scores_clustered)

    if (n_clusters < 2L) {
      stop("Cluster-robust covariance requires at least two clusters.")
    }
    if (n_clusters < 30L) {
      warning(
        "Cluster-robust p-values use normal z-tests; with fewer than 30 ",
        "clusters they may be anti-conservative.",
        call. = FALSE
      )
    }

    # Compute meat matrix from clustered scores
    sandwich_scores <- scores_clustered
    meat <- crossprod(scores_clustered)
  } else {
    n_clusters <- NULL
    # Standard (unclustered) meat matrix
    sandwich_scores <- scores
    meat <- crossprod(scores)
  }

  adjustment_factor <- 1
  if (type == "HC1") {
    if (n <= p) {
      stop("`type = \"HC1\"` requires more observations than parameters.")
    }
    adjustment_factor <- adjustment_factor * ((n - 1) / (n - p))
  }
  if (cadjust) {
    adjustment_clusters <- if (is.null(n_clusters)) n else n_clusters
    if (adjustment_clusters < 2L) {
      stop("`cadjust = TRUE` requires at least two independent clusters.")
    }
    adjustment_factor <- adjustment_factor *
      (adjustment_clusters / (adjustment_clusters - 1))
  }
  meat <- meat * adjustment_factor
  dimnames(meat) <- list(names(coefficients), names(coefficients))
  validate_vglm_matrix(
    meat,
    names(coefficients),
    "sandwich meat",
    positive_semidefinite = TRUE
  )

  # --- 4. Compute the sandwich: V = B * M * B ---
  # The crossproduct form is algebraically equivalent for symmetric bread and
  # remains numerically positive semidefinite for ill-conditioned fits.
  transformed_scores <- sandwich_scores %*% bread
  robust_var <- crossprod(transformed_scores) * adjustment_factor

  # Ensure symmetry (numerical stability)
  robust_var <- (robust_var + t(robust_var)) / 2
  dimnames(robust_var) <- dimnames(model_vcov)
  validate_vglm_matrix(
    robust_var,
    names(coefficients),
    "robust covariance",
    positive_semidefinite = TRUE
  )

  # Extract standard errors
  robust_se <- sqrt(pmax(diag(robust_var), 0))
  names(robust_se) <- names(coefficients)

  # --- 5. Compute z-statistics and p-values ---
  z_stats <- coefficients / robust_se
  pvalues <- 2 * pnorm(-abs(z_stats))
  names(z_stats) <- names(coefficients)
  names(pvalues) <- names(coefficients)

  # --- 6. Extract model information from vglm object ---
  # Use slot() for S4 objects, with tryCatch for optional slots
  safe_slot <- function(obj, name) {
    tryCatch(methods::slot(obj, name), error = function(e) NULL)
  }

  # --- 7. Build result object with all model information ---
  result <- list(
    # Model coefficients and inference
    coefficients = coefficients,
    var = robust_var,
    se = robust_se,
    z = z_stats,
    pvalues = pvalues,

    # Original (model-based) variance
    original_vcov = model_vcov,
    original_se = original_se,

    # Sandwich components
    scores = scores,
    influence_beta_obs = influence_beta_obs,
    bread = bread,
    bread_type = bread_type,
    bread_diagnostics = list(
      jacobian_asymmetry = bread_result$jacobian_asymmetry,
      min_information_eigenvalue = bread_result$min_information_eigenvalue
    ),
    meat = meat,

    # Clustering information
    cluster = cluster,
    n = n,
    n_clusters = n_clusters,
    type = type,
    cadjust = cadjust,
    adjustment_factor = adjustment_factor,
    adjust = cadjust,
    convergence = convergence,

    # Original vglm fit (for prediction)
    vglm_fit = fit,

    # Model information from vglm
    fitted.values = fit@fitted.values,
    residuals = fit@residuals,
    family = safe_slot(fit, "family"),
    predictors = safe_slot(fit, "predictors"),
    prior.weights = safe_slot(fit, "prior.weights"),
    df.residual = safe_slot(fit, "df.residual"),
    df.total = safe_slot(fit, "df.total"),
    rank = safe_slot(fit, "rank"),
    criterion = safe_slot(fit, "criterion"),
    constraints = safe_slot(fit, "constraints"),
    control = safe_slot(fit, "control"),
    extra = safe_slot(fit, "extra"),
    misc = safe_slot(fit, "misc"),
    model = safe_slot(fit, "model"),
    x = safe_slot(fit, "x"),
    y = safe_slot(fit, "y"),
    terms = safe_slot(fit, "terms"),
    assign = safe_slot(fit, "assign"),
    xlevels = safe_slot(fit, "xlevels"),
    na.action = safe_slot(fit, "na.action"),
    call = safe_slot(fit, "call"),
    original_call = match.call()
  )

  class(result) <- "robcov_vglm"
  result <- markov_attach_model_data(
    result,
    markov_model_data(fit),
    markov_model_id_var(fit)
  )
  attr(result, "markov_vglm") <- attr(fit, "markov_vglm", exact = TRUE)
  attr(result, "markov_split_assign") <- attr(
    fit,
    "markov_split_assign",
    exact = TRUE
  )
  attr(result, "markov_basis_terms") <- attr(
    fit,
    "markov_basis_terms",
    exact = TRUE
  )
  return(result)
}


#' Compute Observation-Level Score Contributions for vglm
#'
#' Internal function to compute the score (gradient of log-likelihood)
#' contributions for each observation in a vglm model.
#'
#' @param fit A fitted vglm object.
#'
#' @return An n x p matrix of score contributions, where n is the number of
#'   observations and p is the number of parameters.
#'
#' @details
#' The score for observation i is computed as:
#' \deqn{\psi_i = X_i' \cdot \frac{\partial \ell_i}{\partial \eta}}
#' where \eqn{X_i} is the portion of the VLM model matrix for observation i
#' (an M x p submatrix where M is the number of linear predictors), and
#' \eqn{\frac{\partial \ell_i}{\partial \eta}} is the derivative of the
#' log-likelihood with respect to the linear predictors for observation i.
#'
#' The VLM (vector linear model) structure in VGAM means that for an ordinal
#' model with M cutpoints, each observation contributes M rows to the design
#' matrix. The score for each observation is computed by multiplying the
#' M x p design matrix portion by the M x 1 derivative vector.
#'
#' @keywords internal
compute_scores_vglm <- function(fit) {
  if (!inherits(fit, "vglm")) {
    stop("'fit' must be a vglm object")
  }

  # Get dimensions
  n <- nobs(fit, type = "lm")
  M <- VGAM::npred(fit)
  p <- length(coef(fit))

  # Get working weights with derivatives
  # wt_info$deriv is an n x M matrix containing d(log L)/d(eta_m) for each obs
  wt_info <- VGAM::weights(fit, type = "working", deriv.arg = TRUE)
  if (!is.list(wt_info) || is.null(wt_info$deriv)) {
    stop(
      "Could not extract d logLik / d eta from ",
      "VGAM::weights(..., deriv.arg = TRUE)."
    )
  }
  deriv_eta <- wt_info$deriv

  if (!is.matrix(deriv_eta)) {
    deriv_eta <- matrix(deriv_eta, nrow = n)
  }
  if (nrow(deriv_eta) != n || ncol(deriv_eta) != M) {
    stop(
      "Unexpected derivative dimensions: got [",
      nrow(deriv_eta),
      " x ",
      ncol(deriv_eta),
      "], expected [",
      n,
      " x ",
      M,
      "]."
    )
  }

  # Get the VLM model matrix (n*M x p)
  X_vlm <- stats::model.matrix(fit, type = "vlm")

  if (nrow(X_vlm) != n * M || ncol(X_vlm) != p) {
    stop(
      "Unexpected VLM model matrix dimensions: got [",
      nrow(X_vlm),
      " x ",
      ncol(X_vlm),
      "], expected [",
      n * M,
      " x ",
      p,
      "]."
    )
  }
  if (!identical(colnames(X_vlm), names(coef(fit)))) {
    stop(
      "VLM model-matrix columns must exactly match vglm coefficient names ",
      "and order."
    )
  }

  deriv_vec <- as.vector(t(deriv_eta))
  obs_index <- rep(seq_len(n), each = M)
  scores <- rowsum(X_vlm * deriv_vec, group = obs_index, reorder = FALSE)
  scores <- as.matrix(scores)

  colnames(scores) <- names(coef(fit))
  lm_rownames <- tryCatch(
    rownames(stats::model.matrix(fit, type = "lm")),
    error = function(e) NULL
  )
  if (length(lm_rownames) == n) {
    rownames(scores) <- lm_rownames
  }
  return(scores)
}


resolve_vglm_corrections <- function(
  cluster,
  adjust,
  cadjust,
  adjust_missing,
  cadjust_missing
) {
  has_adjust <- !adjust_missing && !is.null(adjust)
  has_cadjust <- !cadjust_missing && !is.null(cadjust)

  if (has_adjust && has_cadjust) {
    stop("Supply only one of legacy `adjust` and `cadjust`, not both.")
  }

  validate_flag <- function(x, arg) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop("`", arg, "` must be TRUE or FALSE.")
    }
    x
  }

  if (has_adjust) {
    adjust <- validate_flag(adjust, "adjust")
    # The legacy argument never adjusted the unclustered covariance.
    return(list(cadjust = adjust && !is.null(cluster)))
  }
  if (has_cadjust) {
    return(list(cadjust = validate_flag(cadjust, "cadjust")))
  }

  list(cadjust = !is.null(cluster))
}


validate_vglm_matrix <- function(
  x,
  coefficient_names,
  label,
  positive_definite = FALSE,
  positive_semidefinite = FALSE
) {
  p <- length(coefficient_names)
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(p, p))) {
    stop("`", label, "` must be a numeric [", p, " x ", p, "] matrix.")
  }
  if (
    !identical(rownames(x), coefficient_names) ||
      !identical(colnames(x), coefficient_names)
  ) {
    stop(
      "Rows and columns of `",
      label,
      "` must exactly match vglm coefficient names and order."
    )
  }
  if (any(!is.finite(x))) {
    stop("`", label, "` contains non-finite values.")
  }

  scale <- max(1, max(abs(x)))
  symmetry_error <- max(abs(x - t(x))) / scale
  if (!is.finite(symmetry_error) || symmetry_error > 1e-8) {
    stop("`", label, "` is not numerically symmetric.")
  }

  if (positive_definite || positive_semidefinite) {
    eigenvalues <- eigen(
      (x + t(x)) / 2,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    eigen_scale <- max(1, max(abs(eigenvalues))) * p
    positive_definite_tolerance <- 100 * .Machine$double.eps * eigen_scale
    positive_semidefinite_tolerance <- sqrt(.Machine$double.eps) * eigen_scale
    if (
      positive_definite &&
        min(eigenvalues) <= positive_definite_tolerance
    ) {
      stop("`", label, "` is not numerically positive definite.")
    }
    if (
      positive_semidefinite &&
        min(eigenvalues) < -positive_semidefinite_tolerance
    ) {
      stop("`", label, "` is not numerically positive semidefinite.")
    }
  }

  invisible(x)
}


validate_vglm_convergence <- function(fit, scores, model_vcov) {
  iter <- tryCatch(fit@iter, error = function(e) NA_integer_)
  maxit <- tryCatch(fit@control$maxit, error = function(e) NA_integer_)
  if (
    length(iter) == 1L &&
      length(maxit) == 1L &&
      is.finite(iter) &&
      is.finite(maxit) &&
      iter >= maxit
  ) {
    stop(
      "The vglm fit reached its IRLS iteration limit without convergence. ",
      "Refit the model before computing a robust covariance."
    )
  }

  regularity <- tryCatch(fit@misc$RegCondOK, error = function(e) NULL)
  if (
    is.logical(regularity) && length(regularity) == 1L && !is.na(regularity)
  ) {
    if (!regularity) {
      stop(
        "VGAM reports that maximum-likelihood regularity conditions were ",
        "violated; robust covariance estimation is not reliable."
      )
    }
  }

  model_information <- solve(model_vcov)
  information_scale <- sqrt(pmax(diag(model_information), .Machine$double.eps))
  standardized_score <- colSums(scores) / information_scale
  score_norm <- max(abs(standardized_score))
  epsilon <- tryCatch(fit@control$epsilon, error = function(e) NULL)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon)) {
    epsilon <- 1e-7
  }
  score_tolerance <- max(1e-3, sqrt(epsilon))
  if (!is.finite(score_norm) || score_norm > score_tolerance) {
    stop(
      "The standardized norm of the summed vglm scores is ",
      format(score_norm, digits = 4),
      ", exceeding the convergence tolerance ",
      format(score_tolerance, digits = 4),
      ". Refit the model before computing a robust covariance."
    )
  }

  list(
    iter = iter,
    maxit = maxit,
    standardized_score_norm = score_norm,
    score_tolerance = score_tolerance,
    regularity_ok = if (is.null(regularity)) NA else regularity
  )
}


vglm_score_evaluator <- function(fit) {
  n <- stats::nobs(fit, type = "lm")
  M <- VGAM::npred(fit)
  beta_hat <- stats::coef(fit)
  X_vlm <- stats::model.matrix(fit, type = "vlm")

  if (!identical(colnames(X_vlm), names(beta_hat))) {
    stop(
      "VLM model-matrix columns must exactly match vglm coefficient names ",
      "and order."
    )
  }

  fitted_linear_predictors <- matrix(
    drop(X_vlm %*% beta_hat),
    nrow = n,
    ncol = M,
    byrow = TRUE
  )
  predictor_offset <- fit@predictors - fitted_linear_predictors
  if (any(!is.finite(predictor_offset))) {
    stop("The fitted vglm predictor offset contains non-finite values.")
  }

  evaluate <- function(beta) {
    candidate <- fit
    eta <- matrix(
      drop(X_vlm %*% beta),
      nrow = n,
      ncol = M,
      byrow = TRUE
    ) +
      predictor_offset

    mu <- tryCatch(
      candidate@family@linkinv(eta = eta, extra = candidate@extra),
      error = function(e) NULL
    )
    if (is.null(mu) || any(!is.finite(mu))) {
      return(NULL)
    }

    valid_params <- candidate@family@validparams
    if (length(body(valid_params))) {
      params_ok <- tryCatch(
        isTRUE(valid_params(
          eta = eta,
          y = candidate@y,
          extra = candidate@extra
        )),
        error = function(e) FALSE
      )
      if (!params_ok) {
        return(NULL)
      }
    }

    valid_fitted <- candidate@family@validfitted
    if (length(body(valid_fitted))) {
      fitted_ok <- tryCatch(
        isTRUE(valid_fitted(mu = mu, y = candidate@y, extra = candidate@extra)),
        error = function(e) FALSE
      )
      if (!fitted_ok) {
        return(NULL)
      }
    }

    fitted_values <- mu
    fitted_template <- candidate@fitted.values
    if (is.null(dim(fitted_values))) {
      if (length(fitted_values) != length(fitted_template)) {
        return(NULL)
      }
      fitted_values <- matrix(
        fitted_values,
        nrow = nrow(fitted_template),
        ncol = ncol(fitted_template)
      )
    }
    if (!identical(dim(fitted_values), dim(fitted_template))) {
      return(NULL)
    }
    dimnames(fitted_values) <- dimnames(fitted_template)

    candidate@coefficients <- beta
    candidate@predictors <- eta
    candidate@fitted.values <- fitted_values
    derivatives <- tryCatch(
      VGAM::weights(
        candidate,
        type = "working",
        deriv.arg = TRUE
      )$deriv,
      error = function(e) NULL
    )
    if (is.null(derivatives)) {
      return(NULL)
    }
    if (!is.matrix(derivatives)) {
      derivatives <- matrix(derivatives, nrow = n)
    }
    if (
      !identical(dim(derivatives), c(n, M)) ||
        any(!is.finite(derivatives))
    ) {
      return(NULL)
    }

    scores <- rowsum(
      X_vlm * as.vector(t(derivatives)),
      group = rep(seq_len(n), each = M),
      reorder = FALSE
    )
    score_sum <- colSums(scores)
    names(score_sum) <- names(beta_hat)
    score_sum
  }

  list(evaluate = evaluate, beta = beta_hat)
}


compute_observed_bread_vglm <- function(fit) {
  evaluator <- vglm_score_evaluator(fit)
  beta <- evaluator$beta
  p <- length(beta)
  jacobian <- matrix(
    NA_real_,
    nrow = p,
    ncol = p,
    dimnames = list(names(beta), names(beta))
  )

  initial_steps <- .Machine$double.eps^(1 / 3) * pmax(1, abs(beta))
  max_step_reductions <- 12L

  for (j in seq_len(p)) {
    derivative <- NULL
    for (reduction in 0:max_step_reductions) {
      step <- initial_steps[j] / (2^reduction)
      beta_plus <- beta_minus <- beta
      beta_plus_half <- beta_minus_half <- beta
      beta_plus[j] <- beta_plus[j] + step
      beta_minus[j] <- beta_minus[j] - step
      beta_plus_half[j] <- beta_plus_half[j] + step / 2
      beta_minus_half[j] <- beta_minus_half[j] - step / 2

      score_plus <- evaluator$evaluate(beta_plus)
      score_minus <- evaluator$evaluate(beta_minus)
      score_plus_half <- evaluator$evaluate(beta_plus_half)
      score_minus_half <- evaluator$evaluate(beta_minus_half)
      candidates <- list(
        score_plus,
        score_minus,
        score_plus_half,
        score_minus_half
      )
      valid <- vapply(
        candidates,
        function(x) length(x) == p && all(is.finite(x)),
        logical(1)
      )
      if (all(valid)) {
        coarse <- (score_plus - score_minus) / (2 * step)
        fine <- (score_plus_half - score_minus_half) / step
        derivative <- (4 * fine - coarse) / 3
        break
      }
    }

    if (is.null(derivative) || any(!is.finite(derivative))) {
      stop(
        "Could not compute a valid central score derivative for coefficient `",
        names(beta)[j],
        "` after adaptive step reduction. The fitted model may be too close ",
        "to an invalid parameter boundary. Use `bread = \"vglm\"` only if ",
        "working-information inference is acceptable."
      )
    }
    jacobian[, j] <- derivative
  }

  jacobian_scale <- max(1, max(abs(jacobian)))
  jacobian_asymmetry <- max(abs(jacobian - t(jacobian))) / jacobian_scale
  if (!is.finite(jacobian_asymmetry) || jacobian_asymmetry > 1e-5) {
    stop(
      "The numerical score Jacobian is materially asymmetric (relative ",
      "asymmetry ",
      format(jacobian_asymmetry, digits = 4),
      "). Use `bread = \"vglm\"` only if working-information inference is ",
      "acceptable."
    )
  }

  information <- -(jacobian + t(jacobian)) / 2
  dimnames(information) <- list(names(beta), names(beta))
  information_eigenvalues <- eigen(
    information,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  information_tolerance <- 100 *
    .Machine$double.eps *
    max(1, max(abs(information_eigenvalues))) *
    p
  if (
    any(!is.finite(information_eigenvalues)) ||
      min(information_eigenvalues) <= information_tolerance
  ) {
    stop(
      "The observed vglm information is singular or not positive definite. ",
      "Use `bread = \"vglm\"` only if working-information inference is ",
      "acceptable."
    )
  }

  bread <- chol2inv(chol(information))
  dimnames(bread) <- list(names(beta), names(beta))
  list(
    bread = bread,
    information = information,
    jacobian_asymmetry = jacobian_asymmetry,
    min_information_eigenvalue = min(information_eigenvalues)
  )
}


align_cluster_vglm <- function(fit, cluster, n) {
  if (length(cluster) == n) {
    return(cluster)
  }

  na_action <- tryCatch(
    methods::slot(fit, "na.action"),
    error = function(e) NULL
  )
  if (is.list(na_action)) {
    na_action <- unlist(na_action, use.names = FALSE)
  }
  na_action <- as.integer(na_action)
  na_action <- na_action[!is.na(na_action)]

  if (
    length(na_action) > 0L &&
      max(na_action) <= length(cluster) &&
      length(cluster[-na_action]) == n
  ) {
    return(cluster[-na_action])
  }

  stop(
    "Length of 'cluster' must equal nobs(fit, type = \"lm\") = ",
    n,
    ", or must be the original pre-NA vector from which fit@na.action ",
    "can be removed."
  )
}


#' Summary Method for robcov_vglm Objects
#'
#' Prints a comprehensive summary of the robust covariance estimation results,
#' including coefficients, robust standard errors, z-statistics, and p-values.
#'
#' @param object A robcov_vglm object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns a list containing the coefficient table and
#'   summary statistics.
#'
#' @examplesIf rlang::is_installed("VGAM")
#' trial <- sim_actt2_brownian(n_patients = 40, follow_up_time = 6, seed = 1)
#' markov_data <- prepare_markov_data(trial, absorbing_state = 8)
#' fit <- vglm_markov(
#'   ordered(y) ~ time + tx + yprev,
#'   family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
#'   data = markov_data
#' )
#' robust_fit <- robcov_vglm(fit, cluster = markov_data$id)
#' summary(robust_fit)
#'
#' @method summary robcov_vglm
#' @export
summary.robcov_vglm <- function(object, ...) {
  cat("Robust (Sandwich) Covariance Estimation for vglm\n")
  cat("================================================\n\n")

  # Model call
  if (!is.null(object$call)) {
    cat("Call:\n")
    print(object$call)
    cat("\n")
  }

  # Sample information
  cat("Number of observations:", object$n, "\n")
  cat("Bread:", object$bread_type, "\n")
  cat("HC type:", object$type, "\n")
  if (!is.null(object$n_clusters)) {
    cat("Number of clusters:", object$n_clusters, "\n")
  }
  if (isTRUE(object$cadjust)) {
    cat("Cluster adjustment: applied (G/(G-1))\n")
  }
  if (!isTRUE(all.equal(object$adjustment_factor, 1))) {
    cat(
      "Combined adjustment factor:",
      format(object$adjustment_factor, digits = 6),
      "\n"
    )
  }
  cat("\n")

  # Build coefficient table
  coef_table <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = object$se,
    `z value` = object$z,
    `Pr(>|z|)` = object$pvalues,
    check.names = FALSE
  )

  # Add significance stars
  signif_stars <- symnum(
    object$pvalues,
    corr = FALSE,
    na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  coef_table$` ` <- as.character(signif_stars)

  cat("Coefficients (Robust SE):\n")
  print(coef_table, digits = 4)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")

  # Model fit statistics if available
  if (!is.null(object$criterion)) {
    cat("Model fit:\n")
    for (nm in names(object$criterion)) {
      cat(sprintf("  %s: %s\n", nm, format(object$criterion[[nm]], digits = 4)))
    }
    cat("\n")
  }

  # Return coefficient table invisibly
  invisible(list(
    coefficients = coef_table,
    n = object$n,
    n_clusters = object$n_clusters,
    bread_type = object$bread_type,
    type = object$type,
    cadjust = object$cadjust,
    adjustment_factor = object$adjustment_factor
  ))
}


#' Print Method for robcov_vglm Objects
#'
#' Prints a concise summary of the robust covariance estimation results.
#'
#' @param x A robcov_vglm object.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the input object.
#'
#' @examplesIf rlang::is_installed("VGAM")
#' trial <- sim_actt2_brownian(n_patients = 40, follow_up_time = 6, seed = 1)
#' markov_data <- prepare_markov_data(trial, absorbing_state = 8)
#' fit <- vglm_markov(
#'   ordered(y) ~ time + tx + yprev,
#'   family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
#'   data = markov_data
#' )
#' robust_fit <- robcov_vglm(fit, cluster = markov_data$id)
#' print(robust_fit)
#'
#' @method print robcov_vglm
#' @export
print.robcov_vglm <- function(x, ...) {
  cat("Robust (Sandwich) Covariance Estimation for vglm\n")
  cat("================================================\n\n")

  # Model call
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }

  # Sample information
  cat("Number of observations:", x$n, "\n")
  cat("Bread:", x$bread_type, "\n")
  cat("HC type:", x$type, "\n")
  if (!is.null(x$n_clusters)) {
    cat("Number of clusters:", x$n_clusters, "\n")
  }
  cat("\n")

  # Coefficients with robust SE
  cat("Coefficients:\n")
  coef_display <- rbind(
    Estimate = x$coefficients,
    `Robust SE` = x$se
  )
  print(coef_display, digits = 4)
  cat("\n")

  cat("(Use summary() for full coefficient table with z-values and p-values)\n")

  invisible(x)
}


#' Extract Coefficients from robcov_vglm Objects
#'
#' @param object A robcov_vglm object.
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector of coefficients.
#'
#' @examplesIf rlang::is_installed("VGAM")
#' trial <- sim_actt2_brownian(n_patients = 40, follow_up_time = 6, seed = 1)
#' markov_data <- prepare_markov_data(trial, absorbing_state = 8)
#' fit <- vglm_markov(
#'   ordered(y) ~ time + tx + yprev,
#'   family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
#'   data = markov_data
#' )
#' robust_fit <- robcov_vglm(fit, cluster = markov_data$id)
#' coef(robust_fit)
#'
#' @method coef robcov_vglm
#' @export
coef.robcov_vglm <- function(object, ...) {
  object$coefficients
}


#' Extract Variance-Covariance Matrix from robcov_vglm Objects
#'
#' Returns the robust (sandwich) variance-covariance matrix.
#'
#' @param object A robcov_vglm object.
#' @param ... Additional arguments (ignored).
#'
#' @return The robust variance-covariance matrix.
#'
#' @examplesIf rlang::is_installed("VGAM")
#' trial <- sim_actt2_brownian(n_patients = 40, follow_up_time = 6, seed = 1)
#' markov_data <- prepare_markov_data(trial, absorbing_state = 8)
#' fit <- vglm_markov(
#'   ordered(y) ~ time + tx + yprev,
#'   family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
#'   data = markov_data
#' )
#' robust_fit <- robcov_vglm(fit, cluster = markov_data$id)
#' vcov(robust_fit)
#'
#' @method vcov robcov_vglm
#' @export
vcov.robcov_vglm <- function(object, ...) {
  object$var
}


#' Compare Standard Errors from orm and vglm
#'
#' Utility function to compare standard errors between `rms::orm()` and
#' `VGAM::vglm()` fits, including robust (sandwich) versions. This is useful
#' for validating that both approaches give similar results.
#'
#' @param orm_fit A fitted orm object from the rms package.
#' @param vglm_fit A fitted vglm object from the VGAM package. Should be fit
#'   with `cumulative(parallel = TRUE, reverse = TRUE)` to match orm's
#'   parameterization for ordinal outcomes, or `binomialff` for binary outcomes.
#' @param cluster Optional cluster variable for robust standard errors. Must
#'   have the same length as the number of observations in both models.
#'
#' @return A data frame comparing standard errors with columns:
#'   \item{parameter}{Parameter name (intercepts numbered, then regression
#'     coefficients)}
#'   \item{se_orm_model}{Model-based SE from orm}
#'   \item{se_vglm_model}{Model-based SE from vglm}
#'   \item{se_orm_robust}{Robust SE from orm}
#'   \item{se_vglm_robust}{Robust SE from vglm}
#'   \item{ratio_model}{Ratio of vglm to orm model-based SEs}
#'   \item{ratio_robust}{Ratio of vglm to orm robust SEs}
#'
#' @details
#' This function aligns parameters between the two model types:
#' - orm uses names like "y>=2", "y>=3", etc. for intercepts (or "Intercept"
#'   for binary outcomes)
#' - vglm uses names like "(Intercept):1", "(Intercept):2", etc.
#'
#' The comparison focuses on:
#' 1. Model-based standard errors (from the inverse Hessian)
#' 2. Robust (sandwich) standard errors
#'
#' Ratios close to 1.0 indicate good agreement. For equivalent cumulative-logit
#' likelihoods, the default observed-score robust standard errors should agree
#' closely after parameter alignment. Native model-based standard errors can
#' still differ because VGAM reports final IRLS working-information covariance.
#'
#' @examples
#' \dontrun{
#' library(rms)
#' library(VGAM)
#'
#' # Fit models
#' dd <- datadist(mydata)
#' options(datadist = "dd")
#' m_orm <- orm(y ~ x1 + x2, data = mydata, x = TRUE, y = TRUE)
#' m_vglm <- vglm(y ~ x1 + x2,
#'                family = cumulative(parallel = TRUE, reverse = TRUE),
#'                data = mydata)
#'
#' # Compare standard errors
#' comparison <- compare_se_orm_vglm(m_orm, m_vglm)
#' print(comparison)
#'
#' # With clustering
#' comparison_cl <- compare_se_orm_vglm(m_orm, m_vglm, cluster = mydata$id)
#' }
#'
#' @importFrom stats vcov
#'
#' @export
compare_se_orm_vglm <- function(orm_fit, vglm_fit, cluster = NULL) {
  # Model-based SEs
  se_orm_model <- sqrt(diag(rms::robcov(orm_fit)$orig.var))
  se_vglm_model <- sqrt(diag(vcov(vglm_fit)))

  # Robust SEs
  if (!is.null(cluster)) {
    orm_robust <- rms::robcov(orm_fit, cluster = cluster)
    se_orm_robust <- sqrt(diag(orm_robust$var))

    vglm_robust <- robcov_vglm(
      vglm_fit,
      cluster = cluster,
      type = "HC0",
      cadjust = FALSE
    )
    se_vglm_robust <- vglm_robust$se
  } else {
    orm_robust <- rms::robcov(orm_fit)
    se_orm_robust <- sqrt(diag(orm_robust$var))

    vglm_robust <- robcov_vglm(vglm_fit)
    se_vglm_robust <- vglm_robust$se
  }

  # Match parameter names (they differ between packages)
  # orm uses names like "y>=2" for ordinal or "Intercept" for binary
  # vglm uses "(Intercept):1" for ordinal or "(Intercept)" for binary

  # Count intercepts - handle both ordinal (y>=) and binary (Intercept) cases
  intercept_pattern_orm <- "^y>=|^Intercept$"
  intercept_pattern_vglm <- "\\(Intercept\\)"

  n_intercepts_orm <- length(grep(intercept_pattern_orm, names(se_orm_model)))
  n_intercepts_vglm <- length(grep(
    intercept_pattern_vglm,
    names(se_vglm_model)
  ))

  # Ensure same number of intercepts (sanity check)
  if (n_intercepts_orm != n_intercepts_vglm) {
    warning(
      "Number of intercepts differs between models: orm has ",
      n_intercepts_orm,
      ", vglm has ",
      n_intercepts_vglm
    )
  }

  # Extract regression coefficient indices
  reg_idx_orm <- grep(intercept_pattern_orm, names(se_orm_model), invert = TRUE)
  reg_idx_vglm <- grep(
    intercept_pattern_vglm,
    names(se_vglm_model),
    invert = TRUE
  )

  # Use the minimum number of intercepts to avoid mismatches
  n_intercepts <- min(n_intercepts_orm, n_intercepts_vglm)

  # Build comparison data frame
  comparison <- data.frame(
    parameter = c(
      paste0("Intercept_", seq_len(n_intercepts)),
      names(se_orm_model)[reg_idx_orm]
    ),
    se_orm_model = c(
      se_orm_model[seq_len(n_intercepts)],
      se_orm_model[reg_idx_orm]
    ),
    se_vglm_model = c(
      se_vglm_model[seq_len(n_intercepts)],
      se_vglm_model[reg_idx_vglm]
    ),
    se_orm_robust = c(
      se_orm_robust[seq_len(n_intercepts)],
      se_orm_robust[reg_idx_orm]
    ),
    se_vglm_robust = c(
      se_vglm_robust[seq_len(n_intercepts)],
      se_vglm_robust[reg_idx_vglm]
    )
  )

  comparison$ratio_model <- comparison$se_vglm_model / comparison$se_orm_model
  comparison$ratio_robust <- comparison$se_vglm_robust /
    comparison$se_orm_robust

  return(comparison)
}
