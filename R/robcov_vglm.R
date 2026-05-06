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
#' @param adjust Logical. If TRUE, applies small-sample correction factor
#'   G/(G-1) for clustered data, where G is the number of clusters. Default is
#'   FALSE to match `rms::robcov()` behavior.
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
#'   \item{influence_beta_obs}{Optional matrix of observation-level influence
#'     contributions for \eqn{\beta} (n x p), computed via `VGAM::Influence()`
#'     when available}
#'   \item{bread}{The model-based covariance matrix, `vcov(fit)`}
#'   \item{meat}{The crossproduct of raw observation scores, or cluster-summed
#'     scores when `cluster` is supplied}
#'
#'   **Clustering information:**
#'   \item{cluster}{The cluster variable (if provided)}
#'   \item{n}{Number of observations}
#'   \item{n_clusters}{Number of clusters (if clustered)}
#'   \item{adjust}{Whether small-sample adjustment was applied}
#'
#' @details
#' The sandwich estimator has the form:
#' \deqn{V = B \cdot M \cdot B}
#' where B is the model-based covariance matrix returned by `vcov(fit)`, based
#' on VGAM's final fitted information matrix, and M is the score crossproduct.
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
#' When `adjust = TRUE`, the meat matrix is multiplied by G/(G-1) for clustered
#' data, providing a small-sample bias correction. This is sometimes recommended
#' but is not applied by `rms::robcov()`, so the default is FALSE for
#' consistency.
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
#' # Compute cluster-robust standard errors
#' robust_vcov_cl <- robcov_vglm(fit, cluster = mydata$cluster_id)
#' robust_vcov_cl$se
#'
#' # With small-sample correction
#' robust_vcov_adj <- robcov_vglm(fit, cluster = mydata$cluster_id, adjust = TRUE)
#' }
#'
#' @seealso \code{\link{compare_se_orm_vglm}} for comparing standard errors
#'   between orm and vglm fits
#'
#' @importFrom stats vcov nobs coef pnorm fitted residuals symnum
#'
#' @export
robcov_vglm <- function(fit, cluster = NULL, adjust = FALSE) {
  if (!inherits(fit, "vglm")) {
    stop("'fit' must be a vglm object")
  }

  # --- 1. Extract key quantities from the vglm object ---
  n <- nobs(fit, type = "lm")
  p <- length(coef(fit))
  coefficients <- coef(fit)

  # Get the bread matrix: VGAM's model-based covariance matrix.
  model_vcov <- vcov(fit)
  original_se <- sqrt(diag(model_vcov))
  names(original_se) <- names(coefficients)
  bread <- model_vcov

  # --- 2. Compute observation-level score contributions ---
  scores <- compute_scores_vglm(fit)
  influence_beta_obs <- tryCatch(
    VGAM::Influence(fit),
    error = function(e) NULL
  )

  # Validate scores dimensions

  if (nrow(scores) != n) {
    stop("Number of score rows does not match number of observations")
  }
  if (ncol(scores) != p) {
    stop("Number of score columns does not match number of parameters")
  }
  if (!is.null(influence_beta_obs)) {
    if (
      !is.matrix(influence_beta_obs) ||
        nrow(influence_beta_obs) != n ||
        ncol(influence_beta_obs) != p
    ) {
      warning(
        "VGAM::Influence(fit) did not return the expected [n x p] matrix; ",
        "setting influence_beta_obs to NULL.",
        call. = FALSE
      )
      influence_beta_obs <- NULL
    }
  }

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
    meat <- crossprod(scores_clustered)

    # Apply small-sample correction if requested
    if (adjust) {
      meat <- meat * (n_clusters / (n_clusters - 1))
    }
  } else {
    n_clusters <- NULL
    # Standard (unclustered) meat matrix
    meat <- crossprod(scores)
  }

  # --- 4. Compute the sandwich: V = B * M * B ---
  robust_var <- bread %*% meat %*% bread

  # Ensure symmetry (numerical stability)
  robust_var <- (robust_var + t(robust_var)) / 2
  dimnames(robust_var) <- dimnames(model_vcov)

  # Extract standard errors
  robust_se <- sqrt(diag(robust_var))
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
    meat = meat,

    # Clustering information
    cluster = cluster,
    n = n,
    n_clusters = n_clusters,
    adjust = adjust,

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
  if (!is.null(object$n_clusters)) {
    cat("Number of clusters:", object$n_clusters, "\n")
    if (object$adjust) {
      cat("Small-sample adjustment: applied (G/(G-1))\n")
    }
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
    n_clusters = object$n_clusters
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
#' Ratios close to 1.0 indicate good agreement. Small differences (~1-2%) are
#' expected for ordinal outcomes due to different internal implementations.
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

    vglm_robust <- robcov_vglm(vglm_fit, cluster = cluster)
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
