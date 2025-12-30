#' Robust (Sandwich) Covariance Matrix Estimation for vglm Objects
#'
#' Computes the Huber-White sandwich covariance matrix estimator for vglm
#' objects, optionally with cluster-robust standard errors. This provides
#' functionality similar to `rms::robcov()` for vglm models.
#'
#' @param fit A fitted vglm object from the VGAM package.
#' @param cluster Optional vector of cluster identifiers. If provided, computes
#'   cluster-robust standard errors. Must have the same length as the number of
#'   observations.
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
#'   \item{bread}{The "bread" matrix (scaled inverse Hessian)}
#'   \item{meat}{The "meat" matrix}
#'
#'   **Clustering information:**
#'   \item{cluster}{The cluster variable (if provided)}
#'   \item{n}{Number of observations}
#'   \item{n_clusters}{Number of clusters (if clustered)}
#'   \item{adjust}{Whether small-sample adjustment was applied}
#'
#' @details
#' The sandwich estimator has the form:
#' \deqn{V = B \cdot M \cdot B / n}
#' where B is the "bread" (n times the inverse of the negative Hessian, i.e.,
#' n * vcov(fit)) and M is the "meat" (the variance of the score contributions).
#'
#' For unclustered data, the meat matrix is:
#' \deqn{M = \frac{1}{n} \sum_{i=1}^{n} \psi_i \psi_i'}
#' where \eqn{\psi_i} is the score (gradient of log-likelihood) for observation i.
#'
#' For clustered data, scores are summed within clusters before computing the
#' meat matrix:
#' \deqn{M = \frac{1}{n} \sum_{g=1}^{G} \left(\sum_{i \in g} \psi_i\right)
#'       \left(\sum_{i \in g} \psi_i\right)'}
#'
#' When `adjust = TRUE`, the meat matrix is multiplied by G/(G-1) for clustered
#' data, providing a small-sample bias correction. This is sometimes recommended
#' but is not applied by `rms::robcov()`, so the default is FALSE for
#' consistency.
#'
#' **Z-statistics and p-values**: The returned object includes z-statistics
#' computed as coefficients divided by robust standard errors, and two-sided
#' p-values from the standard normal distribution. These can be used for
#' inference that is robust to model misspecification or within-cluster
#' correlation.
#'
#' **Comparison with rms::robcov()**: This function produces nearly identical
#' results to `rms::robcov()` when used with equivalent models:
#' - For binary outcomes (vglm with binomialff vs orm), results match to
#'   numerical precision
#' - For ordinal outcomes (vglm with cumulative vs orm), small differences
#'   (~0.1-1.5%) may occur due to different internal parameterizations, but
#'   the results are practically equivalent
#'
#' @examples
#' \\dontrun{
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

  # Get the bread matrix: n * vcov(fit)
  # vcov(fit) is the inverse of the negative Hessian
  model_vcov <- vcov(fit)
  original_se <- sqrt(diag(model_vcov))
  names(original_se) <- names(coefficients)
  bread <- n * model_vcov

  # --- 2. Compute observation-level score contributions ---
  scores <- compute_scores_vglm(fit)

  # Validate scores dimensions

  if (nrow(scores) != n) {
    stop("Number of score rows does not match number of observations")
  }
  if (ncol(scores) != p) {
    stop("Number of score columns does not match number of parameters")
  }

  # --- 3. Compute the meat matrix ---
  if (!is.null(cluster)) {
    if (length(cluster) != n) {
      stop("Length of 'cluster' must equal number of observations (", n, ")")
    }

    # Convert to factor and aggregate scores by cluster
    cluster <- as.factor(cluster)
    n_clusters <- nlevels(cluster)

    # Sum scores within clusters using efficient rowsum
    scores_clustered <- rowsum(scores, cluster, reorder = FALSE)

    # Compute meat matrix from clustered scores
    meat <- crossprod(scores_clustered) / n

    # Apply small-sample correction if requested
    if (adjust && n_clusters > 1) {
      meat <- meat * (n_clusters / (n_clusters - 1))
    }
  } else {
    n_clusters <- NULL
    # Standard (unclustered) meat matrix
    meat <- crossprod(scores) / n
  }

  # --- 4. Compute the sandwich: V = B * M * B / n ---
  robust_var <- (bread %*% meat %*% bread) / n

  # Ensure symmetry (numerical stability)
  robust_var <- (robust_var + t(robust_var)) / 2

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
    tryCatch(slot(obj, name), error = function(e) NULL)
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
    bread = bread,
    meat = meat,

    # Clustering information
    cluster = cluster,
    n = n,
    n_clusters = n_clusters,
    adjust = adjust,

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
  # Get dimensions
  n <- nobs(fit, type = "lm")
  M <- VGAM::npred(fit)
  p <- length(coef(fit))

  # Get working weights with derivatives
  # wt_info$deriv is an n x M matrix containing d(log L)/d(eta_m) for each obs
  wt_info <- VGAM::weights(fit, type = "working", deriv = TRUE)
  deriv_eta <- wt_info$deriv

  # Get the VLM model matrix (n*M x p)
  X_vlm <- model.matrix(fit, type = "vlm")

  # Initialize score matrix
  scores <- matrix(0, nrow = n, ncol = p)

  # Compute scores for each observation
  # Score_i = X_i' %*% deriv_eta_i
  # where X_i is the M x p submatrix of X_vlm for observation i
  for (i in seq_len(n)) {
    vlm_rows <- ((i - 1) * M + 1):(i * M)
    X_i <- X_vlm[vlm_rows, , drop = FALSE]
    deriv_i <- deriv_eta[i, ]
    scores[i, ] <- as.vector(crossprod(X_i, deriv_i))
  }

  colnames(scores) <- names(coef(fit))
  return(scores)
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
