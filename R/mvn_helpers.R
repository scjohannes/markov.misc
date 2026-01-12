#' Set Coefficients for Model Objects
#'
#' Replaces the coefficients in a fitted model object. This is used for
#' simulation-based inference (MVN simulation) where we draw coefficient
#' vectors from a multivariate normal distribution and need to compute
#' predictions using each draw.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, `lrm`).
#' @param new_coefs A numeric vector of coefficients to set. Must match
#'   the length and order of the original coefficients.
#'
#' @return The model object with updated coefficients.
#'
#' @details
#' This function creates a shallow copy of the model and replaces its

#' coefficients. The modified model can then be passed to `predict()` or
#' `soprob_markov()` to compute predictions under the simulated coefficients.
#'
#' **Supported model classes:**
#' \itemize{
#'   \item `orm`, `lrm`, `rms` (S3): Modifies `$coefficients`
#'   \item `vglm`, `vgam` (S4): Modifies `@coefficients` slot
#'   \item `robcov_vglm`: Modifies `$coefficients`
#' }
#'
#' **Note:** Prediction functions (`predict.orm`, `predict.vglm`) recompute
#' linear predictors from the design matrix and coefficients, so we only need
#' to replace the coefficients - not any stored linear predictors.
#'
#' @examples
#' \dontrun{
#' library(VGAM)
#'
#' # Fit a model
#' fit <- vglm(y ~ x1 + x2,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = mydata)
#'
#' # Draw from MVN
#' library(mvtnorm)
#' new_coefs <- rmvnorm(1, mean = coef(fit), sigma = vcov(fit))
#'
#' # Set new coefficients and predict
#' fit_sim <- set_coef(fit, new_coefs[1, ])
#' pred_sim <- predict(fit_sim, newdata = test_data, type = "response")
#' }
#'
#' @seealso [get_vcov_robust()] for obtaining robust covariance matrices
#'
#' @export
set_coef <- function(model, new_coefs) {
  UseMethod("set_coef")
}


#' @rdname set_coef
#' @export
set_coef.default <- function(model, new_coefs) {
  # Default method for S3 objects with $coefficients
  if (is.null(model$coefficients)) {
    stop(
      "Model does not have a 'coefficients' element. ",
      "set_coef() may not support this model class: ", class(model)[1]
    )
  }

  if (length(new_coefs) != length(model$coefficients)) {
    stop(
      "Length of new_coefs (", length(new_coefs), ") does not match ",
      "number of model coefficients (", length(model$coefficients), ")"
    )
  }

  # Preserve names if they exist
  if (!is.null(names(model$coefficients))) {
    names(new_coefs) <- names(model$coefficients)
  }

  model$coefficients <- new_coefs
  model
}


#' @rdname set_coef
#' @export
set_coef.orm <- function(model, new_coefs) {
  set_coef.default(model, new_coefs)
}


#' @rdname set_coef
#' @export
set_coef.lrm <- function(model, new_coefs) {
  set_coef.default(model, new_coefs)
}


#' @rdname set_coef
#' @export
set_coef.rms <- function(model, new_coefs) {
  set_coef.default(model, new_coefs)
}


#' @rdname set_coef
#' @export
set_coef.robcov_vglm <- function(model, new_coefs) {
  set_coef.default(model, new_coefs)
}


#' @rdname set_coef
#' @export
set_coef.vglm <- function(model, new_coefs) {
  # S4 object - use @ for slot access
  if (length(new_coefs) != length(model@coefficients)) {
    stop(
      "Length of new_coefs (", length(new_coefs), ") does not match ",
      "number of model coefficients (", length(model@coefficients), ")"
    )
  }

  # Preserve names if they exist
  if (!is.null(names(model@coefficients))) {
    names(new_coefs) <- names(model@coefficients)
  }

  model@coefficients <- new_coefs
  model
}


#' @rdname set_coef
#' @export
set_coef.vgam <- function(model, new_coefs) {
  set_coef.vglm(model, new_coefs)
}


#' Get Robust Variance-Covariance Matrix
#'
#' Extracts or computes a (cluster-)robust variance-covariance matrix for
#' use in simulation-based inference. This function provides a unified
#' interface for obtaining robust vcov from different model types.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, `lrm`) or a
#'   `robcov_vglm` object (which already contains robust vcov).
#' @param cluster Optional. Specifies how to compute cluster-robust standard
#'   errors:
#'   \itemize{
#'     \item `NULL`: Use standard (model-based) vcov
#'     \item Formula (e.g., `~id`): Extract cluster variable from `data` and
#'       compute cluster-robust vcov
#'     \item Character vector: Cluster variable values directly
#'   }
#' @param data Data frame. Required when `cluster` is a formula.
#' @param adjust Logical. Apply small-sample correction for clustering?
#'   Default is FALSE to match `rms::robcov()` behavior.
#'
#' @return A variance-covariance matrix.
#'
#' @details
#' This function handles several cases:
#'
#' **1. Pre-computed robust vcov (`robcov_vglm` object):**
#' If the model is already a `robcov_vglm` object (from calling
#' `robcov_vglm()`), the stored robust vcov is returned directly.
#'
#' **2. No clustering requested:**
#' If `cluster = NULL`, returns the standard model-based vcov from `vcov()`.
#'
#' **3. Cluster formula:**
#' If `cluster` is a formula like `~id`, extracts the variable from `data`
#' and computes cluster-robust vcov using:
#' \itemize{
#'   \item `robcov_vglm()` for `vglm` models
#'   \item `rms::robcov()` for `orm`/`lrm` models
#' }
#'
#' @examples
#' \dontrun{
#' library(VGAM)
#'
#' # Fit a model
#' fit <- vglm(y ~ yprev + time + tx,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = mydata)
#'
#' # Standard vcov
#' V1 <- get_vcov_robust(fit)
#'
#' # Cluster-robust vcov using formula
#' V2 <- get_vcov_robust(fit, cluster = ~id, data = mydata)
#'
#' # Pre-computed robust vcov
#' fit_robust <- robcov_vglm(fit, cluster = mydata$id)
#' V3 <- get_vcov_robust(fit_robust)  # Same as fit_robust$var
#' }
#'
#' @seealso [robcov_vglm()], [rms::robcov()], [set_coef()]
#'
#' @importFrom stats vcov
#' @export
get_vcov_robust <- function(model, cluster = NULL, data = NULL,
                            adjust = FALSE) {

  # --- Case 1: Already a robcov_vglm object ---
  if (inherits(model, "robcov_vglm")) {
    return(model$var)
  }

  # --- Case 2: rms model with robcov already applied ---
  # rms::robcov() stores robust var in fit$var, original in fit$orig.var
  if (inherits(model, c("orm", "lrm")) && !is.null(model$orig.var)) {
    return(model$var)
  }

  # --- Case 3: No clustering requested ---
  if (is.null(cluster)) {
    return(vcov(model))
  }

  # --- Case 4: Extract cluster variable from formula ---
  if (inherits(cluster, "formula")) {
    cluster_var <- all.vars(cluster)
    if (length(cluster_var) != 1) {
      stop(
        "cluster formula must specify exactly one variable (e.g., ~id). ",
        "Got: ", deparse(cluster)
      )
    }
    if (is.null(data)) {
      stop("data is required when cluster is a formula")
    }
    if (!cluster_var %in% names(data)) {
      stop("Cluster variable '", cluster_var, "' not found in data")
    }
    cluster <- data[[cluster_var]]
  }

  # --- Case 5: Compute robust vcov based on model class ---
  if (inherits(model, "vglm") || inherits(model, "vgam")) {
    # Use our robcov_vglm function
    robcov_result <- robcov_vglm(model, cluster = cluster, adjust = adjust)
    return(robcov_result$var)

  } else if (inherits(model, c("orm", "lrm"))) {
    # Use rms::robcov
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("Package 'rms' is required for robust vcov with orm/lrm models")
    }
    robcov_result <- rms::robcov(model, cluster = cluster)
    return(robcov_result$var)

  } else {
    stop(
      "Unsupported model class for robust vcov: ", class(model)[1], "\n",
      "Supported classes: vglm, vgam, orm, lrm, robcov_vglm"
    )
  }
}


#' Extract Coefficients from Model Objects (Generic Helper)
#'
#' A helper function that extracts coefficients from various model types
#' in a consistent way.
#'
#' @param model A fitted model object.
#'
#' @return A named numeric vector of coefficients.
#'
#' @keywords internal
get_coef <- function(model) {
  if (inherits(model, "robcov_vglm")) {
    return(model$coefficients)
  }
  stats::coef(model)
}
