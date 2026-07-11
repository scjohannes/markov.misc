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
#'   \item `vglm` (S4): Modifies `@coefficients` slot
#'   \item `robcov_vglm`: Modifies `$coefficients`
#' }
#'
#' **Note:** Prediction functions (`predict.orm`, `predict.vglm`) recompute
#' linear predictors from the design matrix and coefficients, so we only need
#' to replace the coefficients - not any stored linear predictors.
#'
#' @examplesIf rlang::is_installed("rms")
#' trial <- sim_actt2_brownian(n_patients = 40, follow_up_time = 6, seed = 1)
#' markov_data <- prepare_markov_data(trial, absorbing_state = 8)
#' fit <- rms::orm(
#'   y ~ time + tx + yprev,
#'   data = markov_data,
#'   x = TRUE,
#'   y = TRUE,
#'   opt_method = "LM",
#'   scale = TRUE
#' )
#'
#' new_coefs <- coef(fit)
#' new_coefs["tx"] <- new_coefs["tx"] - 0.1
#' fit_shifted <- set_coef(fit, new_coefs)
#' coef(fit_shifted)["tx"]
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
      "set_coef() may not support this model class: ",
      class(model)[1]
    )
  }

  if (length(new_coefs) != length(model$coefficients)) {
    stop(
      "Length of new_coefs (",
      length(new_coefs),
      ") does not match ",
      "number of model coefficients (",
      length(model$coefficients),
      ")"
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
  # 1. Update the top-level coefficients (robcov structure)
  model <- set_coef.default(model, new_coefs)

  # 2. Update the internal vglm_fit coefficients if present
  # This ensures that predict(model$vglm_fit, ...) uses the new coefficients
  if (!is.null(model$vglm_fit)) {
    model$vglm_fit <- set_coef(model$vglm_fit, new_coefs)
  }

  model
}


#' @rdname set_coef
#' @export
set_coef.vglm <- function(model, new_coefs) {
  # S4 object - use @ for slot access
  if (length(new_coefs) != length(model@coefficients)) {
    stop(
      "Length of new_coefs (",
      length(new_coefs),
      ") does not match ",
      "number of model coefficients (",
      length(model@coefficients),
      ")"
    )
  }

  # Preserve names if they exist
  if (!is.null(names(model@coefficients))) {
    names(new_coefs) <- names(model@coefficients)
  }

  model@coefficients <- new_coefs
  model
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
#'     \item `NULL`: Use robust (sandwich) vcov, but without clustering
#'     \item Formula (e.g., `~id`): Extract cluster variable from the model's
#'       original fitting data and compute cluster-robust vcov
#'     \item Character/numeric vector: Cluster variable values directly (must
#'       match the number of observations used to fit the model)
#'   }
#' @param data Data frame. Optional. If `cluster` is a formula and the cluster
#'   variable cannot be found in the model's original data, it will be extracted
#'   from this data frame. **Note:** For cluster-robust SEs, the cluster vector
#'   must have the same length as the model's fitting data, not prediction data.
#' @param adjust Logical. Apply small-sample correction for clustering?
#'   Default is TRUE.
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
#' If `cluster = NULL`, returns the robust (sandwich) vcov. For `vglm` models,
#' this is computed via [robcov_vglm()]; for `orm` models, it uses
#' `rms::robcov()`.
#'
#' **3. Cluster formula:**
#' If `cluster` is a formula like `~id`, the function first tries to extract
#' the cluster variable from the model's original fitting data (by evaluating
#' the `data` argument from the model's call). This ensures the cluster vector
#' has the correct length for computing robust standard errors.
#'
#' **Important:** Cluster-robust standard errors are computed from the model's
#' score contributions, which have length equal to the number of observations
#' used to fit the model. The cluster variable must therefore match this length,
#' NOT the length of any new prediction data.
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
#' # Robust (unclustered) vcov
#' V1 <- get_vcov_robust(fit)
#'
#' # Cluster-robust vcov using formula (extracts 'id' from model's data)
#' V2 <- get_vcov_robust(fit, cluster = ~id)
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
get_vcov_robust <- function(
  model,
  cluster = NULL,
  data = NULL,
  adjust = TRUE
) {
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
    if (inherits(model, "vglm")) {
      return(robcov_vglm(model, cluster = NULL, adjust = adjust)$var)
    } else if (inherits(model, c("orm", "lrm"))) {
      if (!requireNamespace("rms", quietly = TRUE)) {
        stop("Package 'rms' is required for robust vcov with orm/lrm models")
      }
      return(rms::robcov(model)$var)
    }
    # Fallback for other models
    return(vcov(model))
  }

  # --- Case 4: Extract cluster variable from formula ---
  if (inherits(cluster, "formula")) {
    cluster_var <- all.vars(cluster)
    if (length(cluster_var) != 1) {
      stop(
        "cluster formula must specify exactly one variable (e.g., ~id). ",
        "Got: ",
        deparse(cluster)
      )
    }

    # Prefer the explicit `data` argument when supplied. Model calls often store
    # only a symbol such as `data`, and evaluating that symbol from a caller can
    # resolve to an unrelated object in interactive/test wrapper environments.
    model_data <- markov_model_data(model)
    if (!is.null(model_data) && !cluster_var %in% names(model_data)) {
      model_data <- NULL
    }
    if (!is.null(data) && cluster_var %in% names(data)) {
      model_data <- data
    }

    if (is.null(model_data) && inherits(model, "vglm")) {
      # For S4 vglm objects, try to evaluate the data from the call
      if (!is.null(model@call$data)) {
        model_data <- tryCatch(
          eval(model@call$data, envir = parent.frame(2)),
          error = function(e) NULL
        )
      }
    } else if (is.null(model_data) && !is.null(model$call$data)) {
      # For S3 objects (orm, lrm), try to evaluate data from call
      model_data <- tryCatch(
        eval(model$call$data, envir = parent.frame(2)),
        error = function(e) NULL
      )
    }

    model_nobs <- tryCatch(
      if (inherits(model, "vglm")) {
        stats::nobs(model, type = "lm")
      } else {
        stats::nobs(model)
      },
      error = function(e) NULL
    )

    align_cluster <- function(cluster_values, source_data) {
      if (is.null(model_nobs) || length(cluster_values) == model_nobs) {
        return(cluster_values)
      }

      model_vars <- tryCatch(
        all.vars(stats::formula(model)),
        error = function(e) character(0)
      )
      model_vars <- intersect(model_vars, names(source_data))
      if (length(model_vars) == 0) {
        return(cluster_values)
      }

      keep <- stats::complete.cases(source_data[, model_vars, drop = FALSE])
      if (sum(keep) == model_nobs) {
        return(cluster_values[keep])
      }

      cluster_values
    }

    if (!is.null(model_data) && cluster_var %in% names(model_data)) {
      cluster <- align_cluster(model_data[[cluster_var]], model_data)
    } else {
      stop(
        "Cluster variable '",
        cluster_var,
        "' not found.\n",
        "Checked: model's original fitting data",
        if (!is.null(data)) " and user-provided 'data' argument" else "",
        ".\nEnsure the cluster variable exists in the data used to fit the model."
      )
    }
  }

  # --- Case 5: Compute robust vcov based on model class ---
  if (inherits(model, "vglm")) {
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
      "Unsupported model class for robust vcov: ",
      class(model)[1],
      "\n",
      "Supported classes: vglm, orm, lrm, robcov_vglm"
    )
  }
}

validate_coef_vcov <- function(beta, Sigma, arg = "vcov") {
  if (methods::is(Sigma, "Matrix")) {
    Sigma <- as.matrix(Sigma)
  }

  if (!is.matrix(Sigma)) {
    stop("`", arg, "` must be a matrix.")
  }

  if (length(beta) != nrow(Sigma) || length(beta) != ncol(Sigma)) {
    stop(
      "Dimension mismatch: coefficients (",
      length(beta),
      ") vs ",
      arg,
      " matrix (",
      nrow(Sigma),
      " x ",
      ncol(Sigma),
      "). For orm models, use the full covariance matrix from ",
      "rms::robcov(fit)$var or rms::robcov(fit)$orig.var, not stats::vcov(fit)."
    )
  }

  beta_names <- names(beta)
  Sigma_names <- rownames(Sigma)
  if (
    !is.null(beta_names) &&
      !is.null(Sigma_names) &&
      !identical(beta_names, Sigma_names)
  ) {
    stop(
      "Coefficient names do not match row names of `",
      arg,
      "`. Expected: ",
      paste(beta_names, collapse = ", "),
      "; got: ",
      paste(Sigma_names, collapse = ", "),
      "."
    )
  }

  Sigma_colnames <- colnames(Sigma)
  if (
    !is.null(beta_names) &&
      !is.null(Sigma_colnames) &&
      !identical(beta_names, Sigma_colnames)
  ) {
    stop(
      "Coefficient names do not match column names of `",
      arg,
      "`."
    )
  }

  Sigma
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
