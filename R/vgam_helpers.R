#' Calculate Effective Coefficients for Fast Prediction
#'
#' Computes the "effective coefficients" matrix for a fitted ordinal model.
#' For `vglm` objects, this matrix combines the raw coefficients with the
#' constraint matrices. For `orm` objects, threshold coefficients are stored in
#' the intercept column and slope coefficients are replicated across
#' thresholds. In both cases this allows fast linear predictor calculation via
#' simple matrix multiplication with the "small" design matrix.
#'
#' @param model A fitted `vglm` or `orm` object.
#' @param beta Optional coefficient vector to use instead of `coef(model)`.
#'   This is used internally for simulation draws.
#'
#' @return A matrix with dimensions `[M x P]`, where:
#'   \item{M}{Number of linear predictors (e.g., number of states - 1)}
#'   \item{P}{Number of columns in the small design matrix (covariates)}
#'
#' @details
#' For a `vglm` term \eqn{j} with constraint matrix \eqn{C_j} and coefficient
#' vector \eqn{\beta_j}, the contribution to the linear predictor vector
#' \eqn{\eta} is \eqn{x_j \cdot (C_j \times \beta_j)}. This function
#' pre-calculates \eqn{\Gamma_j = C_j \times \beta_j} for all terms and
#' assembles them into a single matrix \eqn{\Gamma}. For `orm` models, only full
#' proportional odds models are supported, so each slope coefficient contributes
#' equally to every threshold-specific linear predictor.
#'
#' Then, prediction for a new observation with design vector \eqn{x} is simply:
#' \eqn{\eta = \Gamma \times x^T}
#'
#' @export
get_effective_coefs <- function(model, beta = NULL) {
  if (inherits(model, "vglm")) {
    return(get_effective_coefs_vglm(model, beta = beta))
  }

  if (inherits(model, "orm")) {
    return(get_effective_coefs_orm(model, beta = beta))
  }

  stop("get_effective_coefs only supports vglm and orm objects.")
}

get_effective_coefs_vglm <- function(model, beta = NULL) {
  if (is.null(beta)) {
    beta <- coef(model)
  }

  C_list <- VGAM::constraints(model)

  M <- nrow(C_list[[1]])
  P <- length(C_list)
  term_names <- names(C_list)

  Gamma <- matrix(0, nrow = M, ncol = P)
  colnames(Gamma) <- term_names
  rownames(Gamma) <- paste0("eta", 1:M)

  # 5. Iterate through terms and consume coefficients
  # vglm stores coefficients in the order of the terms in constraints()
  current_idx <- 1

  for (j in seq_along(C_list)) {
    term <- term_names[j]
    C_j <- C_list[[j]] # Matrix [M x k]

    # How many coefs for this term? (k)
    k <- ncol(C_j)

    # Extract the chunk of beta
    beta_chunk <- beta[current_idx:(current_idx + k - 1)]

    # Compute effective coefficients for this term: C_j %*% beta_j
    # Result is a vector of length M
    gamma_j <- C_j %*% beta_chunk

    # Store in Gamma matrix
    Gamma[, term] <- as.vector(gamma_j)

    # Advance index
    current_idx <- current_idx + k
  }

  # Sanity check
  if (current_idx - 1 != length(beta)) {
    warning(
      "Mismatch in coefficient consumption. ",
      "Expected ",
      length(beta),
      ", consumed ",
      current_idx - 1
    )
  }

  return(Gamma)
}

get_effective_coefs_orm <- function(model, beta = NULL) {
  if (is.null(beta)) {
    beta <- stats::coef(model)
  }

  if (is.null(model$non.slopes) || length(model$non.slopes) != 1) {
    stop("Cannot determine the number of orm threshold coefficients.")
  }

  M <- model$non.slopes
  if (length(beta) < M) {
    stop("Coefficient vector is shorter than the number of orm thresholds.")
  }

  if (length(beta) != length(stats::coef(model))) {
    stop(
      "Length of beta (",
      length(beta),
      ") does not match number of model coefficients (",
      length(stats::coef(model)),
      ")."
    )
  }

  if (!is.null(names(stats::coef(model)))) {
    names(beta) <- names(stats::coef(model))
  }

  intercepts <- beta[seq_len(M)]
  slopes <- beta[-seq_len(M)]
  slope_names <- names(stats::coef(model))[-seq_len(M)]

  Gamma <- matrix(
    rep(slopes, each = M),
    nrow = M,
    ncol = length(slopes),
    dimnames = list(
      paste0("eta", seq_len(M)),
      slope_names
    )
  )

  Gamma <- cbind("(Intercept)" = intercepts, Gamma)
  rownames(Gamma) <- paste0("eta", seq_len(M))
  Gamma
}
