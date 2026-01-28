#' Calculate Effective Coefficients for Fast Prediction
#'
#' Computes the "effective coefficients" matrix for a fitted vglm model.
#' This matrix combines the raw coefficients with the constraint matrices,
#' allowing for fast linear predictor calculation via simple matrix multiplication
#' with the "small" (unexpanded) design matrix.
#'
#' @param model A fitted `vglm` object.
#'
#' @return A matrix with dimensions `[M x P]`, where:
#'   \item{M}{Number of linear predictors (e.g., number of states - 1)}
#'   \item{P}{Number of columns in the small design matrix (covariates)}
#'
#' @details
#' For a term \eqn{j} with constraint matrix \eqn{C_j} and coefficient vector \eqn{\beta_j},
#' the contribution to the linear predictor vector \eqn{\eta} is \eqn{x_j \cdot (C_j \times \beta_j)}.
#' This function pre-calculates \eqn{\Gamma_j = C_j \times \beta_j} for all terms and
#' assembles them into a single matrix \eqn{\Gamma}.
#'
#' Then, prediction for a new observation with design vector \eqn{x} is simply:
#' \eqn{\eta = \Gamma \times x^T}
#'
#' @export
get_effective_coefs <- function(model) {
  if (!inherits(model, "vglm")) {
    stop("get_effective_coefs only supports vglm objects.")
  }

  # 1. Get raw coefficients
  # coef() returns a named vector
  beta <- coef(model)

  # 2. Get constraints list
  # The names match the columns of the "small" design matrix
  C_list <- VGAM::constraints(model)

  # 3. Determine number of linear predictors (M) from the first constraint
  # All constraint matrices must have M rows
  M <- nrow(C_list[[1]])

  # 4. Initialize effective coefficient matrix Gamma
  # Rows = M (linear predictors)
  # Cols = P (number of design matrix columns / constraint terms)
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
