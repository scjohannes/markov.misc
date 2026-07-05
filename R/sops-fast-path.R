# SOP fast-path matrix prediction.

# =============================================================================
# FAST PATH HELPER FUNCTIONS
# =============================================================================
# These internal functions provide the fast path for simulation-based inference
# on VGLM models by pre-computing design matrices and bypassing VGAM's
# predict() function.

#' Build Fast Markov Simulation Components (Internal)
#'
#' Pre-calculates design matrices for fast Markov simulation. This is the
#' "heavy lifting" step that should be done ONCE before running thousands of
#' simulation draws.
#'
#' @param model Fitted vglm model
#' @param data Baseline data frame (one row per patient)
#' @param t_covs Time-dependent covariate lookup (optional)
#' @param times Time points to predict
#' @param ylevels State labels
#' @param absorb Absorbing state(s). These transition matrices are not needed
#'   because absorbing states keep their prior probability mass.
#' @param tvarname Name of time variable
#' @param pvarname Name of previous state variable
#' @param gap Optional gap variable. With factor visit time, numeric gap values
#'   must be supplied in `t_covs`.
#' @param ... Ignored
#'
#' @return A list containing pre-calculated matrices and metadata.
#'
#' @keywords internal
markov_msm_build <- function(
  model,
  data,
  t_covs = NULL,
  times = NULL,
  ylevels = 1:6,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  gap = NULL,
  ...
) {
  n_pat <- nrow(data)
  ylevel_names <- as_state_labels(ylevels)
  absorb_names <- as_state_labels(absorb)
  n_states <- length(ylevel_names)
  M <- n_states - 1
  absorb_idx <- if (length(absorb_names)) {
    which(ylevel_names %in% absorb_names)
  } else {
    integer(0)
  }

  time_res <- resolve_sop_times(
    model,
    data,
    times,
    tvarname,
    t_covs = t_covs,
    default = "fast"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap, t_covs, time_info)

  # Need Gamma structure to know which columns to keep
  Gamma_template <- get_effective_coefs(model)
  common_cols <- colnames(Gamma_template)

  # Prepare data
  if (!pvarname %in% names(data)) {
    data[[pvarname]] <- factor(ylevel_names[1], levels = ylevel_names)
  }
  if (!tvarname %in% names(data)) {
    data <- assign_sop_time(data, tvarname, times[1], time_info)
  }

  is_vglm <- inherits(model, "vglm")
  is_orm <- inherits(model, "orm")

  # Terms object
  if (is_vglm) {
    tt <- stats::terms(model)
    tt <- stats::delete.response(tt)
    contrasts_arg <- if (length(model@contrasts)) {
      model@contrasts
    } else {
      NULL
    }
    xlev_arg <- model@xlevels
  } else if (is_orm) {
    tt <- NULL
    contrasts_arg <- NULL
    xlev_arg <- NULL
  } else {
    stop("markov_msm_build() supports only vglm and orm models.")
  }

  # Helper to get design matrix columns used by the effective coefficient
  # matrix. For inline transforms such as rms::rcs(), model.matrix() uses the
  # fitted terms object's predvars, including stored knot locations.
  get_X <- function(d) {
    if (is_orm) {
      X <- orm_model_matrix(model, d, include_intercept = TRUE)
    } else {
      X <- stats::model.matrix(
        tt,
        data = d,
        contrasts.arg = contrasts_arg,
        xlev = xlev_arg
      )
    }
    common <- intersect(colnames(X), common_cols)
    if (length(common) == 0) {
      return(NULL)
    }
    X[, common, drop = FALSE]
  }

  set_time <- function(d, time_idx) {
    assign_sop_visit(
      d,
      tvarname = tvarname,
      times = times,
      index = time_idx,
      t_covs = t_covs,
      gap = gap,
      time_info = time_info
    )
  }

  # A. Baseline matrix for T=1 uses each patient's observed previous state.
  d_init <- set_time(data, 1)
  d_init[[pvarname]] <- normalize_previous_state_column(
    d_init[[pvarname]],
    data[[pvarname]],
    pvarname
  )
  X_init <- get_X(d_init)
  if (is.null(X_init)) {
    stop("Could not construct a design matrix for the fitted model.")
  }
  col_names <- colnames(X_init)

  align_X <- function(X) {
    missing_cols <- setdiff(col_names, colnames(X))
    if (length(missing_cols) > 0) {
      stop(
        "Design matrix columns missing during fast path build: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    X[, col_names, drop = FALSE]
  }

  # B. Transition matrices: one matrix for each time and previous-state value.
  X_transition <- vector("list", length(times))
  for (time_idx in seq_along(times)) {
    X_transition[[time_idx]] <- vector("list", n_states)
    d_time <- set_time(data, time_idx)

    for (k in seq_len(n_states)) {
      if (k %in% absorb_idx) {
        X_transition[[time_idx]][[k]] <- NULL
        next
      }

      d_k <- d_time
      d_k[[pvarname]] <- make_previous_state_column(
        states = ylevel_names[k],
        prototype = data[[pvarname]],
        n = nrow(d_k),
        pvarname = pvarname
      )
      X_transition[[time_idx]][[k]] <- align_X(get_X(d_k))
    }
  }

  list(
    X_init = align_X(X_init),
    X_transition = X_transition,
    n_pat = n_pat,
    n_states = n_states,
    M = M,
    ylevels = ylevel_names,
    col_names = col_names
  )
}


#' Run Fast Markov Simulation (Internal)
#'
#' Runs the Markov loop using pre-calculated components and a specific
#' coefficient vector.
#'
#' @param components List returned by `markov_msm_build`
#' @param Gamma Effective coefficient matrix with dimensions M x P (M thresholds by P predictors)
#' @param times Vector of time points
#' @param absorb Absorbing state(s). Can be NULL for no absorbing states.
#'
#' @return Array of state probabilities with dimensions n_pat x n_times x n_states
#'
#' @keywords internal
markov_msm_run <- function(components, Gamma, times, absorb = NULL) {
  # Unpack
  X_init <- components$X_init
  X_transition <- components$X_transition
  n_pat <- components$n_pat
  n_states <- components$n_states
  M <- components$M
  ylevels <- components$ylevels
  col_names <- components$col_names
  n_times <- length(times)

  # Align Gamma to components
  Gamma <- Gamma[, col_names, drop = FALSE]
  Gamma_t <- t(Gamma) # [P x M]

  # Helper: X * Gamma_t -> LP [N x M]
  calc_lp <- function(X) {
    X %*% Gamma_t
  }

  # Simulation
  P_out <- array(0, dim = c(n_pat, n_times, n_states))
  dimnames(P_out)[[3]] <- ylevels

  absorb_idx <- if (!is.null(absorb)) {
    which(as.character(ylevels) %in% as.character(absorb))
  } else {
    integer(0)
  }
  non_absorb_idx <- setdiff(1:n_states, absorb_idx)

  # T=1
  P_out[, 1, ] <- lp_to_probs(calc_lp(X_init), M)

  # Loop 2..T
  if (n_times >= 2) {
    for (t_idx in 2:n_times) {
      p_current <- matrix(0, nrow = n_pat, ncol = n_states)

      for (k in non_absorb_idx) {
        p_prev_k <- P_out[, t_idx - 1, k]
        if (max(p_prev_k) < 1e-12) {
          next
        }

        LP_k <- calc_lp(X_transition[[t_idx]][[k]])
        trans_k <- lp_to_probs(LP_k, M)
        p_current <- p_current + (trans_k * p_prev_k)
      }

      for (a in absorb_idx) {
        p_current[, a] <- p_current[, a] + P_out[, t_idx - 1, a]
      }
      P_out[, t_idx, ] <- p_current
    }
  }

  P_out
}


#' Convert Cumulative Log-Odds to Probabilities (Internal)
#'
#' Converts a matrix of cumulative log-odds to category probabilities.
#'
#' @param eta Matrix of linear predictors with dimensions N x M (N observations by M thresholds)
#' @param M Number of thresholds (one less than the number of categories)
#' @return Matrix of probabilities with dimensions N x (M+1)
#'
#' @keywords internal
lp_to_probs <- function(eta, M) {
  cum_probs <- stats::plogis(eta) # [N x M]
  n_rows <- nrow(eta)
  probs <- matrix(0, nrow = n_rows, ncol = M + 1)

  probs[, 1] <- 1 - cum_probs[, 1]

  if (M > 1) {
    probs[, 2:M] <- cum_probs[, 1:(M - 1)] - cum_probs[, 2:M]
  }

  probs[, M + 1] <- cum_probs[, M]
  probs[probs < 0] <- 0

  normalize_probability_rows(probs)
}

normalize_probability_rows <- function(probs) {
  totals <- rowSums(probs, na.rm = TRUE)
  valid <- is.finite(totals) & totals > 0
  probs[valid, ] <- probs[valid, , drop = FALSE] / totals[valid]
  probs
}

normalize_probability_array <- function(probs) {
  totals <- apply(probs, c(1L, 2L), sum, na.rm = TRUE)
  valid <- is.finite(totals) & totals > 0
  for (k in seq_len(dim(probs)[3])) {
    slice <- probs[,, k]
    slice[valid] <- slice[valid] / totals[valid]
    probs[,, k] <- slice
  }
  probs
}


#' Compute Effective Coefficients from Beta Vector (Internal)
#'
#' Transforms a vector of VGLM coefficients into an effective coefficient matrix
#' by applying the constraint matrices.
#'
#' @param beta Coefficient vector from a VGLM fit
#' @param C_list Constraint matrices list from VGAM::constraints()
#' @return Effective coefficient matrix with dimensions M x P (M linear predictors by P terms)
#'
#' @keywords internal
compute_Gamma <- function(beta, C_list) {
  M <- nrow(C_list[[1]])
  P <- length(C_list)
  term_names <- names(C_list)

  G <- matrix(0, nrow = M, ncol = P)
  colnames(G) <- term_names
  rownames(G) <- paste0("eta", 1:M)

  curr <- 1
  for (j in seq_along(C_list)) {
    k <- ncol(C_list[[j]])
    chunk <- beta[curr:(curr + k - 1)]
    G[, j] <- C_list[[j]] %*% chunk
    curr <- curr + k
  }
  G
}
