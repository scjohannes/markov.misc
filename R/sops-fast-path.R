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
#' @param time_covariates Time-dependent covariate lookup (optional)
#' @param times Time points to predict
#' @param y_levels State labels
#' @param absorb Absorbing state(s). These transition matrices are not needed
#'   because absorbing states keep their prior probability mass.
#' @param time_var Name of time variable
#' @param p_var Name of previous state variable
#' @param gap_var Optional gap_var variable. With factor visit time, numeric gap_var values
#'   must be supplied in `time_covariates`.
#' @param ... Ignored
#'
#' @return A list containing pre-calculated matrices and metadata.
#'
#' @keywords internal
markov_msm_build <- function(
  model,
  newdata = NULL,
  data = NULL,
  time_covariates = NULL,
  times = NULL,
  y_levels = 1:6,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  gap_var = NULL,
  ...
) {
  data <- newdata %||% data
  if (is.null(data)) {
    stop("`newdata` must be supplied.")
  }
  n_pat <- nrow(data)
  ylevel_names <- as_state_labels(y_levels)
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
    time_var,
    time_covariates = time_covariates,
    default = "fast"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap_var, time_covariates, time_info)

  is_vglm <- inherits(model, "vglm")
  is_orm <- inherits(model, "orm")
  if (!is_vglm && !is_orm) {
    stop("markov_msm_build() supports only vglm and orm models.")
  }

  # Need Gamma structure to know which columns to keep
  Gamma_template <- get_effective_coefs(model)
  validate_sop_threshold_count(Gamma_template, y_levels)
  common_cols <- colnames(Gamma_template)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)
  budget <- preflight_execution_plan_designs(
    n_pat = n_pat,
    n_times = length(times),
    rows_per_transition = n_pat * length(non_absorb_idx),
    n_cols = ncol(Gamma_template)
  )

  # Prepare data
  if (!p_var %in% names(data)) {
    data[[p_var]] <- factor(ylevel_names[1], levels = ylevel_names)
  }
  if (!time_var %in% names(data)) {
    data <- assign_sop_time(data, time_var, times[1], time_info)
  }

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
      time_var = time_var,
      times = times,
      index = time_idx,
      time_covariates = time_covariates,
      gap_var = gap_var,
      time_info = time_info
    )
  }

  # A. Baseline matrix for T=1 uses each patient's observed previous state.
  d_init <- set_time(data, 1)
  d_init[[p_var]] <- normalize_previous_state_column(
    d_init[[p_var]],
    data[[p_var]],
    p_var
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

  # B. Retained transition designs. Each visit stores one stacked matrix
  # with patient rows nested inside non-absorbing origin-state blocks.
  X_transition <- vector("list", length(times))
  if (length(times) >= 2L) {
    for (time_idx in 2:length(times)) {
      origin_blocks <- vector("list", length(non_absorb_idx))
      d_time <- set_time(data, time_idx)

      for (origin_pos in seq_along(non_absorb_idx)) {
        k <- non_absorb_idx[origin_pos]
        d_k <- d_time
        d_k[[p_var]] <- make_previous_state_column(
          states = ylevel_names[k],
          prototype = data[[p_var]],
          n = nrow(d_k),
          p_var = p_var
        )
        origin_blocks[[origin_pos]] <- align_X(get_X(d_k))
      }
      X_transition[[time_idx]] <- do.call(rbind, origin_blocks)
    }
  }

  design_bytes <- measure_execution_plan_designs(X_init, X_transition)
  validate_execution_plan_designs(design_bytes, budget$workspace_bytes)

  list(
    X_init = align_X(X_init),
    X_transition = X_transition,
    n_pat = n_pat,
    n_states = n_states,
    M = M,
    y_levels = ylevel_names,
    col_names = col_names,
    transition_layout = "stacked",
    transition_origins = non_absorb_idx,
    workspace_bytes = budget$workspace_bytes,
    design_bytes = design_bytes
  )
}

markov_msm_build_batched <- function(
  model,
  newdata,
  time_covariates = NULL,
  times,
  y_levels,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  gap_var = NULL
) {
  n_pat <- nrow(newdata)
  ylevel_names <- as_state_labels(y_levels)
  absorb_idx <- which(ylevel_names %in% as_state_labels(absorb))
  non_absorb_idx <- setdiff(seq_along(ylevel_names), absorb_idx)
  time_res <- resolve_sop_times(
    model,
    newdata,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "fast"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap_var, time_covariates, time_info)

  Gamma <- get_effective_coefs(model)
  validate_sop_threshold_count(Gamma, y_levels)
  budget <- preflight_execution_plan_designs(
    n_pat = n_pat,
    n_times = length(times),
    rows_per_transition = n_pat * length(non_absorb_idx),
    n_cols = ncol(Gamma)
  )

  if (!p_var %in% names(newdata)) {
    stop("Previous-state variable `", p_var, "` not found in `newdata`.")
  }

  rows_per_visit <- n_pat * length(non_absorb_idx)
  if (rows_per_visit > getOption("markov.misc.max_batched_design_rows", 2e6)) {
    return(markov_msm_build(
      model = model,
      newdata = newdata,
      time_covariates = time_covariates,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      gap_var = gap_var
    ))
  }

  set_time <- function(d, time_idx) {
    assign_sop_visit(
      d,
      time_var = time_var,
      times = times,
      index = time_idx,
      time_covariates = time_covariates,
      gap_var = gap_var,
      time_info = time_info
    )
  }

  d_init <- set_time(newdata, 1L)
  d_init[[p_var]] <- normalize_previous_state_column(
    d_init[[p_var]],
    newdata[[p_var]],
    p_var
  )

  tt <- if (inherits(model, "vglm")) {
    stats::delete.response(stats::terms(model))
  } else {
    NULL
  }
  get_X <- function(d) {
    X <- if (inherits(model, "orm")) {
      orm_model_matrix(model, d, include_intercept = TRUE)
    } else {
      stats::model.matrix(
        tt,
        data = d,
        contrasts.arg = if (length(model@contrasts)) model@contrasts else NULL,
        xlev = model@xlevels
      )
    }
    missing_cols <- setdiff(colnames(Gamma), colnames(X))
    if (length(missing_cols) > 0L) {
      stop(
        "Design matrix columns missing during batched build: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    X[, colnames(Gamma), drop = FALSE]
  }

  X_init <- get_X(d_init)
  X_transition <- vector("list", length(times))
  if (length(times) >= 2L) {
    for (time_idx in 2:length(times)) {
      d_time <- set_time(newdata, time_idx)
      d_visit <- d_time[
        rep(seq_len(n_pat), times = length(non_absorb_idx)),
        ,
        drop = FALSE
      ]
      d_visit[[p_var]] <- make_previous_state_column(
        states = ylevel_names[non_absorb_idx],
        prototype = newdata[[p_var]],
        n = n_pat,
        p_var = p_var
      )
      X_transition[[time_idx]] <- get_X(d_visit)
    }
  }

  design_bytes <- measure_execution_plan_designs(X_init, X_transition)
  validate_execution_plan_designs(design_bytes, budget$workspace_bytes)

  list(
    X_init = X_init,
    X_transition = X_transition,
    n_pat = n_pat,
    n_states = length(ylevel_names),
    M = length(ylevel_names) - 1L,
    y_levels = ylevel_names,
    col_names = colnames(X_init),
    transition_layout = "stacked",
    transition_origins = non_absorb_idx,
    workspace_bytes = budget$workspace_bytes,
    design_bytes = design_bytes
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
  y_levels <- components$y_levels
  col_names <- components$col_names
  n_times <- length(times)
  validate_sop_threshold_count(Gamma, y_levels)

  # Align Gamma to components
  Gamma <- Gamma[, col_names, drop = FALSE]
  Gamma_t <- t(Gamma) # [P x M]

  # Helper: X * Gamma_t -> LP [N x M]
  calc_lp <- function(X) {
    X %*% Gamma_t
  }

  absorb_idx <- if (!is.null(absorb)) {
    which(as.character(y_levels) %in% as.character(absorb))
  } else {
    integer(0)
  }
  non_absorb_idx <- setdiff(1:n_states, absorb_idx)

  po <- markov_po_structure(Gamma, col_names, X_init)
  initial <- if (is.null(po)) {
    lp_to_probs(calc_lp(X_init), M)
  } else {
    scalar <- drop(X_init %*% po$beta)
    lp_to_probs(outer(scalar, po$cutpoints, "+"), M)
  }
  check_transition_probabilities(initial, "the first SOP time point")
  P_out <- array(0, dim = c(n_pat, n_times, n_states))
  P_out[, 1L, ] <- initial
  previous <- initial

  if (n_times >= 2L) {
    for (t_idx in 2:n_times) {
      X_visit <- X_transition[[t_idx]]
      if (!identical(components$transition_layout, "stacked")) {
        X_visit <- do.call(rbind, X_visit[non_absorb_idx])
      } else {
        built_origins <- components$transition_origins %||% non_absorb_idx
        origin_positions <- match(non_absorb_idx, built_origins)
        if (anyNA(origin_positions)) {
          stop("Transition design is missing a required origin state.")
        }
        if (!identical(origin_positions, seq_along(built_origins))) {
          block_rows <- unlist(lapply(origin_positions, function(position) {
            start <- (position - 1L) * n_pat + 1L
            start:(start + n_pat - 1L)
          }))
          X_visit <- X_visit[block_rows, , drop = FALSE]
        }
      }
      if (!is.matrix(X_visit)) {
        stop("Missing transition design for SOP visit ", t_idx, ".")
      }
      if (is.null(po)) {
        logits <- calc_lp(X_visit)
        check_transition_probabilities(logits, paste0("SOP time point ", t_idx))
        previous <- markov_update_logits_native(
          previous,
          logits,
          as.integer(non_absorb_idx),
          as.integer(absorb_idx)
        )
      } else {
        scalar <- drop(X_visit %*% po$beta)
        check_transition_probabilities(scalar, paste0("SOP time point ", t_idx))
        previous <- markov_update_po_native(
          previous,
          scalar,
          po$cutpoints,
          as.integer(non_absorb_idx),
          as.integer(absorb_idx)
        )
      }
      P_out[, t_idx, ] <- previous
    }
  }
  dimnames(P_out) <- list(NULL, NULL, y_levels)
  P_out
}

markov_po_structure <- function(Gamma, col_names, X_init) {
  intercept <- grep("Intercept", col_names, fixed = TRUE)
  if (length(intercept) != 1L || any(X_init[, intercept] != 1)) {
    return(NULL)
  }
  slopes <- setdiff(seq_along(col_names), intercept)
  if (
    length(slopes) &&
      any(vapply(
        slopes,
        function(j) !all(Gamma[, j] == Gamma[1L, j]),
        logical(1)
      ))
  ) {
    return(NULL)
  }
  beta <- Gamma[1L, ]
  beta[intercept] <- 0
  list(
    beta = beta,
    cutpoints = as.numeric(Gamma[, intercept]),
    intercept = intercept
  )
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
  normalize_probability_array_native(probs)
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
