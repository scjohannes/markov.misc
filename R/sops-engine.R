# SOP Markov recursion engines.

#' Calculate State Occupation Probabilities for First-Order Markov Models
#'
#' Estimates state occupation probabilities over time by iterating a transition
#' matrix derived from a fitted model object (e.g., `vglm`, `rms`, or `rmsb`).
#' This function supports linear time iteration as well as complex, non-linear
#' time specifications (e.g., splines) via a covariate lookup table.
#'
#' @param object A fitted model object. Supported classes include:
#'   \code{"lrm"}, \code{"orm"} (from package `rms`),
#'   \code{"blrm"} (from package `rmsb`), and
#'   \code{"vglm"}, \code{"vgam"} (from package `VGAM`).
#' @param data A data frame containing the baseline covariates for the prediction.
#'   Rows represent unique patients. Columns must contain baseline covariates and
#'   initial values for time-varying variables.
#' @param times Visit-scale time points to iterate over. For numeric time
#'   variables this is a numeric vector. For factor-valued visit indices,
#'   values are matched to the fitted factor levels; if `NULL`, all fitted
#'   visit levels are used.
#' @param ylevels A character vector defining the names of the outcome levels (states).
#'   These must match the levels used in the fitted `object`.
#' @param absorb (Optional) A character vector of absorbing states (states from which
#'   transitions out are impossible). Defaults to \code{NULL}.
#' @param tvarname A character string specifying the name of the time variable in the
#'   model formula. Defaults to \code{"time"}.
#' @param pvarname A character string specifying the name of the previous state variable
#'   used in the model formula. Defaults to \code{"yprev"}.
#' @param p2varname Optional character string specifying the second previous
#'   state variable. `NULL` uses first-order recursion; any non-`NULL` value
#'   uses second-order recursion with histories `(p2varname, pvarname)`.
#' @param gap (Optional) A character string specifying the name of the variable representing
#'   the time gap (delta time) between observations, if used in the model. Defaults to \code{NULL}.
#' @param t_covs (Optional) A data frame used for explicit time-basis columns
#'   or other time-varying covariates.
#'   \itemize{
#'     \item **Structure:** The number of rows in \code{t_covs} must exactly match the length of \code{times}.
#'     \item **Columns:** Column names must match the specific basis variables used in the model formula (e.g., \code{t1}, \code{t2}).
#'     \item **Usage:** At step \code{i}, the values from the \code{i}-th row of \code{t_covs} are injected into the prediction data.
#'   }
#'   Inline model terms such as \code{rms::rcs(time, 4)} do not require
#'   \code{t_covs}; prediction reuses the fitted model's stored transform
#'   metadata.
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions. This does not
#'   integrate over the random-effects distribution; it uses the fitted effects
#'   for known IDs.
#' @param id_var Character ID column used when `include_re = TRUE`. If `NULL`
#'   for `blrm`, inferred from `model$clusterInfo$name` when available,
#'   otherwise `"id"`.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible `blrm` draw sampling.
#' @param ... Reserved for internal use.
#'
#' @details
#'
#' \strong{1. Data Expansion:}
#' We construct a "long" expansion dataset at every time point. For \eqn{N} patients
#' and \eqn{K} non-absorbing states, the expansion dataset contains \eqn{N \times K} rows.
#' \itemize{
#'   \item Rows \eqn{1 \dots N}: All patients assuming \eqn{y_{prev} = \text{State } 1}
#'   \item Rows \eqn{(N+1) \dots 2N}: All patients assuming \eqn{y_{prev} = \text{State } 2}
#'   \item ... and so on.
#' }
#' This allows a single \code{predict()} call to generate transition probabilities for the entire
#' cohort for all possible previous states in one step.
#'
#' \strong{2. Element-wise weighted sum}
#' We use an element-wise weighted sum approach based on the
#' Law of Total Probability:
#' \deqn{P(S_t = k) = \sum_{j} P(S_{t-1} = j) \times P(S_t = k | S_{t-1} = j)}
#'
#' where:
#' \itemize{
#'   \item \eqn{S_t}: State occupied at time \eqn{t}.
#'   \item \eqn{S_{t-1}}: State occupied at time \eqn{t-1}.
#'   \item \eqn{k}: The target state at time \eqn{t}.
#'   \item \eqn{j}: The origin state at time \eqn{t-1}.
#' }
#'
#' Let \eqn{\mathbf{v}_{prev, j}} be a vector of length \eqn{N} containing the probability that each
#' patient was in state \eqn{j} at time \eqn{t-1}.
#' Let \eqn{\mathbf{M}_{trans, j \to k}} be a vector of length \eqn{N} containing the transition
#' probability from \eqn{j} to \eqn{k} for each patient (derived from the batched prediction).
#'
#' The probability of being in state \eqn{k} at time \eqn{t} is updated as:
#' \deqn{\mathbf{v}_{curr, k} = \sum_{j} (\mathbf{v}_{prev, j} \odot \mathbf{M}_{trans, j \to k})}
#' where \eqn{\odot} denotes element-wise multiplication.
#'
#' \strong{3. Absorbing States:}
#' If an absorbing state \eqn{a} is present, the update logic handles the accumulation of probability mass:
#' \deqn{P(S_t = a) = P(S_{t-1} = a) + \sum_{j \neq a} P(S_{t-1} = j) \times P(S_t = a | S_{t-1} = j)}
#'
#' @return
#' An array of state probabilities.
#' \itemize{
#'   \item **Frequentist fit:** An array of dimension \code{[n_patients x n_times x n_states]}.
#'   \item **Bayesian fit:** An array of dimension \code{[n_draws x n_patients x n_times x n_states]}.
#' }
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
#' baseline <- markov_data[!duplicated(markov_data$id), , drop = FALSE]
#' probabilities <- soprob_markov(
#'   fit,
#'   data = baseline,
#'   times = 1:3,
#'   ylevels = fit$yunique,
#'   absorb = "8"
#' )
#' dim(probabilities)
#'
#' @export
soprob_markov <- function(
  object,
  data,
  times = NULL,
  ylevels,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  p2varname = NULL,
  gap = NULL,
  t_covs = NULL,
  include_re = FALSE,
  id_var = NULL,
  n_draws = 100L,
  seed = NULL,
  ...
) {
  # --- 1. Initial Checks & Setup ---
  dots <- list(...)
  unknown_dots <- setdiff(names(dots), c(".draw_indices", ".gamma_draws"))
  if (length(unknown_dots) > 0) {
    stop("Unused arguments: ", paste(unknown_dots, collapse = ", "))
  }
  draw_indices_arg <- dots$.draw_indices
  gamma_draws_arg <- dots$.gamma_draws

  cl <- if (inherits(object, "blrm")) {
    "blrm"
  } else if (inherits(object, "robcov_vglm")) {
    "robcov_vglm"
  } else if (inherits(object, "orm")) {
    "orm"
  } else if (inherits(object, "vglm")) {
    "vglm"
  } else if (inherits(object, "vgam")) {
    "vgam"
  } else {
    class(object)[1]
  }
  ftypes <- c(
    orm = "rms",
    blrm = "rmsb",
    vglm = "vgam",
    vgam = "vgam",
    robcov_vglm = "robcov"
  )
  ftype <- ftypes[cl]

  if (is.na(ftype)) {
    stop("Object class not supported")
  }

  # Validate model is compatible with Markov simulation
  validate_markov_model(object)

  if (!pvarname %in% names(data)) {
    stop("Previous-state variable `", pvarname, "` not found in `data`.")
  }
  if (!is.null(p2varname) && !p2varname %in% names(data)) {
    stop("Second previous-state variable `", p2varname, "` not found in `data`.")
  }
  time_res <- resolve_sop_times(
    object,
    data,
    times,
    tvarname,
    t_covs = t_covs,
    default = "unique"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap, t_covs, time_info)

  draw_indices <- NULL
  if (ftype == "rmsb") {
    draw_indices <- draw_indices_arg %||% select_posterior_draws(object, n_draws, seed)
  }

  # Define prediction function
  prd <- switch(
    ftype,
    rms = function(obj, d) predict_orm_response_markov(obj, d),
    vgam = function(obj, d) predict_vglm_response_markov(obj, d),
    rmsb = function(obj, d) {
      predict_blrm_response_markov(
        obj,
        d,
        include_re = include_re,
        id_var = id_var,
        draw_indices = draw_indices,
        gamma_draws = gamma_draws_arg
      )
    },
    robcov = function(obj, d) {
      if (is.null(obj$vglm_fit)) {
        stop(
          "robcov_vglm object does not contain the original vglm fit. ",
          "Please re-run robcov_vglm() with the latest version of the package."
        )
      }
      predict_vglm_response_markov(obj$vglm_fit, d)
    }
  )

  # Prepare dimensions
  n_pat <- nrow(data)
  n_times <- length(times)
  ylevel_names <- as_state_labels(ylevels)
  absorb_names <- as_state_labels(absorb)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)
  yna <- ylevel_names[non_absorb_idx] # Non-absorbing states
  n_yna <- length(yna)

  # Check Bayesian draws
  nd <- if (ftype == "rmsb" && length(object$draws)) length(draw_indices) else 0

  # Initialize Output Array
  # Structure: [Patients, Time, States]
  if (nd == 0) {
    P <- array(
      0,
      dim = c(n_pat, n_times, n_states),
      dimnames = list(rownames(data), times, ylevel_names)
    )
  } else {
    P <- array(
      0,
      dim = c(nd, n_pat, n_times, n_states),
      dimnames = list(draw_indices, rownames(data), times, ylevel_names)
    )
  }

  # --- 2. Time 1 Initialization (Start) ---
  # Update time variables for T1
  data <- assign_sop_visit(
    data,
    tvarname = tvarname,
    times = times,
    index = 1L,
    t_covs = t_covs,
    gap = gap,
    time_info = time_info
  )

  # Predict probabilities at T1
  p_t1 <- prd(object, data) # Returns [n_pat x n_states] or [nd x n_pat x n_states]

  if (nd == 0) {
    P[, 1, ] <- p_t1
  } else {
    P[,, 1, ] <- p_t1
  }

  if (!is.null(p2varname)) {
    return(soprob_markov_second_order_run(
      P = P,
      object = object,
      data = data,
      prd = prd,
      nd = nd,
      times = times,
      ylevel_names = ylevel_names,
      absorb_names = absorb_names,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      time_info = time_info
    ))
  }

  if (n_times < 2) {
    return(P)
  }

  # --- 3. Prepare Expansion Data for Transitions ---
  # We create a long dataframe where every patient is repeated for every possible PREVIOUS state.
  # Order: Patient 1 (State A), Patient 2 (State A)... Patient 1 (State B)...
  # This ordering is crucial for the vectorized update logic later.

  edata_base <- data[rep(1:n_pat, times = n_yna), , drop = FALSE]

  # Assign the previous state variable
  # We repeat each state n_pat times

  edata_base[[pvarname]] <- make_previous_state_column(
    states = yna,
    prototype = data[[pvarname]],
    n = n_pat,
    pvarname = pvarname
  )

  # --- 4. Iterate Through Time (Vectorized over Patients) ---
  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      tvarname = tvarname,
      times = times,
      index = it,
      t_covs = t_covs,
      gap = gap,
      time_info = time_info
    )

    # Get Transition Probabilities
    # This returns matrix: [ (n_pat * n_yna) x n_states ]
    # Rows 1:N are transitions from State 1, Rows N+1:2N from State 2, etc.
    trans_probs <- prd(object, edata_base)

    # --- 5. The Update Step (Markov) ---
    # Formula: P(S_t = k) = Sum_over_j [ P(S_t-1 = j) * P(S_t=k | S_t-1=j) ]

    if (nd == 0) {
      # Initialize current time probabilities with 0
      p_current <- matrix(0, nrow = n_pat, ncol = n_states)
      colnames(p_current) <- ylevel_names

      # Add contribution from Non-Absorbing Previous States
      for (i in seq_along(yna)) {
        prev_state_name <- yna[i]

        # Extract prob of being in this previous state for all patients
        # Vector of length n_pat
        prob_prev <- P[, it - 1, prev_state_name]

        # Extract transition probs GIVEN this previous state
        # Rows correspond to the block for this state in edata_base
        row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
        probs_transition <- trans_probs[row_indices, ]

        # Weighted sum: multiply column-wise by the probability of being in prev state
        p_current <- p_current + (probs_transition * prob_prev)
      }

      # Add contribution from Absorbing States (if any)
      # Absorbing states transition to themselves with Prob 1
      if (length(absorb_idx) > 0) {
        for (a_state in absorb_idx) {
          # If you were in absorb state at t-1, you are in absorb state at t
          p_current[, a_state] <- p_current[, a_state] + P[, it - 1, a_state]
        }
      }

      P[, it, ] <- p_current
    } else {
      # trans_probs is [draw x (patient * previous state) x state].
      p_current <- array(0, dim = c(nd, n_pat, n_states))
      dimnames(p_current) <- list(dimnames(P)[[1]], rownames(data), ylevel_names)

      for (i in seq_along(yna)) {
        prev_state_name <- yna[i]
        prob_prev <- P[, , it - 1, prev_state_name]
        row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
        probs_transition <- trans_probs[, row_indices, , drop = FALSE]

        for (k in seq_len(n_states)) {
          p_current[, , k] <- p_current[, , k] +
            probs_transition[, , k] * prob_prev
        }
      }

      if (length(absorb_idx) > 0) {
        for (a_state in absorb_idx) {
          # total dead at t = new deaths + already dead
          p_current[, , a_state] <- p_current[, , a_state] +
            P[, , it - 1, a_state]
        }
      }
      P[, , it, ] <- p_current
    }
  }

  return(P)
}

match_state_indices <- function(values, ylevel_names, varname) {
  idx <- match(as.character(values), ylevel_names)
  if (anyNA(idx)) {
    bad <- unique(as.character(values[is.na(idx)]))
    stop(
      "Values in `", varname, "` are not among `ylevels`: ",
      paste(utils::head(bad, 5), collapse = ", "),
      if (length(bad) > 5) " ..." else ""
    )
  }
  idx
}

soprob_markov_second_order_run <- function(
  P,
  object,
  data,
  prd,
  nd,
  times,
  ylevel_names,
  absorb_names,
  tvarname,
  pvarname,
  p2varname,
  gap,
  t_covs,
  time_info
) {
  n_pat <- nrow(data)
  n_times <- length(times)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)
  prev_idx <- match_state_indices(data[[pvarname]], ylevel_names, pvarname)

  if (nd == 0) {
    joint_prev <- array(
      0,
      dim = c(n_pat, n_states, n_states),
      dimnames = list(rownames(data), ylevel_names, ylevel_names)
    )
    for (i in seq_len(n_pat)) {
      if (prev_idx[i] %in% absorb_idx) {
        joint_prev[i, prev_idx[i], prev_idx[i]] <- 1
        P[i, 1, ] <- 0
        P[i, 1, prev_idx[i]] <- 1
      } else {
        joint_prev[i, prev_idx[i], ] <- P[i, 1, ]
      }
    }
  } else {
    joint_prev <- array(
      0,
      dim = c(nd, n_pat, n_states, n_states),
      dimnames = list(dimnames(P)[[1]], rownames(data), ylevel_names, ylevel_names)
    )
    for (i in seq_len(n_pat)) {
      if (prev_idx[i] %in% absorb_idx) {
        joint_prev[, i, prev_idx[i], prev_idx[i]] <- 1
        P[, i, 1, ] <- 0
        P[, i, 1, prev_idx[i]] <- 1
      } else {
        joint_prev[, i, prev_idx[i], ] <- P[, i, 1, ]
      }
    }
  }

  if (n_times < 2) {
    return(P)
  }

  pair_grid <- expand.grid(
    h = seq_len(n_states),
    j = seq_len(n_states),
    KEEP.OUT.ATTRS = FALSE
  )
  n_pairs <- nrow(pair_grid)
  predictable_pair <- pair_grid$h %in% non_absorb_idx &
    pair_grid$j %in% non_absorb_idx
  block_rows <- lapply(seq_len(n_pairs), function(i) {
    ((i - 1) * n_pat + 1):(i * n_pat)
  })

  edata_base <- data[rep(seq_len(n_pat), times = n_pairs), , drop = FALSE]
  edata_base[[p2varname]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$h],
    prototype = data[[p2varname]],
    n = n_pat,
    pvarname = p2varname
  )
  edata_base[[pvarname]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$j],
    prototype = data[[pvarname]],
    n = n_pat,
    pvarname = pvarname
  )

  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      tvarname = tvarname,
      times = times,
      index = it,
      t_covs = t_covs,
      gap = gap,
      time_info = time_info
    )

    predict_rows <- unlist(block_rows[predictable_pair], use.names = FALSE)
    trans_probs <- prd(object, edata_base[predict_rows, , drop = FALSE])
    cursor <- 1L

    if (nd == 0) {
      joint_current <- array(0, dim = c(n_pat, n_states, n_states))
      for (pair_i in seq_len(n_pairs)) {
        h <- pair_grid$h[pair_i]
        j <- pair_grid$j[pair_i]
        prob_prev <- joint_prev[, h, j]
        if (!any(prob_prev > 0)) {
          if (predictable_pair[pair_i]) {
            cursor <- cursor + n_pat
          }
          next
        }
        if (j %in% absorb_idx) {
          joint_current[, j, j] <- joint_current[, j, j] + prob_prev
        } else if (predictable_pair[pair_i]) {
          rows <- cursor:(cursor + n_pat - 1L)
          transition <- trans_probs[rows, , drop = FALSE]
          for (l in seq_len(n_states)) {
            joint_current[, j, l] <- joint_current[, j, l] + transition[, l] * prob_prev
          }
          cursor <- cursor + n_pat
        }
      }
      P[, it, ] <- apply(joint_current, c(1, 3), sum)
      joint_prev <- joint_current
    } else {
      joint_current <- array(0, dim = c(nd, n_pat, n_states, n_states))
      for (pair_i in seq_len(n_pairs)) {
        h <- pair_grid$h[pair_i]
        j <- pair_grid$j[pair_i]
        if (j %in% absorb_idx) {
          joint_current[, , j, j] <- joint_current[, , j, j] +
            joint_prev[, , h, j]
        } else if (predictable_pair[pair_i]) {
          rows <- cursor:(cursor + n_pat - 1L)
          prob_prev <- joint_prev[, , h, j]
          if (any(prob_prev > 0)) {
            transition <- trans_probs[, rows, , drop = FALSE]
            for (l in seq_len(n_states)) {
              joint_current[, , j, l] <- joint_current[, , j, l] +
                transition[, , l] * prob_prev
            }
          }
          cursor <- cursor + n_pat
        }
      }
      P[, , it, ] <- apply(joint_current, c(1, 2, 4), sum)
      joint_prev <- joint_current
    }
  }

  P
}
