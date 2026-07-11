# Internal transition traces for model-based diagnostic plots.

markov_transition_trace <- function(
  object,
  data,
  times = NULL,
  y_levels,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  include_re = FALSE,
  id_var = NULL,
  n_draws = 100L,
  seed = NULL,
  return_kernels = FALSE,
  ...
) {
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

  validate_markov_model(object)

  if (!p_var %in% names(data)) {
    stop("Previous-state variable `", p_var, "` not found in `data`.")
  }
  if (!is.null(p2_var) && !p2_var %in% names(data)) {
    stop(
      "Second previous-state variable `",
      p2_var,
      "` not found in `data`."
    )
  }

  time_res <- resolve_sop_times(
    object,
    data,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "unique"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap_var, time_covariates, time_info)

  draw_indices <- NULL
  if (ftype == "rmsb") {
    draw_indices <- draw_indices_arg %||%
      select_posterior_draws(object, n_draws, seed)
  }

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

  n_pat <- nrow(data)
  n_times <- length(times)
  ylevel_names <- as_state_labels(y_levels)
  absorb_names <- as_state_labels(absorb)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  nd <- if (ftype == "rmsb" && length(object$draws)) {
    length(draw_indices)
  } else {
    0L
  }

  kernels <- NULL
  if (nd == 0) {
    P <- array(
      0,
      dim = c(n_pat, n_times, n_states),
      dimnames = list(rownames(data), times, ylevel_names)
    )
    transitions <- array(
      0,
      dim = c(n_pat, n_times, n_states, n_states),
      dimnames = list(rownames(data), times, ylevel_names, ylevel_names)
    )
    if (isTRUE(return_kernels)) {
      kernels <- if (is.null(p2_var)) {
        array(
          0,
          dim = c(n_pat, n_times, n_states, n_states),
          dimnames = list(rownames(data), times, ylevel_names, ylevel_names)
        )
      } else {
        array(
          0,
          dim = c(n_pat, n_times, n_states, n_states, n_states),
          dimnames = list(
            rownames(data),
            times,
            ylevel_names,
            ylevel_names,
            ylevel_names
          )
        )
      }
    }
  } else {
    P <- array(
      0,
      dim = c(nd, n_pat, n_times, n_states),
      dimnames = list(draw_indices, rownames(data), times, ylevel_names)
    )
    transitions <- array(
      0,
      dim = c(nd, n_pat, n_times, n_states, n_states),
      dimnames = list(
        draw_indices,
        rownames(data),
        times,
        ylevel_names,
        ylevel_names
      )
    )
    if (isTRUE(return_kernels)) {
      kernels <- if (is.null(p2_var)) {
        array(
          0,
          dim = c(nd, n_pat, n_times, n_states, n_states),
          dimnames = list(
            draw_indices,
            rownames(data),
            times,
            ylevel_names,
            ylevel_names
          )
        )
      } else {
        array(
          0,
          dim = c(nd, n_pat, n_times, n_states, n_states, n_states),
          dimnames = list(
            draw_indices,
            rownames(data),
            times,
            ylevel_names,
            ylevel_names,
            ylevel_names
          )
        )
      }
    }
  }

  data <- assign_sop_visit(
    data,
    time_var = time_var,
    times = times,
    index = 1L,
    time_covariates = time_covariates,
    gap_var = gap_var,
    time_info = time_info
  )

  prev_idx <- match_state_indices(data[[p_var]], ylevel_names, p_var)
  predict_rows <- !prev_idx %in% absorb_idx

  if (nd == 0) {
    if (any(predict_rows)) {
      p_t1 <- prd(object, data[predict_rows, , drop = FALSE])
      check_transition_probabilities(p_t1, "the first transition time point")
      P[predict_rows, 1, ] <- p_t1
      transitions <- trace_record_first_transition(
        transitions = transitions,
        probs = p_t1,
        rows = which(predict_rows),
        prev_idx = prev_idx,
        time_index = 1L
      )
    }
    if (any(!predict_rows)) {
      recorded <- trace_record_absorbing_transition(
        P = P,
        transitions = transitions,
        rows = which(!predict_rows),
        prev_idx = prev_idx,
        time_index = 1L,
        draw = FALSE
      )
      P <- recorded$P
      transitions <- recorded$transitions
    }
  } else {
    if (any(predict_rows)) {
      p_t1 <- prd(object, data[predict_rows, , drop = FALSE])
      check_transition_probabilities(p_t1, "the first transition time point")
      P[, predict_rows, 1, ] <- p_t1
      transitions <- trace_record_first_transition(
        transitions = transitions,
        probs = p_t1,
        rows = which(predict_rows),
        prev_idx = prev_idx,
        time_index = 1L
      )
    }
    if (any(!predict_rows)) {
      recorded <- trace_record_absorbing_transition(
        P = P,
        transitions = transitions,
        rows = which(!predict_rows),
        prev_idx = prev_idx,
        time_index = 1L,
        draw = TRUE
      )
      P <- recorded$P
      transitions <- recorded$transitions
    }
  }

  if (!is.null(p2_var)) {
    return(markov_transition_trace_second_order(
      P = P,
      transitions = transitions,
      kernels = kernels,
      object = object,
      data = data,
      prd = prd,
      nd = nd,
      times = times,
      ylevel_names = ylevel_names,
      absorb_names = absorb_names,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      time_info = time_info,
      draw_indices = draw_indices
    ))
  }

  markov_transition_trace_first_order(
    P = P,
    transitions = transitions,
    kernels = kernels,
    object = object,
    data = data,
    prd = prd,
    nd = nd,
    times = times,
    ylevel_names = ylevel_names,
    absorb_names = absorb_names,
    time_var = time_var,
    p_var = p_var,
    gap_var = gap_var,
    time_covariates = time_covariates,
    time_info = time_info,
    draw_indices = draw_indices
  )
}

trace_record_first_transition <- function(
  transitions,
  probs,
  rows,
  prev_idx,
  time_index
) {
  if (length(dim(probs)) == 2L) {
    for (pos in seq_along(rows)) {
      row <- rows[pos]
      transitions[row, time_index, prev_idx[row], ] <- probs[pos, ]
    }
    return(transitions)
  }

  for (pos in seq_along(rows)) {
    row <- rows[pos]
    transitions[, row, time_index, prev_idx[row], ] <- probs[, pos, ]
  }
  transitions
}

trace_record_absorbing_transition <- function(
  P,
  transitions,
  rows,
  prev_idx,
  time_index,
  draw
) {
  for (row in rows) {
    state <- prev_idx[row]
    if (isTRUE(draw)) {
      P[, row, time_index, state] <- 1
      transitions[, row, time_index, state, state] <- 1
    } else {
      P[row, time_index, state] <- 1
      transitions[row, time_index, state, state] <- 1
    }
  }
  list(P = P, transitions = transitions)
}

markov_transition_trace_result <- function(
  P,
  transitions,
  draw_indices,
  kernels = NULL
) {
  out <- list(
    sops = P,
    transitions = transitions,
    kernels = kernels,
    draw_indices = draw_indices
  )
  class(out) <- "markov_transition_trace"
  out
}

markov_transition_trace_first_order <- function(
  P,
  transitions,
  kernels,
  object,
  data,
  prd,
  nd,
  times,
  ylevel_names,
  absorb_names,
  time_var,
  p_var,
  gap_var,
  time_covariates,
  time_info,
  draw_indices
) {
  n_pat <- nrow(data)
  n_times <- length(times)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)
  yna <- ylevel_names[non_absorb_idx]

  if (n_times < 2) {
    return(markov_transition_trace_result(
      P,
      transitions,
      draw_indices,
      kernels
    ))
  }

  edata_base <- data[rep(seq_len(n_pat), times = length(yna)), , drop = FALSE]
  edata_base[[p_var]] <- make_previous_state_column(
    states = yna,
    prototype = data[[p_var]],
    n = n_pat,
    p_var = p_var
  )

  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      time_var = time_var,
      times = times,
      index = it,
      time_covariates = time_covariates,
      gap_var = gap_var,
      time_info = time_info
    )

    trans_probs <- prd(object, edata_base)
    check_transition_probabilities(
      trans_probs,
      paste0("transition time point ", it)
    )

    if (nd == 0) {
      transition_current <- array(
        0,
        dim = c(n_pat, n_states, n_states),
        dimnames = list(rownames(data), ylevel_names, ylevel_names)
      )
      for (i in seq_along(yna)) {
        prev_state_name <- yna[i]
        prev_state_idx <- match(prev_state_name, ylevel_names)
        prob_prev <- P[, it - 1, prev_state_name]
        row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
        probs_transition <- trans_probs[row_indices, ]
        if (!is.null(kernels)) {
          kernels[, it, prev_state_idx, ] <- probs_transition
        }
        transition_current[, prev_state_idx, ] <-
          probs_transition * prob_prev
      }
      for (a_state in absorb_idx) {
        transition_current[, a_state, a_state] <- P[, it - 1, a_state]
        if (!is.null(kernels)) {
          kernels[, it, a_state, a_state] <- 1
        }
      }
      transitions[, it, , ] <- transition_current
      P[, it, ] <- apply(transition_current, c(1, 3), sum)
    } else {
      transition_current <- array(
        0,
        dim = c(nd, n_pat, n_states, n_states),
        dimnames = list(
          dimnames(P)[[1]],
          rownames(data),
          ylevel_names,
          ylevel_names
        )
      )
      for (i in seq_along(yna)) {
        prev_state_name <- yna[i]
        prev_state_idx <- match(prev_state_name, ylevel_names)
        prob_prev <- P[,, it - 1, prev_state_name]
        row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
        probs_transition <- trans_probs[, row_indices, , drop = FALSE]
        if (!is.null(kernels)) {
          kernels[,, it, prev_state_idx, ] <- probs_transition
        }
        for (k in seq_len(n_states)) {
          transition_current[,, prev_state_idx, k] <-
            probs_transition[,, k] * prob_prev
        }
      }
      for (a_state in absorb_idx) {
        transition_current[,, a_state, a_state] <- P[,, it - 1, a_state]
        if (!is.null(kernels)) {
          kernels[,, it, a_state, a_state] <- 1
        }
      }
      transitions[,, it, , ] <- transition_current
      P[,, it, ] <- apply(transition_current, c(1, 2, 4), sum)
    }
  }

  markov_transition_trace_result(P, transitions, draw_indices, kernels)
}

markov_transition_trace_second_order <- function(
  P,
  transitions,
  kernels,
  object,
  data,
  prd,
  nd,
  times,
  ylevel_names,
  absorb_names,
  time_var,
  p_var,
  p2_var,
  gap_var,
  time_covariates,
  time_info,
  draw_indices
) {
  n_pat <- nrow(data)
  n_times <- length(times)
  n_states <- length(ylevel_names)
  absorb_idx <- which(ylevel_names %in% absorb_names)
  non_absorb_idx <- setdiff(seq_len(n_states), absorb_idx)

  if (nd == 0) {
    joint_prev <- array(
      transitions[, 1, , ],
      dim = c(n_pat, n_states, n_states),
      dimnames = list(rownames(data), ylevel_names, ylevel_names)
    )
  } else {
    joint_prev <- array(
      transitions[,, 1, , ],
      dim = c(nd, n_pat, n_states, n_states),
      dimnames = list(
        dimnames(P)[[1]],
        rownames(data),
        ylevel_names,
        ylevel_names
      )
    )
  }

  if (n_times < 2) {
    return(markov_transition_trace_result(
      P,
      transitions,
      draw_indices,
      kernels
    ))
  }

  pair_grid <- expand.grid(
    h = seq_len(n_states),
    j = seq_len(n_states),
    KEEP.OUT.ATTRS = FALSE
  )
  n_pairs <- nrow(pair_grid)
  predictable_pair <- pair_grid$h %in%
    non_absorb_idx &
    pair_grid$j %in% non_absorb_idx
  block_rows <- lapply(seq_len(n_pairs), function(i) {
    ((i - 1) * n_pat + 1):(i * n_pat)
  })

  edata_base <- data[rep(seq_len(n_pat), times = n_pairs), , drop = FALSE]
  edata_base[[p2_var]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$h],
    prototype = data[[p2_var]],
    n = n_pat,
    p_var = p2_var
  )
  edata_base[[p_var]] <- make_previous_state_column(
    states = ylevel_names[pair_grid$j],
    prototype = data[[p_var]],
    n = n_pat,
    p_var = p_var
  )

  for (it in 2:n_times) {
    edata_base <- assign_sop_visit(
      edata_base,
      time_var = time_var,
      times = times,
      index = it,
      time_covariates = time_covariates,
      gap_var = gap_var,
      time_info = time_info
    )

    predict_rows <- unlist(block_rows[predictable_pair], use.names = FALSE)
    trans_probs <- prd(object, edata_base[predict_rows, , drop = FALSE])
    check_transition_probabilities(
      trans_probs,
      paste0("second-order transition time point ", it)
    )
    cursor <- 1L

    if (nd == 0) {
      joint_current <- array(
        0,
        dim = c(n_pat, n_states, n_states),
        dimnames = list(rownames(data), ylevel_names, ylevel_names)
      )
      for (pair_i in seq_len(n_pairs)) {
        h <- pair_grid$h[pair_i]
        j <- pair_grid$j[pair_i]
        prob_prev <- joint_prev[, h, j]
        if (j %in% absorb_idx) {
          joint_current[, j, j] <- joint_current[, j, j] + prob_prev
          if (!is.null(kernels)) {
            kernels[, it, h, j, j] <- 1
          }
        } else if (predictable_pair[pair_i]) {
          rows <- cursor:(cursor + n_pat - 1L)
          transition <- trans_probs[rows, , drop = FALSE]
          if (!is.null(kernels)) {
            kernels[, it, h, j, ] <- transition
          }
          if (any(prob_prev > 0)) {
            for (l in seq_len(n_states)) {
              joint_current[, j, l] <- joint_current[, j, l] +
                transition[, l] * prob_prev
            }
          }
          cursor <- cursor + n_pat
        }
      }
      transitions[, it, , ] <- joint_current
      P[, it, ] <- apply(joint_current, c(1, 3), sum)
      joint_prev <- joint_current
    } else {
      joint_current <- array(
        0,
        dim = c(nd, n_pat, n_states, n_states),
        dimnames = list(
          dimnames(P)[[1]],
          rownames(data),
          ylevel_names,
          ylevel_names
        )
      )
      for (pair_i in seq_len(n_pairs)) {
        h <- pair_grid$h[pair_i]
        j <- pair_grid$j[pair_i]
        if (j %in% absorb_idx) {
          joint_current[,, j, j] <- joint_current[,, j, j] +
            joint_prev[,, h, j]
          if (!is.null(kernels)) {
            kernels[,, it, h, j, j] <- 1
          }
        } else if (predictable_pair[pair_i]) {
          rows <- cursor:(cursor + n_pat - 1L)
          prob_prev <- joint_prev[,, h, j]
          transition <- trans_probs[, rows, , drop = FALSE]
          if (!is.null(kernels)) {
            kernels[,, it, h, j, ] <- transition
          }
          if (any(prob_prev > 0)) {
            for (l in seq_len(n_states)) {
              joint_current[,, j, l] <- joint_current[,, j, l] +
                transition[,, l] * prob_prev
            }
          }
          cursor <- cursor + n_pat
        }
      }
      transitions[,, it, , ] <- joint_current
      P[,, it, ] <- apply(joint_current, c(1, 2, 4), sum)
      joint_prev <- joint_current
    }
  }

  markov_transition_trace_result(P, transitions, draw_indices, kernels)
}
