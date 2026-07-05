# Internal correlation-summary data builders for diagnostic plots.

plot_correlation_input_data <- function(
  object,
  newdata,
  refit_data,
  times,
  ylevels,
  absorb,
  id_var,
  time_var,
  y_var,
  pvarname,
  p2varname,
  facet_var,
  gap,
  t_covs,
  seed,
  n_draws,
  triangle
) {
  if (markov_supported_model(object)) {
    return(plot_correlation_model_data(
      model = object,
      newdata = newdata,
      refit_data = refit_data,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      id_var = id_var,
      time_var = time_var,
      pvarname = pvarname,
      p2varname = p2varname,
      facet_var = facet_var,
      gap = gap,
      t_covs = t_covs,
      seed = seed,
      n_draws = n_draws,
      triangle = triangle
    ))
  }

  if (!is.data.frame(object)) {
    stop("`object` must be a data frame or a supported Markov model.")
  }
  if (!is.null(newdata) || !is.null(refit_data)) {
    stop("`newdata` and `refit_data` are only used for model-based plots.")
  }
  if (!is.null(p2varname) || !is.null(gap) || !is.null(t_covs)) {
    stop(
      "`p2varname`, `gap`, and `t_covs` are only used for model-based plots."
    )
  }
  if (!is.null(seed)) {
    stop("`seed` is only used for `blrm` model-based plots.")
  }

  data <- object
  plot_validate_columns(data, c(id_var, time_var, y_var), "`data`")
  if (!is.null(times)) {
    keep <- as.character(data[[time_var]]) %in% as.character(times)
    data <- data[keep, , drop = FALSE]
  }

  plot_correlation_data(
    data = data,
    id_var = id_var,
    time_var = time_var,
    y_var = y_var,
    facet_var = facet_var,
    ylevels = ylevels,
    triangle = triangle
  )
}

plot_correlation_model_data <- function(
  model,
  newdata,
  refit_data,
  times,
  ylevels,
  absorb,
  id_var,
  time_var,
  pvarname,
  p2varname,
  facet_var,
  gap,
  t_covs,
  seed,
  n_draws,
  triangle
) {
  setup <- plot_transition_model_setup(
    model = model,
    newdata = newdata,
    refit_data = refit_data,
    variables = NULL,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    time_var = time_var,
    pvarname = pvarname,
    p2varname = p2varname,
    id_var = id_var,
    gap = gap,
    t_covs = t_covs
  )

  plot_correlation_model_summary(
    model = model,
    setup = setup,
    facet_var = facet_var,
    pvarname = pvarname,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs,
    seed = seed,
    n_draws = n_draws,
    triangle = triangle
  )
}

plot_correlation_model_summary <- function(
  model,
  setup,
  facet_var,
  pvarname,
  p2varname,
  gap,
  t_covs,
  seed,
  n_draws,
  triangle
) {
  second_order <- !is.null(p2varname)

  plot_model_trace_summaries(
    model = model,
    setup = setup,
    facet_var = facet_var,
    pvarname = pvarname,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs,
    seed = seed,
    n_draws = n_draws,
    return_kernels = TRUE,
    summarize_trace = function(trace, draw) {
      plot_correlation_trace_summary(
        trace = trace,
        data = setup$data,
        plot_indices = setup$plot_indices,
        plot_times = setup$plot_times,
        ylevels = setup$ylevels,
        facet_var = facet_var,
        draw = draw,
        second_order = second_order,
        triangle = triangle
      )
    },
    summarize_draws = function(data) {
      plot_correlation_summarize_draws(
        data = data,
        facet_var = facet_var,
        time_keys = as.character(setup$plot_times)
      )
    }
  )
}

plot_correlation_trace_summary <- function(
  trace,
  data,
  plot_indices,
  plot_times,
  ylevels,
  facet_var,
  draw,
  second_order,
  triangle
) {
  groups <- plot_facet_groups(data, facet_var)
  out <- vector("list", nrow(groups))
  time_keys <- as.character(plot_times)
  scores <- seq_along(ylevels)

  for (i in seq_len(nrow(groups))) {
    keep <- plot_group_subset(data, groups[i, , drop = FALSE], facet_var)
    rows <- match(rownames(keep), rownames(data))
    group <- if (is.null(facet_var)) NULL else groups[i, , drop = FALSE]
    out[[i]] <- plot_correlation_trace_group_summary(
      trace = trace,
      rows = rows,
      time_indices = plot_indices,
      time_keys = time_keys,
      scores = scores,
      group = group,
      facet_var = facet_var,
      draw = draw,
      second_order = second_order,
      triangle = triangle
    )
  }

  out <- bind_rows_fill(out)
  plot_correlation_finalize(out, facet_var = facet_var, time_keys = time_keys)
}

plot_correlation_trace_group_summary <- function(
  trace,
  rows,
  time_indices,
  time_keys,
  scores,
  group,
  facet_var,
  draw,
  second_order,
  triangle
) {
  if (isTRUE(draw)) {
    draw_ids <- dimnames(trace$sops)[[1]]
    out <- vector("list", length(draw_ids))
    for (i in seq_along(draw_ids)) {
      kernels_i <- if (isTRUE(second_order)) {
        drop_first_array_dim(trace$kernels[i, , , , , , drop = FALSE])
      } else {
        drop_first_array_dim(trace$kernels[i, , , , , drop = FALSE])
      }
      corr <- plot_correlation_matrix_from_trace(
        sops = drop_first_array_dim(trace$sops[i, , , , drop = FALSE]),
        transitions = drop_first_array_dim(
          trace$transitions[i, , , , , drop = FALSE]
        ),
        kernels = kernels_i,
        rows = rows,
        time_indices = time_indices,
        scores = scores,
        second_order = second_order
      )
      out[[i]] <- plot_correlation_matrix_long(
        corr,
        group = group,
        facet_var = facet_var,
        triangle = triangle
      )
      out[[i]]$draw_id <- draw_ids[i]
    }
    return(bind_rows_fill(out))
  }

  corr <- plot_correlation_matrix_from_trace(
    sops = trace$sops,
    transitions = trace$transitions,
    kernels = trace$kernels,
    rows = rows,
    time_indices = time_indices,
    scores = scores,
    second_order = second_order
  )
  plot_correlation_matrix_long(
    corr,
    group = group,
    facet_var = facet_var,
    triangle = triangle
  )
}

drop_first_array_dim <- function(x) {
  array(
    x,
    dim = dim(x)[-1],
    dimnames = dimnames(x)[-1]
  )
}

plot_correlation_matrix_from_trace <- function(
  sops,
  transitions,
  kernels,
  rows,
  time_indices,
  scores,
  second_order
) {
  n_plot_times <- length(time_indices)
  time_keys <- dimnames(sops)[[2]][time_indices]
  corr <- matrix(
    NA_real_,
    nrow = n_plot_times,
    ncol = n_plot_times,
    dimnames = list(time_keys, time_keys)
  )
  if (length(rows) == 0L) {
    return(corr)
  }

  means <- plot_trace_state_means(sops, rows, time_indices, scores)
  variances <- means$mean2 - means$mean^2
  second_order_cross <- if (isTRUE(second_order)) {
    plot_second_order_cross_moment_matrix(
      transitions = transitions,
      kernels = kernels,
      rows = rows,
      time_indices = time_indices,
      scores = scores
    )
  } else {
    NULL
  }

  for (i in seq_len(n_plot_times)) {
    corr[i, i] <- 1
    if (i == n_plot_times) {
      next
    }
    for (j in seq.int(i + 1L, n_plot_times)) {
      cross <- if (isTRUE(second_order)) {
        second_order_cross[i, j]
      } else {
        plot_first_order_cross_moment(
          sops = sops,
          kernels = kernels,
          rows = rows,
          start = time_indices[i],
          end = time_indices[j],
          scores = scores
        )
      }
      denom <- sqrt(variances[i] * variances[j])
      value <- if (is.finite(denom) && denom > 0) {
        (cross - means$mean[i] * means$mean[j]) / denom
      } else {
        NA_real_
      }
      corr[i, j] <- value
      corr[j, i] <- value
    }
  }

  corr
}

plot_trace_state_means <- function(sops, rows, time_indices, scores) {
  mean <- numeric(length(time_indices))
  mean2 <- numeric(length(time_indices))
  scores2 <- scores^2
  for (i in seq_along(time_indices)) {
    prob <- sops[rows, time_indices[i], , drop = FALSE]
    prob <- matrix(prob, nrow = length(rows), ncol = length(scores))
    mean[i] <- mean(as.vector(prob %*% scores), na.rm = TRUE)
    mean2[i] <- mean(as.vector(prob %*% scores2), na.rm = TRUE)
  }
  list(mean = mean, mean2 = mean2)
}

plot_first_order_cross_moment <- function(
  sops,
  kernels,
  rows,
  start,
  end,
  scores
) {
  values <- numeric(length(rows))
  for (i in seq_along(rows)) {
    row <- rows[i]
    future <- scores
    for (time in seq.int(end, start + 1L)) {
      future <- as.numeric(kernels[row, time, , ] %*% future)
    }
    values[i] <- sum(sops[row, start, ] * scores * future)
  }
  mean(values, na.rm = TRUE)
}

plot_second_order_cross_moment_matrix <- function(
  transitions,
  kernels,
  rows,
  time_indices,
  scores
) {
  n_plot_times <- length(time_indices)
  out <- matrix(NA_real_, nrow = n_plot_times, ncol = n_plot_times)
  if (n_plot_times < 2L) {
    return(out)
  }

  for (start_pos in seq_len(n_plot_times - 1L)) {
    start <- time_indices[start_pos]
    end_positions <- seq.int(start_pos + 1L, n_plot_times)
    max_end <- max(time_indices[end_positions])
    anchor_pair <- plot_second_order_anchor_pair_array(
      pair_history = transitions[rows, start, , , drop = FALSE],
      n_states = length(scores)
    )
    for (time in seq.int(start + 1L, max_end)) {
      anchor_pair <- plot_second_order_anchor_step_array(
        anchor_pair = anchor_pair,
        kernel = kernels[rows, time, , , , drop = FALSE]
      )
      end_pos <- match(time, time_indices)
      if (!is.na(end_pos) && end_pos > start_pos) {
        out[start_pos, end_pos] <- plot_second_order_anchor_cross_moment(
          anchor_pair,
          scores
        )
      }
    }
  }

  out
}

plot_second_order_anchor_pair_array <- function(pair_history, n_states) {
  n_rows <- dim(pair_history)[1]
  pair_history <- array(pair_history, dim = c(n_rows, n_states, n_states))
  out <- array(0, dim = c(n_rows, n_states, n_states, n_states))

  for (previous_state in seq_len(n_states)) {
    for (anchor_state in seq_len(n_states)) {
      out[, anchor_state, previous_state, anchor_state] <-
        pair_history[, previous_state, anchor_state]
    }
  }

  out
}

plot_second_order_anchor_step_array <- function(anchor_pair, kernel) {
  n_rows <- dim(anchor_pair)[1]
  n_states <- dim(anchor_pair)[2]
  kernel <- array(kernel, dim = c(n_rows, n_states, n_states, n_states))
  out <- array(0, dim = dim(anchor_pair))

  for (second_previous in seq_len(n_states)) {
    for (previous in seq_len(n_states)) {
      mass <- array(
        anchor_pair[,, second_previous, previous, drop = FALSE],
        dim = c(n_rows, n_states)
      )
      if (!any(mass > 0)) {
        next
      }
      probs <- array(
        kernel[, second_previous, previous, , drop = FALSE],
        dim = c(n_rows, n_states)
      )
      prob_sum <- rowSums(probs)
      reachable <- rowSums(mass) > 0
      if (any(reachable & prob_sum == 0)) {
        stop(
          "Second-order model-implied correlation recursion encountered a ",
          "reachable previous-state pair with no transition probabilities."
        )
      }
      for (state in seq_len(n_states)) {
        out[,, previous, state] <-
          out[,, previous, state] + mass * probs[, state]
      }
    }
  }

  if (abs(sum(out) - sum(anchor_pair)) > 1e-8 * max(1, n_states)) {
    stop(
      "Second-order model-implied correlation recursion lost probability mass."
    )
  }

  out
}

plot_second_order_anchor_cross_moment <- function(anchor_pair, scores) {
  n_rows <- dim(anchor_pair)[1]
  n_states <- length(scores)
  score_product <- outer(scores, scores)
  values <- numeric(n_rows)

  for (anchor_state in seq_len(n_states)) {
    for (state in seq_len(n_states)) {
      joint <- array(
        anchor_pair[, anchor_state, , state, drop = FALSE],
        dim = c(n_rows, n_states)
      )
      values <- values + rowSums(joint) * score_product[anchor_state, state]
    }
  }

  mean(values, na.rm = TRUE)
}

plot_correlation_summarize_draws <- function(data, facet_var, time_keys) {
  data$time_1 <- as.character(data$time_1)
  data$time_2 <- as.character(data$time_2)
  group_cols <- c("time_1", "time_2", facet_var)
  out <- stats::aggregate(
    data["correlation"],
    data[group_cols],
    FUN = function(x) {
      if (all(is.na(x))) {
        NA_real_
      } else {
        mean(x, na.rm = TRUE)
      }
    }
  )
  plot_correlation_finalize(out, facet_var = facet_var, time_keys = time_keys)
}

plot_correlation_data <- function(
  data,
  id_var,
  time_var,
  y_var,
  facet_var,
  ylevels,
  triangle
) {
  plot_validate_columns(data, c(id_var, time_var, y_var), "`data`")
  plot_validate_facets(data, facet_var)

  groups <- plot_facet_groups(data, facet_var)
  out <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    subset <- plot_group_subset(data, groups[i, , drop = FALSE], facet_var)
    mat <- plot_state_time_matrix(
      subset,
      id_var = id_var,
      time_var = time_var,
      y_var = y_var,
      ylevels = ylevels
    )
    corr <- stats::cor(mat, method = "pearson", use = "pairwise.complete.obs")
    out[[i]] <- plot_correlation_matrix_long(
      corr,
      group = groups[i, , drop = FALSE],
      facet_var = facet_var,
      triangle = triangle
    )
  }

  out <- bind_rows_fill(out)
  time_keys <- as.character(plot_ordered_values(data[[time_var]]))
  plot_correlation_finalize(
    out,
    facet_var = facet_var,
    time_keys = time_keys
  )
}

plot_state_time_matrix <- function(
  data,
  id_var,
  time_var,
  y_var,
  ylevels = NULL
) {
  ids <- unique(as.character(data[[id_var]]))
  times <- plot_ordered_values(data[[time_var]])
  time_labels <- as.character(times)
  id_index <- match(as.character(data[[id_var]]), ids)
  time_index <- match(as.character(data[[time_var]]), time_labels)

  key <- paste(id_index, time_index, sep = "\r")
  if (anyDuplicated(key)) {
    stop(
      "`data` must contain at most one row per `",
      id_var,
      "` and `",
      time_var,
      "` combination within each correlation stratum."
    )
  }

  mat <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(time_labels),
    dimnames = list(ids, time_labels)
  )
  mat[cbind(id_index, time_index)] <- plot_state_scores(data[[y_var]], ylevels)
  mat
}

plot_state_scores <- function(x, ylevels = NULL) {
  if (!is.null(ylevels)) {
    ylevel_names <- as_state_labels(ylevels)
    out <- match(as.character(x), ylevel_names)
    bad <- is.na(out) & !is.na(x)
    if (any(bad)) {
      stop("`data` contains state values that are not in `ylevels`.")
    }
    return(out)
  }

  if (is.factor(x)) {
    return(as.integer(x))
  }
  out <- suppressWarnings(as.integer(as.character(x)))
  if (anyNA(out) && any(!is.na(x))) {
    levels <- sort(unique(as.character(x)))
    out <- as.integer(factor(as.character(x), levels = levels))
  }
  out
}

plot_correlation_matrix_long <- function(corr, group, facet_var, triangle) {
  keep <- switch(
    triangle,
    upper = upper.tri(corr, diag = FALSE),
    full = row(corr) != col(corr)
  )
  idx <- which(keep, arr.ind = TRUE)
  out <- data.frame(
    time_1 = rownames(corr)[idx[, 1]],
    time_2 = colnames(corr)[idx[, 2]],
    correlation = corr[idx],
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_var)) {
    for (var in facet_var) {
      out[[var]] <- as.character(group[[var]][1])
    }
  }
  out$time_1 <- factor(out$time_1, levels = colnames(corr))
  out$time_2 <- factor(out$time_2, levels = colnames(corr))
  out
}

plot_correlation_finalize <- function(data, facet_var, time_keys) {
  data$time_1 <- factor(as.character(data$time_1), levels = time_keys)
  data$time_2 <- factor(as.character(data$time_2), levels = time_keys)
  if (length(facet_var) > 0L) {
    data$.panel <- plot_generic_panel(data, facet_var)
    data$.panel <- factor(data$.panel, levels = unique(data$.panel))
  }
  data
}

plot_variogram_data <- function(corr) {
  time_1 <- as.character(corr$time_1)
  time_2 <- as.character(corr$time_2)
  t1 <- suppressWarnings(as.numeric(time_1))
  t2 <- suppressWarnings(as.numeric(time_2))
  if (anyNA(t1) || anyNA(t2)) {
    levels <- unique(c(time_1, time_2))
    t1 <- match(time_1, levels)
    t2 <- match(time_2, levels)
  }
  corr$delta <- abs(t1 - t2)
  corr
}
