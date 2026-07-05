# Average comparison metric reducers.

avg_comparison_from_avg_sops <- function(
  x,
  metric,
  state_sets,
  comparison,
  time_map,
  origin_time,
  xout,
  origin,
  time_unit,
  return_draws
) {
  if (
    metric == "time_in_state" && (!is.null(time_map) || !is.null(origin_time))
  ) {
    if (is.null(time_map)) {
      stop("`time_map` must be supplied for real-time AUC.")
    }
    x <- interpolate_sops(
      x,
      time_map = time_map,
      xout = xout,
      origin_time = origin_time,
      origin = origin
    )
  }

  variables <- attr(x, "avg_args")$variables
  by <- attr(x, "avg_args")$by
  real_time <- inherits(x, "markov_interpolated_sops")
  point <- reduce_avg_sop_metric_df(
    data = as.data.frame(x),
    metric = metric,
    state_sets = state_sets,
    variables = variables,
    by = by,
    comparison = comparison,
    value_col = "estimate",
    real_time = real_time,
    time_unit = time_unit
  )

  draw_attr <- sop_draw_attr_name(x)
  draws <- if (!is.null(draw_attr)) attr(x, draw_attr) else NULL

  if (!is.null(draws)) {
    value_col <- sop_draw_value_col(draws)
    draw_reduced <- reduce_avg_sop_metric_df(
      data = as.data.frame(draws),
      metric = metric,
      state_sets = state_sets,
      variables = variables,
      by = by,
      comparison = comparison,
      value_col = value_col,
      real_time = real_time,
      time_unit = time_unit
    )
    if (value_col != "estimate") {
      names(draw_reduced)[names(draw_reduced) == value_col] <- "estimate"
    }

    group_cols <- setdiff(names(draw_reduced), c("estimate", "draw_id"))
    if (identical(attr(x, "method"), "posterior")) {
      point <- summarize_comparison_draws(
        draw_reduced,
        group_cols = group_cols,
        posterior_summary = attr(x, "posterior_summary") %||% "mean",
        conf_level = attr(x, "conf_level") %||% 0.95
      )
    } else {
      ci <- compute_ci_from_draws(
        draws_df = draw_reduced,
        group_cols = group_cols,
        conf_level = attr(x, "conf_level") %||% 0.95,
        conf_type = attr(x, "conf_type") %||% "perc"
      )
      point <- merge_comparison_ci(point, ci, group_cols)
    }

    if (isTRUE(return_draws)) {
      attr(point, draw_attr) <- draw_reduced
    }
  }

  point
}

reduce_avg_sop_metric_df <- function(
  data,
  metric,
  state_sets,
  variables,
  by,
  comparison,
  value_col,
  real_time,
  time_unit
) {
  if (metric == "sop") {
    out <- reduce_sop_comparison_df(
      data = data,
      state_sets = state_sets,
      variables = variables,
      by = by,
      comparison = comparison,
      value_col = value_col
    )
  } else if (metric == "time_in_state") {
    out <- reduce_time_in_state_comparison_df(
      data = data,
      state_sets = state_sets,
      variables = variables,
      by = by,
      comparison = comparison,
      value_col = value_col,
      real_time = real_time
    )
  } else {
    stop("Unknown comparison metric: ", metric)
  }

  if (!is.null(time_unit)) {
    out$time_unit <- time_unit
  }
  out
}

reduce_sop_comparison_df <- function(
  data,
  state_sets,
  variables,
  by,
  comparison,
  value_col
) {
  varname <- names(variables)[1]
  draw_cols <- intersect("draw_id", names(data))
  pieces <- vector("list", length(state_sets))

  for (i in seq_along(state_sets)) {
    states <- state_sets[[i]]
    keep <- as.character(data$state) %in% as.character(states)
    work <- data[keep, , drop = FALSE]
    group_cols <- unique(c(draw_cols, "time", varname, by))
    agg <- aggregate_value(work, value_col, group_cols, sum)
    comp <- compare_counterfactual_levels(
      data = agg,
      value_col = value_col,
      key_cols = group_cols,
      variables = variables,
      comparison = comparison
    )
    comp$metric <- "sop"
    comp$state_set <- names(state_sets)[i]
    pieces[[i]] <- comp
  }

  out <- bind_rows_fill(pieces)
  reorder_columns(
    out,
    c(
      "metric",
      "time",
      "state_set",
      "variable",
      "reference_level",
      "comparison_level",
      "contrast",
      "comparison",
      by,
      draw_cols,
      value_col
    )
  )
}

reduce_time_in_state_comparison_df <- function(
  data,
  state_sets,
  variables,
  by,
  comparison,
  value_col,
  real_time
) {
  varname <- names(variables)[1]
  draw_cols <- intersect("draw_id", names(data))
  pieces <- vector("list", length(state_sets))

  for (i in seq_along(state_sets)) {
    states <- state_sets[[i]]
    keep <- as.character(data$state) %in% as.character(states)
    work <- data[keep, , drop = FALSE]
    by_time_cols <- unique(c(draw_cols, "time", varname, by))
    by_time <- aggregate_value(work, value_col, by_time_cols, sum)

    reduce_cols <- unique(c(draw_cols, varname, by))
    if (isTRUE(real_time)) {
      totals <- aggregate_auc(by_time, value_col, reduce_cols)
    } else {
      totals <- aggregate_value(by_time, value_col, reduce_cols, sum)
    }

    comp <- compare_counterfactual_levels(
      data = totals,
      value_col = value_col,
      key_cols = reduce_cols,
      variables = variables,
      comparison = comparison
    )
    comp$metric <- "time_in_state"
    comp$state_set <- names(state_sets)[i]
    pieces[[i]] <- comp
  }

  out <- bind_rows_fill(pieces)
  reorder_columns(
    out,
    c(
      "metric",
      "state_set",
      "variable",
      "reference_level",
      "comparison_level",
      "contrast",
      "comparison",
      by,
      draw_cols,
      value_col
    )
  )
}

aggregate_value <- function(data, value_col, group_cols, fun) {
  if (nrow(data) == 0L) {
    stop("No rows are available for the requested state set.")
  }
  if (length(group_cols) == 0L) {
    out <- data.frame(.dummy = 1L)
    out[[value_col]] <- fun(data[[value_col]], na.rm = TRUE)
    out$.dummy <- NULL
    return(out)
  }

  f <- stats::as.formula(
    paste(value_col, "~", paste(group_cols, collapse = " + "))
  )
  stats::aggregate(f, data = data, FUN = fun, na.rm = TRUE)
}

aggregate_auc <- function(data, value_col, group_cols) {
  groups <- split(seq_len(nrow(data)), split_key(data, group_cols), drop = TRUE)
  pieces <- lapply(groups, function(idx) {
    group <- data[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    meta[[value_col]] <- trapezoid_auc(group$time, group[[value_col]])
    meta
  })
  bind_rows_fill(pieces)
}

compare_counterfactual_levels <- function(
  data,
  value_col,
  key_cols,
  variables,
  comparison
) {
  varname <- names(variables)[1]
  values <- variables[[1]]
  ref <- values[1]
  levels <- values[-1]
  base_cols <- setdiff(key_cols, varname)
  use_dummy <- length(base_cols) == 0L

  work <- data
  if (use_dummy) {
    work$.dummy_key <- 1L
    base_cols <- ".dummy_key"
  }

  ref_rows <- as.character(work[[varname]]) == as.character(ref)
  ref_df <- work[ref_rows, c(base_cols, value_col), drop = FALSE]
  if (nrow(ref_df) == 0L) {
    stop(
      "Reference level `",
      as.character(ref),
      "` was not found in comparison data."
    )
  }
  names(ref_df)[names(ref_df) == value_col] <- ".reference"

  pieces <- lapply(levels, function(level) {
    hi_rows <- as.character(work[[varname]]) == as.character(level)
    hi_df <- work[hi_rows, c(base_cols, value_col), drop = FALSE]
    if (nrow(hi_df) == 0L) {
      stop(
        "Comparison level `",
        as.character(level),
        "` was not found in comparison data."
      )
    }
    names(hi_df)[names(hi_df) == value_col] <- ".comparison"

    merged <- merge(hi_df, ref_df, by = base_cols, all.x = TRUE, sort = FALSE)
    out <- merged[, base_cols, drop = FALSE]
    out$variable <- varname
    out$reference_level <- ref
    out$comparison_level <- level
    out$contrast <- comparison_contrast_label(level, ref, comparison)
    out$comparison <- comparison
    out[[value_col]] <- comparison_value(
      hi = merged$.comparison,
      lo = merged$.reference,
      comparison = comparison
    )
    if (use_dummy) {
      out$.dummy_key <- NULL
    }
    out
  })

  bind_rows_fill(pieces)
}

comparison_contrast_label <- function(hi, lo, comparison) {
  operator <- switch(
    comparison,
    difference = " - ",
    ratio = " / ",
    stop("Unknown comparison: ", comparison)
  )
  paste0(as.character(hi), operator, as.character(lo))
}

comparison_value <- function(hi, lo, comparison) {
  if (comparison == "difference") {
    return(hi - lo)
  }
  if (comparison == "ratio") {
    out <- hi / lo
    out[lo == 0] <- NA_real_
    return(out)
  }
  stop("Unknown comparison: ", comparison)
}

avg_comparison_time_benefit_point <- function(
  model,
  setup,
  ylevels,
  absorb,
  tvarname,
  pvarname,
  p2varname,
  gap,
  t_covs,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  comparison,
  time_map,
  origin_time,
  xout,
  origin,
  time_unit,
  ...
) {
  ind <- sops(
    model = model,
    newdata = setup$newdata_pred,
    refit_data = setup$refit_data,
    times = setup$times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    id_var = setup$id_var,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    return_draws = inherits(model, "blrm") || isTRUE(return_draws),
    ...
  )

  if (!is.null(time_map) || !is.null(origin_time)) {
    if (is.null(time_map)) {
      stop("`time_map` must be supplied for real-time AUC.")
    }
    ind <- interpolate_sops(
      ind,
      time_map = time_map,
      xout = xout,
      origin_time = origin_time,
      origin = origin
    )
  }

  setup$times <- sort(unique(ind$time))
  setup$ylevels <- attr(ind, "ylevels")
  point <- reduce_time_benefit_sops_df(
    data = as.data.frame(ind),
    setup = setup,
    comparison = comparison,
    value_col = "estimate",
    real_time = inherits(ind, "markov_interpolated_sops"),
    weights = NULL,
    time_unit = time_unit
  )

  draw_attr <- sop_draw_attr_name(ind)
  draws <- if (!is.null(draw_attr)) attr(ind, draw_attr) else NULL
  if (!is.null(draws)) {
    value_col <- sop_draw_value_col(draws)
    draw_reduced <- reduce_time_benefit_sops_df(
      data = as.data.frame(draws),
      setup = setup,
      comparison = comparison,
      value_col = value_col,
      real_time = inherits(ind, "markov_interpolated_sops"),
      weights = NULL,
      time_unit = time_unit
    )
    if (value_col != "estimate") {
      names(draw_reduced)[names(draw_reduced) == value_col] <- "estimate"
    }

    group_cols <- setdiff(names(draw_reduced), c("estimate", "draw_id"))
    if (inherits(model, "blrm")) {
      point <- summarize_comparison_draws(
        draw_reduced,
        group_cols = group_cols,
        posterior_summary = posterior_summary,
        conf_level = conf_level
      )
    }
    if (isTRUE(return_draws)) {
      attr(point, draw_attr) <- draw_reduced
    }
  }

  setup <- copy_sops_attrs_to_setup(setup, ind)
  list(result = point, setup = setup)
}

copy_sops_attrs_to_setup <- function(setup, x) {
  call_args <- attr(x, "call_args")
  setup$times <- call_args$times %||% setup$times
  setup$ylevels <- attr(x, "ylevels") %||% setup$ylevels
  setup$call_args <- call_args %||% setup$call_args
  setup$tvarname <- attr(x, "tvarname") %||% setup$tvarname
  setup$pvarname <- attr(x, "pvarname") %||% setup$pvarname
  setup$p2varname <- attr(x, "p2varname") %||% setup$p2varname
  setup$absorb <- attr(x, "absorb") %||% setup$absorb
  setup$gap <- attr(x, "gap") %||% setup$gap
  setup$t_covs <- attr(x, "t_covs") %||% setup$t_covs
  setup
}

reduce_time_benefit_sops_df <- function(
  data,
  setup,
  comparison,
  value_col,
  real_time,
  weights,
  time_unit
) {
  draw_ids <- if ("draw_id" %in% names(data)) sort(unique(data$draw_id)) else NA
  pieces <- vector("list", length(draw_ids))

  for (i in seq_along(draw_ids)) {
    draw_id <- draw_ids[i]
    work <- if (is.na(draw_id)) {
      data
    } else {
      data[data$draw_id == draw_id, , drop = FALSE]
    }
    arr <- counterfactual_array_from_tidy(
      data = work,
      value_col = value_col,
      setup = setup
    )
    out <- time_benefit_from_counterfactual_array(
      sops_array = arr,
      setup = setup,
      comparison = comparison,
      weights = weights,
      real_time = real_time
    )
    if (!is.na(draw_id)) {
      out$draw_id <- draw_id
    }
    pieces[[i]] <- out
  }

  out <- bind_rows_fill(pieces)
  if (!is.null(time_unit)) {
    out$time_unit <- time_unit
  }
  out
}

counterfactual_array_from_tidy <- function(data, value_col, setup) {
  states <- as_state_labels(setup$ylevels)
  times <- setup$times %||% sort(unique(data$time))
  n_states <- length(states)
  n_times <- length(times)
  arr <- array(
    NA_real_,
    dim = c(setup$n_cf, setup$n_each, n_times, n_states)
  )

  rowid <- data$rowid
  scenario <- ((rowid - 1L) %/% setup$n_each) + 1L
  profile <- ((rowid - 1L) %% setup$n_each) + 1L
  time_idx <- match(as.character(data$time), as.character(times))
  state_idx <- match(as.character(data$state), states)
  keep <- !is.na(scenario) &
    !is.na(profile) &
    !is.na(time_idx) &
    !is.na(state_idx)

  arr[cbind(
    scenario[keep],
    profile[keep],
    time_idx[keep],
    state_idx[keep]
  )] <- data[[value_col]][keep]

  arr
}

time_benefit_from_counterfactual_array <- function(
  sops_array,
  setup,
  comparison,
  weights,
  real_time
) {
  values <- setup$variables[[1]]
  varname <- names(setup$variables)[1]
  ref_idx <- 1L
  comp_idx <- seq.int(2L, length(values))
  state_seq <- seq_along(as_state_labels(setup$ylevels))
  score <- outer(state_seq, state_seq, function(ref, hi) sign(ref - hi))
  times <- setup$times %||% seq_len(dim(sops_array)[3])
  by <- setup$by
  baseline <- setup$baseline_data
  pieces <- vector("list", length(comp_idx))

  for (i in seq_along(comp_idx)) {
    level_idx <- comp_idx[i]
    ref <- sops_array[ref_idx, , , , drop = FALSE]
    hi <- sops_array[level_idx, , , , drop = FALSE]
    ref <- array(ref, dim = dim(ref)[-1])
    hi <- array(hi, dim = dim(hi)[-1])

    benefit <- matrix(NA_real_, nrow = dim(ref)[1], ncol = dim(ref)[2])
    for (time_i in seq_len(dim(ref)[2])) {
      ref_t <- ref[, time_i, , drop = FALSE]
      hi_t <- hi[, time_i, , drop = FALSE]
      ref_t <- matrix(ref_t, nrow = dim(ref)[1])
      hi_t <- matrix(hi_t, nrow = dim(hi)[1])
      benefit[, time_i] <- rowSums((ref_t %*% score) * hi_t)
    }

    profile_value <- if (isTRUE(real_time)) {
      apply(benefit, 1, function(x) trapezoid_auc(as.numeric(times), x))
    } else {
      rowSums(benefit, na.rm = TRUE)
    }

    summary <- average_profile_values(
      profile_value = profile_value,
      baseline = baseline,
      by = by,
      weights = weights
    )
    summary$metric <- "time_benefit"
    summary$variable <- varname
    summary$reference_level <- values[1]
    summary$comparison_level <- values[level_idx]
    summary$contrast <- comparison_contrast_label(
      values[level_idx],
      values[1],
      comparison
    )
    summary$comparison <- comparison
    pieces[[i]] <- summary
  }

  out <- bind_rows_fill(pieces)
  reorder_columns(
    out,
    c(
      "metric",
      "variable",
      "reference_level",
      "comparison_level",
      "contrast",
      "comparison",
      by,
      "estimate"
    )
  )
}

average_profile_values <- function(
  profile_value,
  baseline,
  by,
  weights = NULL
) {
  if (length(profile_value) != nrow(baseline)) {
    stop("Profile-level comparison values are not aligned with baseline data.")
  }
  weights <- validate_sops_weights(weights, length(profile_value))
  work <- baseline[, by %||% character(), drop = FALSE]
  work$estimate <- profile_value
  if (!is.null(weights)) {
    work$.markov_misc_weight <- weights
  }

  if (is.null(by)) {
    estimate <- if (is.null(weights)) {
      mean(profile_value, na.rm = TRUE)
    } else {
      weighted_sop_mean(profile_value, weights)
    }
    return(data.frame(estimate = estimate))
  }

  aggregate_sops_estimates(
    work,
    group_cols = by,
    weight_col = if (is.null(weights)) NULL else ".markov_misc_weight"
  )
}

summarize_comparison_draws <- function(
  draws_df,
  group_cols,
  posterior_summary,
  conf_level
) {
  conf_level <- validate_conf_level(conf_level)
  alpha <- 1 - conf_level
  groups <- split(
    seq_len(nrow(draws_df)),
    split_key(draws_df, group_cols),
    drop = TRUE
  )
  out <- draws_df[vapply(groups, `[`, integer(1), 1L), group_cols, drop = FALSE]
  values <- lapply(groups, function(idx) draws_df$estimate[idx])
  out$estimate <- vapply(
    values,
    function(x) {
      switch(
        posterior_summary,
        mean = mean(x, na.rm = TRUE),
        median = stats::median(x, na.rm = TRUE)
      )
    },
    numeric(1)
  )
  out$conf.low <- vapply(
    values,
    stats::quantile,
    numeric(1),
    probs = alpha / 2,
    na.rm = TRUE
  )
  out$conf.high <- vapply(
    values,
    stats::quantile,
    numeric(1),
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )
  out$std.error <- vapply(values, stats::sd, numeric(1), na.rm = TRUE)
  rownames(out) <- NULL
  out
}

merge_comparison_ci <- function(point, ci, group_cols) {
  keep_ci <- ci[,
    c(group_cols, "conf.low", "conf.high", "std.error"),
    drop = FALSE
  ]
  point$.comparison_order <- seq_len(nrow(point))
  out <- merge(point, keep_ci, by = group_cols, all.x = TRUE, sort = FALSE)
  out <- out[order(out$.comparison_order), , drop = FALSE]
  out$.comparison_order <- NULL
  rownames(out) <- NULL
  out
}
