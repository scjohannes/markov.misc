# Average SOP comparisons.

#' Average Comparisons From Markov State Occupancy Probabilities
#'
#' Computes average contrasts between counterfactual state occupancy summaries.
#' The first value in `variables[[1]]` is the reference level, and every later
#' value is compared to it. Reverse the order of values to reverse the contrast.
#'
#' @param model A fitted Markov transition model supported by [avg_sops()].
#' @param newdata Optional data frame of standardization profiles. See
#'   [avg_sops()].
#' @param variables A named list with one counterfactual variable. The variable
#'   must contain at least two values, for example `list(tx = c(0, 1))`.
#' @param metric Character scalar. One of `"sop"`, `"time_in_state"`, or
#'   `"time_benefit"`.
#' @param states State selection for `"sop"` and `"time_in_state"`. `NULL`
#'   returns one result per state. An atomic vector is treated as one lumped
#'   state set. A named list returns one result per named state set.
#' @param comparison Character scalar. `"difference"` computes comparison level
#'   minus reference level. `"ratio"` computes comparison level divided by
#'   reference level. `"time_benefit"` only supports `"difference"`.
#' @param by Optional character vector of variables to stratify by after
#'   standardization.
#' @param times Required visit-scale time points. See [avg_sops()].
#' @param ylevels A vector of state levels. See [avg_sops()].
#' @param absorb The absorbing state. See [avg_sops()].
#' @param time_map Optional named numeric vector or data frame mapping visit
#'   labels to real elapsed times. Used by `"time_in_state"` and
#'   `"time_benefit"`.
#' @param origin_time Optional real time for an empirical baseline anchor. See
#'   [interpolate_sops()].
#' @param xout Optional numeric real-time grid when `time_map` is supplied.
#' @param origin Origin handling when `origin_time` is supplied. See
#'   [interpolate_sops()].
#' @param time_unit Optional label stored in output.
#' @param refit_data Optional full longitudinal data used only by refit-bootstrap
#'   inference. It is not used for point estimates. See [avg_sops()].
#' @param id_var Name of the patient ID variable. See [avg_sops()].
#' @param tvarname Name of the time variable in the model. See [avg_sops()].
#' @param pvarname Name of the previous state variable in the model. See
#'   [avg_sops()].
#' @param p2varname Optional second previous-state variable. See [avg_sops()].
#' @param gap Name of the time gap variable, if used. See [avg_sops()].
#' @param t_covs Optional time-varying covariate lookup table for explicit
#'   precomputed time-basis columns. See [avg_sops()].
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws.
#' @param seed Optional random seed for reproducible posterior draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior comparison draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior comparison draws.
#'   Frequentist uncertainty draws are stored by [inferences()].
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A data frame of class `markov_avg_comparisons`.
#'
#' @details
#' `metric = "sop"` and `metric = "time_in_state"` are computed from marginal
#' SOPs. `metric = "time_benefit"` is computed before marginalization from
#' paired patient/profile-level counterfactual state distributions, because it
#' is nonlinear in the two state distributions.
#'
#' Frequentist uncertainty follows the existing package workflow:
#' `avg_comparisons(...) |> inferences(...)`. Bayesian `blrm` models use their
#' posterior SOP draws directly and therefore return uncertainty from
#' `avg_comparisons()`.
#'
#' @seealso [avg_sops()], [sops()], [inferences()]
#'
#' @examples
#' \dontrun{
#' avg_comparisons(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   metric = "time_in_state",
#'   states = "1",
#'   times = 1:30,
#'   ylevels = 1:6,
#'   absorb = 6
#' ) |>
#'   inferences(method = "simulation", n_sim = 500)
#' }
#'
#' @export
avg_comparisons <- function(
  model,
  newdata = NULL,
  variables,
  metric = c("sop", "time_in_state", "time_benefit"),
  states = NULL,
  comparison = c("difference", "ratio"),
  by = NULL,
  times,
  ylevels = NULL,
  absorb = NULL,
  time_map = NULL,
  origin_time = NULL,
  xout = NULL,
  origin = c("empirical_baseline", "none"),
  time_unit = NULL,
  refit_data = NULL,
  id_var = NULL,
  tvarname = "time",
  pvarname = "yprev",
  p2varname = NULL,
  gap = NULL,
  t_covs = NULL,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  metric <- match.arg(metric)
  comparison <- match.arg(comparison)
  origin <- match.arg(origin)
  posterior_summary <- match.arg(posterior_summary)

  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied to `avg_comparisons()`.")
  }
  variables <- validate_avg_comparison_variables(variables)
  validate_avg_comparison_metric(metric, states, comparison)

  if (metric == "time_benefit") {
    setup <- avg_comparison_setup(
      model = model,
      newdata = newdata,
      refit_data = refit_data,
      variables = variables,
      by = by,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      id_var = id_var,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      ...
    )
    tb <- avg_comparison_time_benefit_point(
      model = model,
      setup = setup,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      comparison = comparison,
      time_map = time_map,
      origin_time = origin_time,
      xout = xout,
      origin = origin,
      time_unit = time_unit,
      ...
    )
    result <- tb$result
    setup <- tb$setup
  } else {
    avg <- avg_comparison_replay_avg_sops(
      model = model,
      newdata = newdata,
      refit_data = refit_data,
      variables = variables,
      by = by,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      id_var = id_var,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = inherits(model, "blrm") || isTRUE(return_draws),
      extra_args = list(...)
    )
    state_sets <- normalize_comparison_state_sets(states, attr(avg, "ylevels"))
    result <- avg_comparison_from_avg_sops(
      avg,
      metric = metric,
      state_sets = state_sets,
      comparison = comparison,
      time_map = time_map,
      origin_time = origin_time,
      xout = xout,
      origin = origin,
      time_unit = time_unit,
      return_draws = return_draws
    )
    setup <- avg_comparison_setup_from_sops(avg)
  }

  result <- set_avg_comparison_attrs(
    result = result,
    model = model,
    setup = setup,
    metric = metric,
    states = states,
    comparison = comparison,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    time_map = time_map,
    origin_time = origin_time,
    xout = xout,
    origin = origin,
    time_unit = time_unit,
    extra_args = list(...)
  )

  result
}

validate_avg_comparison_variables <- function(variables) {
  if (missing(variables) || is.null(variables)) {
    stop("`variables` must be a named list with one counterfactual variable.")
  }
  if (!is.list(variables) || is.null(names(variables))) {
    stop("`variables` must be a named list, for example `list(tx = c(0, 1))`.")
  }
  if (length(variables) != 1L || !nzchar(names(variables)[1])) {
    stop("`avg_comparisons()` supports exactly one named variable in v1.")
  }

  values <- variables[[1]]
  if (length(values) < 2L) {
    stop("`variables[[1]]` must contain at least two values.")
  }
  if (anyDuplicated(as.character(values))) {
    stop("`variables[[1]]` must not contain duplicate values.")
  }

  variables
}

validate_avg_comparison_metric <- function(metric, states, comparison) {
  if (metric == "time_benefit") {
    if (!is.null(states)) {
      stop("`states` is not used with `metric = \"time_benefit\"`.")
    }
    if (comparison != "difference") {
      stop(
        "`metric = \"time_benefit\"` only supports `comparison = \"difference\"`."
      )
    }
  }
  invisible(NULL)
}

avg_comparison_setup <- function(
  model,
  newdata,
  refit_data,
  variables,
  by,
  times,
  ylevels,
  absorb,
  tvarname,
  pvarname,
  id_var,
  p2varname,
  gap,
  t_covs,
  include_re,
  ...
) {
  validate_markov_model(model)

  data_res <- resolve_markov_source_data(model, newdata, refit_data)
  newdata_orig <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied

  if (inherits(model, "blrm") && isTRUE(include_re)) {
    id_var <- resolve_blrm_id_var(model, newdata_orig, id_var)
    validate_markov_id_var(id_var, newdata_orig, "newdata")
  } else {
    id_var <- markov_model_id_var(model, id_var) %||% "id"
  }

  if (!is.null(refit_data)) {
    validate_markov_id_var(id_var, refit_data, "refit_data")
  }
  if (!newdata_supplied) {
    validate_markov_id_var(id_var, newdata_orig, "stored model data")
  }

  varname <- names(variables)[1]
  if (!varname %in% names(newdata_orig)) {
    stop("Variables not in data: ", varname)
  }

  baseline_data <- if (newdata_supplied) {
    newdata_orig
  } else {
    resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      tvarname = tvarname
    )
  }
  validate_sops_by(by, baseline_data)
  baseline_data <- ensure_markov_rowid(baseline_data)

  grid <- do.call(expand.grid, variables)
  newdata_pred <- create_counterfactual_data(baseline_data, grid, variables)

  list(
    variables = variables,
    by = by,
    times = times,
    tvarname = tvarname,
    pvarname = pvarname,
    id_var = id_var,
    p2varname = p2varname,
    ylevels = ylevels,
    absorb = absorb,
    gap = gap,
    t_covs = t_covs,
    newdata_orig = newdata_orig,
    refit_data = refit_data,
    newdata_supplied = newdata_supplied,
    baseline_data = baseline_data,
    grid = grid,
    newdata_pred = newdata_pred,
    n_cf = nrow(grid),
    n_each = nrow(baseline_data)
  )
}

avg_comparison_setup_from_sops <- function(x) {
  avg_args <- attr(x, "avg_args")
  variables <- avg_args$variables
  grid <- do.call(expand.grid, variables)
  newdata_pred <- attr(x, "newdata_pred")
  n_each <- if (!is.null(newdata_pred) && nrow(grid) > 0L) {
    nrow(newdata_pred) / nrow(grid)
  } else {
    NA_real_
  }

  list(
    variables = variables,
    by = avg_args$by,
    times = avg_args$times,
    id_var = avg_args$id_var,
    newdata_orig = attr(x, "newdata_orig"),
    refit_data = attr(x, "refit_data"),
    newdata_supplied = isTRUE(attr(x, "newdata_supplied")),
    baseline_data = if (!is.null(newdata_pred) && is.finite(n_each)) {
      newdata_pred[seq_len(n_each), , drop = FALSE]
    } else {
      NULL
    },
    grid = grid,
    newdata_pred = newdata_pred,
    n_cf = nrow(grid),
    n_each = as.integer(n_each),
    call_args = attr(x, "call_args"),
    tvarname = attr(x, "tvarname"),
    pvarname = attr(x, "pvarname"),
    p2varname = attr(x, "p2varname"),
    ylevels = attr(x, "ylevels"),
    absorb = attr(x, "absorb"),
    gap = attr(x, "gap"),
    t_covs = attr(x, "t_covs")
  )
}

avg_comparison_replay_avg_sops <- function(
  model,
  newdata,
  refit_data,
  variables,
  by,
  times,
  ylevels,
  absorb,
  tvarname,
  pvarname,
  id_var,
  p2varname,
  gap,
  t_covs,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  extra_args
) {
  args <- c(
    list(
      model = model,
      newdata = newdata,
      refit_data = refit_data,
      variables = variables,
      by = by,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      id_var = id_var,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws
    ),
    extra_args
  )
  do.call(avg_sops, args)
}

normalize_comparison_state_sets <- function(states, ylevels) {
  state_labels <- as_state_labels(ylevels)
  if (is.null(states)) {
    out <- as.list(state_labels)
    names(out) <- state_labels
    return(out)
  }

  if (is.list(states)) {
    out <- lapply(states, as_state_labels)
    if (is.null(names(out))) {
      names(out) <- vapply(out, state_set_label, character(1))
    }
  } else {
    out <- list(as_state_labels(states))
    names(out) <- state_set_label(out[[1]])
  }

  missing_states <- setdiff(
    unique(unlist(out, use.names = FALSE)),
    state_labels
  )
  if (length(missing_states) > 0L) {
    stop(
      "State(s) not found in SOP output: ",
      paste(missing_states, collapse = ", ")
    )
  }

  out
}

state_set_label <- function(states) {
  paste(as_state_labels(states), collapse = "+")
}

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

set_avg_comparison_attrs <- function(
  result,
  model,
  setup,
  metric,
  states,
  comparison,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  time_map,
  origin_time,
  xout,
  origin,
  time_unit,
  extra_args
) {
  attrs <- list(
    model = model,
    avg_args = list(
      variables = setup$variables,
      by = setup$by,
      times = setup$times,
      id_var = setup$id_var
    ),
    comparison_args = list(
      metric = metric,
      states = states,
      comparison = comparison,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      time_map = time_map,
      origin_time = origin_time,
      xout = xout,
      origin = origin,
      time_unit = time_unit,
      extra_args = extra_args
    ),
    newdata_orig = setup$newdata_orig,
    newdata_pred = setup$newdata_pred,
    refit_data = setup$refit_data,
    newdata_supplied = setup$newdata_supplied,
    comparison_baseline_data = setup$baseline_data,
    comparison_grid = setup$grid,
    comparison_n_each = setup$n_each,
    id_var = setup$id_var
  )

  for (nm in names(attrs)) {
    if (!is.null(attrs[[nm]])) {
      attr(result, nm) <- attrs[[nm]]
    }
  }

  for (nm in c(
    "call_args",
    "tvarname",
    "pvarname",
    "p2varname",
    "ylevels",
    "absorb",
    "gap",
    "t_covs"
  )) {
    if (!is.null(setup[[nm]])) {
      attr(result, nm) <- setup[[nm]]
    }
  }

  if (inherits(model, "blrm")) {
    attr(result, "method") <- "posterior"
    attr(result, "engine") <- "posterior"
    attr(result, "conf_level") <- conf_level
    attr(result, "posterior_summary") <- posterior_summary
  }

  class(result) <- c(
    "markov_avg_comparisons",
    setdiff(class(result), "markov_avg_comparisons")
  )
  result
}

inferences_avg_comparisons <- function(
  object,
  method,
  engine,
  score_weight_dist,
  n_sim,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  update_datadist,
  use_coefstart
) {
  args <- attr(object, "comparison_args")
  metric <- args$metric

  if (metric %in% c("sop", "time_in_state")) {
    return(inferences_avg_comparisons_linear(
      object = object,
      method = method,
      engine = engine,
      score_weight_dist = score_weight_dist,
      n_sim = n_sim,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    ))
  }

  if (method == "simulation") {
    return(inferences_avg_comparisons_time_benefit_simulation(
      object = object,
      engine = engine,
      score_weight_dist = score_weight_dist,
      n_sim = n_sim,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws
    ))
  }

  bootstrap_engine <- if (engine == "mvn") "standard" else engine
  inferences_avg_comparisons_time_benefit_bootstrap(
    object = object,
    engine = bootstrap_engine,
    n_boot = n_sim,
    workers = workers,
    conf_level = conf_level,
    return_draws = return_draws,
    update_datadist = update_datadist,
    use_coefstart = use_coefstart
  )
}

inferences_avg_comparisons_linear <- function(
  object,
  method,
  engine,
  score_weight_dist,
  n_sim,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  update_datadist,
  use_coefstart
) {
  args <- attr(object, "comparison_args")
  avg_args <- attr(object, "avg_args")
  newdata <- if (isTRUE(attr(object, "newdata_supplied"))) {
    attr(object, "newdata_orig")
  } else {
    NULL
  }
  avg <- avg_comparison_replay_avg_sops(
    model = attr(object, "model"),
    newdata = newdata,
    refit_data = attr(object, "refit_data"),
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    ylevels = attr(object, "ylevels"),
    absorb = attr(object, "absorb"),
    tvarname = attr(object, "tvarname") %||% "time",
    pvarname = attr(object, "pvarname") %||% "yprev",
    id_var = avg_args$id_var,
    p2varname = attr(object, "p2varname"),
    gap = attr(object, "gap"),
    t_covs = attr(object, "t_covs"),
    include_re = args$include_re,
    n_draws = args$n_draws,
    seed = args$seed,
    posterior_summary = args$posterior_summary,
    conf_level = conf_level,
    return_draws = FALSE,
    extra_args = args$extra_args
  )
  avg <- inferences(
    avg,
    method = method,
    engine = engine,
    score_weight_dist = score_weight_dist,
    n_sim = n_sim,
    vcov = vcov,
    cluster = cluster,
    workers = workers,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = TRUE,
    update_datadist = update_datadist,
    use_coefstart = use_coefstart
  )

  state_sets <- normalize_comparison_state_sets(
    args$states,
    attr(avg, "ylevels")
  )
  result <- avg_comparison_from_avg_sops(
    avg,
    metric = args$metric,
    state_sets = state_sets,
    comparison = args$comparison,
    time_map = args$time_map,
    origin_time = args$origin_time,
    xout = args$xout,
    origin = args$origin,
    time_unit = args$time_unit,
    return_draws = return_draws
  )
  result <- restore_avg_comparison_inference_attrs(result, object, avg)
  if (!isTRUE(return_draws)) {
    attr(result, "simulation_draws") <- NULL
    attr(result, "bootstrap_draws") <- NULL
  }
  result
}

restore_avg_comparison_inference_attrs <- function(result, object, inferred) {
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  for (a in c(
    "n_sim",
    "n_boot",
    "n_successful",
    "conf_level",
    "conf_type",
    "method",
    "engine",
    "score_weight_dist",
    "fwb_weight_type",
    "fwb_weight_scale",
    "draw_weights_attached",
    "draw_weight_col",
    "draw_weight_omission_reason"
  )) {
    if (!is.null(attr(inferred, a))) {
      attr(result, a) <- attr(inferred, a)
    }
  }
  class(result) <- class(object)
  result
}

inferences_avg_comparisons_time_benefit_simulation <- function(
  object,
  engine,
  score_weight_dist,
  n_sim,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws
) {
  if (!is.null(workers) && workers > 1) {
    warning(
      "`workers` is ignored for `time_benefit` comparison inference in v1.",
      call. = FALSE
    )
  }
  if (engine %in% c("standard", "fwb")) {
    stop(
      "`engine = \"",
      engine,
      "\"` is only used when `method = \"bootstrap\"`."
    )
  }

  setup <- setup_from_avg_comparison_object(object)
  args <- attr(object, "comparison_args")
  model <- attr(object, "model")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))

  if (engine == "mvn") {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      stop(
        "Package 'mvtnorm' is required for MVN simulation inference.\n",
        "Install with: install.packages('mvtnorm')"
      )
    }
    beta_hat <- get_coef(model)
    Sigma <- if (!is.null(vcov) && is.matrix(vcov)) {
      validate_coef_vcov(beta_hat, vcov, arg = "vcov")
    } else {
      validate_coef_vcov(beta_hat, get_vcov_robust(model), arg = "model vcov")
    }
    beta_draws <- mvtnorm::rmvnorm(n_sim, mean = beta_hat, sigma = Sigma)
    baseline_weights_draws <- NULL
  } else if (engine == "score_bootstrap") {
    if (!is.null(vcov)) {
      stop(
        "`vcov` cannot be supplied when `engine = \"score_bootstrap\"` because ",
        "draws are generated from score perturbations."
      )
    }
    return_baseline_weights <- !newdata_supplied
    if (return_baseline_weights) {
      validate_prediction_weight_ids(setup$baseline_data, setup$id_var)
    }
    score_draws <- generate_score_bootstrap_draws(
      model = model,
      baseline_data = if (return_baseline_weights) {
        setup$baseline_data
      } else {
        NULL
      },
      id_var = setup$id_var,
      n_sim = n_sim,
      score_weight_dist = score_weight_dist,
      cluster = cluster,
      return_baseline_weights = return_baseline_weights
    )
    beta_draws <- score_draws$beta_draws
    baseline_weights_draws <- if (return_baseline_weights) {
      score_draws$baseline_weights
    } else {
      warn_fixed_profile_bootstrap_weights("score-bootstrap")
      NULL
    }
  } else {
    stop("Unknown simulation engine: ", engine)
  }

  draw_results <- vector("list", nrow(beta_draws))
  for (i in seq_len(nrow(beta_draws))) {
    model_i <- set_coef(model, beta_draws[i, ])
    sops_array <- tryCatch(
      soprob_markov(
        object = model_i,
        data = setup$newdata_pred,
        times = setup$times,
        ylevels = setup$ylevels,
        absorb = setup$absorb,
        tvarname = setup$tvarname,
        pvarname = setup$pvarname,
        p2varname = setup$p2varname,
        gap = setup$gap,
        t_covs = setup$t_covs
      ),
      error = function(e) {
        warning("soprob_markov failed in draw ", i, ": ", e$message)
        NULL
      }
    )
    if (is.null(sops_array)) {
      next
    }
    baseline_weights <- if (!is.null(baseline_weights_draws)) {
      baseline_weights_draws[i, ]
    } else {
      NULL
    }
    draw_i <- reduce_time_benefit_array_for_setup(
      sops_array = sops_array,
      setup = setup,
      comparison = args$comparison,
      weights = baseline_weights,
      time_unit = args$time_unit
    )
    draw_i$draw_id <- i
    draw_results[[i]] <- draw_i
  }

  draw_results <- Filter(Negate(is.null), draw_results)
  if (length(draw_results) == 0L) {
    stop("All simulation draws failed.")
  }
  draws_df <- bind_rows_fill(draw_results)
  finalize_avg_comparison_draw_inference(
    object = object,
    draws_df = draws_df,
    draw_attr = "simulation_draws",
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = return_draws,
    metadata = list(
      n_sim = n_sim,
      n_successful = length(draw_results),
      conf_level = conf_level,
      conf_type = conf_type,
      method = "simulation",
      engine = engine,
      score_weight_dist = if (engine == "score_bootstrap") {
        score_weight_dist
      } else {
        NULL
      }
    )
  )
}

inferences_avg_comparisons_time_benefit_bootstrap <- function(
  object,
  engine,
  n_boot,
  workers,
  conf_level,
  return_draws,
  update_datadist,
  use_coefstart
) {
  engine <- match.arg(engine, choices = c("standard", "fwb"))
  setup <- setup_from_avg_comparison_object(object)
  args <- attr(object, "comparison_args")
  model <- attr(object, "model")
  refit_data <- attr(object, "refit_data") %||% attr(object, "newdata_orig")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))

  if (is.null(refit_data)) {
    stop("Full refit data not stored. Cannot perform bootstrap.")
  }
  validate_markov_id_var(setup$id_var, refit_data, "refit_data")
  if (newdata_supplied && engine == "fwb") {
    warn_fixed_profile_bootstrap_weights("FWB")
  }

  factor_cols <- c("y", setup$pvarname, setup$p2varname)
  factor_cols <- intersect(factor_cols, names(refit_data))

  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    fit_weights <- if (engine == "fwb") boot_data$fwb_weight else NULL
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      ylevels = setup$ylevels,
      absorb = setup$absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart,
      fit_weights = fit_weights
    )
    m_boot <- boot_result$model
    if (is.null(m_boot)) {
      return(NULL)
    }

    if (newdata_supplied) {
      newdata_cf <- setup$newdata_pred
      baseline_data <- setup$baseline_data
      baseline_weights <- NULL
    } else if (engine == "standard") {
      baseline_data <- resolve_markov_prediction_data(
        boot_result$data,
        id_var = "new_id",
        tvarname = setup$tvarname,
        data_label = "bootstrap data"
      )
      baseline_data <- ensure_markov_rowid(baseline_data)
      newdata_cf <- create_counterfactual_data(
        baseline_data,
        setup$grid,
        setup$variables
      )
      baseline_weights <- NULL
    } else {
      baseline_data <- resolve_markov_prediction_data(
        boot_result$data,
        id_var = setup$id_var,
        tvarname = setup$tvarname,
        data_label = "bootstrap data"
      )
      baseline_data <- ensure_markov_rowid(baseline_data)
      newdata_cf <- create_counterfactual_data(
        baseline_data,
        setup$grid,
        setup$variables
      )
      baseline_weights <- fwb_baseline_weights(
        fwb_weights = fwb_weights,
        baseline_data = baseline_data,
        id_var = setup$id_var
      )
    }

    sops_array <- tryCatch(
      soprob_markov(
        object = m_boot,
        data = newdata_cf,
        times = setup$times,
        ylevels = factor(boot_result$ylevels),
        absorb = boot_result$absorb,
        tvarname = setup$tvarname,
        pvarname = setup$pvarname,
        p2varname = setup$p2varname,
        gap = setup$gap,
        t_covs = setup$t_covs
      ),
      error = function(e) {
        warning("soprob_markov failed: ", e$message)
        NULL
      }
    )
    if (is.null(sops_array)) {
      return(NULL)
    }

    sops_array <- complete_comparison_sops_array_states(
      sops_array = sops_array,
      boot_ylevels = boot_result$ylevels,
      ylevels = setup$ylevels
    )
    boot_setup <- setup
    boot_setup$baseline_data <- baseline_data
    boot_setup$newdata_pred <- newdata_cf
    boot_setup$n_each <- nrow(baseline_data)
    reduce_time_benefit_array_for_setup(
      sops_array = sops_array,
      setup = boot_setup,
      comparison = args$comparison,
      weights = baseline_weights,
      time_unit = args$time_unit
    )
  }

  if (engine == "standard") {
    boot_ids <- fast_group_bootstrap(
      data = refit_data,
      id_var = setup$id_var,
      n_boot = n_boot
    )
    boot_results <- apply_to_bootstrap(
      boot_samples = boot_ids,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = setup$id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats"),
      globals = c(
        "model",
        "setup",
        "args",
        "refit_data",
        "newdata_supplied",
        "factor_cols",
        "engine",
        "update_datadist",
        "use_coefstart",
        "complete_comparison_sops_array_states"
      )
    )
  } else {
    fwb_samples <- generate_fwb_bootstrap_weights(
      data = refit_data,
      id_var = setup$id_var,
      n_boot = n_boot
    )
    boot_results <- apply_to_fwb_bootstrap(
      fwb_samples = fwb_samples,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = setup$id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats"),
      globals = c(
        "model",
        "setup",
        "args",
        "refit_data",
        "newdata_supplied",
        "factor_cols",
        "engine",
        "update_datadist",
        "use_coefstart",
        "complete_comparison_sops_array_states"
      )
    )
  }

  boot_results <- Filter(Negate(is.null), boot_results)
  if (length(boot_results) == 0L) {
    stop("All bootstrap iterations failed.")
  }
  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- bind_rows_fill(boot_results)

  finalize_avg_comparison_draw_inference(
    object = object,
    draws_df = boot_df,
    draw_attr = "bootstrap_draws",
    conf_level = conf_level,
    conf_type = "perc",
    return_draws = return_draws,
    metadata = list(
      n_boot = n_boot,
      n_successful = length(boot_results),
      conf_level = conf_level,
      method = "bootstrap",
      engine = engine,
      fwb_weight_type = if (engine == "fwb") "exponential" else NULL,
      fwb_weight_scale = if (engine == "fwb") "cluster_mean_1" else NULL
    )
  )
}

setup_from_avg_comparison_object <- function(object) {
  avg_args <- attr(object, "avg_args")
  call_args <- attr(object, "call_args")
  list(
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    id_var = avg_args$id_var,
    p2varname = attr(object, "p2varname"),
    newdata_orig = attr(object, "newdata_orig"),
    refit_data = attr(object, "refit_data"),
    newdata_supplied = isTRUE(attr(object, "newdata_supplied")),
    baseline_data = attr(object, "comparison_baseline_data"),
    grid = attr(object, "comparison_grid"),
    newdata_pred = attr(object, "newdata_pred"),
    n_cf = nrow(attr(object, "comparison_grid")),
    n_each = attr(object, "comparison_n_each"),
    call_args = call_args,
    tvarname = attr(object, "tvarname") %||% call_args$tvarname %||% "time",
    pvarname = attr(object, "pvarname") %||% call_args$pvarname %||% "yprev",
    ylevels = attr(object, "ylevels") %||% call_args$ylevels,
    absorb = attr(object, "absorb") %||% call_args$absorb,
    gap = attr(object, "gap") %||% call_args$gap,
    t_covs = attr(object, "t_covs") %||% call_args$t_covs
  )
}

reduce_time_benefit_array_for_setup <- function(
  sops_array,
  setup,
  comparison,
  weights,
  time_unit
) {
  arr <- array(
    NA_real_,
    dim = c(setup$n_cf, setup$n_each, dim(sops_array)[2], dim(sops_array)[3])
  )
  for (cf_i in seq_len(setup$n_cf)) {
    rows <- ((cf_i - 1L) * setup$n_each + 1L):(cf_i * setup$n_each)
    arr[cf_i, , , ] <- sops_array[rows, , , drop = FALSE]
  }
  out <- time_benefit_from_counterfactual_array(
    sops_array = arr,
    setup = setup,
    comparison = comparison,
    weights = weights,
    real_time = FALSE
  )
  if (!is.null(time_unit)) {
    out$time_unit <- time_unit
  }
  out
}

complete_comparison_sops_array_states <- function(
  sops_array,
  boot_ylevels,
  ylevels
) {
  boot_labels <- as_state_labels(boot_ylevels)
  target_labels <- as_state_labels(ylevels)
  if (identical(boot_labels, target_labels)) {
    return(sops_array)
  }

  full <- array(
    0,
    dim = c(dim(sops_array)[1], dim(sops_array)[2], length(target_labels))
  )
  for (i in seq_along(boot_labels)) {
    target <- match(boot_labels[i], target_labels)
    if (!is.na(target)) {
      full[,, target] <- sops_array[,, i]
    }
  }
  full
}

finalize_avg_comparison_draw_inference <- function(
  object,
  draws_df,
  draw_attr,
  conf_level,
  conf_type,
  return_draws,
  metadata
) {
  group_cols <- setdiff(names(draws_df), c("estimate", "draw_id"))
  ci <- compute_ci_from_draws(
    draws_df = draws_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )
  result <- merge_comparison_ci(as.data.frame(object), ci, group_cols)
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  for (nm in names(metadata)) {
    if (!is.null(metadata[[nm]])) {
      attr(result, nm) <- metadata[[nm]]
    }
  }
  if (isTRUE(return_draws)) {
    attr(result, draw_attr) <- draws_df
  }
  class(result) <- class(object)
  result
}
