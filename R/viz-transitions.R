# Transition heatmaps.

#' Plot Empirical or Model-Based Transition Proportions
#'
#' `plot_transitions()` creates heatmaps of population transition proportions
#' from observed trajectory data or from transition probabilities implied by a
#' fitted Markov model. Model-based plots are population averaged over the
#' supplied or stored patient profiles. When `comparison = "difference"`,
#' counterfactual transition summaries are computed for the levels in
#' `variables`, and the heatmap shows the difference in transition proportions
#' relative to the first level.
#'
#' @details
#' Minimum inputs depend on the object type. For trajectory data, `object` must
#' contain the current time, current state, and previous-state columns named by
#' `time_var`, `y_var`, and `pvarname`. For model-based plots, pass a fitted
#' Markov model and `times`; use `newdata` for explicit prediction profiles and
#' pass `ylevels` or `absorb` when they are not inferable or absorbing-state
#' behavior matters. Difference plots also require
#' `comparison = "difference"` and `variables = list(var = c(reference,
#' comparison))`.
#'
#' @param object A trajectory data frame or a fitted Markov model.
#' @param newdata Optional data frame of prediction profiles for model-based
#'   plots. If `NULL`, wrapper-fitted models use their stored data and extract
#'   one prediction row per ID.
#' @param refit_data Optional full longitudinal data used only for stored-data
#'   resolution in model-based plots.
#' @param variables Optional named list with one counterfactual variable. For
#'   `comparison = "difference"`, at least two values are required and the
#'   first value is the reference.
#' @param times Optional current-time values to plot. Required for model-based
#'   plots.
#' @param ylevels Optional state levels. If omitted for model-based plots,
#'   levels are inferred from the model when possible.
#' @param absorb Optional absorbing-state label.
#' @param comparison `"none"` for transition proportions or `"difference"` for
#'   counterfactual differences in transition proportions.
#' @param time_var Character name of the time column in trajectory data.
#' @param y_var Character name of the current-state column in trajectory data.
#' @param pvarname Character name of the previous-state column.
#' @param p2varname Optional second previous-state column for second-order
#'   model recursion.
#' @param id_var Optional ID column used to extract one stored prediction row per
#'   patient.
#' @param facet_var Optional observed grouping variable. The plot always facets
#'   by time; `facet_var` adds grouping to the time panels.
#' @param gap Optional time-gap variable used by the fitted model.
#' @param t_covs Optional time-varying covariate lookup table used by the fitted
#'   model.
#' @param seed Optional random seed for `blrm` posterior draw selection.
#' @param show_values Logical. If `TRUE`, print rounded values in each tile.
#' @param digits Number of digits used for tile labels. Defaults to 2 for
#'   proportions and 3 for differences.
#' @param fill_limits Optional fill-scale limits. Difference plots use symmetric
#'   limits around zero when this is `NULL`.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_transitions(data, time_var = "week", pvarname = "yprev")
#'
#' plot_transitions(
#'   fit,
#'   times = 1:7,
#'   variables = list(tx = c(0, 1)),
#'   comparison = "difference",
#'   ylevels = 1:8,
#'   absorb = 8
#' )
#' }
#' @export
plot_transitions <- function(
  object,
  newdata = NULL,
  refit_data = NULL,
  variables = NULL,
  times = NULL,
  ylevels = NULL,
  absorb = NULL,
  comparison = c("none", "difference"),
  time_var = "time",
  y_var = "y",
  pvarname = "yprev",
  p2varname = NULL,
  id_var = NULL,
  facet_var = NULL,
  gap = NULL,
  t_covs = NULL,
  seed = NULL,
  show_values = TRUE,
  digits = NULL,
  fill_limits = NULL
) {
  comparison <- match.arg(comparison)
  if (
    !is.logical(show_values) || length(show_values) != 1L || is.na(show_values)
  ) {
    stop("`show_values` must be TRUE or FALSE.")
  }
  if (!is.null(digits)) {
    plot_validate_scalar(digits, "digits", lower = 0)
    digits <- floor(digits)
  }
  if (is.null(digits)) {
    digits <- if (identical(comparison, "difference")) 3L else 2L
  }

  if (markov_supported_model(object)) {
    data <- plot_transitions_model_data(
      model = object,
      newdata = newdata,
      refit_data = refit_data,
      variables = variables,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      comparison = comparison,
      time_var = time_var,
      pvarname = pvarname,
      p2varname = p2varname,
      id_var = id_var,
      facet_var = facet_var,
      gap = gap,
      t_covs = t_covs,
      seed = seed
    )
  } else {
    if (!is.data.frame(object)) {
      stop("`object` must be a data frame or a supported Markov model.")
    }
    if (!is.null(newdata) || !is.null(refit_data)) {
      stop("`newdata` and `refit_data` are only used for model-based plots.")
    }
    if (!is.null(variables)) {
      stop("`variables` is only used for model-based plots.")
    }
    if (!identical(comparison, "none")) {
      stop("`comparison = \"difference\"` requires a fitted Markov model.")
    }
    data <- plot_transitions_empirical_data(
      data = object,
      times = times,
      time_var = time_var,
      y_var = y_var,
      pvarname = pvarname,
      facet_var = facet_var,
      ylevels = ylevels
    )
  }

  plot_transitions_heatmap(
    data = data,
    comparison = comparison,
    show_values = show_values,
    digits = digits,
    fill_limits = fill_limits
  )
}

plot_transitions_empirical_data <- function(
  data,
  times,
  time_var,
  y_var,
  pvarname,
  facet_var,
  ylevels
) {
  plot_validate_columns(data, c(time_var, y_var, pvarname), "`data`")
  plot_validate_facets(data, facet_var)
  if (!is.null(times)) {
    keep <- as.character(data[[time_var]]) %in% as.character(times)
    data <- data[keep, , drop = FALSE]
  }

  plot_transition_summary(
    data = data,
    time_var = time_var,
    y_var = y_var,
    pvarname = pvarname,
    facet_var = facet_var,
    ylevels = ylevels
  )
}

plot_transitions_model_data <- function(
  model,
  newdata,
  refit_data,
  variables,
  times,
  ylevels,
  absorb,
  comparison,
  time_var,
  pvarname,
  p2varname,
  id_var,
  facet_var,
  gap,
  t_covs,
  seed
) {
  variables <- validate_plot_counterfactual_variables(
    variables,
    require_two = identical(comparison, "difference")
  )
  if (identical(comparison, "difference") && is.null(variables)) {
    stop("`variables` is required when `comparison = \"difference\"`.")
  }
  setup <- plot_transition_model_setup(
    model = model,
    newdata = newdata,
    refit_data = refit_data,
    variables = variables,
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

  facet_summary <- facet_var
  if (!is.null(variables)) {
    facet_summary <- unique(c(facet_summary, names(variables)))
  }

  summary <- plot_transition_model_summary(
    model = model,
    setup = setup,
    facet_var = facet_summary,
    pvarname = pvarname,
    p2varname = p2varname,
    gap = gap,
    t_covs = t_covs,
    seed = seed
  )

  if (identical(comparison, "none")) {
    return(summary)
  }

  plot_transition_difference(
    summary = summary,
    variable = names(variables)[1],
    values = variables[[1]],
    facet_var = facet_var
  )
}

plot_transition_model_setup <- function(
  model,
  newdata,
  refit_data,
  variables,
  times,
  ylevels,
  absorb,
  time_var,
  pvarname,
  p2varname,
  id_var,
  gap,
  t_covs
) {
  validate_markov_model(model)
  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied for model-based plots.")
  }

  data_res <- resolve_markov_source_data(model, newdata, refit_data)
  source_data <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied
  id_var <- markov_model_id_var(model, id_var) %||% "id"

  if (!is.null(refit_data)) {
    validate_markov_id_var(id_var, refit_data, "refit_data")
  }
  if (!newdata_supplied) {
    validate_markov_id_var(id_var, source_data, "stored model data")
  }

  baseline_data <- if (newdata_supplied) {
    source_data
  } else {
    resolve_markov_prediction_data(
      source_data,
      id_var = id_var,
      tvarname = time_var
    )
  }

  plot_validate_columns(baseline_data, pvarname, "`newdata`")
  if (!is.null(p2varname)) {
    plot_validate_columns(baseline_data, p2varname, "`newdata`")
  }

  if (!is.null(variables)) {
    missing_vars <- setdiff(names(variables), names(baseline_data))
    if (length(missing_vars) > 0) {
      stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
    }
    grid <- do.call(
      expand.grid,
      c(variables, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    )
    baseline_data <- create_counterfactual_data(baseline_data, grid, variables)
  }

  baseline_data <- ensure_markov_rowid(baseline_data)

  time_res <- resolve_sop_times(
    model,
    baseline_data,
    times,
    time_var,
    t_covs = NULL,
    default = "unique"
  )
  plot_times <- time_res$times
  time_info <- time_res$time_info
  recursion_times <- complete_plot_recursion_times(plot_times, time_info)
  if (!is.null(t_covs) && nrow(t_covs) != length(recursion_times)) {
    stop(
      "`t_covs` must have one row per recursion time point for model-based ",
      "plots. Sparse plot times are expanded to the full recursion grid; ",
      "expected ",
      length(recursion_times),
      " rows but got ",
      nrow(t_covs),
      "."
    )
  }
  validate_factor_gap(gap, t_covs, time_info)

  ylevels <- markov_model_ylevels(model, ylevels)
  ylevel_names <- as_state_labels(ylevels)
  plot_indices <- match(as.character(plot_times), as.character(recursion_times))

  list(
    data = baseline_data,
    id_var = id_var,
    times = recursion_times,
    plot_times = plot_times,
    plot_indices = plot_indices,
    ylevels = ylevel_names,
    absorb = absorb,
    time_var = time_var
  )
}

plot_transition_model_summary <- function(
  model,
  setup,
  facet_var,
  pvarname,
  p2varname,
  gap,
  t_covs,
  seed
) {
  plot_validate_facets(setup$data, facet_var)

  if (!inherits(model, "blrm")) {
    trace <- markov_transition_trace(
      object = model,
      data = setup$data,
      times = setup$times,
      ylevels = setup$ylevels,
      absorb = setup$absorb,
      tvarname = setup$time_var,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs
    )
    return(plot_transition_trace_summary(
      transitions = trace$transitions,
      data = setup$data,
      plot_indices = setup$plot_indices,
      plot_times = setup$plot_times,
      ylevels = setup$ylevels,
      facet_var = facet_var,
      draw = FALSE
    ))
  }

  draw_indices <- select_posterior_draws(model, n_draws = 100L, seed = seed)
  chunks <- split_draw_indices(
    draw_indices,
    chunk_size = getOption("markov.misc.blrm_avg_chunk_size", 50L)
  )
  gamma_draws <- cache_blrm_gamma_draws(
    model,
    draw_indices = draw_indices,
    include_re = FALSE
  )
  summaries <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]
    trace <- markov_transition_trace(
      object = model,
      data = setup$data,
      times = setup$times,
      ylevels = setup$ylevels,
      absorb = setup$absorb,
      tvarname = setup$time_var,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      id_var = setup$id_var,
      n_draws = NULL,
      .draw_indices = chunk,
      .gamma_draws = subset_cached_draw_matrix(gamma_draws, draw_indices, chunk)
    )
    summaries[[i]] <- plot_transition_trace_summary(
      transitions = trace$transitions,
      data = setup$data,
      plot_indices = setup$plot_indices,
      plot_times = setup$plot_times,
      ylevels = setup$ylevels,
      facet_var = facet_var,
      draw = TRUE
    )
  }

  plot_transition_summarize_draws(
    data = bind_rows_fill(summaries),
    facet_var = facet_var,
    ylevels = setup$ylevels
  )
}

plot_transition_trace_summary <- function(
  transitions,
  data,
  plot_indices,
  plot_times,
  ylevels,
  facet_var,
  draw
) {
  groups <- plot_facet_groups(data, facet_var)
  out <- vector("list", nrow(groups))
  time_keys <- as.character(plot_times)

  for (i in seq_len(nrow(groups))) {
    keep <- plot_group_subset(data, groups[i, , drop = FALSE], facet_var)
    rows <- match(rownames(keep), rownames(data))
    group <- if (is.null(facet_var)) NULL else groups[i, , drop = FALSE]
    out[[i]] <- plot_transition_trace_group_summary(
      transitions = transitions,
      rows = rows,
      time_indices = plot_indices,
      time_keys = time_keys,
      ylevels = ylevels,
      group = group,
      facet_var = facet_var,
      draw = draw
    )
  }

  out <- bind_rows_fill(out)
  plot_transition_finalize_summary(
    out,
    facet_var = facet_var,
    time_keys = time_keys,
    state_levels = ylevels
  )
}

plot_transition_trace_group_summary <- function(
  transitions,
  rows,
  time_indices,
  time_keys,
  ylevels,
  group,
  facet_var,
  draw
) {
  total <- length(rows)
  if (isTRUE(draw)) {
    mass <- transitions[, rows, time_indices, , , drop = FALSE]
    summed <- apply(mass, c(1, 3, 4, 5), sum)
    draw_ids <- dimnames(transitions)[[1]]
    out <- expand.grid(
      draw_id = draw_ids,
      .time_key = time_keys,
      previous_state = ylevels,
      state = ylevels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    mass <- transitions[rows, time_indices, , , drop = FALSE]
    summed <- apply(mass, c(2, 3, 4), sum)
    out <- expand.grid(
      .time_key = time_keys,
      previous_state = ylevels,
      state = ylevels,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }

  out$n <- as.vector(summed)
  out$total <- total
  out$estimate <- if (total > 0) out$n / total else NA_real_
  if (!is.null(facet_var)) {
    for (var in facet_var) {
      out[[var]] <- as.character(group[[var]][1])
    }
  }
  out
}

plot_transition_summarize_draws <- function(data, facet_var, ylevels) {
  group_cols <- c(".time_key", facet_var, "previous_state", "state")
  agg_formula <- stats::as.formula(
    paste("cbind(estimate, n, total) ~", paste(group_cols, collapse = " + "))
  )
  out <- stats::aggregate(
    agg_formula,
    data = data,
    FUN = mean,
    na.rm = TRUE
  )
  time_keys <- as.character(plot_ordered_values(out$.time_key))
  plot_transition_finalize_summary(
    out,
    facet_var = facet_var,
    time_keys = time_keys,
    state_levels = ylevels
  )
}

plot_transition_finalize_summary <- function(
  data,
  facet_var,
  time_keys,
  state_levels
) {
  data$previous_state <- factor(data$previous_state, levels = state_levels)
  data$state <- factor(data$state, levels = state_levels)
  panel_levels <- plot_transition_panel_levels(
    data = data,
    time_keys = time_keys,
    facet_var = facet_var
  )
  data$.panel <- plot_transition_panel(data, facet_var)
  data$.panel <- factor(data$.panel, levels = panel_levels)
  attr(data, "state_levels") <- state_levels
  data
}

plot_transition_summary <- function(
  data,
  time_var,
  y_var,
  pvarname,
  facet_var,
  ylevels
) {
  state_levels <- plot_transition_state_levels(
    data,
    vars = c(y_var, pvarname),
    ylevels = ylevels
  )
  time_values <- plot_ordered_values(data[[time_var]])
  time_keys <- as.character(time_values)

  input <- data.frame(
    .time_key = as.character(data[[time_var]]),
    previous_state = as.character(data[[pvarname]]),
    state = as.character(data[[y_var]]),
    .n = 1L,
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_var)) {
    for (var in facet_var) {
      input[[var]] <- as.character(data[[var]])
    }
  }

  by_count <- c(".time_key", facet_var, "previous_state", "state")
  counts <- stats::aggregate(
    input[".n"],
    input[by_count],
    base::sum
  )

  grid <- expand.grid(
    .time_key = time_keys,
    previous_state = state_levels,
    state = state_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_var)) {
    facet_grid <- unique(input[facet_var])
    grid <- merge(facet_grid, grid, by = NULL)
  }

  out <- merge(grid, counts, by = by_count, all.x = TRUE, sort = FALSE)
  out$.n[is.na(out$.n)] <- 0L

  by_total <- c(".time_key", facet_var)
  totals <- stats::aggregate(
    out[".n"],
    out[by_total],
    base::sum
  )
  names(totals)[names(totals) == ".n"] <- "total"
  out <- merge(out, totals, by = by_total, all.x = TRUE, sort = FALSE)
  out$estimate <- ifelse(out$total > 0, out$.n / out$total, NA_real_)
  out$n <- out$.n
  out$.n <- NULL

  out$previous_state <- factor(out$previous_state, levels = state_levels)
  out$state <- factor(out$state, levels = state_levels)
  panel_levels <- plot_transition_panel_levels(
    data = input,
    time_keys = time_keys,
    facet_var = facet_var
  )
  out$.panel <- plot_transition_panel(out, facet_var)
  out$.panel <- factor(out$.panel, levels = panel_levels)
  attr(out, "state_levels") <- state_levels
  out
}

plot_transition_difference <- function(summary, variable, values, facet_var) {
  reference <- as.character(values[1])
  comparisons <- as.character(values[-1])
  if (!variable %in% names(summary)) {
    stop("`variables` column not found in transition summary.")
  }

  key <- setdiff(
    names(summary),
    c(
      variable,
      "estimate",
      "n",
      "total",
      ".panel",
      "reference_level",
      "comparison_level",
      "contrast"
    )
  )
  key <- setdiff(key, facet_var[facet_var == variable])

  reference_data <- summary[as.character(summary[[variable]]) == reference, ]
  out <- vector("list", length(comparisons))
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    comparison_data <- summary[
      as.character(summary[[variable]]) == comparison,
      ,
      drop = FALSE
    ]
    merged <- merge(
      comparison_data,
      reference_data,
      by = key,
      suffixes = c(".comparison", ".reference"),
      sort = FALSE
    )
    merged$estimate <- merged$estimate.comparison - merged$estimate.reference
    merged$n <- merged$n.comparison
    merged$total <- merged$total.comparison
    merged$variable <- variable
    merged$reference_level <- reference
    merged$comparison_level <- comparison
    merged$contrast <- paste0(comparison, " - ", reference)
    keep <- c(
      key,
      "estimate",
      "n",
      "total",
      "variable",
      "reference_level",
      "comparison_level",
      "contrast"
    )
    out[[i]] <- merged[keep]
  }

  out <- bind_rows_fill(out)
  time_keys <- as.character(plot_ordered_values(out$.time_key))
  panel_levels <- plot_transition_panel_levels(
    data = out,
    time_keys = time_keys,
    facet_var = c(facet_var, "contrast")
  )
  out$.panel <- plot_transition_panel(out, c(facet_var, "contrast"))
  out$.panel <- factor(out$.panel, levels = panel_levels)
  attr(out, "state_levels") <- attr(summary, "state_levels")
  out
}

plot_transition_state_levels <- function(data, vars, ylevels = NULL) {
  if (!is.null(ylevels)) {
    return(as_state_labels(ylevels))
  }

  state_levels <- character()
  for (var in vars) {
    x <- data[[var]]
    if (is.factor(x)) {
      state_levels <- c(state_levels, as.character(levels(x)))
    } else {
      state_levels <- c(state_levels, as.character(stats::na.omit(unique(x))))
    }
  }
  state_levels <- unique(state_levels)
  numeric_levels <- suppressWarnings(as.numeric(state_levels))
  if (!anyNA(numeric_levels)) {
    state_levels <- state_levels[order(numeric_levels)]
  }
  state_levels
}

plot_ordered_values <- function(x) {
  if (is.factor(x)) {
    values <- levels(x)
  } else {
    values <- unique(x)
  }
  if (is.numeric(values) || is.integer(values)) {
    return(sort(values))
  }
  labels <- as.character(values)
  numeric_labels <- suppressWarnings(as.numeric(labels))
  if (length(labels) > 0L && !anyNA(numeric_labels)) {
    return(labels[order(numeric_labels, seq_along(labels))])
  }
  if (is.factor(x)) {
    return(labels)
  }
  sort(labels)
}

plot_transition_panel <- function(data, facet_var) {
  label <- paste0("Time ", data$.time_key)
  facet_var <- facet_var[!is.na(facet_var)]
  facet_var <- setdiff(facet_var, character())
  facet_var <- facet_var[facet_var %in% names(data)]
  if (length(facet_var) == 0L) {
    return(label)
  }

  facet_label <- do.call(
    paste,
    c(
      lapply(facet_var, function(var) paste0(var, "=", data[[var]])),
      sep = ", "
    )
  )
  paste(label, facet_label, sep = " | ")
}

plot_transition_panel_levels <- function(data, time_keys, facet_var) {
  facet_var <- facet_var[!is.na(facet_var)]
  facet_var <- setdiff(facet_var, character())
  facet_var <- facet_var[facet_var %in% names(data)]
  if (length(facet_var) == 0L) {
    return(paste0("Time ", time_keys))
  }

  facet_grid <- unique(data[facet_var])
  out <- vector("list", nrow(facet_grid))
  for (i in seq_len(nrow(facet_grid))) {
    panel_data <- data.frame(
      .time_key = time_keys,
      stringsAsFactors = FALSE
    )
    for (var in facet_var) {
      panel_data[[var]] <- as.character(facet_grid[[var]][i])
    }
    out[[i]] <- panel_data
  }

  plot_transition_panel(do.call(rbind, out), facet_var)
}

plot_transition_label <- function(x, digits) {
  ifelse(
    is.na(x),
    "",
    formatC(round(x, digits), format = "f", digits = digits)
  )
}

plot_transitions_heatmap <- function(
  data,
  comparison,
  show_values,
  digits,
  fill_limits
) {
  data$.label <- plot_transition_label(data$estimate, digits)
  p <- ggplot2::ggplot(data) +
    ggplot2::aes(
      x = .data$previous_state,
      y = .data$state,
      fill = .data$estimate
    ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::facet_wrap(ggplot2::vars(.data$.panel)) +
    ggplot2::labs(
      x = "Previous State",
      y = "Current State",
      fill = if (identical(comparison, "difference")) {
        "Difference in transition proportion"
      } else {
        "Transition proportion"
      }
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      aspect.ratio = 1
    )

  if (show_values) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$.label),
        color = "black",
        size = 2.5
      )
  }

  if (identical(comparison, "difference")) {
    if (is.null(fill_limits)) {
      max_abs <- max(abs(data$estimate), na.rm = TRUE)
      if (!is.finite(max_abs) || max_abs == 0) {
        max_abs <- 1
      }
      fill_limits <- c(-max_abs, max_abs)
    }
    return(
      p +
        ggplot2::scale_fill_gradient2(
          low = "#d73027",
          mid = "white",
          high = "#2c7bb6",
          midpoint = 0,
          limits = fill_limits,
          na.value = "grey90"
        )
    )
  }

  p +
    ggplot2::scale_fill_viridis_c(
      limits = fill_limits,
      na.value = "grey90"
    )
}
