# Internal transition-summary data builders for diagnostic plots.

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
  seed,
  n_draws
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
    seed = seed,
    n_draws = n_draws
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

plot_transition_model_summary <- function(
  model,
  setup,
  facet_var,
  pvarname,
  p2varname,
  gap,
  t_covs,
  seed,
  n_draws
) {
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
    return_kernels = FALSE,
    summarize_trace = function(trace, draw) {
      plot_transition_trace_summary(
        transitions = trace$transitions,
        data = setup$data,
        plot_indices = setup$plot_indices,
        plot_times = setup$plot_times,
        ylevels = setup$ylevels,
        facet_var = facet_var,
        draw = draw
      )
    },
    summarize_draws = function(data) {
      plot_transition_summarize_draws(
        data = data,
        facet_var = facet_var,
        ylevels = setup$ylevels
      )
    }
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
