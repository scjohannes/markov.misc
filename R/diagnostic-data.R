# Shared data helpers for model-based diagnostic plots.

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

plot_model_trace_summaries <- function(
  model,
  setup,
  facet_var,
  pvarname,
  p2varname,
  gap,
  t_covs,
  seed,
  n_draws,
  return_kernels,
  summarize_trace,
  summarize_draws
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
      t_covs = t_covs,
      return_kernels = return_kernels
    )
    return(summarize_trace(trace, draw = FALSE))
  }

  draw_indices <- select_posterior_draws(model, n_draws = n_draws, seed = seed)
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
      return_kernels = return_kernels,
      .draw_indices = chunk,
      .gamma_draws = subset_cached_draw_matrix(gamma_draws, draw_indices, chunk)
    )
    summaries[[i]] <- summarize_trace(trace, draw = TRUE)
  }

  summarize_draws(bind_rows_fill(summaries))
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

plot_facet_groups <- function(data, facet_var) {
  if (is.null(facet_var)) {
    return(data.frame(.all = 1L))
  }
  unique(data[facet_var])
}

plot_group_subset <- function(data, group, facet_var) {
  if (is.null(facet_var)) {
    return(data)
  }
  keep <- rep(TRUE, nrow(data))
  for (var in facet_var) {
    keep <- keep & as.character(data[[var]]) == as.character(group[[var]][1])
  }
  data[keep, , drop = FALSE]
}

plot_generic_panel <- function(data, facet_var) {
  do.call(
    paste,
    c(
      lapply(facet_var, function(var) paste0(var, "=", data[[var]])),
      sep = ", "
    )
  )
}
