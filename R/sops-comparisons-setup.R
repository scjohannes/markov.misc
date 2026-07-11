# Average comparison setup helpers.

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

validate_avg_comparison_estimand <- function(estimand, state_sets, comparison) {
  if (estimand == "time_benefit") {
    if (!is.null(state_sets)) {
      stop("`state_sets` is not used with `estimand = \"time_benefit\"`.")
    }
    if (comparison != "difference") {
      stop(
        "`estimand = \"time_benefit\"` only supports `comparison = \"difference\"`."
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
  y_levels,
  absorb,
  time_var,
  p_var,
  id_var,
  p2_var,
  gap_var,
  time_covariates,
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
      time_var = time_var
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
    time_var = time_var,
    p_var = p_var,
    id_var = id_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    time_covariates = time_covariates,
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
    time_var = attr(x, "time_var"),
    p_var = attr(x, "p_var"),
    p2_var = attr(x, "p2_var"),
    y_levels = attr(x, "y_levels"),
    absorb = attr(x, "absorb"),
    gap_var = attr(x, "gap_var"),
    time_covariates = attr(x, "time_covariates")
  )
}

avg_comparison_replay_avg_sops <- function(
  model,
  newdata,
  refit_data,
  variables,
  by,
  times,
  y_levels,
  absorb,
  time_var,
  p_var,
  id_var,
  p2_var,
  gap_var,
  time_covariates,
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
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      id_var = id_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
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

normalize_comparison_state_sets <- function(states, y_levels) {
  state_labels <- as_state_labels(y_levels)
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


set_avg_comparison_attrs <- function(
  result,
  model,
  setup,
  estimand,
  state_sets,
  comparison,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  time_map,
  origin_time,
  target_times,
  origin,
  time_unit,
  extra_args
) {
  result <- as.data.frame(result)
  attrs <- list(
    model = model,
    avg_args = list(
      variables = setup$variables,
      by = setup$by,
      times = setup$times,
      id_var = setup$id_var
    ),
    comparison_args = list(
      estimand = estimand,
      state_sets = state_sets,
      comparison = comparison,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      time_map = time_map,
      origin_time = origin_time,
      target_times = target_times,
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
    "time_var",
    "p_var",
    "p2_var",
    "y_levels",
    "absorb",
    "gap_var",
    "time_covariates"
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

  result <- order_estimate_columns(result)
  class(result) <- c("markov_avg_comparisons", "data.frame")
  result
}
