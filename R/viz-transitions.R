# Transition heatmaps.

#' Plot Empirical or Model-Based Transition Proportions
#'
#' `plot_transitions()` creates heatmaps of population transition proportions
#' from observed trajectory data or from trajectories simulated from a fitted
#' Markov model. Model-based plots are population averaged over the supplied or
#' stored patient profiles. When `comparison = "difference"`, counterfactual
#' trajectories are simulated for the levels in `variables`, and the heatmap
#' shows the difference in transition proportions relative to the first level.
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
#' @param p2varname Optional second previous-state column for model simulation.
#' @param id_var Optional ID column used to extract one stored prediction row per
#'   patient and to include fitted random effects for `blrm` models.
#' @param facet_var Optional observed grouping variable. The plot always facets
#'   by time; `facet_var` adds grouping to the time panels.
#' @param gap Optional time-gap variable used by the fitted model.
#' @param t_covs Optional time-varying covariate lookup table used by the fitted
#'   model.
#' @param include_re Logical. For `blrm` fits with `cluster()`, include fitted
#'   random effects for known IDs.
#' @param n_rep Number of simulated paths per patient profile for model-based
#'   plots.
#' @param n_draws Number of posterior draws to average over for `blrm` model
#'   predictions.
#' @param seed Optional random seed for model-based trajectory simulation.
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
  include_re = FALSE,
  n_rep = 20L,
  n_draws = 100L,
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
      include_re = include_re,
      n_rep = n_rep,
      n_draws = n_draws,
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
  include_re,
  n_rep,
  n_draws,
  seed
) {
  variables <- validate_plot_counterfactual_variables(
    variables,
    require_two = identical(comparison, "difference")
  )
  if (identical(comparison, "difference") && is.null(variables)) {
    stop("`variables` is required when `comparison = \"difference\"`.")
  }

  sim <- markov_simulate_fitted_paths(
    model = model,
    newdata = newdata,
    refit_data = refit_data,
    variables = variables,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = time_var,
    pvarname = pvarname,
    p2varname = p2varname,
    id_var = id_var,
    gap = gap,
    t_covs = t_covs,
    include_re = include_re,
    n_rep = n_rep,
    n_draws = n_draws,
    seed = seed
  )

  ylevels <- attr(sim, "ylevels")
  facet_summary <- facet_var
  if (!is.null(variables)) {
    facet_summary <- unique(c(facet_summary, names(variables)))
  }

  summary <- plot_transition_summary(
    data = sim,
    time_var = time_var,
    y_var = "y",
    pvarname = pvarname,
    facet_var = facet_summary,
    ylevels = ylevels
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
  out$.panel <- plot_transition_panel(out, facet_var)
  out$.panel <- factor(out$.panel, levels = unique(out$.panel))
  attr(out, "state_levels") <- state_levels
  out
}

plot_transition_difference <- function(summary, variable, values, facet_var) {
  reference <- as.character(values[1])
  comparisons <- as.character(values[-1])
  if (!variable %in% names(summary)) {
    stop("`variables` column not found in simulated transition summary.")
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
  out$.panel <- plot_transition_panel(out, c(facet_var, "contrast"))
  out$.panel <- factor(out$.panel, levels = unique(out$.panel))
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
