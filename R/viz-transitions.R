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
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
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
  n_draws = 100L,
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
      seed = seed,
      n_draws = n_draws
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
