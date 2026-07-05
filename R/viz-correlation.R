# Correlation heatmaps and variograms.

#' Plot Empirical or Model-Implied Correlations Over Time
#'
#' `plot_correlation()` computes correlations between ordinal state scores at
#' pairs of follow-up times. With a fitted Markov model, correlations are
#' computed from model-implied pairwise moments rather than sampled paths.
#'
#' @details
#' Minimum inputs depend on the object type. For trajectory data, `object` must
#' contain the ID, time, and state columns named by `id_var`, `time_var`, and
#' `y_var`. For model-based plots, pass a fitted Markov model and `times`; use
#' `newdata` for explicit prediction profiles and pass `ylevels` or `absorb`
#' when needed. Model-based correlations are computed from exact first- or
#' second-order Markov recursions. For `blrm` fits, posterior draws are kept
#' separate through the moment calculation and summarized only after
#' correlations are computed for each draw.
#'
#' @param object A trajectory data frame or a fitted Markov model.
#' @param newdata Optional data frame of prediction profiles for model-based
#'   plots. If `NULL`, wrapper-fitted models use their stored data and extract
#'   one prediction row per ID.
#' @param refit_data Optional full longitudinal data used only for stored-data
#'   resolution in model-based plots.
#' @param times Optional time values to include. Required for model-based plots.
#' @param ylevels Optional state levels. If supplied for raw trajectory data,
#'   state scores follow this order.
#' @param absorb Optional absorbing-state label for model-based recursions.
#' @param id_var Character ID column for raw data and stored model-data
#'   extraction.
#' @param time_var Character name of the time column.
#' @param y_var Character name of the state column in raw trajectory data.
#' @param pvarname Character name of the previous-state column for model
#'   recursions.
#' @param p2varname Optional second previous-state column for second-order
#'   model recursions.
#' @param facet_var Optional grouping variable. Correlations are computed
#'   separately within each observed stratum.
#' @param gap Optional time-gap variable used by the fitted model.
#' @param t_covs Optional time-varying covariate lookup table used by the fitted
#'   model.
#' @param seed Optional random seed for `blrm` posterior draw selection.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param triangle `"upper"` to show one triangle or `"full"` to show both
#'   off-diagonal triangles.
#' @param show_values Logical. If `TRUE`, print rounded correlations in tiles.
#' @param digits Number of digits used for tile labels.
#' @param fill_limits Fill-scale limits. Defaults to `c(0, 1)`.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_correlation(data, id_var = "id", time_var = "day")
#' plot_correlation(fit, times = 1:14, ylevels = 1:8, absorb = 8)
#' }
#' @export
plot_correlation <- function(
  object,
  newdata = NULL,
  refit_data = NULL,
  times = NULL,
  ylevels = NULL,
  absorb = NULL,
  id_var = "id",
  time_var = "time",
  y_var = "y",
  pvarname = "yprev",
  p2varname = NULL,
  facet_var = NULL,
  gap = NULL,
  t_covs = NULL,
  seed = NULL,
  n_draws = 100L,
  triangle = c("upper", "full"),
  show_values = TRUE,
  digits = 2,
  fill_limits = c(0, 1)
) {
  triangle <- match.arg(triangle)
  if (
    !is.logical(show_values) || length(show_values) != 1L || is.na(show_values)
  ) {
    stop("`show_values` must be TRUE or FALSE.")
  }
  plot_validate_scalar(digits, "digits", lower = 0)

  corr <- plot_correlation_input_data(
    object = object,
    newdata = newdata,
    refit_data = refit_data,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    id_var = id_var,
    time_var = time_var,
    y_var = y_var,
    pvarname = pvarname,
    p2varname = p2varname,
    facet_var = facet_var,
    gap = gap,
    t_covs = t_covs,
    seed = seed,
    n_draws = n_draws,
    triangle = triangle
  )

  plot_correlation_heatmap(
    corr,
    show_values = show_values,
    digits = floor(digits),
    fill_limits = fill_limits
  )
}

#' Plot an Empirical or Model-Implied Correlation Variogram
#'
#' `plot_variogram()` summarizes the same time-by-time correlation matrix used
#' by [plot_correlation()] as correlations against absolute time differences.
#' Its minimum input requirements are the same as [plot_correlation()].
#'
#' @inheritParams plot_correlation
#' @param smooth Logical. If `TRUE`, add a loess smooth.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_variogram(data, id_var = "id", time_var = "day")
#' }
#' @export
plot_variogram <- function(
  object,
  newdata = NULL,
  refit_data = NULL,
  times = NULL,
  ylevels = NULL,
  absorb = NULL,
  id_var = "id",
  time_var = "time",
  y_var = "y",
  pvarname = "yprev",
  p2varname = NULL,
  facet_var = NULL,
  gap = NULL,
  t_covs = NULL,
  seed = NULL,
  n_draws = 100L,
  smooth = TRUE
) {
  if (!is.logical(smooth) || length(smooth) != 1L || is.na(smooth)) {
    stop("`smooth` must be TRUE or FALSE.")
  }

  corr <- plot_correlation_input_data(
    object = object,
    newdata = newdata,
    refit_data = refit_data,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    id_var = id_var,
    time_var = time_var,
    y_var = y_var,
    pvarname = pvarname,
    p2varname = p2varname,
    facet_var = facet_var,
    gap = gap,
    t_covs = t_covs,
    seed = seed,
    n_draws = n_draws,
    triangle = "upper"
  )
  variogram <- plot_variogram_data(corr)

  p <- ggplot2::ggplot(variogram) +
    ggplot2::aes(x = .data$delta, y = .data$correlation) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::labs(
      x = "Absolute Time Difference",
      y = "Correlation"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme_bw()

  if (smooth) {
    p <- p +
      ggplot2::geom_smooth(
        method = "loess",
        formula = y ~ x,
        se = FALSE
      )
  }
  if (".panel" %in% names(variogram)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$.panel))
  }

  p
}

plot_correlation_label <- function(x, digits) {
  ifelse(
    is.na(x),
    "",
    formatC(round(x, digits), format = "f", digits = digits)
  )
}

plot_correlation_heatmap <- function(
  data,
  show_values,
  digits,
  fill_limits
) {
  data$.label <- plot_correlation_label(data$correlation, digits)
  p <- ggplot2::ggplot(data) +
    ggplot2::aes(x = .data$time_1, y = .data$time_2, fill = .data$correlation) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_viridis_c(limits = fill_limits, na.value = "grey90") +
    ggplot2::labs(x = "Time", y = "Time", fill = "") +
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
  if (".panel" %in% names(data)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$.panel))
  }

  p
}
