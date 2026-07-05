# Average comparison plotting.

#' Plot Average Comparisons
#'
#' Plots output from [avg_comparisons()] on the comparison scale.
#'
#' @param data A `markov_avg_comparisons` object returned by
#'   [avg_comparisons()], optionally after [inferences()].
#' @param x_var Character string or `NULL`. X-axis variable. If `NULL`, uses
#'   `"time"` when available, then `"state_set"` when available, and otherwise
#'   `"contrast"`.
#' @param color_var Character string or `NULL`. Optional variable mapped to
#'   color. If `NULL`, a useful grouping variable is inferred when needed.
#' @param linetype_var Character string or `NULL`. Optional variable mapped to
#'   linetype for line plots.
#' @param facet_var Character vector or `NULL`. Optional faceting variable. Use
#'   a length-two vector for `facet_grid(row ~ column)`.
#' @param geom Character string. `"auto"` uses lines when `x_var` is `"time"`
#'   and points otherwise. Use `"line"` or `"point"` to override.
#' @param estimate_var Character string. Column containing comparison estimates.
#' @param show_uncertainty Logical. If `TRUE`, add confidence ribbons for line
#'   plots or intervals for point plots when `conf.low` and `conf.high` exist.
#' @param ribbon_alpha Numeric alpha level for confidence ribbons.
#' @param point_size Numeric point size for point plots.
#' @param line_width Numeric line width for lines and interval geoms.
#'
#' @return A ggplot object.
#'
#' @details
#' The y axis is the contrast scale: difference contrasts are centered on a
#' dashed reference line at 0, and ratio contrasts are centered on a dashed
#' reference line at 1. A single plot must contain only one comparison scale.
#'
#' @examples
#' \dontrun{
#' cmp <- avg_comparisons(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   metric = "time_in_state",
#'   states = c(1, 2),
#'   times = 1:30
#' ) |>
#'   inferences(method = "simulation", n_sim = 500)
#'
#' plot_comparisons(cmp)
#' }
#'
#' @importFrom ggplot2 geom_errorbar geom_hline geom_point
#' @export
plot_comparisons <- function(
  data,
  x_var = NULL,
  color_var = NULL,
  linetype_var = NULL,
  facet_var = NULL,
  geom = c("auto", "line", "point"),
  estimate_var = "estimate",
  show_uncertainty = TRUE,
  ribbon_alpha = 0.12,
  point_size = 2,
  line_width = 0.7
) {
  geom <- match.arg(geom)
  plot_validate_scalar(ribbon_alpha, "ribbon_alpha", lower = 0, upper = 1)
  plot_validate_scalar(point_size, "point_size", lower = 0)
  plot_validate_scalar(line_width, "line_width", lower = 0)
  if (
    !is.logical(show_uncertainty) ||
      length(show_uncertainty) != 1 ||
      is.na(show_uncertainty)
  ) {
    stop("`show_uncertainty` must be TRUE or FALSE.")
  }

  plot_comparisons_validate_data(data, estimate_var)
  data <- as.data.frame(data)

  if (is.null(x_var)) {
    x_var <- plot_comparisons_default_x_var(data)
  }
  plot_validate_columns(data, x_var, "`x_var`")
  plot_validate_facets(data, facet_var)

  if (is.null(color_var)) {
    color_var <- plot_comparisons_default_color_var(data, x_var)
  } else {
    plot_validate_columns(data, color_var, "`color_var`")
  }

  if (identical(geom, "auto")) {
    geom <- if (identical(x_var, "time")) "line" else "point"
  }

  if (
    is.null(linetype_var) &&
      identical(geom, "line") &&
      plot_comparisons_has_multiple(data, "contrast") &&
      !identical(color_var, "contrast")
  ) {
    linetype_var <- "contrast"
  }
  if (!is.null(linetype_var)) {
    plot_validate_columns(data, linetype_var, "`linetype_var`")
  }

  data <- plot_comparisons_order_state_sets(data)
  comparison <- unique(as.character(data$comparison))
  comparison <- comparison[!is.na(comparison)]
  reference <- if (identical(comparison, "ratio")) 1 else 0
  data <- plot_comparisons_add_group(data, color_var, linetype_var)
  aes_mapping <- ggplot2::aes(
    x = .data[[x_var]],
    y = .data[[estimate_var]],
    group = .data[[".plot_group"]]
  )
  if (!is.null(color_var)) {
    aes_mapping$colour <- ggplot2::aes(
      colour = factor(.data[[color_var]])
    )$colour
  }
  if (!is.null(linetype_var)) {
    aes_mapping$linetype <- ggplot2::aes(
      linetype = factor(.data[[linetype_var]])
    )$linetype
  }

  p <- ggplot2::ggplot(data, aes_mapping) +
    ggplot2::geom_hline(
      yintercept = reference,
      linetype = "dashed",
      color = "grey50",
      linewidth = line_width * 0.6
    )

  has_interval <- show_uncertainty &&
    all(c("conf.low", "conf.high") %in% names(data))

  if (identical(geom, "line")) {
    if (has_interval) {
      ribbon_mapping <- ggplot2::aes(
        ymin = .data[["conf.low"]],
        ymax = .data[["conf.high"]],
        group = .data[[".plot_group"]]
      )
      if (!is.null(color_var)) {
        ribbon_mapping$fill <- ggplot2::aes(
          fill = factor(.data[[color_var]])
        )$fill
        p <- p +
          ggplot2::geom_ribbon(
            ribbon_mapping,
            alpha = ribbon_alpha,
            color = NA
          )
      } else {
        p <- p +
          ggplot2::geom_ribbon(
            ribbon_mapping,
            alpha = ribbon_alpha,
            color = NA,
            fill = "grey70"
          )
      }
    }
    p <- p + ggplot2::geom_line(linewidth = line_width)
  } else {
    dodge <- ggplot2::position_dodge(width = 0.35)
    if (has_interval) {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = .data[["conf.low"]],
            ymax = .data[["conf.high"]]
          ),
          width = 0.12,
          linewidth = line_width * 0.7,
          position = dodge
        )
    }
    p <- p + ggplot2::geom_point(size = point_size, position = dodge)
  }

  p <- p +
    ggplot2::labs(
      x = plot_comparisons_axis_label(x_var),
      y = plot_comparisons_y_label(comparison),
      colour = plot_comparisons_axis_label(color_var),
      fill = plot_comparisons_axis_label(color_var),
      linetype = plot_comparisons_axis_label(linetype_var),
      title = plot_comparisons_title(data)
    )

  p <- plot_add_default_scales(p)
  plot_add_facets(p, facet_var)
}

plot_comparisons_validate_data <- function(data, estimate_var) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame returned by `avg_comparisons()`.")
  }
  required <- c("metric", "comparison", "contrast", estimate_var)
  missing_vars <- setdiff(required, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "`data` is missing required comparison column(s): ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  comparison <- unique(as.character(data$comparison))
  comparison <- comparison[!is.na(comparison)]
  if (length(comparison) != 1 || !comparison %in% c("difference", "ratio")) {
    stop(
      "`plot_comparisons()` requires one comparison scale: ",
      "`difference` or `ratio`.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

plot_comparisons_default_x_var <- function(data) {
  if ("time" %in% names(data)) {
    return("time")
  }
  if ("state_set" %in% names(data)) {
    return("state_set")
  }
  "contrast"
}

plot_comparisons_default_color_var <- function(data, x_var) {
  candidates <- c("state_set", "contrast", "comparison_level")
  for (candidate in candidates) {
    if (
      candidate %in%
        names(data) &&
        !identical(candidate, x_var) &&
        plot_comparisons_has_multiple(data, candidate)
    ) {
      return(candidate)
    }
  }
  NULL
}

plot_comparisons_order_state_sets <- function(data) {
  if (!"state_set" %in% names(data)) {
    return(data)
  }
  levels <- plot_comparisons_state_set_levels(data)
  data$state_set <- factor(as.character(data$state_set), levels = levels)
  data
}

plot_comparisons_state_set_levels <- function(data) {
  values <- unique(as.character(data$state_set))
  if (is.factor(data$state_set)) {
    return(intersect(levels(data$state_set), values))
  }

  parsed <- strsplit(values, "+", fixed = TRUE)
  ylevels <- attr(data, "ylevels", exact = TRUE)
  if (!is.null(ylevels)) {
    state_levels <- as_state_labels(ylevels)
    state_order <- stats::setNames(seq_along(state_levels), state_levels)
    keys <- vapply(
      parsed,
      function(x) {
        if (all(x %in% names(state_order))) {
          return(min(state_order[x]))
        }
        NA_real_
      },
      numeric(1)
    )
    if (any(!is.na(keys))) {
      known <- which(!is.na(keys))
      unknown <- which(is.na(keys))
      known <- known[order(keys[known], lengths(parsed)[known], values[known])]
      return(c(values[known], values[unknown]))
    }
  }

  numeric_values <- suppressWarnings(as.numeric(values))
  if (all(!is.na(numeric_values))) {
    return(values[order(numeric_values)])
  }

  values
}

plot_comparisons_has_multiple <- function(data, var) {
  var %in% names(data) && length(unique(data[[var]])) > 1
}

plot_comparisons_add_group <- function(data, color_var, linetype_var) {
  group_vars <- unique(c(color_var, linetype_var))
  group_vars <- group_vars[!is.na(group_vars)]
  if (length(group_vars) == 0) {
    data$.plot_group <- factor(1)
    return(data)
  }
  data$.plot_group <- interaction(data[, group_vars, drop = FALSE], drop = TRUE)
  data
}

plot_comparisons_axis_label <- function(var) {
  if (is.null(var)) {
    return(NULL)
  }
  labels <- c(
    time = "Time",
    state_set = "State set",
    contrast = "Contrast",
    comparison_level = "Comparison level",
    reference_level = "Reference level"
  )
  if (var %in% names(labels)) {
    return(unname(labels[[var]]))
  }
  var
}

plot_comparisons_y_label <- function(comparison) {
  switch(
    comparison,
    difference = "Difference",
    ratio = "Ratio",
    comparison
  )
}

plot_comparisons_title <- function(data) {
  metric <- unique(as.character(data$metric))
  metric <- metric[!is.na(metric)]
  if (length(metric) != 1) {
    return("Average Comparisons")
  }
  metric_label <- switch(
    metric,
    sop = "SOP",
    time_in_state = "Time-in-State",
    time_benefit = "Time-Benefit",
    metric
  )
  paste("Average", metric_label, "Comparisons")
}
