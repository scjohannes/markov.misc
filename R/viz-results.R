# Operating-characteristic plotting.

#' Plot Operating Characteristics
#'
#' Creates one plot for each requested operating-characteristic column.
#'
#' @param x A data frame of operating-characteristic results.
#' @param outcomes Character vector of columns to plot on the y axes.
#' @param x_var Character name of the x-axis column.
#' @param group_var Character name of the grouping column.
#' @param ggplot_options List of ggplot layers or scales added to every plot.
#' @param combine If `TRUE`, return a combined patchwork; otherwise return a
#'   named list of ggplot objects.
#'
#' @return A patchwork object when `combine = TRUE`, otherwise a named list of
#'   ggplot objects.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @export
plot_operchar <- function(
  x,
  outcomes,
  x_var,
  group_var,
  ggplot_options = list(),
  combine = TRUE
) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame.")
  }
  for (arg in c("outcomes", "x_var", "group_var")) {
    value <- get(arg)
    if (!is.character(value) || length(value) == 0L || any(!nzchar(value))) {
      stop("`", arg, "` must contain character column names.")
    }
  }
  if (length(x_var) != 1L || length(group_var) != 1L) {
    stop("`x_var` and `group_var` must each name one column.")
  }
  missing <- setdiff(c(outcomes, x_var, group_var), names(x))
  if (length(missing) > 0L) {
    stop("Columns not found in `x`: ", paste(missing, collapse = ", "))
  }
  if (!is.factor(x[[group_var]])) {
    stop("`group_var` must identify a factor column.")
  }
  if (!is.list(ggplot_options)) {
    stop("`ggplot_options` must be a list.")
  }
  if (!is.logical(combine) || length(combine) != 1L || is.na(combine)) {
    stop("`combine` must be TRUE or FALSE.")
  }

  plots <- lapply(outcomes, function(outcome) {
    geoms <- if (is.factor(x[[x_var]])) {
      list(
        ggplot2::geom_point(stat = "summary", fun = sum),
        ggplot2::stat_summary(fun = sum, geom = "line")
      )
    } else {
      list(ggplot2::geom_point(), ggplot2::geom_line())
    }
    plot <- ggplot2::ggplot(
      x,
      ggplot2::aes(
        x = .data[[x_var]],
        y = .data[[outcome]],
        color = .data[[group_var]],
        group = .data[[group_var]]
      )
    )
    layers <- c(
      geoms,
      list(ggplot2::theme_minimal(), ggplot2::scale_color_viridis_d()),
      ggplot_options
    )
    for (layer in layers) {
      plot <- plot + layer
    }
    plot + ggplot2::labs(y = outcome)
  })
  names(plots) <- outcomes

  if (isTRUE(combine)) {
    return(patchwork::wrap_plots(plotlist = plots))
  }
  plots
}
