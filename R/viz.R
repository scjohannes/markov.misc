#' Plot State Occupancy Proportions Over Time
#'
#' Creates a stacked bar chart showing the proportion of patients in each state
#' over time, with optional faceting by a grouping variable (typically treatment).
#'
#' @param data A data frame containing trajectory data, typically from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`.
#' @param time_var Character string. Name of the time variable column (default: "time").
#' @param y_var Character string. Name of the state variable column (default: "y").
#' @param facet_var Character string. Name of the variable to facet by (default: "tx").
#'
#' @return A ggplot object showing state occupancy proportions over time.
#'
#' @details
#' This function creates a visualization of how the distribution of patients
#' across states changes over time. Each bar represents a time point, with
#' colors indicating the proportion in each state. Faceting allows comparison
#' between groups (e.g., treatment vs. control).
#'
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, lp_function = my_lp)
#'
#' # Basic plot
#' plot_sops(trajectories)
#'
#' # Customize variable names
#' plot_sops(trajectories, time_var = "day", y_var = "state", facet_var = "treatment")
#' }
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap vars labs
#' @importFrom rlang .data
#'
#' @export
plot_sops <- function(
  data,
  time_var = "time",
  y_var = "y",
  facet_var = "tx"
) {
  ggplot2::ggplot(
    data,
    aes(x = .data[[time_var]], fill = factor(.data[[y_var]]))
  ) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::facet_wrap(vars(.data[[facet_var]])) +
    ggplot2::labs(
      x = "Time",
      y = "Proportion",
      fill = "State",
      title = "State Occupancy Proportions Over Time"
    )
}



#' ggplot wrapper to present simulation results
#'
#' A result data frame with operating chracteristic estimates is plotted in a
#' unified way.
#'
#' @param data A data frame of results.
#' @param ... Operating characteristics to be displayed on the y-axis (in different plots)
#' @param x What to plot on the x axis.
#' @param group Typically different analysis strategies
#' @param ggplot_options List of further options passed to ggplot
#' @param combine Should the plots for different operating characteristics
#                     be combined into a plot grid (using cowplot)?
#'
#'
#' @return A ggplot or cowplot
#'
#' @details
#' The input is typically a data frame where in each row the operating characteristics
#' of a strategy in a certain scenario are contained.
#'
#' ... specifies what to plot on the y-axis, which can be different inputs -
#' typically operating characteristics - and a separate plot for each input is generated.
#'
#' x determines the column of values to be plotted on the x-axis.
#'
#' group will typically be the analysis strategies compared and will be lines in
#' different colours. When combined = TRUE, the single plots and a cowplot::plot_grid
#' with all plots is returned.
#'
#' ggplot_options is integrated into a list of ggplot options and directly passed
#' to the ggplot calls.
#'
#' @import ggplot2
#' @importFrom cowplot plot_grid
#'
#'
#' @examples
#' \dontrun{
#' ...
#' }
#'
#'
#' @export
plot_results <- function(data, ..., x, group, ggplot_options = list(),
                         combine = TRUE) {

  group_name <- as.character(substitute(group))
  x_name <- as.character(substitute(x))

  # Capture y variables
  y_vars <- as.character(substitute(list(...)))[-1]

  # Check
  if (length(group_name) != 1) {
    stop("Variable name of '", group_name, "' not recognizable")
  }
  if (!group_name %in% names(data)) {
    stop("Variable '", group_name, "' not found in data")
  }
  if (!is.factor(data[[group_name]])) {
    stop("Variable '", group_name, "' must be a factor")
  }

  # Check x variable and determine plot type
  if (!x_name %in% names(data)) {
    stop("Variable '", x_name, "' not found in data")
  }

  # Plot defaults
  default_design <- list(
    theme_minimal(),
    scale_color_viridis_d())

  # Based on x type, choose plot function
  if(is.factor(data[[x_name]])) {
    plot_fun <- list(
      geom_point(stat='summary', fun=sum),
      stat_summary(fun=sum, geom="line"))
  } else {
    plot_fun <- list(geom_point(), geom_line())
  }

  default_options <- c(plot_fun, default_design)

  # Update options by user input
  opt <- modifyList(default_options, ggplot_options)

  # Loop over the operating characteristics and plot
  plots <- list()
  for(i in 1:length(y_vars)) {
    y_col <- y_vars[i]
    plots[[i]] <-
      ggplot(data, aes(x = {{ x }}, y = .data[[y_col]],
                       color = {{ group }}, group = {{ group }})) +
      opt +
      labs(y = y_col)
  }

  if(combine == TRUE) {
    grid <- cowplot::plot_grid(plotlist = plots)
  } else {grid <- NULL}

  return(
    list(single_plots = plots,
         grid = grid)
  )
}

