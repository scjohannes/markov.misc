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
