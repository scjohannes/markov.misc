#' Plot Empirical State Occupancy Probabilities Over Time
#'
#' Creates a visualization showing the empirical state occupancy probabilities
#' over time, with optional faceting by a grouping variable (typically treatment).
#' Can display as stacked bars (default) or as lines.
#'
#' @param data A data frame containing trajectory data, typically from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`.
#' @param time_var Character string. Name of the time variable column (default: "time").
#' @param y_var Character string. Name of the state variable column (default: "y").
#' @param facet_var Character string or NULL. Name of the variable to facet by (default: "tx").
#'   Set to NULL to disable faceting.
#' @param linetype_var Character string or NULL. Name of the variable to map to linetype
#'   (default: NULL). Useful for comparing groups within the same plot. When both
#'   `facet_var` and `linetype_var` are specified, you can facet by one variable
#'   (e.g., site) and distinguish groups by linetype (e.g., treatment).
#' @param geom Character string. Type of plot: "bar" for stacked bar chart (default)
#'   or "line" for line plot.
#'
#' @return A ggplot object showing empirical state occupancy probabilities over time.
#'
#' @details
#' This function creates a visualization of how the distribution of patients
#' across states changes over time.
#'
#' When `geom = "bar"` (default), each bar represents a time point with
#' colors indicating the empirical state occupancy probability for each state.
#'
#' When `geom = "line"`, separate lines show the empirical state occupancy
#' probability for each state over time, which can be useful for seeing trends
#' in individual states.
#'
#' Faceting allows comparison between groups (e.g., treatment vs. control).
#' The `linetype_var` parameter allows distinguishing groups by linetype instead
#' of (or in addition to) faceting, enabling more flexible plot layouts.
#'
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, lp_function = my_lp)
#'
#' # Basic stacked bar plot
#' plot_sops(trajectories)
#'
#' # Line plot with faceting by treatment
#' plot_sops(trajectories, geom = "line")
#'
#' # Line plot with linetype by treatment (no faceting)
#' plot_sops(trajectories, geom = "line", facet_var = NULL, linetype_var = "tx")
#'
#' # Facet by one variable, linetype by another
#' plot_sops(trajectories, facet_var = "site", linetype_var = "tx", geom = "line")
#'
#' # Customize variable names
#' plot_sops(trajectories, time_var = "day", y_var = "state", facet_var = "treatment")
#' }
#' @importFrom ggplot2 ggplot aes geom_bar geom_line facet_wrap vars labs
#' @importFrom dplyr group_by summarise mutate ungroup n across all_of
#' @importFrom rlang .data
#'
#' @export
plot_sops <- function(
  data,
  time_var = "time",
  y_var = "y",
  facet_var = "tx",
  linetype_var = NULL,
  geom = c("bar", "line")
) {
  # Match argument
  geom <- match.arg(geom)

  # Base plot
  p <- ggplot2::ggplot(data)

  if (geom == "bar") {
    # Stacked bar chart
    p <- p +
      ggplot2::aes(x = .data[[time_var]], fill = factor(.data[[y_var]])) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(
        x = "Time",
        y = "Probability",
        fill = "State",
        title = "Empirical State Occupancy Probabilities Over Time"
      )
  } else {
    # Line plot - compute empirical probabilities at each time point
    # Determine grouping variables for summarization
    group_vars <- c(time_var, y_var)
    if (!is.null(facet_var)) {
      group_vars <- c(facet_var, group_vars)
    }
    if (!is.null(linetype_var)) {
      group_vars <- c(linetype_var, group_vars)
    }
    group_vars <- unique(group_vars)

    # Compute proportions by group variables and time
    data_summary <- data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop_last") |>
      dplyr::mutate(prob = n / sum(n)) |>
      dplyr::ungroup()

    # Build base aesthetic mapping
    aes_mapping <- ggplot2::aes(
      x = .data[[time_var]],
      y = .data[["prob"]],
      color = factor(.data[[y_var]])
    )

    # Add linetype mapping if specified
    if (!is.null(linetype_var)) {
      aes_mapping$linetype <- ggplot2::aes(
        linetype = factor(.data[[linetype_var]])
      )$linetype
      # Update group to include both state and linetype
      aes_mapping$group <- ggplot2::aes(
        group = interaction(
          factor(.data[[y_var]]),
          factor(.data[[linetype_var]])
        )
      )$group
    } else {
      aes_mapping$group <- ggplot2::aes(group = factor(.data[[y_var]]))$group
    }

    p <- ggplot2::ggplot(data_summary, aes_mapping) +
      ggplot2::geom_line() +
      ggplot2::labs(
        x = "Time",
        y = "Probability",
        color = "State",
        title = "Empirical State Occupancy Probabilities Over Time"
      )

    # Add linetype label if used
    if (!is.null(linetype_var)) {
      p <- p + ggplot2::labs(linetype = linetype_var)
    }
  }

  # Add faceting if specified
  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(vars(.data[[facet_var]]))
  }

  return(p)
}


#' ggplot wrapper to present simulation results
#'
#' A result data frame with operating characteristics estimates is plotted in a
#' unified way.
#'
#' @param data A data frame of results.
#' @param ... Operating characteristics to be displayed on the y-axis (in different plots)
#' @param x What to plot on the x axis.
#' @param group Typically different analysis strategies
#' @param ggplot_options List of further options passed to ggplot
#' @param combine Should the plots for different operating characteristics
#                     be combined into a grid?
#'
#'
#' @return A ggplot list and optionally a patchwork object
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
#' different colours. When combined = TRUE, the single plots and a patchwork::wrap_plots
#' with all plots is returned.
#'
#' ggplot_options is integrated into a list of ggplot options and directly passed
#' to the ggplot calls.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#'
#'
#' @examples
#' \dontrun{
#' ...
#' }
#'
#'
#' @export
plot_results <- function(
  data,
  ...,
  x,
  group,
  ggplot_options = list(),
  combine = TRUE
) {
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
    scale_color_viridis_d()
  )

  # Based on x type, choose plot function
  if (is.factor(data[[x_name]])) {
    plot_fun <- list(
      geom_point(stat = 'summary', fun = sum),
      stat_summary(fun = sum, geom = "line")
    )
  } else {
    plot_fun <- list(geom_point(), geom_line())
  }

  default_options <- c(plot_fun, default_design)

  # Update options by user input
  opt <- modifyList(default_options, ggplot_options)

  # Loop over the operating characteristics and plot
  plots <- list()
  for (i in 1:length(y_vars)) {
    y_col <- y_vars[i]
    plots[[i]] <-
      ggplot(
        data,
        aes(
          x = {{ x }},
          y = .data[[y_col]],
          color = {{ group }},
          group = {{ group }}
        )
      ) +
      opt +
      labs(y = y_col)
  }

  if (combine == TRUE) {
    grid <- patchwork::wrap_plots(plotlist = plots)
  } else {
    grid <- NULL
  }

  return(
    list(single_plots = plots, grid = grid)
  )
}


#' Plot Bootstrap Standardized State Occupancy Probabilities with Confidence Bands
#'
#' Creates a visualization showing standardized state occupancy probabilities
#' from bootstrap results with confidence intervals. Plots median probabilities
#' as lines with shaded confidence bands.
#'
#' @param bootstrap_data A data frame from `bootstrap_standardized_sops()` containing
#'   bootstrap samples with columns for time, treatment, state probabilities, and boot_id.
#' @param time_var Character string. Name of the time variable column (default: "time").
#' @param group_var Character string or NULL. Name of the grouping variable (default: "tx").
#'   Typically treatment or other comparison variable. Set to NULL if no grouping.
#' @param facet_var Character string or NULL. Name of the variable to facet by (default: NULL).
#'   Set to a variable name to create separate panels.
#' @param conf_level Numeric. Confidence level for intervals (default: 0.95).
#' @param title Character string. Plot title (default: auto-generated).
#'
#' @return A ggplot object showing standardized SOPs with confidence bands.
#'
#' @details
#' This function takes the output from `bootstrap_standardized_sops()` and creates
#' a publication-ready plot showing:
#' - Median state occupancy probabilities over time (lines)
#' - Confidence bands (shaded ribbons)
#' - Optional grouping by linetype (e.g., treatment vs. control)
#' - Optional faceting by another variable
#'
#' The function automatically:
#' 1. Reshapes wide-format bootstrap data to long format
#' 2. Computes median and quantile-based confidence intervals
#' 3. Creates lines colored by state
#' 4. Adds ribbons for confidence bands
#' 5. Uses linetype to distinguish groups if `group_var` is specified
#'
#' @examples
#' \dontrun{
#' # After bootstrapping
#' bs_sops <- bootstrap_standardized_sops(
#'   model = m1,
#'   data = data,
#'   times = 1:30,
#'   n_boot = 100,
#'   ylevels = factor(1:6),
#'   absorb = "6"
#' )
#'
#' # Basic plot with treatment grouping
#' plot_bootstrap_sops(bs_sops)
#'
#' # Custom confidence level
#' plot_bootstrap_sops(bs_sops, conf_level = 0.90)
#'
#' # With faceting by another variable
#' plot_bootstrap_sops(bs_sops, facet_var = "site")
#'
#' # Without grouping
#' plot_bootstrap_sops(bs_sops, group_var = NULL)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs facet_wrap vars
#' @importFrom dplyr group_by summarise across all_of
#' @importFrom tidyr pivot_longer starts_with
#' @importFrom rlang .data
#'
#' @export
plot_bootstrap_sops <- function(
  bootstrap_data,
  time_var = "time",
  group_var = "tx",
  facet_var = NULL,
  conf_level = 0.95,
  title = NULL
) {
  # Calculate quantile probabilities
  alpha <- 1 - conf_level
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2

  # Reshape to long format
  sops_long <- bootstrap_data |>
    tidyr::pivot_longer(
      cols = tidyr::starts_with("state_"),
      names_to = "state",
      values_to = "probability",
      names_prefix = "state_"
    )

  # Determine grouping variables for summarization
  group_vars <- c(time_var, "state")
  if (!is.null(group_var)) {
    group_vars <- c(group_var, group_vars)
  }
  if (!is.null(facet_var)) {
    group_vars <- c(facet_var, group_vars)
  }
  group_vars <- unique(group_vars)

  # Compute confidence intervals
  ci_bands <- sops_long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      lower = quantile(.data[["probability"]], lower_q),
      median = median(.data[["probability"]]),
      upper = quantile(.data[["probability"]], upper_q),
      .groups = "drop"
    )

  # Build base aesthetic mapping
  aes_mapping <- ggplot2::aes(
    x = .data[[time_var]],
    y = .data[["median"]],
    color = .data[["state"]]
  )

  # Add grouping and linetype if group_var is specified
  if (!is.null(group_var)) {
    aes_mapping$linetype <- ggplot2::aes(
      linetype = factor(.data[[group_var]])
    )$linetype
    aes_mapping$group <- ggplot2::aes(
      group = interaction(.data[["state"]], factor(.data[[group_var]]))
    )$group
  } else {
    aes_mapping$group <- ggplot2::aes(group = .data[["state"]])$group
  }

  # Create plot
  p <- ggplot2::ggplot(ci_bands, aes_mapping) +
    ggplot2::geom_line()

  # Add ribbon with appropriate grouping
  if (!is.null(group_var)) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["lower"]],
          ymax = .data[["upper"]],
          fill = .data[["state"]],
          group = interaction(.data[["state"]], factor(.data[[group_var]]))
        ),
        alpha = 0.2,
        color = NA
      )
  } else {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["lower"]],
          ymax = .data[["upper"]],
          fill = .data[["state"]]
        ),
        alpha = 0.2,
        color = NA
      )
  }

  # Generate default title if not provided
  if (is.null(title)) {
    ci_pct <- round(conf_level * 100)
    title <- paste0(
      "Standardized State Occupation Probabilities with Bootstrapped ",
      ci_pct,
      "% Confidence Bands"
    )
  }

  # Add labels
  p <- p +
    ggplot2::labs(
      x = "Time",
      y = "State Occupation Probability",
      color = "State",
      fill = "State",
      title = title
    )

  # Add linetype label if used
  if (!is.null(group_var)) {
    p <- p + ggplot2::labs(linetype = group_var)
  }

  # Add faceting if specified
  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]]))
  }

  return(p)
}
