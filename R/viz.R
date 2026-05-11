#' Plot State Occupancy Probabilities Over Time
#'
#' Creates a visualization showing empirical or model-derived state occupancy
#' probabilities over time, with optional faceting by grouping variables
#' (typically treatment). Can display as lines (default) or stacked bars.
#'
#' @param data A data frame containing trajectory data, typically from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`, or a
#'   model-derived SOP summary such as a `markov_avg_sops` object from
#'   `avg_sops()` and `inferences()`.
#' @param time_var Character string. Name of the time variable column (default: "time").
#' @param y_var Character string. Name of the observed state variable column for
#'   empirical trajectory data (default: "y").
#' @param facet_var Character vector or NULL. Name of the variable to facet by
#'   (default: "tx"). Use a length-two vector for `facet_grid(row ~ column)`.
#'   Set to NULL to disable faceting.
#' @param linetype_var Character string or NULL. Name of the variable to map to linetype
#'   (default: NULL). Useful for comparing groups within the same plot. When both
#'   `facet_var` and `linetype_var` are specified, you can facet by one variable
#'   (e.g., site) and distinguish groups by linetype (e.g., treatment).
#' @param geom Character string. Type of plot: "line" for line plot (default)
#'   or "bar" for stacked bar chart.
#' @param prob_var Character string or NULL. Name of the probability column for
#'   model-derived SOP summary data. Defaults to `"estimate"` when present.
#' @param n_draws Integer. Maximum number of stored inference draws to overlay
#'   for `markov_avg_sops` bar plots. Draws are selected deterministically.
#' @param draw_alpha Numeric. Alpha level for draw overlays in model-derived
#'   bar plots.
#' @param ribbon_alpha Numeric. Alpha level for confidence bands in
#'   model-derived line plots.
#' @param show_uncertainty Logical. If TRUE, add confidence ribbons for line
#'   plots and draw overlays for `markov_avg_sops` bar plots when available.
#' @return A ggplot object showing state occupancy probabilities over time.
#'
#' @details
#' This function creates a visualization of how the distribution of patients
#' across states changes over time.
#'
#' For empirical trajectory data, probabilities are computed from observed
#' counts at each time point. For model-derived SOP summaries, the function
#' uses the supplied probability column, usually `estimate`.
#'
#' When `geom = "line"` (default), separate lines show the state occupancy
#' probability for each state over time. Model-derived SOP summaries with
#' `conf.low` and `conf.high` columns also show uncertainty bands.
#'
#' When `geom = "bar"`, each bar represents a time point with colors
#' indicating the state occupancy probability for each state. If `data` is a
#' `markov_avg_sops` object with stored draws from `inferences(...,
#' return_draws = TRUE)`, a deterministic subset of draws is overlaid with low
#' alpha to show uncertainty in the stacked probabilities.
#'
#' Faceting allows comparison between groups (e.g., treatment vs. control).
#' The `linetype_var` parameter allows distinguishing groups by linetype instead
#' of (or in addition to) faceting, enabling more flexible plot layouts.
#'
#' `plot_sops()` adds viridis discrete color and fill scales by default. Add
#' another `scale_color_*()` or `scale_fill_*()` layer to replace them.
#' Summary data must have only one row for each plotted time, state, facet, and
#' linetype combination.
#'
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, lp_function = my_lp)
#'
#' # Basic line plot
#' plot_sops(trajectories)
#'
#' # Stacked bar plot with faceting by treatment
#' plot_sops(trajectories, geom = "bar")
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
#' @importFrom ggplot2 ggplot aes geom_bar geom_col geom_line geom_ribbon
#' @importFrom ggplot2 facet_grid facet_wrap vars labs scale_color_viridis_d
#' @importFrom ggplot2 scale_fill_viridis_d
#' @importFrom rlang .data
#'
#' @export
plot_sops <- function(
  data,
  time_var = "time",
  y_var = "y",
  facet_var = "tx",
  linetype_var = NULL,
  geom = c("line", "bar"),
  prob_var = NULL,
  n_draws = 50,
  draw_alpha = 0.06,
  ribbon_alpha = 0.12,
  show_uncertainty = TRUE
) {
  geom <- match.arg(geom)
  plot_sops_validate_scalar(n_draws, "n_draws", lower = 0)
  plot_sops_validate_scalar(draw_alpha, "draw_alpha", lower = 0, upper = 1)
  plot_sops_validate_scalar(ribbon_alpha, "ribbon_alpha", lower = 0, upper = 1)
  if (
    !is.logical(show_uncertainty) ||
      length(show_uncertainty) != 1 ||
      is.na(show_uncertainty)
  ) {
    stop("`show_uncertainty` must be TRUE or FALSE.")
  }

  if (is.null(prob_var) && "estimate" %in% names(data)) {
    prob_var <- "estimate"
  }
  is_summary <- !is.null(prob_var) &&
    prob_var %in% names(data) &&
    "state" %in% names(data)

  plot_sops_validate_columns(data, time_var, "`time_var`")
  plot_sops_validate_facets(data, facet_var)
  if (!is.null(linetype_var)) {
    plot_sops_validate_columns(data, linetype_var, "`linetype_var`")
  }

  if (is_summary) {
    p <- plot_sops_summary(
      data = data,
      time_var = time_var,
      state_var = "state",
      facet_var = facet_var,
      linetype_var = linetype_var,
      geom = geom,
      prob_var = prob_var,
      n_draws = floor(n_draws),
      draw_alpha = draw_alpha,
      ribbon_alpha = ribbon_alpha,
      show_uncertainty = show_uncertainty
    )
  } else {
    plot_sops_validate_columns(data, y_var, "`y_var`")
    p <- plot_sops_empirical(
      data = data,
      time_var = time_var,
      y_var = y_var,
      facet_var = facet_var,
      linetype_var = linetype_var,
      geom = geom
    )
  }

  p <- plot_sops_add_default_scales(p)
  plot_sops_add_facets(p, facet_var)
}

plot_sops_empirical <- function(
  data,
  time_var,
  y_var,
  facet_var,
  linetype_var,
  geom
) {
  if (geom == "bar") {
    ggplot2::ggplot(data) +
      ggplot2::aes(x = .data[[time_var]], fill = factor(.data[[y_var]])) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::labs(
        x = "Time",
        y = "Probability",
        fill = "State",
        title = "Empirical State Occupancy Probabilities Over Time"
      )
  } else {
    data_summary <- plot_sops_empirical_summary(
      data = data,
      time_var = time_var,
      y_var = y_var,
      facet_var = facet_var,
      linetype_var = linetype_var
    )

    aes_mapping <- ggplot2::aes(
      x = .data[[time_var]],
      y = .data[["prob"]],
      color = factor(.data[[y_var]])
    )

    if (!is.null(linetype_var)) {
      aes_mapping$linetype <- ggplot2::aes(
        linetype = factor(.data[[linetype_var]])
      )$linetype
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

    if (!is.null(linetype_var)) {
      p <- p + ggplot2::labs(linetype = linetype_var)
    }
    p
  }
}

plot_sops_summary <- function(
  data,
  time_var,
  state_var,
  facet_var,
  linetype_var,
  geom,
  prob_var,
  n_draws,
  draw_alpha,
  ribbon_alpha,
  show_uncertainty
) {
  plot_sops_validate_columns(data, state_var, "`state`")
  plot_sops_validate_columns(data, prob_var, "`prob_var`")

  state_levels <- plot_sops_state_levels(data, state_var)
  data[[state_var]] <- factor(as.character(data[[state_var]]), levels = state_levels)
  plot_sops_validate_summary_keys(
    data = data,
    time_var = time_var,
    state_var = state_var,
    facet_var = facet_var,
    linetype_var = linetype_var,
    geom = geom
  )

  if (geom == "bar") {
    draw_data <- NULL
    if (show_uncertainty && n_draws > 0) {
      draw_data <- plot_sops_draw_data(
        data = data,
        time_var = time_var,
        state_var = state_var,
        facet_var = facet_var,
        n_draws = n_draws,
        state_levels = state_levels
      )
    }

    if (
      show_uncertainty &&
        n_draws > 0 &&
        inherits(data, "markov_avg_sops") &&
        is.null(draw_data)
    ) {
      warning(
        "No stored draws found. Run inferences(..., return_draws = TRUE) ",
        "to overlay uncertainty on model-derived stacked bar plots.",
        call. = FALSE
      )
    }

    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = .data[[time_var]],
        y = .data[[prob_var]],
        fill = .data[[state_var]]
      )
    ) +
      ggplot2::geom_col(
        width = 0.85,
        alpha = if (is.null(draw_data)) 1 else 0.65,
        color = NA
      )

    if (!is.null(draw_data)) {
      draw_ids <- sort(unique(draw_data$draw_id))
      for (draw_id in draw_ids) {
        draw_i <- draw_data[draw_data$draw_id == draw_id, , drop = FALSE]
        p <- p +
          ggplot2::geom_col(
            data = draw_i,
            ggplot2::aes(
              x = .data[[time_var]],
              y = .data[[".draw"]],
              fill = .data[[state_var]]
            ),
            inherit.aes = FALSE,
            width = 0.85,
            alpha = draw_alpha,
            color = NA
          )
      }
    }

    p +
      ggplot2::labs(
        x = "Time",
        y = "Probability",
        fill = "State",
        title = "State Occupancy Probabilities Over Time"
      )
  } else {
    aes_mapping <- ggplot2::aes(
      x = .data[[time_var]],
      y = .data[[prob_var]],
      color = .data[[state_var]]
    )

    if (!is.null(linetype_var)) {
      aes_mapping$linetype <- ggplot2::aes(
        linetype = factor(.data[[linetype_var]])
      )$linetype
      aes_mapping$group <- ggplot2::aes(
        group = interaction(.data[[state_var]], factor(.data[[linetype_var]]))
      )$group
    } else {
      aes_mapping$group <- ggplot2::aes(group = .data[[state_var]])$group
    }

    p <- ggplot2::ggplot(data, aes_mapping)

    if (show_uncertainty && all(c("conf.low", "conf.high") %in% names(data))) {
      ribbon_mapping <- ggplot2::aes(
        ymin = .data[["conf.low"]],
        ymax = .data[["conf.high"]],
        fill = .data[[state_var]]
      )
      if (!is.null(linetype_var)) {
        ribbon_mapping$group <- ggplot2::aes(
          group = interaction(.data[[state_var]], factor(.data[[linetype_var]]))
        )$group
      } else {
        ribbon_mapping$group <- ggplot2::aes(group = .data[[state_var]])$group
      }
      p <- p +
        ggplot2::geom_ribbon(
          ribbon_mapping,
          alpha = ribbon_alpha,
          color = NA
        )
    }

    p <- p +
      ggplot2::geom_line() +
      ggplot2::labs(
        x = "Time",
        y = "Probability",
        color = "State",
        fill = "State",
        title = "State Occupancy Probabilities Over Time"
      )

    if (!is.null(linetype_var)) {
      p <- p + ggplot2::labs(linetype = linetype_var)
    }
    p
  }
}

plot_sops_empirical_summary <- function(
  data,
  time_var,
  y_var,
  facet_var,
  linetype_var
) {
  group_vars <- c(time_var, y_var)
  if (!is.null(facet_var)) {
    group_vars <- c(facet_var, group_vars)
  }
  if (!is.null(linetype_var)) {
    group_vars <- c(linetype_var, group_vars)
  }
  group_vars <- unique(group_vars)

  data_summary <- stats::aggregate(
    rep(1L, nrow(data)),
    data[, group_vars, drop = FALSE],
    length
  )
  names(data_summary)[ncol(data_summary)] <- "n"
  denom_vars <- setdiff(group_vars, y_var)
  denom <- stats::aggregate(
    data_summary$n,
    data_summary[, denom_vars, drop = FALSE],
    sum
  )
  names(denom)[ncol(denom)] <- "denom"
  data_summary <- left_join_preserve_order(data_summary, denom, by = denom_vars)
  data_summary$prob <- data_summary$n / data_summary$denom
  data_summary$denom <- NULL
  data_summary
}

plot_sops_draw_data <- function(
  data,
  time_var,
  state_var,
  facet_var,
  n_draws,
  state_levels
) {
  draws <- attr(data, "bootstrap_draws")
  if (is.null(draws)) {
    draws <- attr(data, "simulation_draws")
  }
  if (is.null(draws)) {
    draws <- attr(data, "draws")
  }
  if (is.null(draws) || !"draw_id" %in% names(draws)) {
    return(NULL)
  }

  value_var <- if ("draw" %in% names(draws)) "draw" else "estimate"
  required <- c("draw_id", time_var, state_var, value_var)
  if (!all(required %in% names(draws))) {
    return(NULL)
  }

  facet_cols <- character()
  if (!is.null(facet_var)) {
    if (!all(facet_var %in% names(draws))) {
      return(NULL)
    }
    facet_cols <- facet_var
  }

  keep_ids <- plot_sops_select_draw_ids(draws$draw_id, n_draws)
  if (length(keep_ids) == 0) {
    return(NULL)
  }
  draws <- draws[draws$draw_id %in% keep_ids, , drop = FALSE]
  draws[[state_var]] <- factor(
    as.character(draws[[state_var]]),
    levels = state_levels
  )
  draws$.draw <- draws[[value_var]]

  order_cols <- c("draw_id", facet_cols, time_var, state_var)
  draws[do.call(order, draws[, order_cols, drop = FALSE]), , drop = FALSE]
}

plot_sops_select_draw_ids <- function(draw_ids, n_draws) {
  draw_ids <- sort(unique(draw_ids))
  if (n_draws <= 0 || length(draw_ids) == 0) {
    return(draw_ids[0])
  }
  if (length(draw_ids) <= n_draws) {
    return(draw_ids)
  }
  idx <- unique(round(seq(1, length(draw_ids), length.out = n_draws)))
  draw_ids[idx]
}

plot_sops_state_levels <- function(data, state_var) {
  if (is.factor(data[[state_var]])) {
    return(levels(data[[state_var]]))
  }
  unique(as.character(data[[state_var]]))
}

plot_sops_validate_summary_keys <- function(
  data,
  time_var,
  state_var,
  facet_var,
  linetype_var,
  geom
) {
  key_vars <- c(time_var, state_var, facet_var)
  if (geom == "line") {
    key_vars <- c(key_vars, linetype_var)
  }
  key_vars <- unique(key_vars[!is.na(key_vars)])
  key_data <- data[, key_vars, drop = FALSE]
  duplicate_key <- duplicated(key_data) | duplicated(key_data, fromLast = TRUE)

  if (!any(duplicate_key)) {
    return(invisible(NULL))
  }

  guidance <- if (geom == "line") {
    "Include additional grouping variables in `facet_var` or `linetype_var`, "
  } else {
    "Include additional grouping variables in `facet_var`, "
  }

  stop(
    "Model-derived SOP summary data has multiple rows for the same plotted ",
    "combination of ",
    paste(key_vars, collapse = ", "),
    ". ",
    guidance,
    "aggregate before plotting, or filter to one scenario.",
    call. = FALSE
  )
}

plot_sops_add_facets <- function(p, facet_var) {
  if (is.null(facet_var)) {
    return(p)
  }
  if (length(facet_var) == 1) {
    return(p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]])))
  }
  p + ggplot2::facet_grid(
    rows = ggplot2::vars(.data[[facet_var[1]]]),
    cols = ggplot2::vars(.data[[facet_var[2]]])
  )
}

plot_sops_add_default_scales <- function(p) {
  p +
    ggplot2::scale_color_viridis_d() +
    ggplot2::scale_fill_viridis_d()
}

plot_sops_validate_facets <- function(data, facet_var) {
  if (is.null(facet_var)) {
    return(invisible(NULL))
  }
  if (!is.character(facet_var) || !length(facet_var) %in% c(1, 2)) {
    stop("`facet_var` must be NULL or a character vector of length 1 or 2.")
  }
  plot_sops_validate_columns(data, facet_var, "`facet_var`")
}

plot_sops_validate_columns <- function(data, vars, arg) {
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      arg,
      " column not found in `data`: ",
      paste(missing_vars, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(NULL)
}

plot_sops_validate_scalar <- function(x, arg, lower, upper = Inf) {
  if (
    !is.numeric(x) ||
      length(x) != 1 ||
      is.na(x) ||
      !is.finite(x) ||
      x < lower ||
      x > upper
  ) {
    stop(
      "`",
      arg,
      "` must be a numeric scalar between ",
      lower,
      " and ",
      upper,
      "."
    )
  }
  invisible(NULL)
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
  opt <- utils::modifyList(default_options, ggplot_options)

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
  sops_long <- pivot_state_columns_long(bootstrap_data)

  # Determine grouping variables for summarization
  group_vars <- c(time_var, "state")
  if (!is.null(group_var)) {
    group_vars <- c(group_var, group_vars)
  }
  if (!is.null(facet_var)) {
    group_vars <- c(facet_var, group_vars)
  }
  group_vars <- unique(group_vars)

  # Compute confidence intervals.
  group_key <- do.call(
    interaction,
    c(sops_long[, group_vars, drop = FALSE], drop = TRUE, sep = "\r")
  )
  ci_bands <- bind_rows_fill(lapply(split(seq_len(nrow(sops_long)), group_key), function(idx) {
    group_data <- sops_long[idx, , drop = FALSE]
    out <- group_data[1, group_vars, drop = FALSE]
    out$lower <- as.numeric(stats::quantile(group_data$probability, lower_q))
    out$median <- stats::median(group_data$probability)
    out$upper <- as.numeric(stats::quantile(group_data$probability, upper_q))
    out
  }))

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
