#' Plot State Occupancy Probabilities Over Time
#'
#' Creates a visualization showing empirical or model-derived state occupancy
#' probabilities over time, with optional faceting by grouping variables
#' (typically treatment). Can display as lines (default) or stacked bars.
#'
#' @param x A data frame containing trajectory data, typically from
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
#' @param estimate_var Character string or NULL. Name of the probability column for
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
  x,
  time_var = "time",
  y_var = "y",
  facet_var = "tx",
  linetype_var = NULL,
  geom = c("line", "bar"),
  estimate_var = NULL,
  n_draws = 50,
  draw_alpha = 0.06,
  ribbon_alpha = 0.12,
  show_uncertainty = TRUE
) {
  data <- x
  prob_var <- estimate_var
  geom <- match.arg(geom)
  plot_validate_scalar(n_draws, "n_draws", lower = 0)
  plot_validate_scalar(draw_alpha, "draw_alpha", lower = 0, upper = 1)
  plot_validate_scalar(ribbon_alpha, "ribbon_alpha", lower = 0, upper = 1)
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

  plot_validate_columns(data, time_var, "`time_var`")
  plot_validate_facets(data, facet_var)
  if (!is.null(linetype_var)) {
    plot_validate_columns(data, linetype_var, "`linetype_var`")
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
    plot_validate_columns(data, y_var, "`y_var`")
    p <- plot_sops_empirical(
      data = data,
      time_var = time_var,
      y_var = y_var,
      facet_var = facet_var,
      linetype_var = linetype_var,
      geom = geom
    )
  }

  p <- plot_add_default_scales(p)
  plot_add_facets(p, facet_var)
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
  plot_validate_columns(data, state_var, "`state`")
  plot_validate_columns(data, prob_var, "`prob_var`")

  state_levels <- plot_sops_state_levels(data, state_var)
  data[[state_var]] <- factor(
    as.character(data[[state_var]]),
    levels = state_levels
  )
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
  draws <- attr(data, "draws")
  if (is.null(draws)) {
    draws <- attr(data, "draws")
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
  y_levels <- attr(data, "y_levels", exact = TRUE)
  if (!is.null(y_levels)) {
    return(as_state_labels(y_levels))
  }

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
