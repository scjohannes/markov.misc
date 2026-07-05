# Correlation heatmaps and variograms.

#' Plot Empirical or Model-Based Correlations Over Time
#'
#' `plot_correlation()` computes Spearman correlations between ordinal states
#' at pairs of follow-up times. With a fitted Markov model, trajectories are
#' simulated from the model for each patient in their observed group and the
#' correlations are computed from those simulated paths.
#'
#' @param object A trajectory data frame or a fitted Markov model.
#' @param newdata Optional data frame of prediction profiles for model-based
#'   plots. If `NULL`, wrapper-fitted models use their stored data and extract
#'   one prediction row per ID.
#' @param refit_data Optional full longitudinal data used only for stored-data
#'   resolution in model-based plots.
#' @param times Optional time values to include. Required for model-based plots.
#' @param ylevels Optional state levels for model-based simulation.
#' @param absorb Optional absorbing-state label for model-based simulation.
#' @param id_var Character ID column for raw data and stored model-data
#'   extraction.
#' @param time_var Character name of the time column.
#' @param y_var Character name of the state column in raw trajectory data.
#' @param pvarname Character name of the previous-state column for model
#'   simulation.
#' @param p2varname Optional second previous-state column for model simulation.
#' @param facet_var Optional grouping variable. Correlations are computed
#'   separately within each observed stratum.
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
#' @param method Correlation method passed to [stats::cor()].
#' @param use Missing-data handling passed to [stats::cor()].
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
  include_re = FALSE,
  n_rep = 20L,
  n_draws = 100L,
  seed = NULL,
  method = "spearman",
  use = "pairwise.complete.obs",
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

  data <- plot_correlation_input_data(
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
    gap = gap,
    t_covs = t_covs,
    include_re = include_re,
    n_rep = n_rep,
    n_draws = n_draws,
    seed = seed
  )

  corr <- plot_correlation_data(
    data = data$data,
    id_var = data$id_var,
    time_var = time_var,
    y_var = data$y_var,
    facet_var = facet_var,
    method = method,
    use = use,
    triangle = triangle
  )

  plot_correlation_heatmap(
    corr,
    show_values = show_values,
    digits = floor(digits),
    fill_limits = fill_limits
  )
}

#' Plot an Empirical or Model-Based Correlation Variogram
#'
#' `plot_variogram()` summarizes the same time-by-time correlation matrix used
#' by [plot_correlation()] as correlations against absolute time differences.
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
  include_re = FALSE,
  n_rep = 20L,
  n_draws = 100L,
  seed = NULL,
  method = "spearman",
  use = "pairwise.complete.obs",
  smooth = TRUE
) {
  if (!is.logical(smooth) || length(smooth) != 1L || is.na(smooth)) {
    stop("`smooth` must be TRUE or FALSE.")
  }

  data <- plot_correlation_input_data(
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
    gap = gap,
    t_covs = t_covs,
    include_re = include_re,
    n_rep = n_rep,
    n_draws = n_draws,
    seed = seed
  )

  corr <- plot_correlation_data(
    data = data$data,
    id_var = data$id_var,
    time_var = time_var,
    y_var = data$y_var,
    facet_var = facet_var,
    method = method,
    use = use,
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

plot_correlation_input_data <- function(
  object,
  newdata,
  refit_data,
  times,
  ylevels,
  absorb,
  id_var,
  time_var,
  y_var,
  pvarname,
  p2varname,
  gap,
  t_covs,
  include_re,
  n_rep,
  n_draws,
  seed
) {
  if (markov_supported_model(object)) {
    sim <- markov_simulate_fitted_paths(
      model = object,
      newdata = newdata,
      refit_data = refit_data,
      variables = NULL,
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
    return(list(data = sim, id_var = ".sim_id", y_var = "y"))
  }

  if (!is.data.frame(object)) {
    stop("`object` must be a data frame or a supported Markov model.")
  }
  if (!is.null(newdata) || !is.null(refit_data)) {
    stop("`newdata` and `refit_data` are only used for model-based plots.")
  }

  data <- object
  plot_validate_columns(data, c(id_var, time_var, y_var), "`data`")
  if (!is.null(times)) {
    keep <- as.character(data[[time_var]]) %in% as.character(times)
    data <- data[keep, , drop = FALSE]
  }
  list(data = data, id_var = id_var, y_var = y_var)
}

plot_correlation_data <- function(
  data,
  id_var,
  time_var,
  y_var,
  facet_var,
  method,
  use,
  triangle
) {
  plot_validate_columns(data, c(id_var, time_var, y_var), "`data`")
  plot_validate_facets(data, facet_var)

  groups <- plot_facet_groups(data, facet_var)
  out <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    subset <- plot_group_subset(data, groups[i, , drop = FALSE], facet_var)
    mat <- plot_state_time_matrix(
      subset,
      id_var = id_var,
      time_var = time_var,
      y_var = y_var
    )
    corr <- stats::cor(mat, method = method, use = use)
    out[[i]] <- plot_correlation_matrix_long(
      corr,
      group = groups[i, , drop = FALSE],
      facet_var = facet_var,
      triangle = triangle
    )
  }

  out <- bind_rows_fill(out)
  if (length(facet_var) > 0L) {
    out$.panel <- plot_generic_panel(out, facet_var)
    out$.panel <- factor(out$.panel, levels = unique(out$.panel))
  }
  out
}

plot_facet_groups <- function(data, facet_var) {
  if (is.null(facet_var)) {
    return(data.frame(.all = 1L))
  }
  unique(data[facet_var])
}

plot_group_subset <- function(data, group, facet_var) {
  if (is.null(facet_var)) {
    return(data)
  }
  keep <- rep(TRUE, nrow(data))
  for (var in facet_var) {
    keep <- keep & as.character(data[[var]]) == as.character(group[[var]][1])
  }
  data[keep, , drop = FALSE]
}

plot_state_time_matrix <- function(data, id_var, time_var, y_var) {
  ids <- unique(as.character(data[[id_var]]))
  times <- plot_ordered_values(data[[time_var]])
  time_labels <- as.character(times)
  id_index <- match(as.character(data[[id_var]]), ids)
  time_index <- match(as.character(data[[time_var]]), time_labels)

  key <- paste(id_index, time_index, sep = "\r")
  if (anyDuplicated(key)) {
    stop(
      "`data` must contain at most one row per `",
      id_var,
      "` and `",
      time_var,
      "` combination within each correlation stratum."
    )
  }

  mat <- matrix(
    NA_real_,
    nrow = length(ids),
    ncol = length(time_labels),
    dimnames = list(ids, time_labels)
  )
  mat[cbind(id_index, time_index)] <- plot_state_scores(data[[y_var]])
  mat
}

plot_state_scores <- function(x) {
  if (is.factor(x)) {
    return(as.integer(x))
  }
  out <- suppressWarnings(as.integer(as.character(x)))
  if (anyNA(out) && any(!is.na(x))) {
    levels <- sort(unique(as.character(x)))
    out <- as.integer(factor(as.character(x), levels = levels))
  }
  out
}

plot_correlation_matrix_long <- function(corr, group, facet_var, triangle) {
  keep <- switch(
    triangle,
    upper = upper.tri(corr, diag = FALSE),
    full = row(corr) != col(corr)
  )
  idx <- which(keep, arr.ind = TRUE)
  out <- data.frame(
    time_1 = rownames(corr)[idx[, 1]],
    time_2 = colnames(corr)[idx[, 2]],
    correlation = corr[idx],
    stringsAsFactors = FALSE
  )
  if (!is.null(facet_var)) {
    for (var in facet_var) {
      out[[var]] <- as.character(group[[var]][1])
    }
  }
  out$time_1 <- factor(out$time_1, levels = colnames(corr))
  out$time_2 <- factor(out$time_2, levels = colnames(corr))
  out
}

plot_generic_panel <- function(data, facet_var) {
  do.call(
    paste,
    c(
      lapply(facet_var, function(var) paste0(var, "=", data[[var]])),
      sep = ", "
    )
  )
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

plot_variogram_data <- function(corr) {
  time_1 <- as.character(corr$time_1)
  time_2 <- as.character(corr$time_2)
  t1 <- suppressWarnings(as.numeric(time_1))
  t2 <- suppressWarnings(as.numeric(time_2))
  if (anyNA(t1) || anyNA(t2)) {
    levels <- unique(c(time_1, time_2))
    t1 <- match(time_1, levels)
    t2 <- match(time_2, levels)
  }
  corr$delta <- abs(t1 - t2)
  corr
}
