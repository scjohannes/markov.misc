# Linear-predictor contrast plots.

#' Plot Treatment Differences in the Linear Predictor
#'
#' `plot_lp_difference()` evaluates a fitted ordinal Markov model at one patient
#' profile and plots differences in the model linear predictor between
#' counterfactual levels of one variable over time. The plot is conditional on
#' the supplied profile and faceted by previous state.
#'
#' @details
#' `plot_lp_difference()` is model-based only. At minimum, supply `model`, a
#' one-row `profile`, `variables`, and `times`. The profile should contain the
#' model covariates for the patient or scenario being plotted; `times` and
#' `previous_states` define the plotting grid.
#'
#' @param x A fitted `orm`, `blrm`, `vglm`, or `robcov_vglm` Markov
#'   model.
#' @param profile One-row data frame containing the covariate profile.
#' @param variables Named list with one variable and at least two values. The
#'   first value is the reference.
#' @param times Time values at which to evaluate the linear predictor.
#' @param previous_states Optional previous-state values to facet by. Defaults
#'   to all non-absorbing states in `y_levels`.
#' @param y_levels Optional state levels. If omitted, levels are inferred from
#'   the model when possible.
#' @param absorb Optional absorbing-state label excluded from the default
#'   previous-state grid.
#' @param time_var Character name of the time column.
#' @param p_var Character name of the previous-state column.
#' @param gap_var Optional time-gap_var variable used by the fitted model.
#' @param time_covariates Optional time-varying covariate lookup table used by the fitted
#'   model.
#' @param line_width Line width for the plotted contrasts.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_lp_difference(
#'   fit,
#'   profile = data.frame(tx = 0, age = 59, sex = "Male", yprev = 3),
#'   variables = list(tx = c(0, 1)),
#'   times = 1:28,
#'   y_levels = 1:8,
#'   absorb = 8
#' )
#' }
#' @export
plot_lp_difference <- function(
  x,
  profile,
  variables,
  times,
  previous_states = NULL,
  y_levels = NULL,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  gap_var = NULL,
  time_covariates = NULL,
  line_width = 0.7
) {
  model <- x
  if (!inherits(model, c("orm", "blrm", "vglm", "robcov_vglm"))) {
    stop(
      "`model` must be an `orm`, `blrm`, `vglm`, or ",
      "`robcov_vglm` object."
    )
  }
  validate_markov_model(model)
  if (!is.data.frame(profile) || nrow(profile) != 1L) {
    stop("`profile` must be a one-row data frame.")
  }
  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied.")
  }
  variables <- validate_avg_comparison_variables(variables)
  variable <- names(variables)[1]
  if (!variable %in% names(profile)) {
    stop("Variables not in profile: ", variable)
  }
  plot_validate_scalar(line_width, "line_width", lower = 0)

  y_levels <- markov_model_ylevels(model, y_levels)
  ylevel_names <- as_state_labels(y_levels)
  previous_states <- previous_states %||%
    setdiff(ylevel_names, as_state_labels(absorb))

  data <- plot_lp_grid(
    model = model,
    profile = profile,
    variables = variables,
    times = times,
    previous_states = previous_states,
    time_var = time_var,
    p_var = p_var,
    gap_var = gap_var,
    time_covariates = time_covariates
  )
  eta <- markov_linear_predictor_matrix(model, data)
  lp_data <- plot_lp_difference_data(
    data = data,
    eta = eta,
    variable = variable,
    values = variables[[1]],
    time_var = time_var,
    p_var = p_var
  )

  p <- ggplot2::ggplot(lp_data) +
    ggplot2::aes(
      x = .data[[time_var]],
      y = .data$estimate,
      color = .data$threshold,
      group = interaction(.data$threshold, .data$contrast)
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::facet_wrap(ggplot2::vars(.data$previous_state)) +
    ggplot2::labs(
      x = "Time",
      y = "Difference in linear predictor",
      color = "Threshold"
    ) +
    ggplot2::theme_bw()

  if (length(unique(lp_data$contrast)) > 1L) {
    p <- p +
      ggplot2::aes(linetype = .data$contrast) +
      ggplot2::labs(linetype = "Contrast")
  }

  p + ggplot2::scale_color_viridis_d()
}

plot_lp_grid <- function(
  model,
  profile,
  variables,
  times,
  previous_states,
  time_var,
  p_var,
  gap_var,
  time_covariates
) {
  variable <- names(variables)[1]
  time_res <- resolve_sop_times(
    model,
    profile,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "unique"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap_var, time_covariates, time_info)

  rows <- vector(
    "list",
    length(previous_states) * length(times) * length(variables[[1]])
  )
  cursor <- 1L
  for (prev in previous_states) {
    for (i in seq_along(times)) {
      for (j in seq_along(variables[[1]])) {
        value <- variables[[1]][j]
        row <- profile
        row[[variable]] <- value
        row[[p_var]] <- plot_lp_previous_state(
          value = prev,
          model = model,
          profile = profile,
          p_var = p_var
        )
        row <- assign_sop_visit(
          row,
          time_var = time_var,
          times = times,
          index = i,
          time_covariates = time_covariates,
          gap_var = gap_var,
          time_info = time_info
        )
        row$previous_state <- as.character(prev)
        rows[[cursor]] <- row
        cursor <- cursor + 1L
      }
    }
  }

  out <- bind_rows_fill(rows)
  rownames(out) <- NULL
  out
}

plot_lp_previous_state <- function(value, model, profile, p_var) {
  if (p_var %in% names(profile)) {
    return(coerce_previous_state_values(value, profile[[p_var]], p_var))
  }

  model_levels <- get_model_factor_levels(model, p_var)
  if (!is.null(model_levels)) {
    return(factor(as.character(value), levels = model_levels))
  }

  value
}

markov_linear_predictor_matrix <- function(model, newdata) {
  if (inherits(model, "robcov_vglm")) {
    if (is.null(model$vglm_fit)) {
      stop(
        "robcov_vglm object does not contain the original vglm fit. ",
        "Please re-run robcov_vglm() with the latest version of the package."
      )
    }
    return(markov_vglm_link_matrix(model$vglm_fit, newdata))
  }
  if (inherits(model, "vglm")) {
    return(markov_vglm_link_matrix(model, newdata))
  }
  if (inherits(model, "blrm")) {
    return(plot_blrm_lp_median_matrix(model, newdata))
  }
  if (inherits(model, "orm")) {
    gamma <- get_effective_coefs(model)
    x <- orm_model_matrix(model, newdata, include_intercept = TRUE)
    x <- x[, colnames(gamma), drop = FALSE]
    out <- x %*% t(gamma)
    colnames(out) <- rownames(gamma) %||% paste0("eta", seq_len(ncol(out)))
    return(out)
  }

  stop("Unsupported model class.")
}

plot_blrm_lp_median_matrix <- function(model, newdata) {
  y_levels <- markov_model_ylevels(model)
  n_thresholds <- length(y_levels) - 1L
  if (n_thresholds < 1L) {
    stop("`blrm` linear predictor plots require at least two state levels.")
  }

  out <- vapply(
    seq_len(n_thresholds),
    function(kint) {
      as.numeric(plot_blrm_predict(
        model,
        newdata = newdata,
        type = "lp",
        kint = kint,
        posterior.summary = "median",
        cint = FALSE
      ))
    },
    numeric(nrow(newdata))
  )
  out <- as.matrix(out)
  colnames(out) <- paste0("eta", seq_len(ncol(out)))
  rownames(out) <- rownames(newdata)
  out
}

plot_blrm_predict <- function(model, newdata, type, ...) {
  if (!requireNamespace("rmsb", quietly = TRUE)) {
    stop("Package 'rmsb' is required for `blrm` plotting.")
  }

  stats::predict(model, newdata = newdata, type = type, ...)
}

markov_vglm_link_matrix <- function(model, newdata) {
  if (
    !isTRUE(attr(model, "markov_vglm")) &&
      !isTRUE(attr(model, "markov_split_assign"))
  ) {
    out <- VGAM::predict(model, newdata, type = "link")
    return(as.matrix(out))
  }

  tt <- stats::delete.response(stats::terms(model))
  x <- stats::model.matrix(
    tt,
    newdata,
    contrasts.arg = if (length(model@contrasts)) model@contrasts else NULL,
    xlev = model@xlevels
  )
  attr(x, "assign") <- model@assign

  lm2vlm <- utils::getFromNamespace("lm2vlm.model.matrix", "VGAM")
  x_vlm <- lm2vlm(
    x,
    Hlist = model@constraints,
    M = model@misc$M,
    xij = model@control$xij,
    Xm2 = NULL
  )
  eta <- matrix(
    x_vlm %*% VGAM::coefvlm(model),
    nrow = nrow(x),
    ncol = model@misc$M,
    byrow = TRUE
  )
  colnames(eta) <- paste0("eta", seq_len(ncol(eta)))
  eta
}

plot_lp_difference_data <- function(
  data,
  eta,
  variable,
  values,
  time_var,
  p_var
) {
  eta <- as.matrix(eta)
  if (is.null(colnames(eta))) {
    colnames(eta) <- paste0("eta", seq_len(ncol(eta)))
  }

  long <- data[rep(seq_len(nrow(data)), times = ncol(eta)), , drop = FALSE]
  long$threshold <- rep(colnames(eta), each = nrow(data))
  long$lp <- as.vector(eta)
  rownames(long) <- NULL

  reference <- as.character(values[1])
  comparisons <- as.character(values[-1])
  keys <- c(time_var, p_var, "previous_state", "threshold")
  reference_data <- long[as.character(long[[variable]]) == reference, ]

  out <- vector("list", length(comparisons))
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    comparison_data <- long[
      as.character(long[[variable]]) == comparison,
      ,
      drop = FALSE
    ]
    merged <- merge(
      comparison_data,
      reference_data,
      by = keys,
      suffixes = c(".comparison", ".reference"),
      sort = FALSE
    )
    merged$estimate <- merged$lp.comparison - merged$lp.reference
    merged$variable <- variable
    merged$reference_level <- reference
    merged$comparison_level <- comparison
    merged$contrast <- paste0(comparison, " - ", reference)
    out[[i]] <- merged[
      c(
        keys,
        "estimate",
        "variable",
        "reference_level",
        "comparison_level",
        "contrast"
      )
    ]
  }

  out <- bind_rows_fill(out)
  out$previous_state <- factor(
    out$previous_state,
    levels = unique(data$previous_state)
  )
  out
}
