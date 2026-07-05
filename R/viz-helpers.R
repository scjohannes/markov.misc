# Shared ggplot helper utilities.

plot_add_facets <- function(p, facet_var) {
  if (is.null(facet_var)) {
    return(p)
  }
  if (length(facet_var) == 1) {
    return(p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]])))
  }
  p +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data[[facet_var[1]]]),
      cols = ggplot2::vars(.data[[facet_var[2]]])
    )
}

plot_add_default_scales <- function(p) {
  p +
    ggplot2::scale_color_viridis_d() +
    ggplot2::scale_fill_viridis_d()
}

plot_validate_facets <- function(data, facet_var) {
  if (is.null(facet_var)) {
    return(invisible(NULL))
  }
  if (!is.character(facet_var) || !length(facet_var) %in% c(1, 2)) {
    stop("`facet_var` must be NULL or a character vector of length 1 or 2.")
  }
  plot_validate_columns(data, facet_var, "`facet_var`")
}

plot_validate_columns <- function(data, vars, arg) {
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

plot_validate_scalar <- function(x, arg, lower, upper = Inf) {
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
