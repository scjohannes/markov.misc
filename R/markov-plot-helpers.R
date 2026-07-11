# Internal fitted-model plotting helpers.

markov_supported_model <- function(object) {
  inherits(object, c("orm", "blrm", "vglm", "vgam", "robcov_vglm"))
}

markov_model_ylevels <- function(model, y_levels = NULL) {
  if (!is.null(y_levels)) {
    return(y_levels)
  }

  if (inherits(model, "vglm")) {
    out <- tryCatch(model@extra$colnames.y, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  if (inherits(model, "robcov_vglm")) {
    out <- model$extra$colnames.y %||%
      tryCatch(model$vglm_fit@extra$colnames.y, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  if (inherits(model, "blrm")) {
    out <- model$ylevels
    if (!is.null(out)) {
      return(out)
    }
  }
  if (inherits(model, "orm")) {
    out <- model$yunique
    if (!is.null(out)) {
      return(out)
    }
  }

  stop("`y_levels` cannot be NULL.")
}

validate_plot_counterfactual_variables <- function(
  variables,
  require_two = FALSE
) {
  if (is.null(variables)) {
    return(NULL)
  }
  if (!is.list(variables) || is.null(names(variables))) {
    stop("`variables` must be a named list, for example `list(tx = c(0, 1))`.")
  }
  if (length(variables) != 1L || !nzchar(names(variables)[1])) {
    stop("`variables` must contain exactly one named variable.")
  }
  if (length(variables[[1]]) < 1L) {
    stop("`variables[[1]]` must contain at least one value.")
  }
  if (require_two && length(variables[[1]]) < 2L) {
    stop("`variables[[1]]` must contain at least two values.")
  }
  if (anyDuplicated(as.character(variables[[1]]))) {
    stop("`variables[[1]]` must not contain duplicate values.")
  }
  variables
}

complete_plot_recursion_times <- function(times, time_info) {
  if (isTRUE(time_info$is_factor)) {
    labels <- as.character(times)
    idx <- match(labels, time_info$levels)
    if (anyNA(idx)) {
      stop("Could not align plot times with fitted visit levels.")
    }
    return(factor(
      time_info$levels[seq_len(max(idx))],
      levels = time_info$levels,
      ordered = isTRUE(time_info$ordered)
    ))
  }

  if (is.numeric(times) || is.integer(times)) {
    numeric_times <- as.numeric(times)
    whole_number_times <- all(
      is.finite(numeric_times) &
        abs(numeric_times - round(numeric_times)) < sqrt(.Machine$double.eps)
    )
    if (whole_number_times) {
      start <- min(c(1, numeric_times), na.rm = TRUE)
      end <- max(numeric_times, na.rm = TRUE)
      return(seq.int(start, end))
    }
  }

  sort(unique(times))
}
