# SOP result shaping helpers.

validate_sops_by <- function(by, data, data_arg = "newdata") {
  if (is.null(by)) {
    return(invisible(NULL))
  }
  if (!is.character(by)) {
    stop("'by' must be a character vector of variable names.")
  }
  missing_vars <- setdiff(by, names(data))
  if (length(missing_vars) > 0) {
    stop(
      "Variables specified in 'by' not found in ",
      data_arg,
      ": ",
      paste(missing_vars, collapse = ", ")
    )
  }
  invisible(NULL)
}

sops_call_args <- function(
  times,
  ylevels,
  absorb,
  tvarname,
  pvarname,
  p2varname,
  gap,
  t_covs,
  by = NULL,
  ...
) {
  c(
    list(
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      by = by
    ),
    list(...)
  )
}

set_sops_attrs <- function(
  result,
  class_name,
  model,
  call_args,
  tvarname,
  pvarname,
  p2varname,
  ylevels,
  absorb,
  gap,
  t_covs,
  by = NULL,
  newdata_orig = NULL,
  avg_args = NULL,
  extra_attrs = list()
) {
  attr(result, "model") <- model
  attr(result, "call_args") <- call_args
  attr(result, "tvarname") <- tvarname
  attr(result, "pvarname") <- pvarname
  attr(result, "p2varname") <- p2varname
  attr(result, "ylevels") <- ylevels
  attr(result, "absorb") <- absorb
  attr(result, "gap") <- gap
  attr(result, "t_covs") <- t_covs
  attr(result, "by") <- by
  attr(result, "newdata_orig") <- newdata_orig
  attr(result, "avg_args") <- avg_args

  for (nm in names(extra_attrs)) {
    if (!is.null(extra_attrs[[nm]])) {
      attr(result, nm) <- extra_attrs[[nm]]
    }
  }

  class(result) <- c(class_name, setdiff(class(result), class_name))
  result
}

restore_sops_attrs <- function(result, object) {
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  class(result) <- class(object)
  result
}

#' Create Counterfactual Datasets for G-Computation
#'
#' Creates copies of baseline data with treatment variable set to each level.
#'
#' @param baseline_data Data frame with one row per patient.
#' @param grid Data frame of variable combinations.
#' @param variables Named list of variable values.
#'
#' @return Data frame with counterfactual data stacked.
#'
#' @keywords internal
create_counterfactual_data <- function(baseline_data, grid, variables) {
  cf_data_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt_copy <- baseline_data
    for (v in names(grid)) {
      dt_copy[[v]] <- grid[i, v]
    }
    cf_data_list[[i]] <- dt_copy
  }
  do.call(rbind, cf_data_list)
}

validate_sops_weights <- function(weights, n, arg = "weights") {
  if (is.null(weights)) {
    return(NULL)
  }
  if (
    !is.numeric(weights) ||
      length(weights) != n ||
      anyNA(weights) ||
      any(weights < 0)
  ) {
    stop(
      "`",
      arg,
      "` must be a non-negative numeric vector of length ",
      n,
      "."
    )
  }
  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    stop("`", arg, "` must have a positive finite sum.")
  }
  weights
}

weighted_sop_mean <- function(x, weights) {
  keep <- !is.na(x)
  if (!any(keep)) {
    return(NA_real_)
  }

  x <- x[keep]
  weights <- weights[keep]
  weight_sum <- sum(weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    return(NA_real_)
  }

  sum(x * weights) / weight_sum
}

aggregate_sops_estimates <- function(result, group_cols, weight_col = NULL) {
  missing_vars <- setdiff(group_cols, names(result))
  if (length(missing_vars) > 0) {
    stop("Grouping variables missing: ", paste(missing_vars, collapse = ", "))
  }

  if (is.null(weight_col)) {
    agg_formula <- stats::as.formula(
      paste("estimate ~", paste(group_cols, collapse = " + "))
    )
    return(stats::aggregate(
      agg_formula,
      data = result,
      FUN = mean,
      na.rm = TRUE
    ))
  }

  if (!weight_col %in% names(result)) {
    stop("Weight column '", weight_col, "' not found in result.")
  }

  split_key <- interaction(result[, group_cols, drop = FALSE], drop = TRUE)
  groups <- split(seq_len(nrow(result)), split_key, drop = TRUE)
  first_rows <- vapply(groups, `[`, integer(1), 1L)
  out <- result[first_rows, group_cols, drop = FALSE]
  out$estimate <- vapply(
    groups,
    function(idx) {
      weighted_sop_mean(result$estimate[idx], result[[weight_col]][idx])
    },
    numeric(1)
  )
  rownames(out) <- NULL
  out
}

validate_prediction_weight_ids <- function(prediction_data, id_var) {
  if (is.null(id_var) || !id_var %in% names(prediction_data)) {
    stop(
      "Prediction data must contain unique `id_var` values before bootstrap ",
      "draw weights can be attached."
    )
  }

  ids <- as.character(prediction_data[[id_var]])
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop(
      "Prediction data must contain exactly one non-missing row per `id_var` ",
      "before bootstrap draw weights can be attached."
    )
  }

  invisible(ids)
}

validate_draw_weight_column_available <- function(newdata, weight_col) {
  if (is.null(weight_col)) {
    return(invisible(NULL))
  }
  if (weight_col %in% names(newdata)) {
    stop(
      "Prediction data already contains a column named `",
      weight_col,
      "`, which is reserved for bootstrap draw weights. Rename this column ",
      "before calling inferences().",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Marginalize SOPs Array Over Patients
#'
#' Averages individual-level SOPs to get population-average (marginal) SOPs.
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param grid Data frame of variable combinations.
#' @param times Vector of time points.
#' @param ylevels Vector of state levels.
#' @param variables Named list of variables.
#' @param n_cf Number of counterfactual scenarios.
#' @param n_each Number of patients per scenario.
#' @param weights Optional patient-level weights for one empirical cohort.
#' @param by Optional character vector of variables to aggregate by.
#' @param newdata Counterfactual prediction data, required when `by` is used.
#'
#' @return Data frame with marginalized SOPs.
#'
#' @keywords internal
marginalize_sops_array <- function(
  sops_array,
  grid,
  times,
  ylevels,
  variables,
  n_cf,
  n_each,
  weights = NULL,
  by = NULL,
  newdata = NULL
) {
  n_times <- length(times)
  n_states <- length(ylevels)
  weights <- validate_sops_weights(weights, n_each)

  if (!is.null(by)) {
    if (is.null(newdata) || !is.data.frame(newdata)) {
      stop("`newdata` must be supplied when `by` is used.")
    }
    if (nrow(newdata) != n_cf * n_each) {
      stop("`newdata` is not aligned with the counterfactual SOP array.")
    }
    validate_sops_by(by, newdata)

    expanded_weights <- if (is.null(weights)) {
      NULL
    } else {
      rep(weights, times = n_cf)
    }
    weight_col <- if (is.null(expanded_weights)) NULL else ".markov_misc_weight"

    result <- array_to_df_individual(
      sops_array = sops_array,
      times = times,
      ylevels = ylevels,
      newdata = newdata,
      by = NULL,
      weights = expanded_weights,
      weight_col = weight_col
    )
    group_cols <- unique(c("time", "state", names(variables), by))
    return(aggregate_sops_estimates(
      result,
      group_cols,
      weight_col = weight_col
    ))
  }

  # Average within each counterfactual group
  avg_sops_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    start_idx <- (cf_i - 1) * n_each + 1
    end_idx <- cf_i * n_each

    sops_cf <- sops_array[start_idx:end_idx, , , drop = FALSE]
    if (is.null(weights)) {
      avg_sops_mat <- apply(sops_cf, c(2, 3), mean)
    } else {
      weights_norm <- weights / sum(weights)
      sops_cf_mat <- matrix(
        aperm(sops_cf, c(2, 3, 1)),
        nrow = n_times * n_states,
        ncol = n_each
      )
      avg_sops_mat <- matrix(
        as.vector(sops_cf_mat %*% weights_norm),
        nrow = n_times,
        ncol = n_states
      )
    }

    avg_sops_list[[cf_i]] <- avg_sops_mat
  }

  # Format as data frame
  result_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    avg_sops_mat <- avg_sops_list[[cf_i]]

    df <- expand.grid(time = times, state = ylevels)
    df$estimate <- as.vector(avg_sops_mat)

    for (v in names(grid)) {
      df[[v]] <- grid[cf_i, v]
    }

    result_list[[cf_i]] <- df
  }

  bind_rows_fill(result_list)
}


#' Convert SOPs Array to Individual-Level Data Frame
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param times Vector of time points.
#' @param ylevels Vector of state levels.
#' @param newdata Original data with rowid.
#' @param by Optional character vector of variables to aggregate by.
#' @param weights Optional patient-level weights aligned to `newdata`.
#' @param weight_col Optional output column name for draw-specific weights.
#'
#' @return Data frame with individual SOPs (or aggregated if by is specified).
#'
#' @keywords internal
array_to_df_individual <- function(
  sops_array,
  times,
  ylevels,
  newdata,
  by = NULL,
  weights = NULL,
  weight_col = NULL
) {
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]
  weights <- validate_sops_weights(weights, n_pat)

  # Flatten array
  probs_flat <- as.vector(sops_array)

  # Construct indices
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)

  # Build result by repeating newdata rows
  result <- newdata[idx_pat, , drop = FALSE]

  # Add rowid if it exists
  if ("rowid" %in% names(newdata)) {
    result$rowid <- newdata$rowid[idx_pat]
  } else {
    result$rowid <- idx_pat
  }

  result$time <- idx_time
  result$state <- idx_state
  result$estimate <- probs_flat
  if (!is.null(weights)) {
    weight_col <- weight_col %||% ".markov_misc_weight"
    validate_draw_weight_column_available(newdata, weight_col)
    result[[weight_col]] <- weights[idx_pat]
  }
  rownames(result) <- NULL

  # Apply stratified aggregation if 'by' is specified
  if (!is.null(by)) {
    # Validate that by variables exist
    missing_vars <- setdiff(by, names(result))
    if (length(missing_vars) > 0) {
      warning(
        "Variables in 'by' not found in result: ",
        paste(missing_vars, collapse = ", "),
        ". Skipping aggregation."
      )
      return(result)
    }

    # Group by time, state, and stratification variables, then average
    group_cols <- unique(c("time", "state", by))

    result <- aggregate_sops_estimates(
      result,
      group_cols,
      weight_col = if (is.null(weights)) NULL else weight_col
    )
  }

  result
}
