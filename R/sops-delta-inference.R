# Analytical delta-method inference for SOP objects.

delta_model_for_plan <- function(model) {
  if (inherits(model, "robcov_vglm")) model$vglm_fit else model
}

delta_validate_target <- function(object, target) {
  if (inherits(object, "markov_sops")) {
    target <- target %||% "fixed"
    if (!identical(target, "fixed")) {
      stop(
        "`sops()` delta inference supports only `target = \"fixed\"`. ",
        "Empirical and superpopulation targets apply to averaged SOP objects.",
        call. = FALSE
      )
    }
    return(target)
  }

  target <- target %||% "empirical"
  if (!target %in% c("empirical", "superpopulation")) {
    stop(
      "Averaged delta inference supports `target = \"empirical\"` or ",
      "`target = \"superpopulation\"`.",
      call. = FALSE
    )
  }
  if (
    identical(target, "superpopulation") &&
      isTRUE(attr(object, "newdata_supplied"))
  ) {
    stop(
      "`target = \"superpopulation\"` requires the stored fitted-patient ",
      "cohort. User-supplied `newdata` is treated as a fixed cohort.",
      call. = FALSE
    )
  }
  target
}

delta_validate_sop_scope <- function(object, model, by) {
  if (!is.null(by)) {
    stop(
      "Analytical delta inference with `by` is not yet supported.",
      call. = FALSE
    )
  }
  if (!is.null(attr(object, "p2_var"))) {
    stop(
      "Analytical delta inference for second-order Markov models is not yet ",
      "supported.",
      call. = FALSE
    )
  }
  model_plan <- delta_model_for_plan(model)
  if (!inherits(model_plan, c("orm", "vglm")) || inherits(model_plan, "blrm")) {
    stop(
      "Analytical delta inference currently supports frequentist `orm`, ",
      "`vglm`, and `robcov_vglm` models.",
      call. = FALSE
    )
  }
  invisible(model_plan)
}

delta_compile_plan <- function(object, model_plan, newdata, times, output) {
  plan <- compile_sop_execution_plan(
    model = model_plan,
    newdata = newdata,
    times = times,
    y_levels = attr(object, "y_levels"),
    absorb = attr(object, "absorb"),
    time_var = attr(object, "time_var") %||% "time",
    p_var = attr(object, "p_var") %||% "yprev",
    p2_var = NULL,
    gap_var = attr(object, "gap_var"),
    time_covariates = attr(object, "time_covariates"),
    builder = "streamed",
    output = output
  )

  Gamma <- get_effective_coefs(model_plan)
  po <- markov_po_structure(
    Gamma,
    plan$components$col_names,
    plan$components$X_init
  )
  if (is.null(po)) {
    stop(
      "Analytical delta inference currently requires a full ",
      "proportional-odds model; partial proportional odds are not yet ",
      "supported.",
      call. = FALSE
    )
  }
  plan
}

delta_validate_plan_result <- function(result, n_profiles, n_times, n_states) {
  if (
    !is.list(result) ||
      is.null(result$probabilities) ||
      is.null(result$jacobian)
  ) {
    stop("The analytical SOP engine returned an invalid result.", call. = FALSE)
  }
  expected_prob <- c(n_profiles, n_times, n_states)
  if (!identical(dim(result$probabilities), expected_prob)) {
    stop("The analytical SOP probabilities have unexpected dimensions.")
  }
  jac_dim <- dim(result$jacobian)
  if (
    length(jac_dim) != 4L ||
      !identical(jac_dim[seq_len(3L)], expected_prob)
  ) {
    stop("The analytical SOP Jacobian has unexpected dimensions.")
  }

  coefficient_names <- result$coef_names %||% names(result$coef)
  if (is.null(coefficient_names) || length(coefficient_names) != jac_dim[4L]) {
    stop("The analytical SOP Jacobian is missing raw coefficient names.")
  }
  if (anyDuplicated(coefficient_names)) {
    stop("The analytical SOP Jacobian coefficient names must be unique.")
  }
  list(
    probabilities = result$probabilities,
    jacobian = result$jacobian,
    coefficient_names = coefficient_names,
    individual_probabilities = result$individual_probabilities
  )
}

delta_flatten_jacobian <- function(jacobian, coefficient_names) {
  out <- matrix(
    as.numeric(jacobian),
    nrow = prod(dim(jacobian)[seq_len(3L)]),
    ncol = dim(jacobian)[4L]
  )
  colnames(out) <- coefficient_names
  out
}

delta_match_cells <- function(object, calculated, key_cols) {
  object_key <- sop_draw_cell_key(object, key_cols)
  calculated_key <- sop_draw_cell_key(calculated, key_cols)
  index <- match(object_key, calculated_key)
  if (anyNA(index)) {
    stop("Analytical SOP cells do not align with the original point estimate.")
  }
  object_estimate <- object$estimate
  calculated_estimate <- calculated$estimate[index]
  if (
    !is.numeric(object_estimate) ||
      !is.numeric(calculated_estimate) ||
      length(object_estimate) != length(calculated_estimate) ||
      any(!is.finite(object_estimate)) ||
      any(!is.finite(calculated_estimate)) ||
      any(abs(object_estimate - calculated_estimate) >= 1e-12)
  ) {
    stop(
      "Analytical SOP replay did not reproduce the stored point estimates.",
      call. = FALSE
    )
  }
  index
}

delta_align_coefficient_vcov <- function(result, coefficient_names) {
  V <- result$vcov
  if (!is.matrix(V)) {
    stop("The analytical coefficient covariance must be a matrix.")
  }
  if (
    is.null(rownames(V)) ||
      is.null(colnames(V)) ||
      anyDuplicated(rownames(V)) ||
      anyDuplicated(colnames(V)) ||
      !setequal(coefficient_names, rownames(V)) ||
      !setequal(coefficient_names, colnames(V))
  ) {
    stop(
      "The analytical coefficient covariance names must match every raw ",
      "model coefficient.",
      call. = FALSE
    )
  }
  V <- V[coefficient_names, coefficient_names, drop = FALSE]
  if (any(!is.finite(V))) {
    stop("The analytical coefficient covariance must contain finite values.")
  }
  (V + t(V)) / 2
}

delta_coefficient_state <- function(model, cluster, vcov, jacobian) {
  covariance <- get_delta_cluster_vcov(
    model = model,
    cluster = cluster,
    vcov = vcov
  )
  V <- delta_align_coefficient_vcov(covariance, colnames(jacobian))
  variance <- delta_jacobian_variance(jacobian, V)
  list(
    std.error = sqrt(variance),
    coefficient_vcov = V,
    metadata = covariance$metadata %||% list(source = "coefficient_vcov")
  )
}

delta_jacobian_variance <- function(jacobian, coefficient_vcov) {
  working_bytes <- as.double(nrow(jacobian)) * ncol(jacobian) * 8
  delta_assert_bytes(working_bytes, "Analytical standard-error calculation")
  JV <- jacobian %*% coefficient_vcov
  variance <- rowSums(JV * jacobian)
  tolerance <- sqrt(.Machine$double.eps) * pmax(1, rowSums(abs(JV * jacobian)))
  variance[variance < 0 & abs(variance) <= tolerance] <- 0
  if (any(variance < 0)) {
    stop("Analytical propagation produced a negative variance.", call. = FALSE)
  }
  variance
}

delta_superpopulation_cluster <- function(cluster) {
  cluster
}

delta_profile_ids <- function(baseline_data, id_var) {
  if (is.null(id_var) || !id_var %in% names(baseline_data)) {
    stop(
      "Superpopulation delta inference requires the stored starting-profile ",
      "ID column.",
      call. = FALSE
    )
  }
  ids <- as.character(baseline_data[[id_var]])
  if (anyNA(ids) || anyDuplicated(ids)) {
    stop(
      "Superpopulation delta inference requires exactly one non-missing ",
      "starting profile per fitted patient.",
      call. = FALSE
    )
  }
  ids
}

delta_superpopulation_state <- function(
  individual_values,
  average_jacobian,
  profile_ids,
  score_components
) {
  cell_names <- paste0("cell_", seq_len(ncol(individual_values)))
  colnames(individual_values) <- cell_names
  rownames(average_jacobian) <- cell_names
  stacked <- delta_stacked_influence(
    individual_values = individual_values,
    average_jacobian = average_jacobian,
    profile_ids = profile_ids,
    score_components = score_components
  )
  influence <- as.matrix(stacked$influence)
  if (
    nrow(influence) != length(profile_ids) ||
      ncol(influence) != nrow(average_jacobian)
  ) {
    stop("The stacked SOP influence values have unexpected dimensions.")
  }
  colnames(influence) <- NULL
  standard_error <- stacked$std.error
  if (is.null(standard_error) || length(standard_error) != ncol(influence)) {
    if (nrow(influence) < 2L) {
      stop("Superpopulation delta inference requires at least two patients.")
    }
    standard_error <- sqrt(diag(stats::cov(influence) / nrow(influence)))
  }
  metadata <- stacked$metadata %||%
    score_components$metadata %||%
    list()
  metadata$source <- metadata$source %||% "stacked_influence"
  list(
    std.error = as.numeric(standard_error),
    influence = influence,
    average_jacobian = average_jacobian,
    metadata = metadata
  )
}

delta_interval_bounds <- function(
  estimate,
  standard_error,
  conf_level,
  conf_type
) {
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  if (identical(conf_type, "wald")) {
    return(list(
      conf.low = estimate - critical * standard_error,
      conf.high = estimate + critical * standard_error
    ))
  }
  if (!identical(conf_type, "logit")) {
    stop(
      "Delta SOP inference supports `conf_type = \"logit\"` or `\"wald\"`.",
      call. = FALSE
    )
  }

  lower <- upper <- rep(NA_real_, length(estimate))
  interior <- is.finite(estimate) & estimate > 0 & estimate < 1
  if (any(interior)) {
    logit_se <- standard_error[interior] /
      (estimate[interior] * (1 - estimate[interior]))
    center <- stats::qlogis(estimate[interior])
    lower[interior] <- stats::plogis(center - critical * logit_se)
    upper[interior] <- stats::plogis(center + critical * logit_se)
  }

  zero_se <- standard_error <= sqrt(.Machine$double.eps)
  structural_zero <- estimate == 0 & zero_se
  structural_one <- estimate == 1 & zero_se
  lower[structural_zero] <- upper[structural_zero] <- 0
  lower[structural_one] <- upper[structural_one] <- 1

  undefined <- !interior & !structural_zero & !structural_one
  if (any(undefined)) {
    warning(
      "Logit-delta limits are undefined for nonstructural boundary SOP ",
      "estimates; their confidence limits were set to NA.",
      call. = FALSE
    )
  }
  list(conf.low = lower, conf.high = upper)
}

delta_copy_result_attrs <- function(result, object) {
  for (name in names(attributes(object))) {
    if (!name %in% c("names", "row.names", "class")) {
      attr(result, name) <- attr(object, name)
    }
  }
  class(result) <- class(object)
  result
}

delta_finalize_result <- function(
  object,
  standard_error,
  conf_level,
  conf_type,
  target,
  analytical
) {
  # Check the retained analytical state before copying the result data frame or
  # attaching any of its potentially large attributes.
  analytical <- delta_attach_byte_accounting(analytical)
  bounds <- delta_interval_bounds(
    estimate = object$estimate,
    standard_error = standard_error,
    conf_level = conf_level,
    conf_type = conf_type
  )
  result <- as.data.frame(object)
  old_inference <- setdiff(inference_columns(), "estimate")
  result <- result[, setdiff(names(result), old_inference), drop = FALSE]
  result$conf.low <- bounds$conf.low
  result$conf.high <- bounds$conf.high
  result$std.error <- standard_error
  result <- delta_copy_result_attrs(result, object)

  for (name in c(
    "draws",
    "n_draws",
    "n_boot",
    "n_successful",
    "engine",
    "score_weight_dist",
    "draw_weights_attached",
    "draw_weight_col",
    "draw_weight_omission_reason"
  )) {
    attr(result, name) <- NULL
  }

  keys <- delta_row_key_frame(result)
  analytical$version <- 1L
  analytical$target <- target
  analytical$interval <- conf_type
  analytical$row_keys <- keys
  analytical$row_key <- delta_row_key(keys)
  if (anyDuplicated(analytical$row_key)) {
    stop("Analytical result-row keys are not unique.")
  }

  attr(result, "analytical") <- analytical
  attr(result, "method") <- "delta"
  attr(result, "target") <- target
  attr(result, "conf_level") <- conf_level
  attr(result, "conf_type") <- conf_type
  attr(
    result,
    "covariance_source"
  ) <- analytical$covariance_metadata$source %||%
    analytical$representation
  attr(result, "covariance_metadata") <- analytical$covariance_metadata
  result <- order_estimate_columns(result)
  class(result) <- class(object)
  result
}

delta_individual_cells <- function(
  probabilities,
  times,
  y_levels,
  newdata
) {
  array_to_df_individual(
    sops_array = probabilities,
    times = times,
    y_levels = y_levels,
    newdata = newdata,
    by = NULL
  )
}

delta_average_arrays <- function(
  probabilities,
  jacobian,
  grid,
  times,
  y_levels,
  variables,
  n_each,
  individual_probabilities = NULL,
  grouped = FALSE
) {
  n_cf <- nrow(grid)
  q <- dim(jacobian)[4L]
  if (isTRUE(grouped)) {
    if (
      !identical(
        dim(probabilities),
        c(n_cf, length(times), length(y_levels))
      ) ||
        !identical(
          dim(jacobian)[seq_len(3L)],
          c(n_cf, length(times), length(y_levels))
        )
    ) {
      stop("Grouped analytical SOP arrays have unexpected dimensions.")
    }
    marginal_probabilities <- probabilities
    marginal_n_each <- 1L
  } else {
    marginal_probabilities <- probabilities
    marginal_n_each <- n_each
    individual_probabilities <- probabilities
  }
  cells <- marginalize_sops_array(
    sops_array = marginal_probabilities,
    grid = grid,
    times = times,
    y_levels = y_levels,
    variables = variables,
    n_cf = n_cf,
    n_each = marginal_n_each
  )
  average_jacobian <- matrix(NA_real_, nrow = nrow(cells), ncol = q)
  individual_values <- if (is.null(individual_probabilities)) {
    NULL
  } else {
    matrix(NA_real_, nrow = n_each, ncol = nrow(cells))
  }
  cursor <- 1L
  for (scenario in seq_len(n_cf)) {
    if (isTRUE(grouped)) {
      scenario_jac <- jacobian[scenario, , , , drop = FALSE]
    } else {
      rows <- ((scenario - 1L) * n_each + 1L):(scenario * n_each)
      scenario_jac <- jacobian[rows, , , , drop = FALSE]
    }
    n_cells <- length(times) * length(y_levels)
    output_rows <- cursor:(cursor + n_cells - 1L)
    average_jacobian[output_rows, ] <- if (isTRUE(grouped)) {
      matrix(scenario_jac, nrow = n_cells, ncol = q)
    } else {
      matrix(
        apply(scenario_jac, c(2L, 3L, 4L), mean),
        nrow = n_cells,
        ncol = q
      )
    }
    if (!is.null(individual_values)) {
      rows <- ((scenario - 1L) * n_each + 1L):(scenario * n_each)
      scenario_prob <- individual_probabilities[rows, , , drop = FALSE]
      individual_values[, output_rows] <- matrix(
        scenario_prob,
        nrow = n_each,
        ncol = n_cells
      )
    }
    cursor <- cursor + n_cells
  }
  list(
    cells = cells,
    jacobian = average_jacobian,
    individual_values = individual_values
  )
}

inferences_delta_sops <- function(
  object,
  target,
  vcov,
  cluster,
  conf_level,
  conf_type
) {
  is_average <- inherits(object, "markov_avg_sops")
  target <- delta_validate_target(object, target)
  if (identical(target, "superpopulation") && !is.null(vcov)) {
    stop(
      "`vcov` cannot be supplied with `target = \"superpopulation\"`; ",
      "superpopulation inference uses fitted-model score components and the ",
      "stacked patient influence function.",
      call. = FALSE
    )
  }
  model <- attr(object, "model")
  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform delta inference.")
  }
  avg_args <- attr(object, "avg_args")
  call_args <- attr(object, "call_args")
  by <- if (is_average) avg_args$by else call_args$by %||% attr(object, "by")
  model_plan <- delta_validate_sop_scope(object, model, by)
  score_components <- if (identical(target, "superpopulation")) {
    get_delta_score_components(
      model = model,
      cluster = delta_superpopulation_cluster(cluster)
    )
  } else {
    NULL
  }
  newdata <- attr(object, "newdata_pred")
  if (is.null(newdata) || !is.data.frame(newdata)) {
    stop("Stored SOP prediction data are required for delta inference.")
  }
  times <- if (is_average) avg_args$times else call_args$times
  y_levels <- attr(object, "y_levels")
  if (is_average) {
    variables <- avg_args$variables
    grid <- do.call(expand.grid, variables)
    n_cf <- nrow(grid)
    if (nrow(newdata) %% n_cf != 0L) {
      stop(
        "Stored prediction data are not aligned with the counterfactual grid."
      )
    }
    n_each <- as.integer(nrow(newdata) / n_cf)
  } else {
    variables <- grid <- NULL
    n_cf <- n_each <- NULL
  }
  plan <- delta_compile_plan(
    object = object,
    model_plan = model_plan,
    newdata = newdata,
    times = times,
    output = if (is_average) "average" else "individual"
  )
  engine_result <- delta_validate_plan_result(
    run_sop_delta_plan(
      plan,
      model_plan,
      average_group_size = if (is_average) n_each else NULL,
      retain_individual_probabilities = identical(target, "superpopulation")
    ),
    n_profiles = if (is_average) n_cf else nrow(newdata),
    n_times = length(times),
    n_states = length(y_levels)
  )

  if (!is_average) {
    calculated <- delta_individual_cells(
      engine_result$probabilities,
      times,
      y_levels,
      newdata
    )
    index <- delta_match_cells(object, calculated, c("rowid", "time", "state"))
    jacobian <- delta_flatten_jacobian(
      engine_result$jacobian,
      engine_result$coefficient_names
    )[index, , drop = FALSE]
    coefficient <- delta_coefficient_state(model, cluster, vcov, jacobian)
    analytical <- list(
      representation = "coefficient",
      jacobian = jacobian,
      coefficient_vcov = coefficient$coefficient_vcov,
      coefficient_names = colnames(jacobian),
      covariance_metadata = coefficient$metadata
    )
    return(delta_finalize_result(
      object,
      standard_error = coefficient$std.error,
      conf_level = conf_level,
      conf_type = conf_type,
      target = target,
      analytical = analytical
    ))
  }

  reduced <- delta_average_arrays(
    probabilities = engine_result$probabilities,
    jacobian = engine_result$jacobian,
    grid = grid,
    times = times,
    y_levels = y_levels,
    variables = variables,
    n_each = n_each,
    individual_probabilities = engine_result$individual_probabilities,
    grouped = TRUE
  )
  key_cols <- c("time", "state", names(variables))
  index <- delta_match_cells(object, reduced$cells, key_cols)
  average_jacobian <- reduced$jacobian[index, , drop = FALSE]
  colnames(average_jacobian) <- engine_result$coefficient_names

  if (identical(target, "superpopulation")) {
    if (is.null(reduced$individual_values)) {
      stop(
        "Superpopulation analytical SOP inference lost individual ",
        "probabilities."
      )
    }
    individual_values <- reduced$individual_values[, index, drop = FALSE]
    baseline_data <- newdata[seq_len(n_each), , drop = FALSE]
    id_var <- avg_args$id_var %||% attr(object, "id_var")
    profile_ids <- delta_profile_ids(baseline_data, id_var)
    superpopulation <- delta_superpopulation_state(
      individual_values = individual_values,
      average_jacobian = average_jacobian,
      profile_ids = profile_ids,
      score_components = score_components
    )
    analytical <- list(
      representation = "influence",
      average_jacobian = superpopulation$average_jacobian,
      influence = superpopulation$influence,
      profile_ids = profile_ids,
      coefficient_names = colnames(average_jacobian),
      covariance_metadata = superpopulation$metadata
    )
    standard_error <- superpopulation$std.error
  } else {
    coefficient <- delta_coefficient_state(
      model,
      cluster,
      vcov,
      average_jacobian
    )
    analytical <- list(
      representation = "coefficient",
      jacobian = average_jacobian,
      coefficient_vcov = coefficient$coefficient_vcov,
      coefficient_names = colnames(average_jacobian),
      covariance_metadata = coefficient$metadata
    )
    standard_error <- coefficient$std.error
  }

  delta_finalize_result(
    object,
    standard_error = standard_error,
    conf_level = conf_level,
    conf_type = conf_type,
    target = target,
    analytical = analytical
  )
}
