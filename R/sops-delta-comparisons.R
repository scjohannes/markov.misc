# Analytical propagation for average SOP comparisons.

delta_validate_comparison_scope <- function(
  object,
  target,
  vcov,
  conf_type
) {
  args <- attr(object, "comparison_args")
  avg_args <- attr(object, "avg_args")
  if (is.null(args) || is.null(avg_args)) {
    stop(
      "Stored average-comparison arguments are required for delta inference."
    )
  }
  if (!identical(conf_type, "wald")) {
    stop(
      "Analytical average comparisons require `conf_type = \"wald\"`.",
      call. = FALSE
    )
  }
  if (!identical(args$comparison, "difference")) {
    stop(
      "Analytical average comparisons currently support differences only; ",
      "ratios are not yet supported.",
      call. = FALSE
    )
  }
  if (!args$estimand %in% c("sop", "time_in_state")) {
    stop(
      "Analytical average comparisons currently support `estimand = \"sop\"` ",
      "or `\"time_in_state\"`; time benefit is not yet supported.",
      call. = FALSE
    )
  }
  if (length(avg_args$by %||% character()) > 0L) {
    stop(
      "Analytical delta inference with `by` is not yet supported.",
      call. = FALSE
    )
  }

  target <- delta_validate_target(object, target)
  if (identical(target, "superpopulation") && !is.null(vcov)) {
    stop(
      "`vcov` cannot be supplied with `target = \"superpopulation\"`; ",
      "superpopulation inference uses fitted-model score components and the ",
      "stacked patient influence function.",
      call. = FALSE
    )
  }
  list(args = args, avg_args = avg_args, target = target)
}

delta_replay_average_comparison <- function(
  object,
  args,
  avg_args,
  conf_level
) {
  newdata <- if (isTRUE(attr(object, "newdata_supplied"))) {
    attr(object, "newdata_orig")
  } else {
    NULL
  }
  avg_comparison_replay_avg_sops(
    model = attr(object, "model"),
    newdata = newdata,
    refit_data = attr(object, "refit_data"),
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    y_levels = attr(object, "y_levels"),
    absorb = attr(object, "absorb"),
    time_var = attr(object, "time_var") %||% "time",
    p_var = attr(object, "p_var") %||% "yprev",
    id_var = avg_args$id_var,
    p2_var = attr(object, "p2_var"),
    gap_var = attr(object, "gap_var"),
    time_covariates = attr(object, "time_covariates"),
    include_re = args$include_re,
    n_draws = args$n_draws,
    seed = args$seed,
    posterior_summary = args$posterior_summary,
    conf_level = conf_level,
    return_draws = FALSE,
    extra_args = args$extra_args
  )
}

delta_trapezoid_weights <- function(times) {
  n <- length(times)
  if (n < 2L) {
    return(numeric(n))
  }
  steps <- diff(times)
  if (n == 2L) {
    return(rep(steps / 2, 2L))
  }
  c(
    steps[1L] / 2,
    (utils::head(steps, -1L) + utils::tail(steps, -1L)) / 2,
    steps[length(steps)] / 2
  )
}

delta_interpolation_operator <- function(source_times, target_times) {
  n_source <- length(source_times)
  if (n_source == 0L) {
    return(matrix(0, nrow = length(target_times), ncol = 0L))
  }
  if (n_source == 1L) {
    out <- matrix(0, nrow = length(target_times), ncol = 1L)
    out[target_times == source_times, 1L] <- 1
    return(out)
  }

  plan <- compile_linear_interpolation_plan(source_times, target_times)
  if (is.null(plan)) {
    stop("Could not compile the analytical real-time interpolation operator.")
  }
  out <- t(apply_linear_interpolation_plan(diag(n_source), plan))
  out[is.na(out)] <- 0
  out
}

delta_real_time_visit_weights <- function(avg, args) {
  time_map <- standardize_time_map(args$time_map)
  source_labels <- unique(as.character(avg$time))
  source_real <- map_sop_time_values(source_labels, time_map)
  origin <- args$origin %||% "empirical_baseline"
  use_origin <- !is.null(args$origin_time) &&
    identical(origin, "empirical_baseline")

  if (!is.null(args$origin_time)) {
    if (
      !is.numeric(args$origin_time) ||
        length(args$origin_time) != 1L ||
        is.na(args$origin_time) ||
        !is.finite(args$origin_time)
    ) {
      stop("`origin_time` must be a single finite numeric value.")
    }
  }
  mapped_range <- range(time_map$real_time)
  lower <- if (use_origin) args$origin_time else mapped_range[1L]
  upper <- mapped_range[2L]
  if (lower > upper) {
    stop("`origin_time` must not be greater than the largest mapped time.")
  }
  target_times <- args$target_times
  if (is.null(target_times)) {
    target_times <- sort(unique(c(
      if (use_origin) args$origin_time else NULL,
      source_real
    )))
  }
  target_times <- validate_sop_xout(target_times, lower, upper)

  augmented_times <- c(if (use_origin) args$origin_time else NULL, source_real)
  node_times <- sort(unique(augmented_times))
  node_index <- match(source_real, node_times)
  counts <- tabulate(
    match(augmented_times, node_times),
    nbins = length(node_times)
  )
  collapse <- matrix(
    0,
    nrow = length(node_times),
    ncol = length(source_real)
  )
  collapse[cbind(node_index, seq_along(source_real))] <- 1 / counts[node_index]

  interpolate <- delta_interpolation_operator(node_times, target_times)
  target_weights <- delta_trapezoid_weights(target_times)
  visit_weights <- drop(target_weights %*% interpolate %*% collapse)
  stats::setNames(visit_weights, source_labels)
}

delta_comparison_operator <- function(object, avg, args, avg_args) {
  source <- as.data.frame(avg)
  source <- source[,
    setdiff(names(source), setdiff(inference_columns(), "estimate")),
    drop = FALSE
  ]
  variables <- avg_args$variables
  varname <- names(variables)[1L]
  required_source <- c("time", "state", varname, "estimate")
  if (!all(required_source %in% names(source))) {
    stop("The replayed average SOP cells are missing comparison identifiers.")
  }

  state_sets <- normalize_comparison_state_sets(
    args$state_sets,
    attr(avg, "y_levels")
  )
  if (anyDuplicated(names(state_sets))) {
    stop("Analytical comparison state-set names must be unique.")
  }
  state_index <- match(as.character(object$state_set), names(state_sets))
  if (anyNA(state_index)) {
    stop("Comparison state sets do not align with the replayed average SOPs.")
  }

  real_time <- identical(args$estimand, "time_in_state") &&
    (!is.null(args$time_map) || !is.null(args$origin_time))
  visit_weights <- if (real_time) {
    if (is.null(args$time_map)) {
      stop("`time_map` must be supplied for real-time AUC.")
    }
    delta_real_time_visit_weights(avg, args)
  } else if (identical(args$estimand, "time_in_state")) {
    stats::setNames(
      rep(1, length(unique(as.character(source$time)))),
      unique(as.character(source$time))
    )
  } else {
    NULL
  }

  n_result <- nrow(object)
  n_source <- nrow(source)
  delta_assert_bytes(
    as.double(n_result) * n_source * 8,
    "The analytical comparison operator"
  )
  operator <- matrix(0, nrow = n_result, ncol = n_source)
  source_level <- as.character(source[[varname]])
  source_state <- as.character(source$state)
  source_time <- as.character(source$time)

  for (i in seq_len(n_result)) {
    states <- as.character(state_sets[[state_index[i]]])
    level_sign <-
      as.numeric(source_level == as.character(object$comparison_level[i])) -
      as.numeric(source_level == as.character(object$reference_level[i]))
    state_weight <- as.numeric(source_state %in% states)
    time_weight <- if (identical(args$estimand, "sop")) {
      as.numeric(source_time == as.character(object$time[i]))
    } else {
      unname(visit_weights[match(source_time, names(visit_weights))])
    }
    if (anyNA(time_weight)) {
      stop("Comparison times do not align with the replayed average SOPs.")
    }
    operator[i, ] <- level_sign * state_weight * time_weight
  }

  propagated <- drop(operator %*% source$estimate)
  if (
    any(!is.finite(propagated)) ||
      any(abs(propagated - object$estimate) >= 1e-12)
  ) {
    stop(
      "The analytical comparison operator did not reproduce the stored point ",
      "estimates.",
      call. = FALSE
    )
  }
  operator
}

delta_propagate_comparison_state <- function(avg, operator) {
  source <- delta_analytical(avg)
  if (identical(source$representation, "coefficient")) {
    jacobian <- operator %*% source$jacobian
    colnames(jacobian) <- colnames(source$jacobian)
    variance <- delta_jacobian_variance(jacobian, source$coefficient_vcov)
    return(list(
      std.error = sqrt(variance),
      analytical = list(
        representation = "coefficient",
        jacobian = jacobian,
        coefficient_vcov = source$coefficient_vcov,
        coefficient_names = source$coefficient_names %||% colnames(jacobian),
        covariance_metadata = source$covariance_metadata
      )
    ))
  }

  if (!identical(source$representation, "influence")) {
    stop("Unknown analytical representation on replayed average SOPs.")
  }
  influence <- source$influence %*% t(operator)
  n <- nrow(influence)
  if (n < 2L) {
    stop("Superpopulation delta inference requires at least two patients.")
  }
  centered <- sweep(influence, 2L, colMeans(influence), "-")
  variance <- colSums(centered^2) / ((n - 1) * n)
  average_jacobian <- operator %*% source$average_jacobian
  colnames(average_jacobian) <- colnames(source$average_jacobian)
  list(
    std.error = sqrt(variance),
    analytical = list(
      representation = "influence",
      average_jacobian = average_jacobian,
      influence = influence,
      profile_ids = source$profile_ids,
      coefficient_names = source$coefficient_names %||%
        colnames(average_jacobian),
      covariance_metadata = source$covariance_metadata
    )
  )
}

inferences_delta_comparisons <- function(
  object,
  target,
  vcov,
  cluster,
  conf_level,
  conf_type
) {
  scope <- delta_validate_comparison_scope(
    object = object,
    target = target,
    vcov = vcov,
    conf_type = conf_type
  )
  avg <- delta_replay_average_comparison(
    object = object,
    args = scope$args,
    avg_args = scope$avg_args,
    conf_level = conf_level
  )
  avg <- inferences(
    avg,
    method = "delta",
    vcov = vcov,
    cluster = cluster,
    conf_level = conf_level,
    conf_type = "wald",
    target = scope$target,
    return_draws = FALSE
  )
  operator <- delta_comparison_operator(
    object = object,
    avg = avg,
    args = scope$args,
    avg_args = scope$avg_args
  )
  propagated <- delta_propagate_comparison_state(avg, operator)
  delta_finalize_result(
    object = object,
    standard_error = propagated$std.error,
    conf_level = conf_level,
    conf_type = "wald",
    target = scope$target,
    analytical = propagated$analytical
  )
}
