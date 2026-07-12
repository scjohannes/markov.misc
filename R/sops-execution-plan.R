# Internal serializable execution plans for frequentist SOPs.

validate_sop_threshold_count <- function(Gamma, y_levels) {
  expected <- length(y_levels) - 1L
  actual <- nrow(Gamma)
  if (!identical(actual, expected)) {
    stop(
      "`y_levels` defines ",
      length(y_levels),
      " states, but the fitted model defines ",
      actual + 1L,
      " states through ",
      actual,
      " threshold coefficients.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

execution_plan_max_bytes <- function() {
  value <- getOption("markov.misc.execution_plan_max_bytes", 256 * 1024^2)
  if (
    length(value) != 1L ||
      !is.numeric(value) ||
      is.na(value) ||
      value <= 0
  ) {
    stop(
      "Option `markov.misc.execution_plan_max_bytes` must be a positive number.",
      call. = FALSE
    )
  }
  as.double(value)
}

stop_execution_plan_too_large <- function(required_bytes, max_bytes) {
  message <- paste0(
    "Compiled SOP designs require ",
    format(required_bytes, scientific = FALSE),
    " bytes, above the execution-plan limit of ",
    format(max_bytes, scientific = FALSE),
    " bytes."
  )
  stop(structure(
    list(
      message = message,
      call = NULL,
      required_bytes = required_bytes,
      max_bytes = max_bytes
    ),
    class = c(
      "markov_misc_execution_plan_too_large",
      "error",
      "condition"
    )
  ))
}

preflight_execution_plan_designs <- function(
  n_pat,
  n_times,
  rows_per_transition,
  n_cols
) {
  initial_cells <- as.double(n_pat) * as.double(n_cols)
  transition_cells <-
    as.double(max(n_times - 1L, 0L)) *
    as.double(rows_per_transition) *
    as.double(n_cols)
  design_bytes <- (initial_cells + transition_cells) * 8
  workspace_bytes <- execution_plan_max_bytes()
  if (!is.finite(design_bytes) || design_bytes > workspace_bytes) {
    stop_execution_plan_too_large(design_bytes, workspace_bytes)
  }
  list(
    design_bytes = design_bytes,
    workspace_bytes = workspace_bytes
  )
}

measure_execution_plan_designs <- function(X_init, X_transition) {
  matrices <- c(list(X_init), Filter(Negate(is.null), X_transition))
  sum(vapply(
    matrices,
    function(x) as.double(length(x)) * 8,
    numeric(1)
  ))
}

validate_execution_plan_designs <- function(design_bytes, workspace_bytes) {
  if (!is.finite(design_bytes) || design_bytes > workspace_bytes) {
    stop_execution_plan_too_large(design_bytes, workspace_bytes)
  }
  invisible(NULL)
}

compile_sop_execution_plan <- function(
  model,
  newdata,
  times,
  y_levels,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  builder = c("batched", "streamed"),
  output = "array"
) {
  builder <- match.arg(builder)
  validate_markov_model(model)
  Gamma_template <- get_effective_coefs(model)
  validate_sop_threshold_count(Gamma_template, y_levels)
  if (!is.null(p2_var)) {
    return(compile_second_order_execution_plan(
      model = model,
      newdata = newdata,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      output = output,
      Gamma = Gamma_template
    ))
  }
  resolved_times <- resolve_sop_times(
    model,
    newdata,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "fast"
  )$times
  build <- if (identical(builder, "batched")) {
    markov_msm_build_batched
  } else {
    markov_msm_build
  }
  components <- build(
    model = model,
    newdata = newdata,
    time_covariates = time_covariates,
    times = resolved_times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    gap_var = gap_var
  )
  structure(
    list(
      version = 1L,
      model_class = class(model)[1L],
      link = "logit",
      recursion_order = 1L,
      times = resolved_times,
      y_levels = as_state_labels(y_levels),
      absorb = as_state_labels(absorb),
      time_var = time_var,
      p_var = p_var,
      p2_var = NULL,
      gap_var = gap_var,
      output = output,
      workspace_bytes = components$workspace_bytes,
      design_bytes = components$design_bytes,
      basis_terms = attr(model, "markov_basis_terms", exact = TRUE),
      components = components
    ),
    class = "markov_sop_exec_plan"
  )
}

run_sop_execution_plan <- function(plan, Gamma) {
  if (!inherits(plan, "markov_sop_exec_plan") || plan$version != 1L) {
    stop("Invalid or unsupported SOP execution plan.")
  }
  validate_sop_threshold_count(Gamma, plan$y_levels)
  if (identical(plan$recursion_order, 2L)) {
    return(run_second_order_execution_plan(plan, Gamma))
  }
  markov_msm_run(
    plan$components,
    Gamma,
    times = plan$times,
    absorb = plan$absorb
  )
}

compile_second_order_execution_plan <- function(
  model,
  newdata,
  times,
  y_levels,
  absorb,
  time_var,
  p_var,
  p2_var,
  gap_var,
  time_covariates,
  output,
  Gamma
) {
  if (!inherits(model, c("orm", "vglm")) || inherits(model, "blrm")) {
    stop("Compiled second-order plans require a frequentist orm or vglm model.")
  }
  if (!all(c(p_var, p2_var) %in% names(newdata))) {
    stop("Second-order prediction variables are missing from `newdata`.")
  }
  time_res <- resolve_sop_times(
    model,
    newdata,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "fast"
  )
  times <- time_res$times
  time_info <- time_res$time_info
  validate_factor_gap(gap_var, time_covariates, time_info)
  ylevel_names <- as_state_labels(y_levels)
  absorb_idx <- which(ylevel_names %in% as_state_labels(absorb))
  non_absorb_idx <- setdiff(seq_along(ylevel_names), absorb_idx)
  pair_grid <- expand.grid(
    older = non_absorb_idx,
    current = non_absorb_idx,
    KEEP.OUT.ATTRS = FALSE
  )
  validate_sop_threshold_count(Gamma, y_levels)
  budget <- preflight_execution_plan_designs(
    n_pat = nrow(newdata),
    n_times = length(times),
    rows_per_transition = nrow(newdata) * nrow(pair_grid),
    n_cols = ncol(Gamma)
  )
  terms <- if (inherits(model, "vglm")) {
    stats::delete.response(stats::terms(model))
  } else {
    NULL
  }
  get_X <- function(data) {
    X <- if (inherits(model, "orm")) {
      orm_model_matrix(model, data, include_intercept = TRUE)
    } else {
      stats::model.matrix(
        terms,
        data = data,
        contrasts.arg = if (length(model@contrasts)) model@contrasts else NULL,
        xlev = model@xlevels
      )
    }
    missing <- setdiff(colnames(Gamma), colnames(X))
    if (length(missing) > 0L) {
      stop(
        "Second-order design is missing fitted columns: ",
        paste(missing, collapse = ", ")
      )
    }
    X[, colnames(Gamma), drop = FALSE]
  }
  set_visit <- function(data, index) {
    assign_sop_visit(
      data,
      time_var = time_var,
      times = times,
      index = index,
      time_covariates = time_covariates,
      gap_var = gap_var,
      time_info = time_info
    )
  }
  initial_data <- set_visit(newdata, 1L)
  initial_data[[p_var]] <- normalize_previous_state_column(
    initial_data[[p_var]],
    newdata[[p_var]],
    p_var
  )
  initial_data[[p2_var]] <- normalize_previous_state_column(
    initial_data[[p2_var]],
    newdata[[p2_var]],
    p2_var
  )
  X_init <- get_X(initial_data)
  designs <- vector("list", length(times))
  if (length(times) >= 2L) {
    for (visit in 2:length(times)) {
      data_visit <- set_visit(newdata, visit)
      expanded <- data_visit[
        rep(seq_len(nrow(newdata)), times = nrow(pair_grid)),
        ,
        drop = FALSE
      ]
      expanded[[p2_var]] <- make_previous_state_column(
        ylevel_names[pair_grid$older],
        newdata[[p2_var]],
        nrow(newdata),
        p2_var
      )
      expanded[[p_var]] <- make_previous_state_column(
        ylevel_names[pair_grid$current],
        newdata[[p_var]],
        nrow(newdata),
        p_var
      )
      designs[[visit]] <- get_X(expanded)
    }
  }
  design_bytes <- measure_execution_plan_designs(X_init, designs)
  validate_execution_plan_designs(design_bytes, budget$workspace_bytes)
  structure(
    list(
      version = 1L,
      model_class = class(model)[1L],
      link = "logit",
      recursion_order = 2L,
      times = times,
      y_levels = ylevel_names,
      absorb = as_state_labels(absorb),
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      output = output,
      workspace_bytes = budget$workspace_bytes,
      design_bytes = design_bytes,
      components = list(
        X_init = X_init,
        X_transition = designs,
        previous = match_state_indices(newdata[[p_var]], ylevel_names, p_var),
        older = pair_grid$older,
        current = pair_grid$current,
        absorb = absorb_idx,
        non_absorb = non_absorb_idx,
        n_pat = nrow(newdata),
        n_states = length(ylevel_names),
        col_names = colnames(X_init)
      )
    ),
    class = "markov_sop_exec_plan"
  )
}

run_second_order_execution_plan <- function(plan, Gamma) {
  components <- plan$components
  Gamma <- Gamma[, components$col_names, drop = FALSE]
  po <- markov_po_structure(Gamma, components$col_names, components$X_init)
  probabilities <- function(X) {
    if (is.null(po)) {
      return(lp_to_probs(X %*% t(Gamma), nrow(Gamma)))
    }
    lp_to_probs(
      outer(drop(X %*% po$beta), po$cutpoints, "+"),
      length(po$cutpoints)
    )
  }
  initial <- probabilities(components$X_init)
  check_transition_probabilities(initial, "the first SOP time point")
  output <- array(
    0,
    dim = c(components$n_pat, length(plan$times), components$n_states)
  )
  output[, 1L, ] <- initial
  joint <- array(
    0,
    dim = c(components$n_pat, components$n_states, components$n_states)
  )
  for (patient in seq_len(components$n_pat)) {
    previous <- components$previous[patient]
    if (previous %in% components$absorb) {
      joint[patient, previous, previous] <- 1
      output[patient, 1L, ] <- 0
      output[patient, 1L, previous] <- 1
    } else {
      joint[patient, previous, ] <- initial[patient, ]
    }
  }
  if (length(plan$times) >= 2L) {
    for (visit in 2:length(plan$times)) {
      if (is.null(po)) {
        transition <- probabilities(components$X_transition[[visit]])
        joint <- markov_update_second_order_native(
          joint,
          transition,
          components$older,
          components$current,
          components$absorb
        )
      } else {
        scalar <- drop(components$X_transition[[visit]] %*% po$beta)
        joint <- markov_update_second_order_po_native(
          joint,
          scalar,
          po$cutpoints,
          components$older,
          components$current,
          components$absorb
        )
      }
      output[, visit, ] <- apply(joint, c(1L, 3L), sum)
    }
  }
  dimnames(output) <- list(NULL, plan$times, plan$y_levels)
  output
}
