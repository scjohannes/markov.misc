# Internal fitted-model trajectory simulation helpers.

markov_supported_model <- function(object) {
  inherits(object, c("orm", "blrm", "vglm", "vgam", "robcov_vglm"))
}

markov_model_backend <- function(object) {
  if (inherits(object, "blrm")) {
    return("blrm")
  }
  if (inherits(object, "robcov_vglm")) {
    return("robcov_vglm")
  }
  if (inherits(object, "orm")) {
    return("orm")
  }
  if (inherits(object, "vglm")) {
    return("vglm")
  }
  if (inherits(object, "vgam")) {
    return("vgam")
  }
  stop("Object class not supported.")
}

markov_model_ylevels <- function(model, ylevels = NULL) {
  if (!is.null(ylevels)) {
    return(ylevels)
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

  stop("`ylevels` cannot be NULL.")
}

markov_prediction_function <- function(
  model,
  include_re,
  id_var,
  n_draws
) {
  backend <- markov_model_backend(model)
  draw_indices <- NULL
  gamma_draws <- NULL

  if (identical(backend, "blrm")) {
    draw_indices <- select_posterior_draws(
      model,
      n_draws = n_draws,
      seed = NULL
    )
    gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  }

  switch(
    backend,
    orm = function(data) predict_orm_response_markov(model, data),
    vglm = function(data) predict_vglm_response_markov(model, data),
    vgam = function(data) predict_vglm_response_markov(model, data),
    robcov_vglm = function(data) {
      if (is.null(model$vglm_fit)) {
        stop(
          "robcov_vglm object does not contain the original vglm fit. ",
          "Please re-run robcov_vglm() with the latest version of the package."
        )
      }
      predict_vglm_response_markov(model$vglm_fit, data)
    },
    blrm = function(data) {
      predict_blrm_response_markov(
        model,
        data,
        include_re = include_re,
        id_var = id_var,
        draw_indices = draw_indices,
        gamma_draws = gamma_draws
      )
    }
  )
}

markov_prediction_matrix <- function(p, n_states, context) {
  check_transition_probabilities(p, context)

  if (length(dim(p)) == 3L) {
    p <- apply(p, c(2, 3), mean)
  }
  p <- as.matrix(p)
  if (ncol(p) != n_states) {
    stop(
      "Model prediction returned ",
      ncol(p),
      " state probabilities, but `ylevels` has length ",
      n_states,
      "."
    )
  }

  p[p < 0 & p > -sqrt(.Machine$double.eps)] <- 0
  rowsum <- rowSums(p)
  bad <- !is.finite(rowsum) | rowsum <= 0
  if (any(bad)) {
    stop(
      "Model prediction returned non-positive transition-probability sums ",
      "during ",
      context,
      "."
    )
  }
  p / rowsum
}

markov_sample_states <- function(probabilities, ylevel_names) {
  idx <- vapply(
    seq_len(nrow(probabilities)),
    function(i) {
      sample.int(length(ylevel_names), size = 1L, prob = probabilities[i, ])
    },
    integer(1)
  )
  ylevel_names[idx]
}

markov_previous_state_prototype <- function(x, ylevel_names) {
  if (!is.factor(x)) {
    return(x)
  }

  factor(
    as.character(x),
    levels = unique(c(levels(x), ylevel_names)),
    ordered = is.ordered(x)
  )
}

validate_markov_n_rep <- function(n_rep) {
  if (
    !is.numeric(n_rep) ||
      length(n_rep) != 1L ||
      is.na(n_rep) ||
      !is.finite(n_rep) ||
      n_rep < 1
  ) {
    stop("`n_rep` must be a positive integer.")
  }
  as.integer(n_rep)
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

markov_simulation_seed_scope <- function(seed) {
  if (is.null(seed)) {
    return(function() NULL)
  }

  old_seed_exists <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )
  old_seed <- if (old_seed_exists) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  set.seed(seed)

  function() {
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    NULL
  }
}

complete_plot_simulation_times <- function(times, time_info) {
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

markov_simulate_fitted_paths <- function(
  model,
  newdata = NULL,
  refit_data = NULL,
  variables = NULL,
  times,
  ylevels = NULL,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  p2varname = NULL,
  id_var = NULL,
  gap = NULL,
  t_covs = NULL,
  include_re = FALSE,
  n_rep = 20L,
  n_draws = 100L,
  seed = NULL
) {
  validate_markov_model(model)
  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied for model-based plots.")
  }
  n_rep <- validate_markov_n_rep(n_rep)

  restore_seed <- markov_simulation_seed_scope(seed)
  on.exit(restore_seed(), add = TRUE)

  data_res <- resolve_markov_source_data(model, newdata, refit_data)
  source_data <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied

  if (inherits(model, "blrm") && isTRUE(include_re)) {
    id_var <- resolve_blrm_id_var(model, source_data, id_var)
    validate_markov_id_var(id_var, source_data, "newdata")
  } else {
    id_var <- markov_model_id_var(model, id_var) %||% "id"
  }

  if (!is.null(refit_data)) {
    validate_markov_id_var(id_var, refit_data, "refit_data")
  }
  if (!newdata_supplied) {
    validate_markov_id_var(id_var, source_data, "stored model data")
  }

  baseline_data <- if (newdata_supplied) {
    source_data
  } else {
    resolve_markov_prediction_data(
      source_data,
      id_var = id_var,
      tvarname = tvarname
    )
  }

  if (!pvarname %in% names(baseline_data)) {
    stop("Previous-state variable `", pvarname, "` not found in `newdata`.")
  }
  if (!is.null(p2varname) && !p2varname %in% names(baseline_data)) {
    stop(
      "Second previous-state variable `",
      p2varname,
      "` not found in `newdata`."
    )
  }

  variables <- validate_plot_counterfactual_variables(variables)
  if (!is.null(variables)) {
    missing_vars <- setdiff(names(variables), names(baseline_data))
    if (length(missing_vars) > 0) {
      stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
    }
    grid <- do.call(
      expand.grid,
      c(variables, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    )
    baseline_data <- create_counterfactual_data(baseline_data, grid, variables)
  }

  baseline_data <- ensure_markov_rowid(baseline_data)
  sim_data <- baseline_data[
    rep(seq_len(nrow(baseline_data)), each = n_rep),
    ,
    drop = FALSE
  ]
  rownames(sim_data) <- NULL
  sim_data$.rep <- rep(seq_len(n_rep), times = nrow(baseline_data))
  sim_data$.sim_id <- seq_len(nrow(sim_data))

  time_res <- resolve_sop_times(
    model,
    baseline_data,
    times,
    tvarname,
    t_covs = NULL,
    default = "unique"
  )
  plot_times <- time_res$times
  time_info <- time_res$time_info
  times <- complete_plot_simulation_times(plot_times, time_info)
  if (!is.null(t_covs) && nrow(t_covs) != length(times)) {
    stop(
      "`t_covs` must have one row per simulated time point for model-based ",
      "plots. Sparse plot times are expanded to the full simulation grid; ",
      "expected ",
      length(times),
      " rows but got ",
      nrow(t_covs),
      "."
    )
  }
  validate_factor_gap(gap, t_covs, time_info)

  ylevels <- markov_model_ylevels(model, ylevels)
  ylevel_names <- as_state_labels(ylevels)
  absorb_names <- as_state_labels(absorb)
  n_states <- length(ylevel_names)

  sim_data[[pvarname]] <- markov_previous_state_prototype(
    sim_data[[pvarname]],
    ylevel_names
  )
  if (!is.null(p2varname)) {
    sim_data[[p2varname]] <- markov_previous_state_prototype(
      sim_data[[p2varname]],
      ylevel_names
    )
  }

  prd <- markov_prediction_function(
    model,
    include_re = include_re,
    id_var = id_var,
    n_draws = n_draws
  )

  prev <- normalize_previous_state_column(
    sim_data[[pvarname]],
    sim_data[[pvarname]],
    pvarname
  )
  prev2 <- if (!is.null(p2varname)) {
    normalize_previous_state_column(
      sim_data[[p2varname]],
      sim_data[[p2varname]],
      p2varname
    )
  } else {
    NULL
  }

  output <- vector("list", length(times))
  for (i in seq_along(times)) {
    step_data <- sim_data
    step_data[[pvarname]] <- normalize_previous_state_column(
      prev,
      sim_data[[pvarname]],
      pvarname
    )
    if (!is.null(p2varname)) {
      step_data[[p2varname]] <- normalize_previous_state_column(
        prev2,
        sim_data[[p2varname]],
        p2varname
      )
    }
    step_data <- assign_sop_visit(
      step_data,
      tvarname = tvarname,
      times = times,
      index = i,
      t_covs = t_covs,
      gap = gap,
      time_info = time_info
    )

    current <- as.character(prev)
    predict_rows <- !as.character(prev) %in% absorb_names
    if (any(predict_rows)) {
      probs <- prd(step_data[predict_rows, , drop = FALSE])
      probs <- markov_prediction_matrix(
        probs,
        n_states = n_states,
        context = paste0("simulated time point ", i)
      )
      current[predict_rows] <- markov_sample_states(probs, ylevel_names)
    }

    out <- step_data
    out$y <- factor(current, levels = ylevel_names)
    output[[i]] <- out

    if (!is.null(p2varname)) {
      prev2 <- prev
    }
    prev <- normalize_previous_state_column(
      current,
      sim_data[[pvarname]],
      pvarname
    )
  }

  plot_indices <- match(as.character(plot_times), as.character(times))
  output <- output[plot_indices]
  out <- do.call(rbind, output)
  rownames(out) <- NULL
  attr(out, "ylevels") <- ylevel_names
  attr(out, "sim_times") <- times
  attr(out, "plot_times") <- plot_times
  out
}
