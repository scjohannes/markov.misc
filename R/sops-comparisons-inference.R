# Average comparison inference helpers.

inferences_avg_comparisons <- function(
  object,
  method,
  engine,
  score_weight_dist,
  n_sim,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  update_datadist,
  use_coefstart
) {
  args <- attr(object, "comparison_args")
  metric <- args$metric

  if (metric %in% c("sop", "time_in_state")) {
    return(inferences_avg_comparisons_linear(
      object = object,
      method = method,
      engine = engine,
      score_weight_dist = score_weight_dist,
      n_sim = n_sim,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    ))
  }

  if (method == "simulation") {
    return(inferences_avg_comparisons_time_benefit_simulation(
      object = object,
      engine = engine,
      score_weight_dist = score_weight_dist,
      n_sim = n_sim,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws
    ))
  }

  bootstrap_engine <- if (engine == "mvn") "standard" else engine
  inferences_avg_comparisons_time_benefit_bootstrap(
    object = object,
    engine = bootstrap_engine,
    n_boot = n_sim,
    workers = workers,
    conf_level = conf_level,
    return_draws = return_draws,
    update_datadist = update_datadist,
    use_coefstart = use_coefstart
  )
}

inferences_avg_comparisons_linear <- function(
  object,
  method,
  engine,
  score_weight_dist,
  n_sim,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws,
  update_datadist,
  use_coefstart
) {
  args <- attr(object, "comparison_args")
  avg_args <- attr(object, "avg_args")
  newdata <- if (isTRUE(attr(object, "newdata_supplied"))) {
    attr(object, "newdata_orig")
  } else {
    NULL
  }
  avg <- avg_comparison_replay_avg_sops(
    model = attr(object, "model"),
    newdata = newdata,
    refit_data = attr(object, "refit_data"),
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    ylevels = attr(object, "ylevels"),
    absorb = attr(object, "absorb"),
    tvarname = attr(object, "tvarname") %||% "time",
    pvarname = attr(object, "pvarname") %||% "yprev",
    id_var = avg_args$id_var,
    p2varname = attr(object, "p2varname"),
    gap = attr(object, "gap"),
    t_covs = attr(object, "t_covs"),
    include_re = args$include_re,
    n_draws = args$n_draws,
    seed = args$seed,
    posterior_summary = args$posterior_summary,
    conf_level = conf_level,
    return_draws = FALSE,
    extra_args = args$extra_args
  )
  avg <- inferences(
    avg,
    method = method,
    engine = engine,
    score_weight_dist = score_weight_dist,
    n_sim = n_sim,
    vcov = vcov,
    cluster = cluster,
    workers = workers,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = TRUE,
    update_datadist = update_datadist,
    use_coefstart = use_coefstart
  )

  state_sets <- normalize_comparison_state_sets(
    args$states,
    attr(avg, "ylevels")
  )
  result <- avg_comparison_from_avg_sops(
    avg,
    metric = args$metric,
    state_sets = state_sets,
    comparison = args$comparison,
    time_map = args$time_map,
    origin_time = args$origin_time,
    xout = args$xout,
    origin = args$origin,
    time_unit = args$time_unit,
    return_draws = return_draws
  )
  result <- restore_avg_comparison_inference_attrs(result, object, avg)
  if (!isTRUE(return_draws)) {
    attr(result, "simulation_draws") <- NULL
    attr(result, "bootstrap_draws") <- NULL
  }
  result
}

restore_avg_comparison_inference_attrs <- function(result, object, inferred) {
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  for (a in c(
    "n_sim",
    "n_boot",
    "n_successful",
    "conf_level",
    "conf_type",
    "method",
    "engine",
    "score_weight_dist",
    "fwb_weight_type",
    "fwb_weight_scale",
    "draw_weights_attached",
    "draw_weight_col",
    "draw_weight_omission_reason"
  )) {
    if (!is.null(attr(inferred, a))) {
      attr(result, a) <- attr(inferred, a)
    }
  }
  class(result) <- class(object)
  result
}

inferences_avg_comparisons_time_benefit_simulation <- function(
  object,
  engine,
  score_weight_dist,
  n_sim,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  return_draws
) {
  if (engine %in% c("standard", "fwb")) {
    stop(
      "`engine = \"",
      engine,
      "\"` is only used when `method = \"bootstrap\"`."
    )
  }

  setup <- setup_from_avg_comparison_object(object)
  args <- attr(object, "comparison_args")
  model <- attr(object, "model")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))

  simulation_draws <- generate_sop_coefficient_draws(
    model = model,
    engine = engine,
    score_weight_dist = score_weight_dist,
    n_sim = n_sim,
    vcov = vcov,
    cluster = cluster,
    baseline_data = setup$baseline_data,
    id_var = setup$id_var,
    newdata_supplied = newdata_supplied,
    is_avg = TRUE,
    by = setup$by,
    return_draws = return_draws,
    prediction_data = setup$newdata_pred,
    attach_individual_weights = FALSE
  )
  beta_draws <- simulation_draws$beta_draws
  baseline_weights_draws <- simulation_draws$baseline_weights_draws

  analysis_fn <- function(i) {
    model_i <- set_coef(model, beta_draws[i, ])
    sops_array <- tryCatch(
      soprob_markov(
        object = model_i,
        data = setup$newdata_pred,
        times = setup$times,
        ylevels = setup$ylevels,
        absorb = setup$absorb,
        tvarname = setup$tvarname,
        pvarname = setup$pvarname,
        p2varname = setup$p2varname,
        gap = setup$gap,
        t_covs = setup$t_covs
      ),
      error = function(e) {
        warning("soprob_markov failed in draw ", i, ": ", e$message)
        NULL
      }
    )
    if (is.null(sops_array)) {
      return(NULL)
    }
    baseline_weights <- if (!is.null(baseline_weights_draws)) {
      baseline_weights_draws[i, ]
    } else {
      NULL
    }
    draw_i <- reduce_time_benefit_array_for_setup(
      sops_array = sops_array,
      setup = setup,
      comparison = args$comparison,
      weights = baseline_weights,
      time_unit = args$time_unit,
      time_map = args$time_map,
      origin_time = args$origin_time,
      xout = args$xout,
      origin = args$origin
    )
    draw_i$draw_id <- i
    draw_i
  }

  draw_results <- apply_sop_simulation_draws(
    n_draws = nrow(beta_draws),
    analysis_fn = analysis_fn,
    workers = workers,
    globals = c(
      "model",
      "beta_draws",
      "setup",
      "args",
      "baseline_weights_draws"
    )
  )

  draw_results <- Filter(Negate(is.null), draw_results)
  if (length(draw_results) == 0L) {
    stop("All simulation draws failed.")
  }
  draws_df <- bind_rows_fill(draw_results)
  finalize_avg_comparison_draw_inference(
    object = object,
    draws_df = draws_df,
    draw_attr = "simulation_draws",
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = return_draws,
    metadata = list(
      n_sim = n_sim,
      n_successful = length(draw_results),
      conf_level = conf_level,
      conf_type = conf_type,
      method = "simulation",
      engine = engine,
      score_weight_dist = if (engine == "score_bootstrap") {
        score_weight_dist
      } else {
        NULL
      },
      draw_weight_omission_reason = simulation_draws$draw_weight_omission_reason
    )
  )
}

inferences_avg_comparisons_time_benefit_bootstrap <- function(
  object,
  engine,
  n_boot,
  workers,
  conf_level,
  return_draws,
  update_datadist,
  use_coefstart
) {
  engine <- match.arg(engine, choices = c("standard", "fwb"))
  setup <- setup_from_avg_comparison_object(object)
  args <- attr(object, "comparison_args")
  model <- attr(object, "model")
  refit_data <- attr(object, "refit_data") %||% attr(object, "newdata_orig")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))

  if (is.null(refit_data)) {
    stop("Full refit data not stored. Cannot perform bootstrap.")
  }
  validate_refit_bootstrap_data(refit_data, setup$id_var, setup$tvarname)
  if (newdata_supplied && engine == "fwb") {
    warn_fixed_profile_bootstrap_weights("FWB")
  }

  factor_cols <- c("y", setup$pvarname, setup$p2varname)
  factor_cols <- intersect(factor_cols, names(refit_data))

  analysis_fn <- function(boot_data, fwb_weights = NULL) {
    fit_weights <- if (engine == "fwb") boot_data$fwb_weight else NULL
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = refit_data,
      ylevels = setup$ylevels,
      absorb = setup$absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart,
      fit_weights = fit_weights
    )
    m_boot <- boot_result$model
    if (is.null(m_boot)) {
      return(NULL)
    }

    if (newdata_supplied) {
      newdata_cf <- setup$newdata_pred
      baseline_data <- setup$baseline_data
      baseline_weights <- NULL
    } else if (engine == "standard") {
      baseline_data <- resolve_markov_prediction_data(
        boot_result$data,
        id_var = "new_id",
        tvarname = setup$tvarname,
        data_label = "bootstrap data"
      )
      baseline_data <- ensure_markov_rowid(baseline_data)
      newdata_cf <- create_counterfactual_data(
        baseline_data,
        setup$grid,
        setup$variables
      )
      baseline_weights <- NULL
    } else {
      baseline_data <- resolve_markov_prediction_data(
        boot_result$data,
        id_var = setup$id_var,
        tvarname = setup$tvarname,
        data_label = "bootstrap data"
      )
      baseline_data <- ensure_markov_rowid(baseline_data)
      newdata_cf <- create_counterfactual_data(
        baseline_data,
        setup$grid,
        setup$variables
      )
      baseline_weights <- fwb_baseline_weights(
        fwb_weights = fwb_weights,
        baseline_data = baseline_data,
        id_var = setup$id_var
      )
    }

    sops_array <- tryCatch(
      soprob_markov(
        object = m_boot,
        data = newdata_cf,
        times = setup$times,
        ylevels = factor(boot_result$ylevels),
        absorb = boot_result$absorb,
        tvarname = setup$tvarname,
        pvarname = setup$pvarname,
        p2varname = setup$p2varname,
        gap = setup$gap,
        t_covs = setup$t_covs
      ),
      error = function(e) {
        warning("soprob_markov failed: ", e$message)
        NULL
      }
    )
    if (is.null(sops_array)) {
      return(NULL)
    }

    sops_array <- complete_comparison_sops_array_states(
      sops_array = sops_array,
      boot_ylevels = boot_result$ylevels,
      ylevels = setup$ylevels
    )
    boot_setup <- setup
    boot_setup$baseline_data <- baseline_data
    boot_setup$newdata_pred <- newdata_cf
    boot_setup$n_each <- nrow(baseline_data)
    reduce_time_benefit_array_for_setup(
      sops_array = sops_array,
      setup = boot_setup,
      comparison = args$comparison,
      weights = baseline_weights,
      time_unit = args$time_unit,
      time_map = args$time_map,
      origin_time = args$origin_time,
      xout = args$xout,
      origin = args$origin
    )
  }

  if (engine == "standard") {
    boot_ids <- fast_group_bootstrap(
      data = refit_data,
      id_var = setup$id_var,
      n_boot = n_boot
    )
    boot_results <- apply_to_bootstrap(
      boot_samples = boot_ids,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = setup$id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats"),
      globals = c(
        "model",
        "setup",
        "args",
        "refit_data",
        "newdata_supplied",
        "factor_cols",
        "engine",
        "update_datadist",
        "use_coefstart",
        "complete_comparison_sops_array_states",
        "reduce_time_benefit_array_for_setup",
        "time_benefit_sops_from_array"
      )
    )
  } else {
    fwb_samples <- generate_fwb_bootstrap_weights(
      data = refit_data,
      id_var = setup$id_var,
      n_boot = n_boot
    )
    boot_results <- apply_to_fwb_bootstrap(
      fwb_samples = fwb_samples,
      analysis_fn = analysis_fn,
      data = refit_data,
      id_var = setup$id_var,
      workers = workers,
      packages = c("rms", "VGAM", "Hmisc", "stats"),
      globals = c(
        "model",
        "setup",
        "args",
        "refit_data",
        "newdata_supplied",
        "factor_cols",
        "engine",
        "update_datadist",
        "use_coefstart",
        "complete_comparison_sops_array_states",
        "reduce_time_benefit_array_for_setup",
        "time_benefit_sops_from_array"
      )
    )
  }

  boot_results <- Filter(Negate(is.null), boot_results)
  if (length(boot_results) == 0L) {
    stop("All bootstrap iterations failed.")
  }
  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- bind_rows_fill(boot_results)

  finalize_avg_comparison_draw_inference(
    object = object,
    draws_df = boot_df,
    draw_attr = "bootstrap_draws",
    conf_level = conf_level,
    conf_type = "perc",
    return_draws = return_draws,
    metadata = list(
      n_boot = n_boot,
      n_successful = length(boot_results),
      conf_level = conf_level,
      method = "bootstrap",
      engine = engine,
      fwb_weight_type = if (engine == "fwb") "exponential" else NULL,
      fwb_weight_scale = if (engine == "fwb") "cluster_mean_1" else NULL
    )
  )
}

setup_from_avg_comparison_object <- function(object) {
  avg_args <- attr(object, "avg_args")
  call_args <- attr(object, "call_args")
  list(
    variables = avg_args$variables,
    by = avg_args$by,
    times = avg_args$times,
    id_var = avg_args$id_var,
    p2varname = attr(object, "p2varname"),
    newdata_orig = attr(object, "newdata_orig"),
    refit_data = attr(object, "refit_data"),
    newdata_supplied = isTRUE(attr(object, "newdata_supplied")),
    baseline_data = attr(object, "comparison_baseline_data"),
    grid = attr(object, "comparison_grid"),
    newdata_pred = attr(object, "newdata_pred"),
    n_cf = nrow(attr(object, "comparison_grid")),
    n_each = attr(object, "comparison_n_each"),
    call_args = call_args,
    tvarname = attr(object, "tvarname") %||% call_args$tvarname %||% "time",
    pvarname = attr(object, "pvarname") %||% call_args$pvarname %||% "yprev",
    ylevels = attr(object, "ylevels") %||% call_args$ylevels,
    absorb = attr(object, "absorb") %||% call_args$absorb,
    gap = attr(object, "gap") %||% call_args$gap,
    t_covs = attr(object, "t_covs") %||% call_args$t_covs
  )
}

reduce_time_benefit_array_for_setup <- function(
  sops_array,
  setup,
  comparison,
  weights,
  time_unit,
  time_map = NULL,
  origin_time = NULL,
  xout = NULL,
  origin = "empirical_baseline"
) {
  if (!is.null(time_map) || !is.null(origin_time)) {
    if (is.null(time_map)) {
      stop("`time_map` must be supplied for real-time AUC.")
    }
    ind <- time_benefit_sops_from_array(sops_array, setup)
    ind <- interpolate_sops(
      ind,
      time_map = time_map,
      xout = xout,
      origin_time = origin_time,
      origin = origin
    )
    real_setup <- setup
    real_setup$times <- sort(unique(ind$time))
    real_setup$ylevels <- attr(ind, "ylevels")
    out <- reduce_time_benefit_sops_df(
      data = as.data.frame(ind),
      setup = real_setup,
      comparison = comparison,
      value_col = "estimate",
      real_time = TRUE,
      weights = weights,
      time_unit = time_unit
    )
    return(out)
  }

  arr <- array(
    NA_real_,
    dim = c(setup$n_cf, setup$n_each, dim(sops_array)[2], dim(sops_array)[3])
  )
  for (cf_i in seq_len(setup$n_cf)) {
    rows <- ((cf_i - 1L) * setup$n_each + 1L):(cf_i * setup$n_each)
    arr[cf_i, , , ] <- sops_array[rows, , , drop = FALSE]
  }
  out <- time_benefit_from_counterfactual_array(
    sops_array = arr,
    setup = setup,
    comparison = comparison,
    weights = weights,
    real_time = FALSE
  )
  if (!is.null(time_unit)) {
    out$time_unit <- time_unit
  }
  out
}

time_benefit_sops_from_array <- function(sops_array, setup) {
  newdata <- setup$newdata_pred
  if (is.null(newdata)) {
    stop(
      "`setup$newdata_pred` is required for real-time time-benefit inference."
    )
  }
  newdata$rowid <- seq_len(nrow(newdata))

  out <- array_to_df_individual(
    sops_array = sops_array,
    times = setup$times,
    ylevels = setup$ylevels,
    newdata = newdata,
    by = NULL
  )
  attr(out, "ylevels") <- setup$ylevels
  attr(out, "call_args") <- setup$call_args %||%
    list(
      times = setup$times,
      ylevels = setup$ylevels,
      absorb = setup$absorb,
      tvarname = setup$tvarname,
      pvarname = setup$pvarname,
      p2varname = setup$p2varname,
      gap = setup$gap,
      t_covs = setup$t_covs
    )
  attr(out, "newdata_orig") <- setup$newdata_orig %||% setup$baseline_data
  attr(out, "newdata_pred") <- newdata
  attr(out, "newdata_supplied") <- setup$newdata_supplied
  attr(out, "id_var") <- setup$id_var
  attr(out, "tvarname") <- setup$tvarname
  attr(out, "pvarname") <- setup$pvarname
  attr(out, "p2varname") <- setup$p2varname
  attr(out, "absorb") <- setup$absorb
  attr(out, "gap") <- setup$gap
  attr(out, "t_covs") <- setup$t_covs
  class(out) <- c("markov_sops", class(out))
  out
}

complete_comparison_sops_array_states <- function(
  sops_array,
  boot_ylevels,
  ylevels
) {
  boot_labels <- as_state_labels(boot_ylevels)
  target_labels <- as_state_labels(ylevels)
  if (identical(boot_labels, target_labels)) {
    return(sops_array)
  }

  full <- array(
    0,
    dim = c(dim(sops_array)[1], dim(sops_array)[2], length(target_labels))
  )
  for (i in seq_along(boot_labels)) {
    target <- match(boot_labels[i], target_labels)
    if (!is.na(target)) {
      full[,, target] <- sops_array[,, i]
    }
  }
  full
}

finalize_avg_comparison_draw_inference <- function(
  object,
  draws_df,
  draw_attr,
  conf_level,
  conf_type,
  return_draws,
  metadata
) {
  group_cols <- setdiff(names(draws_df), c("estimate", "draw_id"))
  ci <- compute_ci_from_draws(
    draws_df = draws_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )
  result <- merge_comparison_ci(as.data.frame(object), ci, group_cols)
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(object, a)
    }
  }
  for (nm in names(metadata)) {
    if (!is.null(metadata[[nm]])) {
      attr(result, nm) <- metadata[[nm]]
    }
  }
  if (isTRUE(return_draws)) {
    attr(result, draw_attr) <- draws_df
  }
  class(result) <- class(object)
  result
}
