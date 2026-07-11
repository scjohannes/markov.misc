# Shared coefficient-draw helpers for simulation inference.

generate_sop_coefficient_draws <- function(
  model,
  engine,
  score_weight_dist,
  n_draws,
  vcov,
  cluster,
  baseline_data,
  id_var,
  newdata_supplied,
  is_avg,
  by,
  return_draws,
  prediction_data = baseline_data,
  attach_individual_weights = TRUE
) {
  out <- list(
    beta_draws = NULL,
    baseline_weights_draws = NULL,
    draw_weight_col = NULL,
    draw_weights_attached = FALSE,
    draw_weight_omission_reason = NULL
  )

  if (engine == "mvn") {
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      stop(
        "Package 'mvtnorm' is required for MVN simulation inference.\n",
        "Install with: install.packages('mvtnorm')"
      )
    }

    beta_hat <- get_coef(model)
    Sigma <- if (!is.null(vcov)) {
      validate_coef_vcov(beta_hat, vcov, arg = "vcov")
    } else {
      validate_coef_vcov(beta_hat, get_vcov_robust(model), arg = "model vcov")
    }
    out$beta_draws <- mvtnorm::rmvnorm(n_draws, mean = beta_hat, sigma = Sigma)
    return(out)
  }

  if (engine != "score_bootstrap") {
    stop("Unknown simulation engine: ", engine)
  }
  if (!is.null(vcov)) {
    stop(
      "`vcov` cannot be supplied when `engine = \"score_bootstrap\"` because ",
      "draws are generated from score perturbations."
    )
  }

  return_baseline_weights <- !newdata_supplied
  if (return_baseline_weights) {
    validate_prediction_weight_ids(baseline_data, id_var)
  }

  score_draws <- generate_score_bootstrap_draws(
    model = model,
    baseline_data = if (return_baseline_weights) baseline_data else NULL,
    id_var = id_var,
    n_draws = n_draws,
    score_weight_dist = score_weight_dist,
    cluster = cluster,
    return_baseline_weights = return_baseline_weights
  )
  out$beta_draws <- score_draws$beta_draws

  out$baseline_weights_draws <- if (return_baseline_weights) {
    score_draws$baseline_weights
  } else {
    if (is_avg || !is.null(by)) {
      warn_fixed_profile_bootstrap_weights("score-bootstrap")
    } else if (isTRUE(return_draws)) {
      warn_fixed_profile_draw_weights("score-bootstrap")
    }
    out$draw_weight_omission_reason <- "user_supplied_newdata"
    NULL
  }

  if (
    attach_individual_weights &&
      !is_avg &&
      is.null(by) &&
      !is.null(out$baseline_weights_draws)
  ) {
    out$draw_weight_col <- "score_weight"
    validate_draw_weight_column_available(prediction_data, out$draw_weight_col)
    out$draw_weights_attached <- isTRUE(return_draws)
  } else if (!is_avg && !is.null(by)) {
    out$draw_weight_omission_reason <- "grouped_sops"
  }
  if (is_avg && is.null(out$draw_weight_omission_reason)) {
    out$draw_weight_omission_reason <- "averaged_sops"
  }

  out
}

apply_sop_simulation_draws <- function(
  n_draws,
  analysis_fn,
  workers,
  globals,
  packages = c("rms", "VGAM", "stats", "markov.misc")
) {
  use_parallel <- !is.null(workers) && workers > 1
  draw_ids <- seq_len(n_draws)

  if (!use_parallel) {
    return(lapply(draw_ids, analysis_fn))
  }
  if (!requireNamespace("furrr", quietly = TRUE)) {
    stop("Package 'furrr' is required for parallel processing")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future.callr::callr, workers = workers)

  furrr::future_map(
    draw_ids,
    analysis_fn,
    .options = furrr::furrr_options(
      seed = TRUE,
      globals = globals,
      packages = packages
    ),
    .progress = FALSE
  )
}
