# SOP inference dispatcher and simulation engine.

#' Inference for State Occupation Probabilities
#'
#' Adds confidence intervals to SOP objects using simulation-based or bootstrap
#' methods. The default method is multivariate-normal coefficient simulation.
#' For objects produced from
#' `rmsb::blrm()` models, posterior uncertainty is already computed by
#' `sops()`/`avg_sops()`, so `inferences()` ignores `method` and returns the
#' object unchanged.
#'
#' @param x A `markov_avg_sops` object from `avg_sops()`, a
#'   `markov_sops` object from `sops()`, or a `markov_avg_comparisons` object
#'   from [avg_comparisons()].
#' @param method Character. Inference method:
#'   \itemize{
#'     \item `"mvn"`: Multivariate-normal coefficient draws.
#'     \item `"score_bootstrap"`: One-step score perturbation using fixed
#'       exponential cluster weights.
#'     \item `"bootstrap"`: Ordinary patient-level refit bootstrap.
#'     \item `"fwb"`: Fractional weighted refits using fixed exponential
#'       patient weights.
#'   }
#' @param n_draws Number of simulation draws (for simulation) or bootstrap
#'   iterations (for bootstrap). Default is 1000. For `blrm` SOP objects this
#'   argument is ignored; rerun `sops()`/`avg_sops()` with `n_draws` and `seed`
#'   to change posterior draws.
#' @param vcov Optional custom variance-covariance matrix. If provided,
#'   overrides the vcov extracted from the model.
#' @param cluster Optional row-level cluster vector for score bootstrap with
#'   `orm` models. Values should match
#'   `id_var` in the baseline data used by `avg_sops()`.
#' @param workers Number of parallel workers. If NULL (default) or 1, uses
#'   sequential processing. If > 1, uses parallel processing with that many
#'   workers.
#' @param conf_level Confidence level for intervals (default 0.95). For `blrm`
#'   SOP objects this argument is ignored; rerun `sops()`/`avg_sops()` with
#'   `conf_level` to change posterior intervals.
#' @param seed Optional integer seed. The caller's complete RNG state is
#'   restored on exit.
#' @param conf_type Type of frequentist confidence interval:
#'   \itemize{
#'     \item `"perc"` (default): Percentile-based intervals from the simulation
#'       distribution.
#'     \item `"wald"`: Uses simulation standard errors with normal quantiles.
#'   }
#' @param null Optional single finite numeric null value. Supplying it adds
#'   Wald `statistic`, `p.value`, and `s.value` columns. A zero null is rejected
#'   for known ratio comparisons.
#' @param return_draws Logical. If TRUE, stores all individual simulation/bootstrap
#'   draws as an attribute. Extract with `get_draws()`. Default is TRUE
#' @param update_datadist Logical. Whether to update datadist for rms models
#'   during bootstrap. Default is TRUE.
#' @param use_coefstart Logical. Use original coefficients as starting values
#'   for bootstrap refitting. Default is FALSE.
#'
#' @return The input object with added columns:
#'   \item{conf.low}{Lower confidence bound}
#'   \item{conf.high}{Upper confidence bound}
#'   \item{std.error}{Standard error from simulation/bootstrap}
#'
#'   If `return_draws = TRUE`, the object also has a `"draws"` attribute
#'   containing all individual draws. Extract
#'   with `get_draws()`. For ungrouped `sops()` objects evaluated on the stored
#'   empirical prediction cohort, score-bootstrap and FWB draws include the
#'   draw-specific `score_weight` or `fwb_weight` column.
#'
#' @details
#' ## Simulation Method
#'
#' The simulation method works as follows:
#' 1. Extract coefficient vector beta_hat and (robust) variance-covariance Sigma
#' 2. Generate n_draws coefficient vectors via the selected engine:
#'    \itemize{
#'      \item `engine = "mvn"`: MVN draws from `N(beta_hat, Sigma)`.
#'      \item `engine = "score_bootstrap"`: one-step score perturbation with
#'        cluster-level exponential multipliers.
#'    }
#' 3. For each draw, replace model coefficients and compute SOPs
#' 4. Compute confidence intervals from the empirical distribution
#'
#' - Works for both individual-level (`sops()`) and averaged (`avg_sops()`) SOPs
#'   for `engine = "mvn"`
#' - `engine = "score_bootstrap"` supports `avg_sops()` and `sops()` with
#'   `robcov_vglm` models and with `orm` models when `cluster` is supplied.
#'   When the prediction rows are the stored empirical cohort, the same
#'   cluster-level weights used for the score perturbation are used for every
#'   empirical averaging step. With `by`, weights are normalized within each
#'   subgroup; without `by`, they are normalized over the full empirical cohort.
#'
#' ## Bootstrap Method
#'
#' The bootstrap method refits the model on bootstrap-weighted data:
#' 1. Sample patient IDs with replacement (`engine = "standard"`) or draw
#'    exponential patient weights normalized to mean 1 (`engine = "fwb"`)
#' 2. Refit model on bootstrap sample (handles missing states through releveling)
#' 3. Compute SOPs using G-computation
#' 4. Compute percentile-based confidence intervals
#' Bootstrap requires the full longitudinal dataset (all time points) either in
#' the original `newdata`/`refit_data` passed to `sops()` or `avg_sops()`, or
#' stored on the model by [orm_markov()], [blrm_markov()], or [vglm_markov()].
#' Standard bootstrap is intentionally limited to marginal `avg_sops()` objects
#' because ordinary resampling can drop outcome-state support needed by fixed
#' individual prediction rows. Fractional weighted bootstrap keeps every row in
#' the refit data with positive patient-level weights, so it is available for
#' individual `sops()` objects. For FWB, the same patient weights used for
#' refitting are used for marginal or subgroup averaging. Ungrouped stored-data
#' `sops()` draws expose the patient weight as `fwb_weight`; grouped `sops()` and
#' `avg_sops()` use the weights internally and return already averaged draws.
#' Prediction data for such ungrouped draw output must not already contain a
#' column named `fwb_weight` or `score_weight`, because those names are reserved
#' for bootstrap draw weights.
#'
#' For marginal `avg_sops()` objects built from user-supplied `newdata`, the
#' supplied rows are fixed standardization profiles. The same fixed-profile rule
#' applies to `sops(newdata = ...)`. Score bootstrap and FWB use the
#' original/refit data for coefficient or refit uncertainty, but do not attach or
#' apply draw weights to the supplied prediction profiles because those rows
#' cannot be assumed to align with the bootstrap clusters.
#'
#'
#' This design ensures consistency: the same vcov is used for both point
#' estimates and inference, regardless of how `inferences()` is called.
#'
#' @seealso [avg_sops()], [sops()], [get_draws()], [robcov_vglm()],
#'   [set_coef()]
#'
#' @examples
#' \dontrun{
#' # Step 1: Fit model with stored data and cluster-robust vcov
#' fit_robust <- vglm_markov(
#'   ordered(y) ~ time + tx + yprev,
#'   family = cumulative(reverse = TRUE, parallel = TRUE),
#'   data = data,
#'   id_var = "id"
#' )
#'
#' # Step 2a: Simulation inference
#' result_sim <- avg_sops(
#'   model = fit_robust,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   y_levels = factor(1:6),
#'   absorb = "6"
#' ) |>
#'   inferences(method = "mvn", n_draws = 1000)
#'
#' # Step 2a-bis: Score bootstrap simulation (requires robcov_vglm)
#' result_score <- avg_sops(
#'   model = fit_robust,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   y_levels = factor(1:6),
#'   absorb = "6"
#' ) |>
#'   inferences(method = "score_bootstrap", n_draws = 1000)
#'
#' # orm_markov() uses rms::robcov() for full robust covariance matrices.
#' dd <- rms::datadist(data)
#' options(datadist = "dd")
#' fit_orm <- orm_markov(y ~ time + tx + yprev, data = data, id_var = "id")
#'
#' avg_orm <- avg_sops(
#'   model = fit_orm,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   y_levels = factor(1:6),
#'   absorb = "6"
#' )
#' result_orm <- avg_orm |>
#'   inferences(method = "mvn", n_draws = 1000)
#'
#' result_orm_score <- avg_orm |>
#'   inferences(
#'     method = "score_bootstrap",
#'     cluster = data$id,
#'     n_draws = 1000
#'   )
#'
#' # Step 2b: Bootstrap inference reuses stored full data
#' result_boot <- avg_sops(
#'   model = fit_robust,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   y_levels = factor(1:6),
#'   absorb = "6"
#' ) |>
#'   inferences(method = "bootstrap", n_draws = 500)
#'
#' # Individual-level SOPs with simulation inference
#' ind_result <- sops(
#'   model = fit_robust,
#'   times = 1:60,
#'   y_levels = factor(1:6),
#'   absorb = "6") |>
#'   inferences(n_draws = 500)
#'
#' # Individual-level refit bootstrap uses FWB.
#' ind_boot <- sops(
#'   model = fit_robust,
#'   times = 1:60,
#'   y_levels = factor(1:6),
#'   absorb = "6") |>
#'   inferences(method = "fwb", n_draws = 500)
#'
#' # Extract draws for custom analyses
#' get_draws(ind_result)
#'
#' # Compute treatment effect on time in state
#' library(dplyr)
#' library(tidyr)
#' treatment_effect <- draws |>
#'   filter(state == 1) |>
#'   group_by(draw_id, tx) |>
#'   summarise(total_time = sum(draw), .groups = "drop") |>
#'   pivot_wider(names_from = tx, values_from = total_time) |>
#'   mutate(effect = `1` - `0`)
#' quantile(treatment_effect$effect, c(0.025, 0.5, 0.975))
#' }
#'
#' @seealso [get_draws()] to extract individual draws from any inference method
#'
#' @export
inferences <- function(
  x,
  method = "mvn",
  n_draws = 1000,
  vcov = NULL,
  cluster = NULL,
  workers = NULL,
  seed = NULL,
  conf_level = 0.95,
  conf_type = "perc",
  null = NULL,
  return_draws = TRUE,
  update_datadist = TRUE,
  use_coefstart = FALSE
) {
  with_local_seed(seed, {
    inferences_impl(
      x = x,
      method = method,
      n_draws = n_draws,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      null = null,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )
  })
}

inferences_impl <- function(
  x,
  method,
  n_draws,
  vcov,
  cluster,
  workers,
  conf_level,
  conf_type,
  null,
  return_draws,
  update_datadist,
  use_coefstart
) {
  # --- Input Validation ---
  if (
    !inherits(
      x,
      c(
        "markov_avg_sops",
        "markov_sops",
        "markov_avg_comparisons"
      )
    )
  ) {
    stop(
      "inferences() requires a 'markov_avg_sops', 'markov_sops', or ",
      "'markov_avg_comparisons' object. ",
      "Got: ",
      paste(class(x), collapse = ", ")
    )
  }

  conf_level <- validate_conf_level(conf_level)
  conf_type <- match.arg(conf_type, choices = c("perc", "wald"))

  method <- match.arg(
    method,
    choices = c("mvn", "score_bootstrap", "bootstrap", "fwb")
  )

  if (inherits(attr(x, "model"), "blrm")) {
    if (!is.null(null)) {
      warning(
        "Wald null tests are not computed for Bayesian posterior outputs.",
        call. = FALSE
      )
    }
    return(x)
  }

  if (
    inherits(x, "markov_avg_comparisons") &&
      identical(unique(x$comparison), "ratio") &&
      identical(null, 0)
  ) {
    stop("Ratio comparisons cannot be tested against `null = 0`.")
  }

  internal_method <- if (method %in% c("mvn", "score_bootstrap")) {
    "simulation"
  } else {
    "bootstrap"
  }
  engine <- switch(
    method,
    mvn = "mvn",
    score_bootstrap = "score_bootstrap",
    bootstrap = "standard",
    fwb = "fwb"
  )

  if (inherits(x, "markov_avg_comparisons")) {
    result <- inferences_avg_comparisons(
      object = x,
      method = internal_method,
      engine = engine,
      score_weight_dist = "exponential",
      n_draws = n_draws,
      vcov = vcov,
      cluster = cluster,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )
  } else if (internal_method == "simulation") {
    result <- inferences_simulation(
      object = x,
      engine = engine,
      score_weight_dist = "exponential",
      n_draws = n_draws,
      vcov = vcov,
      cluster = cluster,
      conf_level = conf_level,
      conf_type = conf_type,
      workers = workers,
      return_draws = return_draws
    )
  } else {
    result <- inferences_bootstrap(
      object = x,
      engine = engine,
      n_boot = n_draws,
      workers = workers,
      conf_level = conf_level,
      conf_type = conf_type,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )
  }

  result <- normalize_inference_result(
    result,
    method = method,
    n_draws = n_draws,
    conf_level = conf_level,
    conf_type = conf_type,
    return_draws = return_draws
  )
  add_null_test(result, null)
}

normalize_inference_result <- function(
  x,
  method,
  n_draws,
  conf_level,
  conf_type,
  return_draws
) {
  primary_class <- class(x)[1]
  draws <- attr(x, "draws")
  attr(x, "draws") <- NULL
  attr(x, "engine") <- NULL
  attr(x, "score_weight_dist") <- NULL
  attr(x, "method") <- method
  attr(x, "n_draws") <- n_draws
  attr(x, "conf_level") <- conf_level
  attr(x, "conf_type") <- conf_type
  attr(x, "draws") <- if (isTRUE(return_draws)) draws else NULL
  x <- order_estimate_columns(x)
  class(x) <- c(primary_class, "data.frame")
  x
}

# =============================================================================
# SIMULATION-BASED INFERENCE (MVN)
# =============================================================================

#' Simulation-Based Inference Using Multivariate Normal
#'
#' Internal function that implements MVN simulation-based inference for SOPs.
#' Draws coefficient vectors from a multivariate normal distribution centered
#' at the point estimates with the (robust) variance-covariance matrix.
#'
#' @param object A `markov_avg_sops` or `markov_sops` object.
#' @param engine Simulation engine (`"mvn"` or `"score_bootstrap"`).
#' @param score_weight_dist Cluster weight distribution for score bootstrap.
#' @param n_draws Number of simulation draws.
#' @param vcov Custom variance-covariance matrix (optional, overrides model vcov).
#' @param cluster Row-level cluster vector for orm score bootstrap.
#' @param conf_level Confidence level.
#' @param conf_type Type of confidence interval ("perc" or "wald").
#' @param workers Number of parallel workers. If NULL or 1, uses sequential
#'   processing. If > 1, uses parallel processing.
#' @param return_draws Store individual draws?
#'
#' @return The input object with added confidence intervals.
#'
#' @keywords internal
inferences_simulation <- function(
  object,
  engine,
  score_weight_dist,
  n_draws,
  vcov,
  cluster,
  conf_level,
  conf_type,
  workers,
  return_draws
) {
  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  newdata_pred_stored <- attr(object, "newdata_pred")
  newdata_supplied <- isTRUE(attr(object, "newdata_supplied"))
  call_args <- attr(object, "call_args")
  avg_args <- attr(object, "avg_args")

  time_var <- attr(object, "time_var")
  p_var <- attr(object, "p_var")
  p2_var <- attr(object, "p2_var")
  y_levels <- attr(object, "y_levels")
  absorb <- attr(object, "absorb")
  gap_var <- attr(object, "gap_var")
  time_covariates <- attr(object, "time_covariates")

  # For avg_sops objects
  is_avg <- inherits(object, "markov_avg_sops")

  if (is_avg) {
    variables <- avg_args$variables
    by <- avg_args$by
    times <- avg_args$times
    id_var <- avg_args$id_var
  } else {
    # For individual sops objects
    times <- call_args$times
    variables <- NULL
    by <- call_args$by %||% attr(object, "by")
    id_var <- attr(object, "id_var") %||% markov_model_id_var(model)
  }

  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform simulation inference.")
  }

  # --- 2. Prepare Prediction Data (COMPUTED ONCE) ---
  if (is_avg) {
    # For avg_sops: create counterfactual datasets
    grid <- do.call(expand.grid, variables)
    n_cf <- nrow(grid)
    if (!is.null(newdata_pred_stored)) {
      newdata_pred <- newdata_pred_stored
      if (nrow(newdata_pred) %% n_cf != 0L) {
        stop("Stored prediction data is not aligned with counterfactual grid.")
      }
      n_each <- nrow(newdata_pred) / n_cf
      baseline_data <- newdata_pred[seq_len(n_each), , drop = FALSE]
    } else {
      baseline_data <- resolve_markov_prediction_data(
        newdata_orig,
        id_var = id_var,
        time_var = time_var
      )
      newdata_pred <- create_counterfactual_data(baseline_data, grid, variables)
      n_each <- nrow(baseline_data)
    }
  } else {
    # For individual sops: use data directly
    newdata_pred <- newdata_pred_stored %||% newdata_orig
    grid <- NULL
    n_cf <- 1
    n_each <- nrow(newdata_pred)
    baseline_data <- if (newdata_supplied) NULL else newdata_pred
  }

  # --- 3. Generate Coefficient Draws ---
  simulation_draws <- generate_sop_coefficient_draws(
    model = model,
    engine = engine,
    score_weight_dist = score_weight_dist,
    n_draws = n_draws,
    vcov = vcov,
    cluster = cluster,
    baseline_data = baseline_data,
    id_var = id_var,
    newdata_supplied = newdata_supplied,
    is_avg = is_avg,
    by = by,
    return_draws = return_draws,
    prediction_data = newdata_pred
  )
  beta_draws <- simulation_draws$beta_draws
  baseline_weights_draws <- simulation_draws$baseline_weights_draws
  draw_weight_col <- simulation_draws$draw_weight_col
  draw_weights_attached <- simulation_draws$draw_weights_attached
  draw_weight_omission_reason <-
    simulation_draws$draw_weight_omission_reason

  n_times <- length(times)
  n_states <- length(y_levels)

  if (is_avg) {
    group_cols <- c("time", "state", names(variables))
    if (!is.null(by)) {
      group_cols <- unique(c(group_cols, by))
    }
  } else if (is.null(by)) {
    group_cols <- c("rowid", "time", "state")
  } else {
    group_cols <- unique(c("time", "state", by))
  }
  object_keys <- sop_draw_cell_key(object, group_cols)

  # --- 4. Detect Fast Path Eligibility ---
  # The fast path uses pre-computed design matrix decompositions for supported
  # ordinal model backends.
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  use_fast_path <- inherits(model_chk, c("vglm", "orm")) && is.null(p2_var)

  if (use_fast_path) {
    # --- FAST PATH: Pre-build components once, then run efficient Markov loop ---
    execution_plan <- tryCatch(
      compile_sop_execution_plan(
        model = model_chk,
        newdata = newdata_pred,
        time_covariates = time_covariates,
        times = times,
        y_levels = y_levels,
        absorb = absorb,
        time_var = time_var,
        p_var = p_var,
        gap_var = gap_var,
        builder = "streamed",
        output = if (is_avg) "average" else "individual"
      ),
      error = function(e) {
        warning(
          "Fast path build failed, falling back to slow path: ",
          e$message
        )
        NULL
      }
    )

    if (!is.null(execution_plan)) {
      # Define fast analysis function
      analysis_fn <- function(i) {
        baseline_weights <- NULL
        if (!is.null(baseline_weights_draws)) {
          baseline_weights <- baseline_weights_draws[i, ]
        }

        # Compute effective coefficients from beta draw
        Gamma_i <- get_effective_coefs(model_chk, beta = beta_draws[i, ])

        # Run fast Markov simulation
        sops_array <- tryCatch(
          run_sop_execution_plan(execution_plan, Gamma_i),
          error = function(e) {
            warning("markov_msm_run failed in draw ", i, ": ", e$message)
            return(NULL)
          }
        )

        if (is.null(sops_array)) {
          return(NULL)
        }

        # Marginalize if needed (for avg_sops)
        if (is_avg) {
          result <- marginalize_sops_array(
            sops_array = sops_array,
            grid = grid,
            times = times,
            y_levels = y_levels,
            variables = variables,
            n_cf = n_cf,
            n_each = n_each,
            weights = baseline_weights,
            by = by,
            newdata = newdata_pred
          )
        } else {
          # Individual-level: convert array to data frame
          result <- array_to_df_individual(
            sops_array,
            times,
            y_levels,
            newdata_pred,
            by = by,
            weights = baseline_weights,
            weight_col = draw_weight_col
          )
        }

        pack_sop_draw_result(
          result,
          group_cols,
          object_keys,
          draw_weight_col
        )
      }
    } else {
      # Fast path build failed, fall back to slow path
      use_fast_path <- FALSE
    }
  }

  if (!use_fast_path) {
    # --- SLOW PATH: Use soprob_markov with coefficient replacement ---
    analysis_fn <- function(i) {
      baseline_weights <- NULL
      if (!is.null(baseline_weights_draws)) {
        baseline_weights <- baseline_weights_draws[i, ]
      }

      # Replace coefficients
      model_i <- set_coef(model, beta_draws[i, ])

      # Compute SOPs
      sops_array <- tryCatch(
        soprob_markov(
          model = model_i,
          newdata = newdata_pred,
          times = times,
          y_levels = y_levels,
          absorb = absorb,
          time_var = time_var,
          p_var = p_var,
          p2_var = p2_var,
          gap_var = gap_var,
          time_covariates = time_covariates
        ),
        error = function(e) {
          warning("soprob_markov failed in draw ", i, ": ", e$message)
          return(NULL)
        }
      )

      if (is.null(sops_array)) {
        return(NULL)
      }

      # Marginalize if needed (for avg_sops)
      if (is_avg) {
        result <- marginalize_sops_array(
          sops_array = sops_array,
          grid = grid,
          times = times,
          y_levels = y_levels,
          variables = variables,
          n_cf = n_cf,
          n_each = n_each,
          weights = baseline_weights,
          by = by,
          newdata = newdata_pred
        )
      } else {
        # Individual-level: convert array to data frame
        result <- array_to_df_individual(
          sops_array,
          times,
          y_levels,
          newdata_pred,
          by = by,
          weights = baseline_weights,
          weight_col = draw_weight_col
        )
      }

      pack_sop_draw_result(
        result,
        group_cols,
        object_keys,
        draw_weight_col
      )
    }
  }

  # --- 6. Apply Across All Draws ---
  globals_list <- c(
    "model",
    "model_chk",
    "beta_draws",
    "newdata_pred",
    "times",
    "y_levels",
    "absorb",
    "time_var",
    "p_var",
    "p2_var",
    "gap_var",
    "time_covariates",
    "is_avg",
    "grid",
    "variables",
    "n_cf",
    "n_each",
    "baseline_weights_draws",
    "draw_weight_col",
    "group_cols",
    "object_keys",
    "use_fast_path",
    "by"
  )
  if (use_fast_path) {
    globals_list <- c(globals_list, "execution_plan")
  }

  sim_results <- apply_sop_simulation_draws(
    n_draws = nrow(beta_draws),
    analysis_fn = analysis_fn,
    workers = workers,
    globals = globals_list
  )

  # --- 6. Combine fixed-order draw cells ---
  successful_ids <- which(!vapply(sim_results, is.null, logical(1)))
  sim_results <- sim_results[successful_ids]

  if (length(sim_results) == 0) {
    stop("All simulation draws failed.")
  }

  draw_values <- do.call(rbind, lapply(sim_results, `[[`, "estimate"))
  weight_values <- NULL
  has_weights <- vapply(sim_results, function(x) !is.null(x$weight), logical(1))
  if (any(has_weights) && !all(has_weights)) {
    stop("Simulation draw weights are not aligned across successful draws.")
  }
  if (all(has_weights)) {
    weight_values <- do.call(rbind, lapply(sim_results, `[[`, "weight"))
  }

  # --- 7. Compute Confidence Intervals ---
  summary_stats <- summarize_sop_draw_matrix(
    draw_values = draw_values,
    conf_level = conf_level,
    conf_type = conf_type,
    point_estimates = object$estimate
  )

  # --- 8. Attach summaries by the already-known cell order ---
  final_result <- object
  final_result$conf.low <- summary_stats$conf.low
  final_result$conf.high <- summary_stats$conf.high
  final_result$std.error <- summary_stats$std.error

  final_result <- restore_sops_attrs(final_result, object)
  # Add metadata
  attr(final_result, "n_draws") <- n_draws
  attr(final_result, "n_successful") <- length(sim_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "conf_type") <- conf_type
  attr(final_result, "method") <- engine
  if (engine == "score_bootstrap") {
    attr(final_result, "score_weight_dist") <- score_weight_dist
    attr(final_result, "draw_weights_attached") <- draw_weights_attached
    attr(final_result, "draw_weight_col") <- draw_weight_col
    if (!is.null(draw_weight_omission_reason)) {
      attr(final_result, "draw_weight_omission_reason") <-
        draw_weight_omission_reason
    }
  }

  if (return_draws) {
    attr(final_result, "draws") <- sop_draw_matrix_to_df(
      draw_values,
      object,
      successful_ids,
      draw_weight_col,
      weight_values
    )
  }

  final_result
}

sop_draw_cell_key <- function(data, group_cols) {
  missing <- setdiff(group_cols, names(data))
  if (length(missing)) {
    stop("SOP draw cells are missing: ", paste(missing, collapse = ", "), ".")
  }
  key <- do.call(
    paste,
    c(lapply(data[, group_cols, drop = FALSE], as.character), sep = "\r")
  )
  if (anyDuplicated(key)) {
    stop("SOP draw cell identifiers are not unique.")
  }
  key
}

pack_sop_draw_result <- function(result, group_cols, object_keys, weight_col) {
  result_keys <- sop_draw_cell_key(result, group_cols)
  index <- match(object_keys, result_keys)
  if (anyNA(index)) {
    stop("A simulation draw did not return every requested SOP cell.")
  }
  list(
    estimate = result$estimate[index],
    weight = if (!is.null(weight_col) && weight_col %in% names(result)) {
      result[[weight_col]][index]
    } else {
      NULL
    }
  )
}

summarize_sop_draw_matrix <- function(
  draw_values,
  conf_level,
  conf_type,
  point_estimates
) {
  alpha <- 1 - validate_conf_level(conf_level)
  standard_error <- apply(draw_values, 2L, stats::sd, na.rm = TRUE)
  if (identical(conf_type, "perc")) {
    lower <- apply(
      draw_values,
      2L,
      stats::quantile,
      probs = alpha / 2,
      na.rm = TRUE,
      names = FALSE,
      type = 7L
    )
    upper <- apply(
      draw_values,
      2L,
      stats::quantile,
      probs = 1 - alpha / 2,
      na.rm = TRUE,
      names = FALSE,
      type = 7L
    )
  } else {
    critical <- abs(stats::qnorm(alpha / 2))
    lower <- point_estimates - critical * standard_error
    upper <- point_estimates + critical * standard_error
  }
  data.frame(
    conf.low = as.numeric(lower),
    conf.high = as.numeric(upper),
    std.error = as.numeric(standard_error)
  )
}

sop_draw_matrix_to_df <- function(
  draw_values,
  object,
  draw_ids,
  weight_col = NULL,
  weight_values = NULL
) {
  metadata_cols <- setdiff(names(object), inference_columns())
  metadata <- as.data.frame(object)[, metadata_cols, drop = FALSE]
  n_cells <- nrow(metadata)
  out <- metadata[
    rep(seq_len(n_cells), times = length(draw_ids)),
    ,
    drop = FALSE
  ]
  out$estimate <- as.vector(t(draw_values))
  if (!is.null(weight_col) && !is.null(weight_values)) {
    out[[weight_col]] <- as.vector(t(weight_values))
  }
  out$draw_id <- rep(draw_ids, each = n_cells)
  rownames(out) <- NULL
  out
}

warn_fixed_profile_bootstrap_weights <- function(engine) {
  warning(
    "User-supplied `newdata` is treated as fixed standardization profiles; ",
    engine,
    " baseline weights have been set to NULL and profiles will be averaged ",
    "equally.",
    call. = FALSE
  )
}

warn_fixed_profile_draw_weights <- function(engine) {
  warning(
    "User-supplied `newdata` is treated as fixed prediction profiles; ",
    engine,
    " draw weights were not attached because prediction rows cannot be ",
    "assumed to align with bootstrap clusters.",
    call. = FALSE
  )
}
