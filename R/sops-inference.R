# SOP inference dispatcher and simulation engine.

#' Inference for State Occupation Probabilities
#'
#' Adds confidence intervals to SOP objects using simulation-based or bootstrap
#' methods. The default method is simulation. For objects produced from
#' `rmsb::blrm()` models, posterior uncertainty is already computed by
#' `sops()`/`avg_sops()`, so `inferences()` ignores `method` and returns the
#' object unchanged.
#'
#' @param object A `markov_avg_sops` object from `avg_sops()`, a
#'   `markov_sops` object from `sops()`, or a `markov_avg_comparisons` object
#'   from [avg_comparisons()].
#' @param method Character. Inference method:
#'   \itemize{
#'     \item `"simulation"` (default): Uses simulation engines that do not
#'       refit models. Works for both individual and averaged SOPs.
#'     \item `"bootstrap"`: Refits the model on resampled or weighted patient
#'       data. Standard resampling works for `markov_avg_sops`; fractional
#'       weighted bootstrap also works for individual `markov_sops`.
#'   }
#' @param engine Character. Inference engine. For `method = "simulation"`:
#'   \itemize{
#'     \item `"mvn"` (default): Draws coefficients from MVN(beta_hat, Sigma).
#'     \item `"score_bootstrap"`: Uses one-step score perturbation with
#'       cluster-level exponential multipliers. Requires a `robcov_vglm` model
#'       or an `orm` model with `cluster`.
#'   }
#'   For `method = "bootstrap"`, use `"standard"` for ordinary patient
#'   resampling with replacement or `"fwb"` for fractional weighted bootstrap
#'   refits using exponential patient-level weights. `engine = "fwb"` is the
#'   only refit-bootstrap engine supported for individual `sops()` objects. The
#'   legacy default `"mvn"` maps to `"standard"` when
#'   `method = "bootstrap"`.
#' @param score_weight_dist Character. Cluster weight distribution for
#'   `engine = "score_bootstrap"`. Currently only `"exponential"` is supported.
#' @param n_sim Number of simulation draws (for simulation) or bootstrap
#'   iterations (for bootstrap). Default is 1000. For `blrm` SOP objects this
#'   argument is ignored; rerun `sops()`/`avg_sops()` with `n_draws` and `seed`
#'   to change posterior draws.
#' @param vcov Optional custom variance-covariance matrix. If provided,
#'   overrides the vcov extracted from the model.
#' @param cluster Optional row-level cluster vector for
#'   `engine = "score_bootstrap"` with `orm` models. Values should match
#'   `id_var` in the baseline data used by `avg_sops()`.
#' @param workers Number of parallel workers. If NULL (default) or 1, uses
#'   sequential processing. If > 1, uses parallel processing with that many
#'   workers.
#' @param conf_level Confidence level for intervals (default 0.95). For `blrm`
#'   SOP objects this argument is ignored; rerun `sops()`/`avg_sops()` with
#'   `conf_level` to change posterior intervals.
#' @param conf_type Type of confidence interval (simulation method only):
#'   \itemize{
#'     \item `"perc"` (default): Percentile-based intervals from the simulation
#'       distribution.
#'     \item `"wald"`: Uses simulation standard errors with normal quantiles.
#'   }
#' @param return_draws Logical. If TRUE, stores all individual simulation/bootstrap
#'   draws as an attribute. Extract with `get_draws()`. Default is TRUE
#' @param update_datadist Logical. Whether to update datadist for rms models
#'   during bootstrap. Default is TRUE.
#' @param use_coefstart Logical. Use original coefficients as starting values
#'   for bootstrap refitting. Default is FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object with added columns:
#'   \item{conf.low}{Lower confidence bound}
#'   \item{conf.high}{Upper confidence bound}
#'   \item{std.error}{Standard error from simulation/bootstrap}
#'
#'   If `return_draws = TRUE`, the object also has a `"simulation_draws"` or
#'   `"bootstrap_draws"` attribute containing all individual draws. Extract
#'   with `get_draws()`. For ungrouped `sops()` objects evaluated on the stored
#'   empirical prediction cohort, score-bootstrap and FWB draws include the
#'   draw-specific `score_weight` or `fwb_weight` column.
#'
#' @details
#' ## Simulation Method
#'
#' The simulation method works as follows:
#' 1. Extract coefficient vector beta_hat and (robust) variance-covariance Sigma
#' 2. Generate n_sim coefficient vectors via the selected engine:
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
#'   ylevels = factor(1:6),
#'   absorb = "6"
#' ) |>
#'   inferences(method = "simulation", n_sim = 1000)
#'
#' # Step 2a-bis: Score bootstrap simulation (requires robcov_vglm)
#' result_score <- avg_sops(
#'   model = fit_robust,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6"
#' ) |>
#'   inferences(method = "simulation", engine = "score_bootstrap", n_sim = 1000)
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
#'   ylevels = factor(1:6),
#'   absorb = "6"
#' )
#' result_orm <- avg_orm |>
#'   inferences(method = "simulation", engine = "mvn", n_sim = 1000)
#'
#' result_orm_score <- avg_orm |>
#'   inferences(
#'     method = "simulation",
#'     engine = "score_bootstrap",
#'     cluster = data$id,
#'     n_sim = 1000
#'   )
#'
#' # Step 2b: Bootstrap inference reuses stored full data
#' result_boot <- avg_sops(
#'   model = fit_robust,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6"
#' ) |>
#'   inferences(method = "bootstrap", n_sim = 500)
#'
#' # Individual-level SOPs with simulation inference
#' ind_result <- sops(
#'   model = fit_robust,
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6") |>
#'   inferences(n_sim = 500)
#'
#' # Individual-level refit bootstrap uses FWB.
#' ind_boot <- sops(
#'   model = fit_robust,
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6") |>
#'   inferences(method = "bootstrap", engine = "fwb", n_sim = 500)
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
  object,
  method = "simulation",
  engine = "mvn",
  score_weight_dist = "exponential",
  n_sim = 1000,
  vcov = NULL,
  cluster = NULL,
  workers = NULL,
  conf_level = 0.95,
  conf_type = "perc",
  return_draws = TRUE,
  update_datadist = TRUE,
  use_coefstart = FALSE,
  ...
) {
  # --- Input Validation ---
  if (
    !inherits(
      object,
      c("markov_avg_sops", "markov_sops", "markov_avg_comparisons")
    )
  ) {
    stop(
      "inferences() requires a 'markov_avg_sops', 'markov_sops', or ",
      "'markov_avg_comparisons' object. ",
      "Got: ",
      paste(class(object), collapse = ", ")
    )
  }

  if (inherits(attr(object, "model"), "blrm")) {
    return(object)
  }

  method <- match.arg(method, choices = c("simulation", "bootstrap"))
  engine <- match.arg(
    engine,
    choices = c("mvn", "score_bootstrap", "standard", "fwb")
  )

  if (method == "simulation" && engine %in% c("standard", "fwb")) {
    stop(
      "`engine = \"",
      engine,
      "\"` is only used when `method = \"bootstrap\"`."
    )
  }
  if (method == "bootstrap" && engine == "score_bootstrap") {
    stop(
      "`engine = \"score_bootstrap\"` is only used when ",
      "`method = \"simulation\"`."
    )
  }

  if (inherits(object, "markov_avg_comparisons")) {
    return(inferences_avg_comparisons(
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

  # --- Dispatch to Method-Specific Implementation ---
  if (method == "simulation") {
    inferences_simulation(
      object = object,
      engine = engine,
      score_weight_dist = score_weight_dist,
      n_sim = n_sim,
      vcov = vcov,
      cluster = cluster,
      conf_level = conf_level,
      conf_type = conf_type,
      workers = workers,
      return_draws = return_draws
    )
  } else if (method == "bootstrap") {
    bootstrap_engine <- if (engine == "mvn") "standard" else engine
    inferences_bootstrap(
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
#' @param n_sim Number of simulation draws.
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
  n_sim,
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

  tvarname <- attr(object, "tvarname")
  pvarname <- attr(object, "pvarname")
  p2varname <- attr(object, "p2varname")
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
  gap <- attr(object, "gap")
  t_covs <- attr(object, "t_covs")

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
        tvarname = tvarname
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
  baseline_weights_draws <- NULL
  draw_weight_col <- NULL
  draw_weights_attached <- FALSE
  draw_weight_omission_reason <- NULL
  if (engine == "mvn") {
    # Check for mvtnorm package
    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
      stop(
        "Package 'mvtnorm' is required for MVN simulation inference.\n",
        "Install with: install.packages('mvtnorm')"
      )
    }

    # The vcov is extracted directly from the model object:
    # - For robcov_vglm: uses the stored cluster-robust vcov
    # - For orm with rms::robcov(): uses the stored robust vcov
    # - For plain vglm/orm: uses model-based vcov
    beta_hat <- get_coef(model)

    if (!is.null(vcov) && is.matrix(vcov)) {
      Sigma <- vcov
      Sigma <- validate_coef_vcov(beta_hat, Sigma, arg = "vcov")
    } else {
      Sigma <- get_vcov_robust(model)
      Sigma <- validate_coef_vcov(beta_hat, Sigma, arg = "model vcov")
    }

    beta_draws <- mvtnorm::rmvnorm(n_sim, mean = beta_hat, sigma = Sigma)
  } else if (engine == "score_bootstrap") {
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
      n_sim = n_sim,
      score_weight_dist = score_weight_dist,
      cluster = cluster,
      return_baseline_weights = return_baseline_weights
    )
    beta_draws <- score_draws$beta_draws
    baseline_weights_draws <- if (return_baseline_weights) {
      score_draws$baseline_weights
    } else {
      if (is_avg || !is.null(by)) {
        warn_fixed_profile_bootstrap_weights("score-bootstrap")
      } else if (isTRUE(return_draws)) {
        warn_fixed_profile_draw_weights("score-bootstrap")
      }
      draw_weight_omission_reason <- "user_supplied_newdata"
      NULL
    }

    if (!is_avg && is.null(by) && !is.null(baseline_weights_draws)) {
      draw_weight_col <- "score_weight"
      validate_draw_weight_column_available(newdata_pred, draw_weight_col)
      draw_weights_attached <- isTRUE(return_draws)
    } else if (!is_avg && !is.null(by)) {
      draw_weight_omission_reason <- "grouped_sops"
    }
    if (is_avg && is.null(draw_weight_omission_reason)) {
      draw_weight_omission_reason <- "averaged_sops"
    }
  } else {
    stop("Unknown simulation engine: ", engine)
  }

  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 4. Detect Fast Path Eligibility ---
  # The fast path uses pre-computed design matrix decompositions for supported
  # ordinal model backends.
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  use_fast_path <- inherits(model_chk, c("vglm", "orm")) && is.null(p2varname)

  if (use_fast_path) {
    # --- FAST PATH: Pre-build components once, then run efficient Markov loop ---
    components <- tryCatch(
      markov_msm_build(
        model = model_chk,
        data = newdata_pred,
        t_covs = t_covs,
        times = times,
        ylevels = ylevels,
        absorb = absorb,
        tvarname = tvarname,
        pvarname = pvarname,
        gap = gap
      ),
      error = function(e) {
        warning(
          "Fast path build failed, falling back to slow path: ",
          e$message
        )
        NULL
      }
    )

    if (!is.null(components)) {
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
          markov_msm_run(components, Gamma_i, times, absorb),
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
            ylevels = ylevels,
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
            ylevels,
            newdata_pred,
            by = by,
            weights = baseline_weights,
            weight_col = draw_weight_col
          )
        }

        result
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
          object = model_i,
          data = newdata_pred,
          times = times,
          ylevels = ylevels,
          absorb = absorb,
          tvarname = tvarname,
          pvarname = pvarname,
          p2varname = p2varname,
          gap = gap,
          t_covs = t_covs
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
          ylevels = ylevels,
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
          ylevels,
          newdata_pred,
          by = by,
          weights = baseline_weights,
          weight_col = draw_weight_col
        )
      }

      result
    }
  }

  # --- 6. Apply Across All Draws ---
  use_parallel <- !is.null(workers) && workers > 1

  if (use_parallel) {
    # Setup parallel
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("Package 'furrr' is required for parallel processing")
    }

    future::plan(future.callr::callr, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    # Build globals list dynamically based on fast/slow path
    globals_list <- c(
      "model",
      "model_chk",
      "beta_draws",
      "newdata_pred",
      "times",
      "ylevels",
      "absorb",
      "tvarname",
      "pvarname",
      "p2varname",
      "gap",
      "t_covs",
      "is_avg",
      "grid",
      "variables",
      "n_cf",
      "n_each",
      "baseline_weights_draws",
      "draw_weight_col",
      "use_fast_path",
      "by"
    )
    if (use_fast_path) {
      globals_list <- c(globals_list, "components")
    }

    sim_results <- furrr::future_map(
      seq_len(n_sim),
      analysis_fn,
      .options = furrr::furrr_options(
        seed = TRUE,
        globals = globals_list,
        packages = c("rms", "VGAM", "stats", "markov.misc")
      ),
      .progress = FALSE
    )
  } else {
    sim_results <- lapply(seq_len(n_sim), analysis_fn)
  }

  # --- 6. Combine Results ---
  sim_results <- Filter(Negate(is.null), sim_results)

  if (length(sim_results) == 0) {
    stop("All simulation draws failed.")
  }

  # Add draw_id
  for (i in seq_along(sim_results)) {
    sim_results[[i]]$draw_id <- i
  }
  draws_df <- bind_rows_fill(sim_results)

  # --- 7. Compute Confidence Intervals ---
  if (is_avg) {
    group_cols <- c("time", "state", names(variables))
    if (!is.null(by)) {
      group_cols <- unique(c(group_cols, by))
    }
  } else {
    # For individual SOPs: include rowid unless 'by' is specified
    if (is.null(by)) {
      group_cols <- c("rowid", "time", "state")
    } else {
      group_cols <- unique(c("time", "state", by))
    }
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = draws_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )

  # --- 8. Merge with Original Object ---
  final_result <- merge(object, summary_stats, by = group_cols, all.x = TRUE)

  final_result <- restore_sops_attrs(final_result, object)
  # Add metadata
  attr(final_result, "n_sim") <- n_sim
  attr(final_result, "n_successful") <- length(sim_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "conf_type") <- conf_type
  attr(final_result, "method") <- "simulation"
  attr(final_result, "engine") <- engine
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
    attr(final_result, "simulation_draws") <- draws_df
  }

  final_result
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
