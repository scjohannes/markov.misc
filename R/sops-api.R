# SOP public APIs and BLRM summaries.

#' Calculate Individual State Occupation Probabilities
#'
#' Computes individual-level state occupation probabilities (SOPs) for each row
#' in a dataset. Optionally aggregates results within strata defined by grouping
#' variables.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, or `blrm`). For
#'   `vglm` models, the family must be cumulative-logit with `reverse = TRUE`.
#'   For `orm` models, the family must be logistic.
#' @param newdata Optional. A data frame of prediction profiles. When supplied,
#'   every row is treated as a separate fixed profile and the internal
#'   `rowid` column is regenerated. If `NULL`, uses the data stored by
#'   [orm_markov()], [blrm_markov()], or [vglm_markov()] and requires exactly
#'   one complete designated starting profile per fitted patient. `refit_data`
#'   is never used as a prediction-profile fallback.
#' @param times Required visit-scale time points to estimate. For numeric time
#'   variables this is usually a numeric vector. For factor-valued visit
#'   indices, values are matched to fitted visit levels.
#' @param y_levels A vector of state levels in fitted outcome order. Its length
#'   must agree with the fitted model's threshold count. If `NULL`, attempts to
#'   infer the levels from `model`.
#' @param absorb The absorbing state.
#' @param by Optional character vector of variable names to stratify by. When
#'   provided, the function aggregates (averages) SOPs within each stratum
#'   defined by combinations of these variables. E.g., `by = "ecog"` aggregates
#'   within ECOG levels; `by = c("ecog", "age_group")` aggregates within each
#'   combination of ECOG and age group. Rows missing any grouping value are
#'   omitted; an error is raised if none remain. This is simple aggregation
#'   within observed strata, not G-computation standardization (use
#'   `avg_sops()` for that).
#' @param refit_data Optional full longitudinal data used only by refit-bootstrap
#'   inference. It is not used for point estimates. Defaults to data stored on
#'   wrapper-fitted models.
#' @param id_var Character ID column used when `include_re = TRUE`. For `blrm`,
#'   `NULL` is inferred from wrapper metadata, then `model$clusterInfo$name`
#'   when available, otherwise `"id"`. For frequentist models, `id_var` is used
#'   for stored starting-profile validation and refit bootstrap metadata, not
#'   for ordinary user-supplied prediction rows.
#' @param time_var Name of the time variable in the model.
#' @param p_var Name of the previous state variable in the model.
#' @param p2_var Optional second previous-state variable. `NULL` uses a
#'   first-order Markov recursion; a non-`NULL` column name uses a second-order
#'   recursion.
#' @param gap_var Name of the time gap_var variable (if used).
#' @param time_covariates Optional time-varying covariate lookup table for
#'   explicit basis columns. Registered inline terms such as
#'   `rms::rcs(time, 4)` and `rms::lsp(time, c(3, 7))` can be used without
#'   supplying `time_covariates`.
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible random draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior SOP draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior SOP draws for
#'   extraction with `get_draws()`. Large draw objects are guarded to protect
#'   memory.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class `markov_sops` containing:
#'   \item{rowid}{Row identifier from newdata (omitted if `by` is specified)}
#'   \item{time}{Time point}
#'   \item{state}{State name}
#'   \item{estimate}{Probability of being in the state (individual-level or stratum average)}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     uncertainty summaries computed directly from SOP draws}
#'   \item{(by variables)}{Stratification variables if `by` is specified}
#'   Plus all columns from `newdata` (for individual-level results).
#'
#' @details
#' This function wraps `soprob_markov()` and converts its array output to a tidy
#' data frame. The output contains one row per patient-time-state combination
#' (or stratum-time-state if `by` is used).
#'
#' Automatic prediction from a wrapper-fitted model uses the cohort-wide
#' designated `start_time` retained before response-driven model-frame omission.
#' The transition response is not required on that row, but ID, modeled
#' predictors, time information, and previous state must be complete. Each
#' included patient must contribute at least one usable fitted transition
#' somewhere; later rows are not substituted for a missing starting profile.
#'
#' For `rmsb::blrm()` models, SOPs are computed on sampled posterior draws and
#' then summarized. The point estimate is the requested posterior summary of the
#' SOP draws, not a plug-in calculation from summarized model parameters.
#' State-wise medians and interval bounds are not constrained to sum to one
#' across states; use `posterior_summary = "mean"` when the displayed estimates
#' themselves need to preserve total probability. Draw-level probabilities
#' stored with `return_draws = TRUE` remain normalized within each draw.
#'
#' **Model Requirements:**
#'
#' For `vglm` models, only the cumulative-logit family with `reverse = TRUE` is
#' supported. For `orm` models, only the logistic family is supported. This is
#' because the package's Markov simulation logic expects
#' higher-numbered states to represent worse outcomes, and uses reverse cumulative
#' probabilities to model the probability of being in state k or worse.
#'
#' **Stratification:**
#'
#' When `by` is specified, the function:
#' 1. Computes individual-level SOPs for all patients
#' 2. Groups results by time, state, and the specified variables
#' 3. Averages the SOPs within each group
#'
#' This aggregation preserves heterogeneity across observed strata without
#' creating counterfactual scenarios (unlike G-computation in `avg_sops()`).
#'
#' For computing marginal/standardized SOPs (G-computation), use `avg_sops()`
#' instead, which creates counterfactual datasets and averages over individuals.
#'
#' @seealso [avg_sops()] for marginal SOPs, [soprob_markov()] for the underlying
#'   computation.
#'
#' @examples
#' \dontrun{
#' # Individual-level SOPs
#' sops_ind <- sops(fit, newdata = baseline_data, times = 1:30, y_levels = 1:6)
#'
#' # Stratified by ECOG performance status
#' sops_ecog <- sops(fit, newdata = baseline_data, times = 1:30, y_levels = 1:6,
#'                   by = "ecog")
#'
#' # Stratified by multiple variables
#' sops_strat <- sops(fit, newdata = baseline_data, times = 1:30, y_levels = 1:6,
#'                    by = c("ecog", "age_group"))
#'
#' # Compatible with inferences()
#' sops_strat |> inferences(method = "mvn", n_draws = 1000)
#' }
#'
#' @importFrom stats model.frame
#' @export
sops <- function(
  model,
  newdata = NULL,
  times,
  y_levels = NULL,
  absorb = NULL,
  by = NULL,
  refit_data = NULL,
  id_var = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  # --- 1. Setup & Defaults ---
  # Validate model compatibility
  validate_markov_model(model)
  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied to `sops()`.")
  }
  conf_level <- validate_conf_level(conf_level)

  data_res <- resolve_markov_source_data(
    model,
    newdata,
    refit_data,
    time_var = time_var,
    p_var = p_var
  )
  newdata_orig <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied
  id_var <- markov_model_id_var(model, id_var)

  if (inherits(model, "blrm") && isTRUE(include_re)) {
    id_var <- resolve_blrm_id_var(model, newdata_orig, id_var)
    validate_markov_id_var(id_var, newdata_orig, "newdata")
  }

  if (!is.null(id_var)) {
    if (!is.null(refit_data)) {
      validate_markov_id_var(id_var, refit_data, "refit_data")
    }
  }

  if (!newdata_supplied) {
    if (is.null(id_var)) {
      stop(
        "Automatic SOP prediction from stored model data requires `id_var`. ",
        "Fit with `orm_markov(id_var = ...)`, `blrm_markov(id_var = ...)`, ",
        "`vglm_markov(id_var = ...)`, or supply explicit `newdata`."
      )
    }
    validate_markov_id_var(id_var, newdata_orig, "stored model data")
  }

  # Validate and coerce supplied times before extracting one prediction row per ID.
  time_res <- resolve_sop_times(
    model,
    newdata_orig,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "unique"
  )
  times <- time_res$times
  validate_factor_gap(gap_var, time_covariates, time_res$time_info)

  newdata <- if (newdata_supplied) {
    newdata_orig
  } else {
    resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      time_var = time_var
    )
  }
  newdata <- ensure_markov_rowid(newdata)

  if (is.null(y_levels)) {
    # Try to infer from model
    if (inherits(model, "vglm")) {
      # VGAM stores response levels
      y_levels <- model@extra$colnames.y
      if (is.null(y_levels)) stop("`y_levels` cannot be NULL")
    } else if (inherits(model, "robcov_vglm")) {
      # robcov_vglm stores extra slot in list
      y_levels <- model$extra$colnames.y
      if (is.null(y_levels)) stop("`y_levels` cannot be NULL")
    } else if (inherits(model, "blrm")) {
      y_levels <- model$ylevels
      if (is.null(y_levels)) stop("`y_levels` cannot be NULL")
    } else if (inherits(model, "orm")) {
      # rms orm stores levels
      y_levels <- model$yunique
      if (is.null(y_levels)) stop("`y_levels` cannot be NULL")
    } else {
      stop("`y_levels` cannot be NULL")
    }
  }

  if (inherits(model, "blrm")) {
    posterior_summary <- match.arg(posterior_summary)
    return(sops_blrm(
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
      by = by,
      include_re = include_re,
      id_var = id_var,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      newdata_orig = newdata_orig,
      refit_data = refit_data,
      newdata_supplied = newdata_supplied
    ))
  }

  # --- 2. Compute SOPs (Vectorized) ---
  sops_array <- soprob_markov(
    model = model,
    newdata = newdata,
    times = times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    gap_var = gap_var,
    time_covariates = time_covariates
  )

  # --- 3. Validate 'by' Parameter ---
  validate_sops_by(by, newdata)

  # --- 4. Tidy the Output ---
  result <- array_to_df_individual(
    sops_array,
    times,
    y_levels,
    newdata,
    by = by
  )

  # --- 6. Store Attributes for Downstream Use ---
  set_sops_attrs(
    result,
    class_name = "markov_sops",
    model = model,
    call_args = sops_call_args(
      times,
      y_levels,
      absorb,
      time_var,
      p_var,
      p2_var,
      gap_var,
      time_covariates,
      by = by
    ),
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    time_covariates = time_covariates,
    by = by,
    newdata_orig = newdata_orig,
    extra_attrs = list(
      newdata_pred = newdata,
      refit_data = refit_data,
      id_var = id_var,
      newdata_supplied = newdata_supplied
    )
  )
}

sops_blrm <- function(
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
  by,
  include_re,
  id_var,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  newdata_orig,
  refit_data,
  newdata_supplied
) {
  validate_sops_by(by, newdata)

  group_setup <- if (is.null(by)) NULL else posterior_group_setup(newdata, by)
  prediction_data <- if (is.null(group_setup)) {
    newdata
  } else {
    newdata[group_setup$rows, , drop = FALSE]
  }
  draw_indices <- select_posterior_draws(model, n_draws, seed)
  n_pat <- nrow(prediction_data)
  n_times <- length(times)
  n_states <- length(y_levels)
  n_cells <- n_pat * n_times * n_states
  n_requested <- length(draw_indices)
  n_output_cells <- if (is.null(group_setup)) {
    n_cells
  } else {
    group_setup$n_groups * n_times * n_states
  }
  # `sops()` keeps individual-level draw summaries in memory. The default
  # guard is 50 million numeric cells, about 400 MB before data-frame overhead.
  max_draw_cells <- getOption("markov.misc.max_sops_draw_cells", 5e7)

  if (n_output_cells * n_requested > max_draw_cells) {
    stop(
      "`sops()` for `blrm` would require ",
      format(n_output_cells * n_requested, scientific = FALSE),
      " posterior draw cells, above the current limit of ",
      format(max_draw_cells, scientific = FALSE),
      ". This limit protects memory because individual-level `sops()` stores ",
      "draw x patient x time x state values. Lower `n_draws`, use ",
      "`avg_sops()` for marginal summaries, or increase option ",
      "`markov.misc.max_sops_draw_cells` if this allocation is intentional."
    )
  }

  draw_values <- matrix(NA_real_, nrow = n_requested, ncol = n_output_cells)
  rownames(draw_values) <- draw_indices
  chunks <- split_draw_indices(
    draw_indices,
    chunk_size = getOption("markov.misc.blrm_chunk_size") %||%
      posterior_chunk_size(n_cells, n_requested)
  )
  gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  prediction_cache <- if (is.null(p2_var)) {
    cache <- new.env(parent = emptyenv())
    cache$cursor <- 0L
    cache$X <- list()
    cache$Z <- list()
    cache
  } else {
    NULL
  }
  row_cursor <- 1L
  for (chunk in chunks) {
    if (!is.null(prediction_cache)) {
      prediction_cache$cursor <- 0L
    }
    arr <- soprob_markov(
      model = model,
      newdata = prediction_data,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      include_re = include_re,
      id_var = id_var,
      n_draws = NULL,
      .draw_indices = chunk,
      .gamma_draws = subset_cached_draw_matrix(
        gamma_draws,
        draw_indices,
        chunk
      ),
      .prediction_cache = prediction_cache
    )
    if (is.null(group_setup)) {
      for (i in seq_along(chunk)) {
        draw_values[row_cursor, ] <- as.vector(arr[i, , , ])
        row_cursor <- row_cursor + 1L
      }
    } else {
      reduced <- reduce_sops_draw_array(
        arr,
        group_setup$groups,
        group_setup$n_groups
      )
      reduced <- aperm(reduced, c(1L, 3L, 4L, 2L))
      rows <- row_cursor:(row_cursor + length(chunk) - 1L)
      draw_values[rows, ] <- matrix(
        as.numeric(reduced),
        nrow = length(chunk)
      )
      row_cursor <- row_cursor + length(chunk)
    }
  }

  if (is.null(by)) {
    stats <- summarize_posterior_draw_matrix(
      draw_values,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    )
    idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
    result <- newdata[idx_pat, , drop = FALSE]
    result$time <- rep(rep(times, each = n_pat), times = n_states)
    result$state <- rep(y_levels, each = n_pat * n_times)
    result$estimate <- stats$estimate
    result$conf.low <- stats$conf.low
    result$conf.high <- stats$conf.high
    result$std.error <- stats$std.error
    rownames(result) <- NULL

    draws_df <- if (return_draws) {
      sops_draw_matrix_to_df(draw_values, result, draw_indices)
    } else {
      NULL
    }
  } else {
    cells <- posterior_group_cells(
      group_setup$metadata,
      times,
      y_levels
    )
    stats <- summarize_posterior_draw_matrix(
      draw_values,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    )
    result <- cbind(cells, stats)
    draws_df <- if (return_draws) {
      posterior_draw_matrix_to_df(draw_values, cells, draw_indices)
    } else {
      NULL
    }
  }

  set_sops_attrs(
    result,
    class_name = "markov_sops",
    model = model,
    call_args = sops_call_args(
      times,
      y_levels,
      absorb,
      time_var,
      p_var,
      p2_var,
      gap_var,
      time_covariates,
      by = by,
      include_re = include_re,
      id_var = id_var,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    ),
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    time_covariates = time_covariates,
    by = by,
    newdata_orig = newdata_orig,
    extra_attrs = list(
      method = "posterior",
      engine = "posterior",
      n_draws = n_requested,
      draw_ids = draw_indices,
      conf_level = conf_level,
      posterior_summary = posterior_summary,
      draws = draws_df,
      newdata_pred = newdata,
      refit_data = refit_data,
      id_var = id_var,
      newdata_supplied = newdata_supplied
    )
  )
}

posterior_group_setup <- function(data, group_cols) {
  rows <- which(stats::complete.cases(data[, group_cols, drop = FALSE]))
  if (!length(rows)) {
    stop(
      "No complete prediction rows remain after omitting missing grouping values.",
      call. = FALSE
    )
  }
  complete_data <- data[rows, , drop = FALSE]
  keys <- interaction(
    complete_data[, group_cols, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- as.integer(keys)
  n_groups <- nlevels(keys)
  first <- vapply(
    seq_len(n_groups),
    function(i) which(groups == i)[1L],
    integer(1)
  )
  metadata <- complete_data[first, group_cols, drop = FALSE]
  rownames(metadata) <- NULL
  list(
    rows = rows,
    groups = groups,
    n_groups = n_groups,
    metadata = metadata
  )
}

posterior_group_cells <- function(metadata, times, y_levels) {
  pieces <- vector("list", nrow(metadata))
  for (i in seq_len(nrow(metadata))) {
    piece <- expand.grid(
      time = times,
      state = y_levels,
      KEEP.OUT.ATTRS = FALSE
    )
    for (nm in names(metadata)) {
      piece[[nm]] <- metadata[[nm]][i]
    }
    pieces[[i]] <- piece
  }
  bind_rows_fill(pieces)
}

posterior_draw_matrix_to_df <- function(draw_values, cells, draw_indices) {
  n_draws <- nrow(draw_values)
  n_cells <- ncol(draw_values)
  out <- cells[rep(seq_len(n_cells), times = n_draws), , drop = FALSE]
  out$draw_id <- rep(draw_indices, each = n_cells)
  out$estimate <- as.vector(t(draw_values))
  rownames(out) <- NULL
  out
}

split_draw_indices <- function(draw_indices, chunk_size = NULL) {
  if (is.null(chunk_size)) {
    chunk_size <- getOption("markov.misc.blrm_chunk_size", 10L)
  }
  chunk_size <- max(1L, as.integer(chunk_size))
  split(draw_indices, ceiling(seq_along(draw_indices) / chunk_size))
}

summarize_posterior_draw_matrix <- function(
  draw_values,
  posterior_summary,
  conf_level
) {
  alpha <- 1 - conf_level
  estimate <- switch(
    posterior_summary,
    mean = colMeans(draw_values, na.rm = TRUE),
    median = apply(draw_values, 2, stats::median, na.rm = TRUE)
  )
  conf.low <- apply(
    draw_values,
    2,
    stats::quantile,
    probs = alpha / 2,
    na.rm = TRUE
  )
  conf.high <- apply(
    draw_values,
    2,
    stats::quantile,
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )
  std.error <- apply(draw_values, 2, stats::sd, na.rm = TRUE)
  data.frame(
    estimate = as.numeric(estimate),
    conf.low = as.numeric(conf.low),
    conf.high = as.numeric(conf.high),
    std.error = as.numeric(std.error)
  )
}

summarize_posterior_draws_df <- function(
  draws_df,
  group_cols,
  posterior_summary,
  conf_level
) {
  alpha <- 1 - conf_level
  split_key <- interaction(draws_df[, group_cols, drop = FALSE], drop = TRUE)
  groups <- split(seq_len(nrow(draws_df)), split_key)
  out <- draws_df[vapply(groups, `[`, integer(1), 1L), group_cols, drop = FALSE]
  values <- lapply(groups, function(idx) draws_df$estimate[idx])
  out$estimate <- vapply(
    values,
    function(x) {
      switch(
        posterior_summary,
        mean = mean(x, na.rm = TRUE),
        median = stats::median(x, na.rm = TRUE)
      )
    },
    numeric(1)
  )
  out$conf.low <- vapply(
    values,
    stats::quantile,
    numeric(1),
    probs = alpha / 2,
    na.rm = TRUE
  )
  out$conf.high <- vapply(
    values,
    stats::quantile,
    numeric(1),
    probs = 1 - alpha / 2,
    na.rm = TRUE
  )
  out$std.error <- vapply(values, stats::sd, numeric(1), na.rm = TRUE)
  rownames(out) <- NULL
  out
}

sops_draw_matrix_to_df <- function(draw_values, result, draw_indices) {
  meta_cols <- setdiff(
    names(result),
    c("estimate", "conf.low", "conf.high", "std.error")
  )
  n_draws <- nrow(draw_values)
  n_cells <- ncol(draw_values)
  out <- result[rep(seq_len(n_cells), times = n_draws), meta_cols, drop = FALSE]
  out$draw_id <- rep(draw_indices, each = n_cells)
  out$estimate <- as.vector(t(draw_values))
  rownames(out) <- NULL
  out
}


#' Calculate Averaged State Occupation Probabilities (Marginal Effects)
#'
#' Computes standardized (marginal) state occupation probabilities using
#' G-computation. Creates counterfactual cohorts by setting all individuals
#' to each level of the treatment variable and averaging over the covariate
#' distribution.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`, or `blrm`). For
#'   `vglm` models, the family must be cumulative-logit with `reverse = TRUE`.
#'   For `orm` models, the family must be logistic.
#' @param newdata Optional data frame of standardization profiles. When
#'   supplied, every row is treated as a separate fixed profile and the
#'   internal `rowid` column is regenerated. If `NULL`, uses data stored by
#'   [orm_markov()], [blrm_markov()], or [vglm_markov()] and requires exactly
#'   one complete designated starting profile per fitted patient. `refit_data`
#'   is never used as a prediction-profile fallback.
#' @param variables A named list specifying the variable(s) to standardize over.
#'   E.g., `list(tx = c(0, 1))` creates counterfactual datasets for treatment
#'   and control.
#' @param by Optional character vector of additional variables to group by
#'   after standardization. Rows missing any grouping value are omitted; an
#'   error is raised if none remain.
#' @param times Required visit-scale time points. Numeric time variables use
#'   numeric values; factor-valued visit indices use fitted visit levels.
#' @param y_levels A vector of state levels in fitted outcome order. Its length
#'   must agree with the fitted model's threshold count. If `NULL`, attempts to
#'   infer the levels from `model`.
#' @param absorb The absorbing state.
#' @param refit_data Optional full longitudinal data used only by refit-bootstrap
#'   inference. It is not used for point estimates. Defaults to data stored on
#'   wrapper-fitted models.
#' @param id_var Name of the patient ID variable. Required for bootstrap
#'   inference and for `blrm` random-effect prediction. If `NULL`, defaults to
#'   `"id"`; for wrapper-fitted models, it is inferred from stored metadata.
#'   For user-supplied `newdata`, `id_var` is required only for `blrm` random
#'   effects. Refit bootstrap inference uses `id_var` from `refit_data` or
#'   wrapper-stored model data. For `blrm` models with `include_re = TRUE`,
#'   `model$clusterInfo$name` is used before the final `"id"` fallback when
#'   wrapper metadata is absent.
#' @param time_var Name of the time variable in the model.
#' @param p_var Name of the previous state variable in the model.
#' @param p2_var Optional second previous-state variable. `NULL` uses a
#'   first-order Markov recursion; a non-`NULL` column name uses a second-order
#'   recursion.
#' @param gap_var Name of the time gap_var variable, if used.
#' @param time_covariates Optional time-varying covariate lookup table for
#'   explicit precomputed time-basis columns. Registered inline terms such as
#'   `rms::rcs(time, 4)` and `rms::lsp(time, c(3, 7))` can be used without it.
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws. Defaults to 100.
#' @param seed Optional random seed for reproducible random draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior SOP draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior SOP draws for
#'   extraction with `get_draws()`.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A data frame of class `markov_avg_sops` with columns:
#'   \item{time}{Time point}
#'   \item{state}{State level}
#'   \item{(variables)}{Value of standardization variable (e.g., tx)}
#'   \item{estimate}{Average probability across individuals}
#'   \item{conf.low, conf.high, std.error}{For `blrm` fits, posterior
#'     uncertainty summaries computed directly from marginalized SOP draws}
#'
#' @details
#' This function implements G-computation (standardization) for Markov SOPs:
#'
#' 1. **Counterfactual Creation**: For each value in `variables`, creates a
#'    copy of `newdata` with that variable set to the specified value.
#'
#' 2. **Individual Prediction**: Calls `sops()` to compute individual-level
#'    SOPs for each patient under each counterfactual scenario.
#'
#' 3. **Marginalization**: Averages individual SOPs across patients within
#'    each time-state-treatment combination, yielding population-average
#'    (marginal) probabilities.
#'
#' The result represents the expected SOP if the entire population received
#' treatment vs. control, averaged over the observed covariate distribution.
#' For automatic wrapper-stored standardization, that distribution contains
#' exactly the fitted patients with a complete designated starting profile and
#' at least one usable likelihood transition. A missing response at the
#' designated row alone does not exclude a patient who contributes a later
#' transition. Patients with no usable fitted transition are excluded, and a
#' fitted patient without a complete designated profile is an error.
#' For `rmsb::blrm()` models, patient-level posterior arrays are chunked and
#' marginalized by draw before summarizing to reduce memory use.
#' State-wise medians and interval bounds are not constrained to sum to one
#' across states; use `posterior_summary = "mean"` when the displayed estimates
#' themselves need to preserve total probability. Draw-level probabilities
#' stored with `return_draws = TRUE` remain normalized within each draw.
#'
#' @seealso [sops()] for individual-level SOPs, [inferences()] for bootstrap
#'   uncertainty and [avg_comparisons()] for specialized estimands.
#'
#' @examples
#' \dontrun{
#' fit <- vglm_markov(
#'   ordered(y) ~ rms::rcs(time, 4) + tx + yprev,
#'   family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
#'   data = data,
#'   id_var = "id"
#' )
#'
#' # Wrapper-fitted models can use their stored designated starting profiles.
#' result <- avg_sops(
#'   model = fit,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   y_levels = 1:6,
#'   absorb = 6
#' ) |> inferences(method = "mvn", n_draws = 500)
#'
#' # Refit bootstrap can reuse the stored full longitudinal data.
#' result_boot <- avg_sops(
#'   model = fit,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   y_levels = 1:6,
#'   absorb = 6
#' ) |> inferences(method = "fwb", n_draws = 500)
#' }
#'
#' @export
avg_sops <- function(
  model,
  newdata = NULL,
  variables = NULL,
  by = NULL,
  times,
  y_levels = NULL,
  absorb = NULL,
  refit_data = NULL,
  id_var = NULL,
  time_var = "time",
  p_var = "yprev",
  p2_var = NULL,
  gap_var = NULL,
  time_covariates = NULL,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  # --- 1. Input Validation ---
  # Validate model compatibility
  validate_markov_model(model)
  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied to `avg_sops()`.")
  }
  conf_level <- validate_conf_level(conf_level)

  if (is.null(variables)) {
    stop(
      "`variables` is required for G-computation. ",
      "Specify the treatment variable, e.g., `variables = list(tx = c(0, 1))`."
    )
  }

  data_res <- resolve_markov_source_data(
    model,
    newdata,
    refit_data,
    time_var = time_var,
    p_var = p_var
  )
  newdata_orig <- data_res$source_data
  refit_data <- data_res$refit_data
  newdata_supplied <- data_res$newdata_supplied

  if (inherits(model, "blrm") && isTRUE(include_re)) {
    id_var <- resolve_blrm_id_var(model, newdata_orig, id_var)
    validate_markov_id_var(id_var, newdata_orig, "newdata")
  } else {
    id_var <- markov_model_id_var(model, id_var) %||% "id"
  }

  if (!is.null(refit_data)) {
    validate_markov_id_var(id_var, refit_data, "refit_data")
  }

  if (!newdata_supplied) {
    validate_markov_id_var(id_var, newdata_orig, "stored model data")
  }

  # Validate variables exist in data
  missing_vars <- setdiff(names(variables), names(newdata_orig))
  if (length(missing_vars) > 0) {
    stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }

  # --- 2. Extract Baseline Data (One Row Per Patient) ---
  # For standardization, we need unique patient baseline covariates
  # Average the counterfactual SOP predictions over the target population.
  baseline_data <- if (newdata_supplied) {
    newdata_orig
  } else {
    resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      time_var = time_var
    )
  }
  baseline_data <- ensure_markov_rowid(baseline_data)

  # --- 3. Create Counterfactual Datasets ---
  # For each combination in variables, create a copy of baseline_data with
  # the variable(s) set to that value

  if (!is.list(variables)) {
    var_list <- list()
    for (i in seq_along(variables)) {
      var_list[[variables[i]]] <- unique(baseline_data[[variables[i]]])
    }
  } else {
    var_list <- variables
  }

  grid <- do.call(expand.grid, var_list)

  newdata_expanded <- create_counterfactual_data(baseline_data, grid, var_list)

  direct_backend <- inherits(model, c("vglm", "robcov_vglm", "orm", "blrm"))
  if (!direct_backend) {
    sops_ind <- sops(
      model,
      newdata = newdata_expanded,
      refit_data = refit_data,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      id_var = id_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      ...
    )
    resolved_times <- attr(sops_ind, "call_args")$times
    group_cols <- unique(c("time", "state", names(var_list), by))
    missing_groups <- setdiff(group_cols, names(sops_ind))
    if (length(missing_groups) > 0L) {
      stop(
        "Grouping variables missing: ",
        paste(missing_groups, collapse = ", ")
      )
    }
    result <- aggregate_sops_estimates(sops_ind, group_cols)
    return(set_sops_attrs(
      result,
      class_name = "markov_avg_sops",
      model = attr(sops_ind, "model"),
      call_args = attr(sops_ind, "call_args"),
      time_var = attr(sops_ind, "time_var"),
      p_var = attr(sops_ind, "p_var"),
      p2_var = attr(sops_ind, "p2_var"),
      y_levels = attr(sops_ind, "y_levels"),
      absorb = attr(sops_ind, "absorb"),
      gap_var = attr(sops_ind, "gap_var"),
      time_covariates = attr(sops_ind, "time_covariates"),
      newdata_orig = newdata_orig,
      avg_args = list(
        variables = var_list,
        by = by,
        times = resolved_times,
        id_var = id_var
      ),
      extra_attrs = list(
        newdata_pred = newdata_expanded,
        refit_data = refit_data,
        id_var = id_var,
        newdata_supplied = newdata_supplied
      )
    ))
  }

  if (is.null(y_levels)) {
    y_levels <- if (inherits(model, "robcov_vglm")) {
      model$extra$colnames.y
    } else if (inherits(model, "vglm")) {
      model@extra$colnames.y
    } else if (inherits(model, c("orm", "blrm"))) {
      model$yunique %||% model$ylevels
    } else {
      NULL
    }
    if (is.null(y_levels)) stop("`y_levels` cannot be NULL")
  }

  if (inherits(model, "blrm")) {
    posterior_summary <- match.arg(posterior_summary)
    result <- avg_sops_blrm(
      model = model,
      newdata_expanded = newdata_expanded,
      grid = grid,
      variables = var_list,
      by = by,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      id_var = id_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      ...
    )
    attr(result, "newdata_orig") <- newdata_orig
    attr(result, "newdata_pred") <- newdata_expanded
    attr(result, "refit_data") <- refit_data
    attr(result, "id_var") <- id_var
    attr(result, "newdata_supplied") <- newdata_supplied
    return(result)
  }

  # --- 4. Compute and marginalize without materializing individual rows ---
  time_res <- resolve_sop_times(
    model,
    newdata_expanded,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "unique"
  )
  resolved_times <- time_res$times
  validate_factor_gap(gap_var, time_covariates, time_res$time_info)
  sops_array <- soprob_markov(
    model = model,
    newdata = newdata_expanded,
    times = resolved_times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    gap_var = gap_var,
    time_covariates = time_covariates,
    ...
  )
  result <- marginalize_sops_array(
    sops_array = sops_array,
    grid = grid,
    times = resolved_times,
    y_levels = y_levels,
    variables = var_list,
    n_cf = nrow(grid),
    n_each = nrow(baseline_data),
    by = by,
    newdata = newdata_expanded
  )

  set_sops_attrs(
    result,
    class_name = "markov_avg_sops",
    model = model,
    call_args = sops_call_args(
      resolved_times,
      y_levels,
      absorb,
      time_var,
      p_var,
      p2_var,
      gap_var,
      time_covariates
    ),
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    time_covariates = time_covariates,
    newdata_orig = newdata_orig,
    avg_args = list(
      variables = var_list,
      by = by,
      times = resolved_times,
      id_var = id_var
    ),
    extra_attrs = list(
      newdata_pred = newdata_expanded,
      refit_data = refit_data,
      id_var = id_var,
      newdata_supplied = newdata_supplied
    )
  )
}

avg_sops_blrm <- function(
  model,
  newdata_expanded,
  grid,
  variables,
  by,
  times,
  y_levels,
  absorb,
  time_var,
  p_var,
  id_var,
  p2_var,
  gap_var,
  time_covariates,
  include_re,
  n_draws,
  seed,
  posterior_summary,
  conf_level,
  return_draws,
  ...
) {
  y_levels <- y_levels %||% model$ylevels

  time_res <- resolve_sop_times(
    model,
    newdata_expanded,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "unique"
  )
  times <- time_res$times
  validate_factor_gap(gap_var, time_covariates, time_res$time_info)

  validate_sops_by(by, newdata_expanded, data_arg = "expanded prediction data")

  n_cf <- nrow(grid)
  n_each <- nrow(newdata_expanded) / n_cf
  group_cols <- unique(c(names(variables), by))
  group_setup <- if (is.null(by)) {
    list(
      rows = seq_len(nrow(newdata_expanded)),
      groups = rep(seq_len(n_cf), each = n_each),
      n_groups = n_cf,
      metadata = grid[, names(variables), drop = FALSE]
    )
  } else {
    posterior_group_setup(newdata_expanded, group_cols)
  }
  prediction_data <- newdata_expanded[group_setup$rows, , drop = FALSE]
  draw_indices <- select_posterior_draws(model, n_draws, seed)
  chunks <- split_draw_indices(
    draw_indices,
    chunk_size = getOption("markov.misc.blrm_avg_chunk_size") %||%
      posterior_chunk_size(
        nrow(prediction_data) * length(times) * length(y_levels),
        length(draw_indices)
      )
  )
  gamma_draws <- cache_blrm_gamma_draws(model, draw_indices, include_re)
  n_cells <- group_setup$n_groups * length(times) * length(y_levels)
  draw_values <- matrix(
    NA_real_,
    nrow = length(draw_indices),
    ncol = n_cells,
    dimnames = list(draw_indices, NULL)
  )
  prediction_cache <- if (is.null(p2_var)) {
    cache <- new.env(parent = emptyenv())
    cache$cursor <- 0L
    cache$X <- list()
    cache$Z <- list()
    cache
  } else {
    NULL
  }
  row_cursor <- 1L

  for (chunk in chunks) {
    if (!is.null(prediction_cache)) {
      prediction_cache$cursor <- 0L
    }
    arr <- soprob_markov(
      model = model,
      newdata = prediction_data,
      times = times,
      y_levels = y_levels,
      absorb = absorb,
      time_var = time_var,
      p_var = p_var,
      p2_var = p2_var,
      gap_var = gap_var,
      time_covariates = time_covariates,
      include_re = include_re,
      id_var = id_var,
      n_draws = NULL,
      .draw_indices = chunk,
      .gamma_draws = subset_cached_draw_matrix(
        gamma_draws,
        draw_indices,
        chunk
      ),
      .prediction_cache = prediction_cache
    )
    reduced <- reduce_sops_draw_array(
      arr,
      group_setup$groups,
      group_setup$n_groups
    )
    reduced <- aperm(reduced, c(1L, 3L, 4L, 2L))
    rows <- row_cursor:(row_cursor + length(chunk) - 1L)
    draw_values[rows, ] <- matrix(as.numeric(reduced), nrow = length(chunk))
    row_cursor <- row_cursor + length(chunk)
  }

  cells <- posterior_group_cells(group_setup$metadata, times, y_levels)
  stats <- summarize_posterior_draw_matrix(
    draw_values,
    posterior_summary = posterior_summary,
    conf_level = conf_level
  )
  result <- cbind(cells, stats)
  draws_df <- if (return_draws) {
    posterior_draw_matrix_to_df(draw_values, cells, draw_indices)
  } else {
    NULL
  }

  set_sops_attrs(
    result,
    class_name = "markov_avg_sops",
    model = model,
    call_args = sops_call_args(
      times,
      y_levels,
      absorb,
      time_var,
      p_var,
      p2_var,
      gap_var,
      time_covariates,
      include_re = include_re,
      id_var = id_var,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level
    ),
    time_var = time_var,
    p_var = p_var,
    p2_var = p2_var,
    y_levels = y_levels,
    absorb = absorb,
    gap_var = gap_var,
    time_covariates = time_covariates,
    by = by,
    avg_args = list(
      variables = variables,
      by = by,
      times = times,
      id_var = id_var
    ),
    extra_attrs = list(
      method = "posterior",
      engine = "posterior",
      n_draws = length(draw_indices),
      draw_ids = draw_indices,
      conf_level = conf_level,
      posterior_summary = posterior_summary,
      draws = draws_df
    )
  )
}

posterior_chunk_size <- function(cells_per_draw, n_draws) {
  budget <- getOption("markov.misc.posterior_working_bytes", 256 * 1024^2)
  bytes_per_draw <- max(1, cells_per_draw * 8 * 3)
  max(1L, min(as.integer(n_draws), floor(budget / bytes_per_draw)))
}
