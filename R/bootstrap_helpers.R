#' Fast Group Bootstrap ID Sampler
#'
#' Generates bootstrap samples by resampling group IDs (e.g., patient IDs) with
#' replacement. Returns only the sampled IDs, not the full data, for memory
#' efficiency. This is much faster and more memory-efficient than
#' rsample::group_bootstraps() for datasets with many groups.
#'
#' @param data A data frame containing the patient data
#' @param id_var Character string specifying the name of the ID variable for
#'   group bootstrap (default "id")
#' @param n_boot Number of bootstrap samples to generate
#'
#' @return A list of length \code{n_boot}, where each element is a data frame
#'   (ID lookup table) with columns:
#'   \itemize{
#'     \item original_id: The sampled group ID from the original data
#'     \item new_id: Unique identifier for resampled groups (e.g., "1_1", "1_2")
#'     \item boot_id: Bootstrap iteration number
#'   }
#'
#' @details
#' The ID lookup tables are later joined with the original data on-demand using
#' \code{\link{materialize_bootstrap_sample}}, so each parallel worker only
#' materializes the data it needs for analysis.
#'
#' The new_id column ensures that if a group is sampled multiple times in
#' the same bootstrap iteration, each instance gets a unique identifier
#' (e.g., "1_1", "1_2"), which is necessary for refitting models that expect
#' unique group identifiers.
#'
#' This approach is dramatically faster than rsample::group_bootstraps() for
#' datasets with many groups or unbalanced group sizes, because rsample uses
#' an inefficient "oversampling then trimming" strategy. See GitHub issue
#' tidymodels/rsample#357 for details.
#'
#' @examples
#' \dontrun{
#' # Generate 2000 bootstrap ID samples
#' boot_ids <- fast_group_bootstrap(my_data, id_var = "id", n_boot = 2000)
#'
#' # Each element contains only the sampled IDs (very small)
#' boot_ids[[1]]
#'
#' # To get a full bootstrap sample, materialize it:
#' boot_sample <- materialize_bootstrap_sample(boot_ids[[1]], my_data, "id")
#'
#' # Use with parallel processing - data is materialized on each worker
#' library(furrr)
#' plan(callr)
#' results <- future_map(boot_ids, function(ids) {
#'   boot_data <- materialize_bootstrap_sample(ids, my_data, "id")
#'   analyze_bootstrap(boot_data)
#' })
#' }
#'
#' @importFrom stats ave
#' @export

fast_group_bootstrap <- function(data, id_var = "id", n_boot) {
  # Input validation
  if (!id_var %in% names(data)) {
    stop("id_var '", id_var, "' not found in data")
  }

  # Get unique group IDs
  unique_ids <- unique(data[[id_var]])
  n_groups <- length(unique_ids)

  # Generate bootstrap ID samples (IDs only, not full data)
  bootstrap_id_samples <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    # Sample group IDs with replacement
    sampled_ids <- sample(unique_ids, size = n_groups, replace = TRUE)

    # Create a lookup table with the sampled IDs
    # Keep track of how many times each ID appears
    id_table <- data.frame(
      original_id = sampled_ids
    )

    # Add occurrence count for duplicate IDs
    id_table$occurrence <- ave(
      id_table$original_id,
      id_table$original_id,
      FUN = seq_along
    )

    # Create unique new_id combining original ID and occurrence
    id_table$new_id <- paste(
      id_table$original_id,
      id_table$occurrence,
      sep = "_"
    )

    # Add bootstrap iteration ID
    id_table$boot_id <- i

    # Remove temporary column
    id_table$occurrence <- NULL

    bootstrap_id_samples[[i]] <- id_table
  }

  return(bootstrap_id_samples)
}


#' Materialize Bootstrap Sample from ID Lookup
#'
#' Joins bootstrap ID lookup table with the original data to create a full
#' bootstrap dataset. This is used internally by bootstrap functions to
#' materialize data on-demand for memory efficiency.
#'
#' @param boot_ids A data frame with columns: original_id, new_id, boot_id
#'   (output from \code{\link{fast_group_bootstrap}})
#' @param data The original data frame
#' @param id_var Name of the grouping variable in data (e.g., "id")
#'
#' @return A data frame containing the bootstrap sample with columns from
#'   the original data plus new_id and boot_id
#'
#' @details
#' This function performs a left join between the bootstrap ID lookup table
#' and the original data, preserving the sampling order and creating unique
#' patient identifiers for groups that were sampled multiple times.
#'
#' The join uses relationship = "many-to-many" because:
#' \itemize{
#'   \item Each ID in boot_ids can match multiple rows in data (longitudinal)
#'   \item Each ID can appear multiple times in boot_ids (bootstrap resampling)
#' }
#'
#' @examples
#' \dontrun{
#' # Generate bootstrap IDs
#' boot_ids <- fast_group_bootstrap(my_data, "id", n_boot = 10)
#'
#' # Materialize first bootstrap sample
#' boot_sample_1 <- materialize_bootstrap_sample(boot_ids[[1]], my_data, "id")
#'
#' # This is done automatically by apply_to_bootstrap
#' }
#'
#' @export

materialize_bootstrap_sample <- function(boot_ids, data, id_var) {
  # Rename original_id to match id_var for joining
  boot_ids_renamed <- boot_ids
  names(boot_ids_renamed)[names(boot_ids_renamed) == "original_id"] <- id_var

  # Join with original data, preserving bootstrap sampling order.
  boot_sample <- left_join_preserve_order(boot_ids_renamed, data, by = id_var)

  return(boot_sample)
}

generate_fwb_bootstrap_weights <- function(
  data,
  id_var = "id",
  n_boot,
  weight_dist = "exponential"
) {
  if (!id_var %in% names(data)) {
    stop("id_var '", id_var, "' not found in data")
  }

  weight_dist <- match.arg(weight_dist, choices = "exponential")
  cluster_ids <- unique(data[[id_var]])
  n_clusters <- length(cluster_ids)
  fwb_samples <- vector("list", n_boot)

  for (i in seq_len(n_boot)) {
    raw_weights <- stats::rexp(n_clusters, rate = 1)
    weight_mean <- mean(raw_weights)
    if (!is.finite(weight_mean) || weight_mean <= 0) {
      stop("Invalid fractional bootstrap weights generated.")
    }

    fwb_samples[[i]] <- data.frame(
      original_id = cluster_ids,
      fwb_weight = raw_weights / weight_mean,
      boot_id = i
    )
  }

  fwb_samples
}

materialize_fwb_bootstrap_sample <- function(fwb_weights, data, id_var) {
  if (!all(c("original_id", "fwb_weight", "boot_id") %in% names(fwb_weights))) {
    stop(
      "`fwb_weights` must contain original_id, fwb_weight, and boot_id columns."
    )
  }

  idx <- match(
    as.character(data[[id_var]]),
    as.character(fwb_weights$original_id)
  )
  if (anyNA(idx)) {
    missing_ids <- unique(as.character(data[[id_var]])[is.na(idx)])
    stop(
      "Some data IDs were not found in the fractional bootstrap weights: ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) " ..." else "",
      "."
    )
  }

  boot_data <- data
  boot_data$fwb_weight <- fwb_weights$fwb_weight[idx]
  boot_data$boot_id <- fwb_weights$boot_id[idx]
  boot_data
}

fwb_baseline_weights <- function(fwb_weights, baseline_data, id_var) {
  if (!id_var %in% names(baseline_data)) {
    stop("ID variable '", id_var, "' not found in baseline data.")
  }

  idx <- match(
    as.character(baseline_data[[id_var]]),
    as.character(fwb_weights$original_id)
  )
  if (anyNA(idx)) {
    missing_ids <- unique(as.character(baseline_data[[id_var]])[is.na(idx)])
    stop(
      "Some baseline IDs were not found in the fractional bootstrap weights: ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) " ..." else "",
      "."
    )
  }

  fwb_weights$fwb_weight[idx]
}


#' Apply Function to Bootstrap Samples with Parallelization
#'
#' Helper function to apply an analysis function to bootstrap samples with
#' optional parallel processing. Handles just-in-time materialization of
#' bootstrap data for memory efficiency.
#'
#' @param boot_samples List of bootstrap ID lookups from
#'   \code{\link{fast_group_bootstrap}}. Each element should be a data frame
#'   with columns: original_id, new_id, boot_id.
#' @param analysis_fn Function to apply to each bootstrap sample. Should accept
#'   a data frame and return a result (can be any type).
#' @param data The original data frame (before bootstrap sampling). This will
#'   be joined with boot_samples on each worker to materialize the data.
#' @param id_var Name of the grouping variable in data (e.g., "id")
#' @param workers Number of workers for parallel processing. If NULL (default)
#'   or 1, sequential processing is used. If > 1, uses parallel processing with
#'   future.callr::callr strategy.
#' @param parallel (Deprecated) Logical indicating whether to use parallel
#'   processing. Use the \code{workers} parameter instead.
#' @param packages Character vector of packages needed by \code{analysis_fn}.
#'   Default is c("rms", "VGAM", "Hmisc", "stats", "dplyr")
#' @param globals Character vector of global variables/functions needed by
#'   \code{analysis_fn}. These will be exported to each worker.
#'
#' @return A list of length \code{length(boot_samples)} containing the results
#'   from applying \code{analysis_fn} to each bootstrap sample.
#'
#' @details
#' This function implements a memory-efficient just-in-time bootstrap approach:
#' \enumerate{
#'   \item Bootstrap ID lookups are distributed to workers (minimal memory)
#'   \item Each worker receives a copy of the original data
#'   \item Worker materializes bootstrap data by joining IDs with original data
#'   \item Worker runs analysis and returns only the result
#'   \item Full bootstrap datasets are never stored in memory simultaneously
#' }
#'
#' If \code{workers} is NULL or 1, uses sequential processing with \code{lapply()}.
#' If \code{workers > 1}, uses future.callr::callr strategy for parallel
#' processing, which provides isolated R processes for each worker. This makes
#' it safe to modify the global environment (e.g., for datadist).
#'
#' @examples
#' \dontrun{
#' # Define analysis function (receives full bootstrap data)
#' analyze_boot <- function(boot_data) {
#'   # Update datadist
#'   dd <- rms::datadist(boot_data)
#'   assign("dd", dd, envir = .GlobalEnv)
#'   options(datadist = "dd")
#'
#'   # Refit model
#'   m <- update(original_model, data = boot_data)
#'
#'   # Return coefficient
#'   coef(m)["tx"]
#' }
#'
#' # Generate bootstrap IDs (not full data)
#' boot_ids <- fast_group_bootstrap(my_data, id_var = "id", n_boot = 2000)
#'
#' # Apply to bootstrap samples with just-in-time materialization
#' results <- apply_to_bootstrap(
#'   boot_samples = boot_ids,
#'   analysis_fn = analyze_boot,
#'   data = my_data,
#'   id_var = "id",
#'   workers = 8,
#'   globals = c("original_model")
#' )
#' }
#'
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan
#' @importFrom future.callr callr
#' @export

apply_to_bootstrap <- function(
  boot_samples,
  analysis_fn,
  data,
  id_var,
  workers = NULL,
  parallel = NULL,
  packages = c("rms", "VGAM", "Hmisc", "stats"),
  globals = character(0)
) {
  # Handle deprecated parallel parameter
  if (!is.null(parallel)) {
    warning(
      "The 'parallel' argument is deprecated. ",
      "Please use 'workers' instead.\n",
      "  - For sequential processing: workers = NULL or workers = 1\n",
      "  - For parallel processing: workers = N (e.g., workers = 8)"
    )
    if (parallel && is.null(workers)) {
      workers <- parallel::detectCores() - 1
    }
  }

  # Determine whether to parallelize based on workers
  use_parallel <- !is.null(workers) && workers > 1

  # Wrapper function that materializes data then runs analysis
  wrapper_fn <- function(boot_ids) {
    # Materialize bootstrap sample from ID lookup
    boot_data <- materialize_bootstrap_sample(boot_ids, data, id_var)

    # Run user's analysis function
    analysis_fn(boot_data)
  }

  # Apply function to bootstrap samples
  if (use_parallel) {
    old_plan <- plan()
    on.exit(plan(old_plan), add = TRUE)
    plan(callr, workers = workers)

    results <- furrr::future_map(
      boot_samples,
      wrapper_fn,
      .options = furrr::furrr_options(
        packages = c(packages),
        globals = c(globals, "data", "id_var", "analysis_fn")
      )
    )
  } else {
    results <- lapply(boot_samples, wrapper_fn)
  }

  return(results)
}

apply_to_fwb_bootstrap <- function(
  fwb_samples,
  analysis_fn,
  data,
  id_var,
  workers = NULL,
  packages = c("rms", "VGAM", "Hmisc", "stats"),
  globals = character(0)
) {
  use_parallel <- !is.null(workers) && workers > 1
  materialize_fwb <- materialize_fwb_bootstrap_sample

  wrapper_fn <- function(fwb_weights) {
    boot_data <- materialize_fwb(fwb_weights, data, id_var)
    analysis_fn(boot_data, fwb_weights)
  }

  if (use_parallel) {
    old_plan <- plan()
    on.exit(plan(old_plan), add = TRUE)
    plan(callr, workers = workers)

    results <- furrr::future_map(
      fwb_samples,
      wrapper_fn,
      .options = furrr::furrr_options(
        packages = c(packages),
        globals = c(globals, "data", "id_var", "analysis_fn", "materialize_fwb")
      )
    )
  } else {
    results <- lapply(fwb_samples, wrapper_fn)
  }

  results
}


#' Bootstrap Analysis Wrapper
#'
#' Common bootstrap analysis steps used across multiple bootstrap functions.
#' Handles releveling factors, updating datadist, and refitting models.
#'
#' @param boot_data Bootstrap sample data frame
#' @param model Original fitted model to update with bootstrap data
#' @param factor_cols Character vector of factor columns to relevel (e.g.,
#'   c("y", "yprev"))
#' @param original_data Original data frame before bootstrap sampling
#' @param y_levels Original state levels (e.g., 1:6)
#' @param absorb Absorbing state in original levels
#' @param update_datadist Logical indicating whether to update datadist
#'   (only needed for rms::orm models)
#' @param use_coefstart Logical indicating whether to use starting coefficients
#'   from the original model when refitting (only for vglm models).
#' @param fit_weights Optional non-negative row weights for the bootstrap fit.
#'
#' @return A list with components:
#'   \itemize{
#'     \item model: Refitted model on bootstrap data (or NULL if fitting failed)
#'     \item data: Bootstrap data with releveled factors
#'     \item y_levels: Updated y_levels (or NULL if not provided)
#'     \item absorb: Updated absorb (or NULL if not provided)
#'     \item missing_states: Character vector of missing state levels
#'   }
#'
#' @details
#' This function wraps common bootstrap operations:
#' 1. Relevel factors to consecutive integers if states are missing
#' 2. Update datadist for rms models
#' 3. Refit the model on the bootstrap data
#'
#' @examples
#' original_data <- data.frame(y = 1:5, x = 1:5)
#' fit <- lm(y ~ x, data = original_data)
#' boot_data <- original_data[c(1, 1, 2, 3, 4), ]
#'
#' result <- bootstrap_analysis_wrapper(
#'   boot_data = boot_data,
#'   model = fit,
#'   factor_cols = character(0),
#'   original_data = original_data,
#'   update_datadist = FALSE
#' )
#' names(result)
#'
#' @seealso [relevel_factors_consecutive()], [fast_group_bootstrap()],
#'   [apply_to_bootstrap()]
#'
#' @export
bootstrap_analysis_wrapper <- function(
  boot_data,
  model,
  factor_cols,
  original_data,
  y_levels = NULL,
  absorb = NULL,
  update_datadist = TRUE,
  use_coefstart = FALSE,
  fit_weights = NULL
) {
  fit_model <- bootstrap_refit_model(model)

  if (!is.null(fit_weights)) {
    if (
      !is.numeric(fit_weights) ||
        length(fit_weights) != nrow(boot_data) ||
        any(!is.finite(fit_weights)) ||
        any(fit_weights < 0)
    ) {
      stop("`fit_weights` must be non-negative finite row weights.")
    }

    original_weights <- bootstrap_original_fit_weights(fit_model, boot_data)
    boot_data$.markov_misc_fit_weight <- original_weights * fit_weights
  }

  # Relevel factors to handle missing states
  releveled <- relevel_factors_consecutive(
    data = boot_data,
    factor_cols = factor_cols,
    original_data = original_data,
    y_levels = y_levels,
    absorb = absorb
  )

  boot_data <- releveled$data
  boot_ylevels <- releveled$ylevels
  boot_absorb <- releveled$absorb
  missing_states <- releveled$missing_levels

  # Update datadist if needed (only for orm models)
  if (update_datadist && inherits(fit_model, "orm")) {
    dd_env <- globalenv()
    old_datadist <- getOption("datadist")
    old_dd_exists <- exists("dd", envir = dd_env, inherits = FALSE)
    old_dd <- if (old_dd_exists) get("dd", envir = dd_env) else NULL
    on.exit(
      {
        options(datadist = old_datadist)
        if (old_dd_exists) {
          assign("dd", old_dd, envir = dd_env)
        } else if (exists("dd", envir = dd_env, inherits = FALSE)) {
          remove(list = "dd", envir = dd_env)
        }
      },
      add = TRUE
    )

    dd <- rms::datadist(bootstrap_datadist_data(boot_data))
    assign("dd", dd, envir = dd_env)
    options(datadist = "dd")
  }

  # Refit model
  m_boot <- tryCatch(
    {
      # Attempt with coefstart if conditions met
      if (
        use_coefstart &&
          inherits(fit_model, "vglm") &&
          length(missing_states) == 0
      ) {
        tryCatch(
          {
            update_bootstrap_model(
              fit_model,
              boot_data,
              fit_weights = fit_weights,
              coefstart = stats::coef(fit_model)
            )
          },
          error = function(e) {
            # If coefstart fails (e.g. non-conformable due to dropped predictor levels),
            # fall back to standard update
            update_bootstrap_model(
              fit_model,
              boot_data,
              fit_weights = fit_weights
            )
          }
        )
      } else {
        # Standard update without coefstart
        update_bootstrap_model(
          fit_model,
          boot_data,
          fit_weights = fit_weights
        )
      }
    },
    error = function(e) {
      warning("Bootstrap model fitting failed: ", e$message)
      return(NULL)
    }
  )

  return(list(
    model = m_boot,
    data = boot_data,
    ylevels = boot_ylevels,
    absorb = boot_absorb,
    missing_states = missing_states
  ))
}

bootstrap_datadist_data <- function(data) {
  metadata_cols <- c(
    "boot_id",
    "new_id",
    "fwb_weight",
    ".markov_misc_fit_weight"
  )
  data[, setdiff(names(data), metadata_cols), drop = FALSE]
}

bootstrap_refit_model <- function(model) {
  if (inherits(model, "robcov_vglm")) {
    if (is.null(model$vglm_fit)) {
      stop(
        "The supplied robcov_vglm object does not contain the original vglm fit."
      )
    }
    return(model$vglm_fit)
  }

  model
}

bootstrap_original_fit_weights <- function(model, data) {
  stored_weights <- bootstrap_stored_fit_weights(model, data)
  if (!is.null(stored_weights)) {
    return(stored_weights)
  }

  model_call <- if (inherits(model, "vglm")) {
    tryCatch(methods::slot(model, "call"), error = function(e) NULL)
  } else {
    model$call
  }

  if (is.null(model_call) || is.null(model_call$weights)) {
    return(rep(1, nrow(data)))
  }

  weights <- tryCatch(
    eval(model_call$weights, data, parent.frame()),
    error = function(e) NULL
  )

  if (is.null(weights)) {
    warning(
      "Could not recover original model weights; using fractional ",
      "bootstrap weights only.",
      call. = FALSE
    )
    return(rep(1, nrow(data)))
  }

  if (
    length(weights) != nrow(data) ||
      any(!is.finite(weights)) ||
      any(weights < 0)
  ) {
    stop("Original model weights must be non-negative finite row weights.")
  }

  as.numeric(weights)
}

bootstrap_stored_fit_weights <- function(model, data) {
  candidates <- list(
    tryCatch(
      stats::model.weights(stats::model.frame(model)),
      error = function(e) NULL
    ),
    bootstrap_model_component(model, "prior.weights"),
    bootstrap_model_component(model, "weights"),
    if (inherits(model, "vglm")) {
      tryCatch(methods::slot(model, "prior.weights"), error = function(e) NULL)
    } else {
      NULL
    }
  )

  for (weights in candidates) {
    weights <- bootstrap_as_row_weights(weights, nrow(data))
    if (!is.null(weights)) {
      return(weights)
    }
  }

  NULL
}

bootstrap_model_component <- function(model, name) {
  if (methods::is(model, "S4")) {
    return(NULL)
  }

  tryCatch(model[[name]], error = function(e) NULL)
}

bootstrap_as_row_weights <- function(weights, n_rows) {
  if (is.null(weights)) {
    return(NULL)
  }

  if (is.matrix(weights) || is.data.frame(weights)) {
    if (nrow(weights) != n_rows) {
      return(NULL)
    }
    if (ncol(weights) == 1) {
      weights <- weights[, 1]
    } else if (all(weights == weights[, 1])) {
      weights <- weights[, 1]
    } else {
      return(NULL)
    }
  }

  if (!is.numeric(weights) || length(weights) != n_rows) {
    return(NULL)
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("Original model weights must be non-negative finite row weights.")
  }

  as.numeric(weights)
}

update_bootstrap_model <- function(
  model,
  boot_data,
  fit_weights = NULL,
  coefstart = NULL
) {
  has_fit_weights <- !is.null(fit_weights)

  if (has_fit_weights && !is.null(coefstart)) {
    return(suppress_orm_bootstrap_weight_warning(
      model,
      stats::update(
        model,
        data = boot_data,
        weights = .markov_misc_fit_weight,
        coefstart = coefstart
      )
    ))
  }
  if (has_fit_weights) {
    return(suppress_orm_bootstrap_weight_warning(
      model,
      stats::update(
        model,
        data = boot_data,
        weights = .markov_misc_fit_weight
      )
    ))
  }
  if (!is.null(coefstart)) {
    return(suppress_orm_bootstrap_weight_warning(
      model,
      stats::update(
        model,
        data = boot_data,
        coefstart = coefstart
      )
    ))
  }

  suppress_orm_bootstrap_weight_warning(
    model,
    stats::update(model, data = boot_data)
  )
}

suppress_orm_bootstrap_weight_warning <- function(model, expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (
        inherits(model, "vglm") &&
          grepl(
            "`id_var` was not supplied to vglm_markov()",
            conditionMessage(w),
            fixed = TRUE
          )
      ) {
        invokeRestart("muffleWarning")
      }
      if (
        inherits(model, "orm") &&
          grepl(
            "currently weights are ignored in model validation and bootstrapping orm fits",
            conditionMessage(w),
            fixed = TRUE
          )
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
