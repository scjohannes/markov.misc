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
#' @importFrom dplyr left_join
#' @export

materialize_bootstrap_sample <- function(boot_ids, data, id_var) {
  # Rename original_id to match id_var for joining
  boot_ids_renamed <- boot_ids
  names(boot_ids_renamed)[names(boot_ids_renamed) == "original_id"] <- id_var

  # Join with original data
  boot_sample <- dplyr::left_join(
    boot_ids_renamed,
    data,
    by = id_var,
    relationship = "many-to-many"
  )

  return(boot_sample)
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
#' @param parallel Logical indicating whether to use parallel processing
#'   (default FALSE)
#' @param workers Number of workers for parallel processing. If NULL (default),
#'   uses parallel::detectCores() - 1
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
#' If \code{parallel = TRUE}, uses future.callr::callr strategy for parallel
#' processing, which provides isolated R processes for each worker. This makes
#' it safe to modify the global environment (e.g., for datadist).
#'
#' If \code{parallel = FALSE}, uses sequential processing with purrr::map().
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
#'   parallel = TRUE,
#'   globals = c("original_model")
#' )
#' }
#'
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan
#' @importFrom future.callr callr
#' @importFrom purrr map
#' @export

apply_to_bootstrap <- function(
  boot_samples,
  analysis_fn,
  data,
  id_var,
  parallel = FALSE,
  workers = NULL,
  packages = c("rms", "VGAM", "Hmisc", "stats", "dplyr"),
  globals = character(0)
) {
  # Set up workers
  if (is.null(workers)) {
    workers <- parallel::detectCores() - 1
  }

  # Wrapper function that materializes data then runs analysis
  wrapper_fn <- function(boot_ids) {
    # Materialize bootstrap sample from ID lookup
    boot_data <- materialize_bootstrap_sample(boot_ids, data, id_var)

    # Run user's analysis function
    analysis_fn(boot_data)
  }

  # Apply function to bootstrap samples
  if (parallel) {
    plan(callr, workers = workers)

    results <- furrr::future_map(
      boot_samples,
      wrapper_fn,
      .options = furrr::furrr_options(
        packages = c(packages),
        globals = c(globals, "data", "id_var", "analysis_fn")
      )
    )

    plan("sequential")
  } else {
    results <- purrr::map(boot_samples, wrapper_fn)
  }

  return(results)
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
#' @param ylevels Original state levels (e.g., 1:6)
#' @param absorb Absorbing state in original levels
#' @param update_datadist Logical indicating whether to update datadist
#'   (only needed for rms::orm models)
#'
#' @return A list with components:
#'   - model: Refitted model on bootstrap data (or NULL if fitting failed)
#'   - data: Bootstrap data with releveled factors
#'   - ylevels: Updated ylevels (or NULL if not provided)
#'   - absorb: Updated absorb (or NULL if not provided)
#'   - missing_states: Character vector of missing state levels
#'
#' @keywords internal
#' @noRd

bootstrap_analysis_wrapper <- function(
  boot_data,
  model,
  factor_cols,
  original_data,
  ylevels = NULL,
  absorb = NULL,
  update_datadist = TRUE,
  use_coefstart = FALSE
) {
  # Relevel factors to handle missing states
  releveled <- relevel_factors_consecutive(
    data = boot_data,
    factor_cols = factor_cols,
    original_data = original_data,
    ylevels = ylevels,
    absorb = absorb
  )

  boot_data <- releveled$data
  boot_ylevels <- releveled$ylevels
  boot_absorb <- releveled$absorb
  missing_states <- releveled$missing_levels

  # Update datadist if needed (only for orm models)
  if (update_datadist && inherits(model, "orm")) {
    dd <- rms::datadist(boot_data)
    assign("dd", dd, envir = .GlobalEnv)
    options(datadist = "dd")
  }

  # Refit model
  m_boot <- tryCatch(
    {
      # Use coefstart only if requested, valid for model type, and safe (no missing states)
      if (
        use_coefstart && inherits(model, "vglm") && length(missing_states) == 0
      ) {
        stats::update(model, data = boot_data, coefstart = stats::coef(model))
      } else {
        stats::update(model, data = boot_data)
      }
    },
    error = function(e) {
      warning("Bootstrap iteration failed: ", e$message)
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
