# Coefficient bootstrap summaries.

#' Bootstrap confidence intervals for model coefficients
#'
#' Performs fast group bootstrap resampling to compute confidence intervals for
#' all model parameters (intercepts and coefficients). This is a general-purpose
#' bootstrap function that works with any \code{rms::orm} or \code{VGAM::vglm}
#' model. Uses a custom fast bootstrap implementation that is dramatically
#' faster than rsample::group_bootstraps() when working with many groups.
#'
#' @param model A fitted model object from \code{rms::orm} or \code{VGAM::vglm}.
#'   For \code{orm}, should be fitted with \code{x = TRUE, y = TRUE} to enable
#'   model updating with bootstrap samples.
#' @param data A data frame containing the patient data. This is required (cannot be NULL)
#'   because the ID variable is needed for group bootstrap but is typically not included
#'   in the model formula.
#' @param n_boot Number of bootstrap samples
#' @param workers Number of workers for parallel processing. If NULL (default)
#'   or 1, sequential processing is used. If > 1, uses parallel processing.
#' @param parallel (Deprecated) Logical indicating whether to use parallel
#'   processing. Use the \code{workers} parameter instead.
#' @param id_var Name of the ID variable for group bootstrap (default "id")
#' @param use_coefstart Logical. If TRUE, uses the original model coefficients as
#'   starting values when refitting. This can speed up convergence but may affect
#'   results if bootstrap samples differ substantially from the original data.
#'   Default is FALSE.
#'
#' @return A tibble with bootstrap results. Each row represents one bootstrap
#'   iteration and contains:
#'   \itemize{
#'     \item boot_id: Bootstrap iteration number
#'     \item All model coefficients (both intercepts and slope parameters)
#'   }
#'
#' @details
#' Uses fast group bootstrap resampling (resampling by patient ID) to preserve
#' the within-patient correlation structure. This implementation uses a
#' memory-efficient just-in-time (JIT) approach that is dramatically faster and
#' more scalable than rsample::group_bootstraps() for datasets with many groups.
#' Memory usage scales with number of groups (not rows), making it feasible to
#' bootstrap large datasets. See GitHub issue tidymodels/rsample#357 for details.
#'
#' For each bootstrap iteration:
#' \enumerate{
#'   \item Extracts the bootstrap sample and creates unique IDs
#'   \item Relevels factor variables to handle missing levels (converts to consecutive integers)
#'   \item Updates the datadist in the global environment (safe with future.callr)
#'   \item Refits the model using \code{update()} with the bootstrap sample
#'   \item Extracts all coefficients from the refitted model
#' }
#'
#' \strong{Handling missing factor levels:} When certain factor levels are absent
#' from a bootstrap sample, the function automatically relevels ordered factors
#' (like state variables) to consecutive integers. This prevents model fitting
#' failures while maintaining the ordinal structure.
#'
#' Parallelization is handled via future.callr::callr strategy, which provides
#' isolated R processes for each worker, making it safe to modify the global
#' environment (for datadist) without conflicts.
#'
#' The returned tibble can be used to compute confidence intervals using
#' quantile-based methods (percentile or BCa intervals).
#'
#' \strong{Previous-state type:} \code{yprev} may be a factor for categorical
#' previous-state effects or numeric for linear/spline effects. Bootstrap
#' releveling preserves the column type used to fit the original model.
#'
#' @keywords bootstrap coefficients confidence intervals
#'
#' @importFrom stats update coef
#'
#' @examples
#' \dontrun{
#' # Fit initial model (note: id is NOT in the formula but is needed for bootstrap)
#' d <- rms::datadist(my_data)
#' options(datadist = "d")
#' m <- orm(y ~ tx + yprev + time, data = my_data, x = TRUE, y = TRUE)
#'
#' # Bootstrap all coefficients (must provide data with id variable)
#' bs_coefs <- bootstrap_model_coefs(
#'   model = m,
#'   data = my_data,  # Required - must contain id variable
#'   n_boot = 1000
#' )
#'
#' # Compute 95% confidence intervals
#' library(dplyr)
#' ci_results <- bs_coefs |>
#'   summarise(across(-boot_id, list(
#'     lower = ~quantile(.x, 0.025, na.rm = TRUE),
#'     upper = ~quantile(.x, 0.975, na.rm = TRUE)
#'   )))
#' }
#' @export

bootstrap_model_coefs <- function(
  model,
  data = NULL,
  n_boot,
  workers = NULL,
  parallel = NULL,
  id_var = "id",
  use_coefstart = FALSE
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

  # Check model class
  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  # Extract data from model if not provided
  if (is.null(data)) {
    stop(
      "No data provided. The 'data' argument is required for bootstrap_model_coefs() ",
      "because the ID variable is typically not included in the model matrix. ",
      "Please provide the original data frame used to fit the model."
    )
  }

  # Check that id_var exists in data
  if (!id_var %in% names(data)) {
    stop("id_var '", id_var, "' not found in data")
  }

  # Identify state columns that may need releveling in bootstrap samples.
  # Numeric yprev is preserved for linear/spline previous-state models.
  factor_cols <- unique(c(names(data)[sapply(data, is.factor)], "yprev"))
  factor_cols <- intersect(factor_cols, names(data))

  # Generate bootstrap ID samples using fast helper (memory-efficient JIT approach)
  boot_ids <- fast_group_bootstrap(
    data = data,
    id_var = id_var,
    n_boot = n_boot
  )

  # Define analysis function
  analysis_fn <- function(boot_data) {
    # Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = data,
      ylevels = NULL,
      absorb = NULL,
      update_datadist = inherits(model, "orm"),
      use_coefstart = use_coefstart
    )

    m_boot <- boot_result$model

    # Extract coefficients
    if (!is.null(m_boot)) {
      coefs <- coef(m_boot)
      return(as.list(coefs))
    } else {
      return(NULL)
    }
  }

  # Apply analysis function to bootstrap samples with JIT materialization
  bs_coefs_list <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = data,
    id_var = id_var,
    workers = workers,
    packages = c("rms", "VGAM", "stats"),
    globals = c(
      "model",
      "factor_cols"
    )
  )

  # Convert to a data frame with boot_id and one column per coefficient.
  result <- named_list_to_wide(bs_coefs_list, id = seq_len(n_boot))

  return(result)
}
