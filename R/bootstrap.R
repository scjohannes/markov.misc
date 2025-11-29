#' Use bootstrap to compute confidence intervals for the time in target state(s)
#'
#' Performs fast group bootstrap resampling and computes confidence intervals for
#' the time spent in specified target state(s). Uses a custom fast bootstrap
#' implementation that is dramatically faster than rsample::group_bootstraps()
#' when working with many groups (e.g., >500 patients).
#'
#' @param model A fitted model object from \code{rms::orm} or \code{VGAM::vglm}.
#'   For \code{orm}, should be fitted with \code{x = TRUE, y = TRUE} to enable
#'   model updating with bootstrap samples.
#' @param data A data frame containing the patient data
#' @param n_boot Number of bootstrap samples
#' @param workers Number of workers used for parallelization. Default is
#'   parallel::detectCores() - 1
#' @param parallel Whether parallelization should be used (default FALSE)
#' @param ylevels States in the data
#' @param absorb Absorbing state
#' @param times Time points in the data. Default is \code{1:max(data[["time"]])}
#' @param target_states Integer vector of target state(s) to calculate time in.
#'   Default is 1 (typically "home" or baseline state). Can specify multiple
#'   states, e.g., c(1, 2) to sum time across states 1 and 2.
#' @param varnames List of variable names in the data
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item id: Bootstrap iteration identifier (e.g., "Bootstrap01", "Bootstrap02", ...)
#'     \item models: List column containing the bootstrap estimates of time in target state(s)
#'   }
#'
#' @details
#' Uses fast group bootstrap resampling (resampling by patient ID) to preserve
#' the within-patient correlation structure. This implementation uses a
#' memory-efficient just-in-time (JIT) approach that is dramatically faster and
#' more scalable than rsample::group_bootstraps() for datasets with many groups.
#' Memory usage scales with number of groups (not rows), making it feasible to
#' bootstrap large datasets (e.g., 250k rows, 1000 groups, 2000 bootstraps).
#' See GitHub issue tidymodels/rsample#357 for details on the rsample performance issues.
#'
#' For each bootstrap iteration:
#' \enumerate{
#'   \item Samples patient IDs with replacement using \code{\link{fast_group_bootstrap}}
#'   \item Creates unique IDs for resampled patients (e.g., "P001_1", "P001_2")
#'   \item Checks which states are present and relevels factors to consecutive integers
#'   \item Adjusts \code{ylevels} and \code{absorb} parameters if states are missing
#'   \item Updates the datadist in the global environment (safe with future.callr)
#'   \item Refits the model using \code{update()} with the bootstrap sample
#'   \item Passes the refitted model to \code{\link{taooh}} to calculate the estimand
#' }
#'
#' \strong{Handling missing states:} When certain states are absent from a
#' bootstrap sample (e.g., if state 3 is missing), the function automatically:
#' \itemize{
#'   \item Relevels \code{y} and \code{yprev} to consecutive integers (1, 2, 4, 5, 6 becomes 1, 2, 3, 4, 5)
#'   \item Updates \code{ylevels} to the new consecutive sequence
#'   \item Adjusts the \code{absorb} parameter to the new position of the absorbing state
#' }
#' This prevents model fitting failures while maintaining the ordinal structure.
#'
#' Parallelization is handled via future.callr::callr strategy, which provides
#' isolated R processes for each worker, making it safe to modify the global
#' environment (for datadist) without conflicts.
#'
#' \strong{Important:} The \code{yprev} variable should be a factor before
#' fitting the model for proper handling of state levels.
#'
#' @keywords bootstrap time alive and out of hospital
#'
#' @importFrom tibble tibble
#' @importFrom stats update
#'
#' @examples
#' \dontrun{
#' # Fit initial model
#' d <- rms::datadist(my_data)
#' options(datadist = "d")
#' m <- orm(y ~ tx + yprev + time, data = my_data, x = TRUE, y = TRUE)
#'
#' # Bootstrap for time at home (state 1)
#' bs_results <- taooh_bootstrap(
#'   model = m,
#'   data = my_data,
#'   n_boot = 1000,
#'   target_states = 1
#' )
#'
#' # Bootstrap for time in states 1 or 2
#' bs_results <- taooh_bootstrap(
#'   model = m,
#'   data = my_data,
#'   n_boot = 1000,
#'   target_states = c(1, 2)
#' )
#' }
#' @export

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Handle missing states in the data (change ylevels and absorb depending on the data)
# - Expand to support vgam models in addition to orm

taooh_bootstrap <- function(
  model,
  data,
  n_boot,
  workers = NULL,
  parallel = FALSE,
  ylevels = as.character(1:6),
  absorb = 6,
  times = NULL,
  target_states = 1,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Check model class
  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  # Set times if not provided
  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  # Generate bootstrap ID samples using fast helper (memory-efficient JIT approach)
  boot_ids <- fast_group_bootstrap(
    data = data,
    id_var = varnames$id,
    n_boot = n_boot
  )

  # Define analysis function
  analysis_fn <- function(boot_data) {
    # Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = c("y", varnames$pvarname),
      original_data = data,
      ylevels = ylevels,
      absorb = absorb,
      update_datadist = TRUE
    )

    m_boot <- boot_result$model
    boot_data <- boot_result$data
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb

    # Calculate taooh with refitted model
    if (!is.null(m_boot)) {
      tryCatch(
        taooh(
          model = m_boot,
          data = boot_data,
          times = times,
          ylevels = boot_ylevels,
          absorb = boot_absorb,
          target_states = target_states,
          varnames = list(
            tvarname = varnames$tvarname,
            pvarname = varnames$pvarname,
            id = "new_id",
            tx = varnames$tx
          )
        )[1],
        error = function(e) {
          warning("taooh calculation failed: ", e$message)
          return(NA_real_)
        }
      )
    } else {
      return(NA_real_)
    }
  }

  # Apply analysis function to bootstrap samples with JIT materialization
  bs_estimates <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = data,
    id_var = varnames$id,
    parallel = parallel,
    workers = workers,
    packages = c("VGAM", "rms", "Hmisc", "stats", "dplyr"),
    globals = c(
      "model",
      "times",
      "ylevels",
      "absorb",
      "target_states",
      "varnames"
    )
  )

  # Create result tibble matching original format
  result <- tibble::tibble(
    id = paste0("Bootstrap", sprintf("%02d", seq_len(n_boot))),
    models = bs_estimates
  )

  return(result)
}


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
#' @param workers Number of workers used for parallelization. Default is
#'   parallel::detectCores() - 1
#' @param parallel Whether parallelization should be used (default FALSE)
#' @param id_var Name of the ID variable for group bootstrap (default "id")
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
#' \strong{Important:} Factor variables (like \code{yprev}) should be factors
#' before fitting the model for proper handling of factor levels.
#'
#' @keywords bootstrap coefficients confidence intervals
#'
#' @importFrom tibble tibble
#' @importFrom stats update coef
#' @importFrom tidyr unnest_wider
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
  parallel = FALSE,
  id_var = "id"
) {
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

  # Check for yprev variable and warn if not a factor
  if ("yprev" %in% names(data) && !is.factor(data$yprev)) {
    warning(
      "Variable 'yprev' should be a factor but is not. "
    )
  }

  # Identify factor columns
  factor_cols <- names(data)[sapply(data, is.factor)]

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
      update_datadist = inherits(model, "orm")
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
    parallel = parallel,
    workers = workers,
    packages = c("rms", "VGAM", "stats", "dplyr"),
    globals = c(
      "model",
      "factor_cols"
    )
  )

  # Convert to tibble with boot_id and unnest coefficients
  result <- tibble::tibble(
    boot_id = seq_len(n_boot),
    coefs = bs_coefs_list
  ) |>
    tidyr::unnest_wider(coefs)

  return(result)
}


#' Bootstrap confidence intervals for standardized state occupancy probabilities
#'
#' Performs fast group bootstrap resampling to compute confidence intervals for
#' standardized state occupancy probabilities (SOPs) across all states and time
#' points. Uses g-computation to estimate marginal (population-averaged) SOPs
#' under treatment and control. Uses a custom fast bootstrap implementation that
#' is dramatically faster than rsample::group_bootstraps() when working with
#' many groups.
#'
#' @param model A fitted model object from \code{rms::orm} or \code{VGAM::vglm}.
#'   For \code{orm}, should be fitted with \code{x = TRUE, y = TRUE} to enable
#'   model updating with bootstrap samples.
#' @param data A data frame containing the patient trajectory data
#' @param n_boot Number of bootstrap samples
#' @param workers Number of workers used for parallelization. Default is
#'   NULL
#' @param parallel Whether parallelization should be used (default FALSE).
#'   Set to TRUE only if this function is being called in isolation. If you're
#'   running multiple simulations in parallel at a higher level, keep this FALSE
#'   to avoid nested parallelization and resource contention.
#' @param ylevels States in the data (e.g., 1:6)
#' @param absorb Absorbing state (e.g., 6 for death)
#' @param times Time points in the data. Default is \code{1:max(data[["time"]])}
#' @param varnames List of variable names in the data with components:
#'   \itemize{
#'     \item tvarname: time variable name (default "time")
#'     \item pvarname: previous state variable name (default "yprev")
#'     \item id: patient identifier (default "id")
#'     \item tx: treatment indicator (default "tx")
#'   }
#'
#' @return A tibble with columns:
#'   \itemize{
#'     \item boot_id: Bootstrap iteration number
#'     \item time: Time point
#'     \item tx: Treatment indicator (0 = control, 1 = treatment)
#'     \item state_1, state_2, ..., state_K: Probability of being in each state
#'       at the given time point under the given treatment
#'   }
#'   Each bootstrap iteration contributes 2 * length(times) rows (one set for
#'   treatment, one set for control).
#'
#' @details
#' For each bootstrap iteration:
#' \enumerate{
#'   \item Extracts the bootstrap sample and creates unique IDs
#'   \item Uses \code{\link{relevel_factors_consecutive}} to handle missing states
#'     by releveling to consecutive integers and updating ylevels/absorb
#'   \item Refits the model using the bootstrap sample
#'   \item Computes standardized SOPs using \code{\link{standardize_sops}}
#'   \item Expands results back to original state space, padding with zeros for missing states
#' }
#'
#' \strong{Handling missing states:}
#' When a state (e.g., state 3) is missing from a bootstrap sample:
#' \itemize{
#'   \item The model is fit on releveled data (states 1,2,4,5,6 become 1,2,3,4,5)
#'   \item SOPs are predicted for the releveled states
#'   \item Results are mapped back to original state numbering (1,2,3,4,5,6)
#'   \item Missing state 3 is assigned probability 0 at all time points
#' }
#'
#' This ensures that all bootstrap iterations return SOPs for all original states,
#' which is required for computing valid confidence intervals.
#'
#' \strong{Output format:} The function returns a single long-format tibble
#' containing both treatment and control SOPs across all bootstrap iterations.
#' Each row represents one time point for one treatment group in one bootstrap
#' iteration. This format makes it easy to:
#' \itemize{
#'   \item Filter by treatment group (tx == 0 or tx == 1)
#'   \item Compute confidence bands using group_by(time, tx) and quantiles
#'   \item Plot results with ggplot2
#'   \item Calculate treatment effects by comparing tx == 1 vs tx == 0
#' }
#'
#' Parallelization is handled via future.callr::callr strategy, which provides
#' isolated R processes for each worker.
#'
#' @keywords bootstrap standardization state occupancy g-computation
#'
#' @importFrom tibble as_tibble
#' @importFrom stats update
#' @importFrom dplyr select bind_rows
#'
#' @examples
#' \dontrun{
#' # Fit initial model
#' d <- rms::datadist(my_data)
#' options(datadist = "d")
#' m <- orm(y ~ tx + yprev + time, data = my_data, x = TRUE, y = TRUE)
#'
#' # Bootstrap standardized SOPs
#' bs_sops <- bootstrap_standardized_sops(
#'   model = m,
#'   data = my_data,
#'   times = 1:30,
#'   n_boot = 1000,
#'   ylevels = 1:6,
#'   absorb = 6
#' )
#'
#' # Compute 95% confidence bands for treatment SOPs
#' library(dplyr)
#' library(tidyr)
#'
#' # Filter to treatment group and reshape to long format
#' sops_tx_long <- bs_sops |>
#'   filter(tx == 1) |>
#'   pivot_longer(
#'     cols = starts_with("state_"),
#'     names_to = "state",
#'     values_to = "probability",
#'     names_prefix = "state_"
#'   )
#'
#' # Compute confidence intervals
#' ci_bands <- sops_tx_long |>
#'   group_by(time, state) |>
#'   summarise(
#'     lower = quantile(probability, 0.025),
#'     median = median(probability),
#'     upper = quantile(probability, 0.975),
#'     .groups = "drop"
#'   )
#'
#' # Plot with confidence bands
#' library(ggplot2)
#' ggplot(ci_bands, aes(x = time, y = median, color = state)) +
#'   geom_line() +
#'   geom_ribbon(aes(ymin = lower, ymax = upper, fill = state), alpha = 0.2) +
#'   facet_wrap(~state) +
#'   labs(title = "Standardized SOPs under Treatment with 95% CI")
#' }
#' @export

bootstrap_standardized_sops <- function(
  model,
  data,
  n_boot,
  workers = NULL,
  parallel = FALSE,
  ylevels = 1:6,
  absorb = 6,
  times = NULL,
  update_datadist = TRUE,
  use_coefstart = TRUE,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Check model class
  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  # Set times if not provided
  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  n_times <- length(times)
  n_states <- length(ylevels)

  # Generate bootstrap ID samples using fast helper (memory-efficient JIT approach)
  boot_ids <- fast_group_bootstrap(
    data = data,
    id_var = varnames$id,
    n_boot = n_boot
  )

  # Define analysis function
  analysis_fn <- function(boot_data) {
    # Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = c("y", varnames$pvarname),
      original_data = data,
      ylevels = ylevels,
      absorb = absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )

    m_boot <- boot_result$model
    boot_data <- boot_result$data
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }

    # Get states present in original numbering (for mapping back later)
    states_present <- setdiff(ylevels, as.numeric(missing_states))

    # Compute standardized SOPs on releveled states
    sop_result <- tryCatch(
      standardize_sops(
        model = m_boot,
        data = boot_data,
        times = times,
        ylevels = boot_ylevels,
        absorb = boot_absorb,
        varnames = list(
          tvarname = varnames$tvarname,
          pvarname = varnames$pvarname,
          id = "new_id",
          tx = varnames$tx
        )
      ),
      error = function(e) {
        warning("standardize_sops failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(sop_result)) {
      return(NULL)
    }

    # Expand back to original state space with zeros for missing states
    if (length(missing_states) > 0) {
      # Initialize full matrices with zeros
      sop_tx_full <- matrix(0, nrow = n_times, ncol = n_states)
      sop_ctrl_full <- matrix(0, nrow = n_times, ncol = n_states)
      colnames(sop_tx_full) <- as.character(ylevels)
      colnames(sop_ctrl_full) <- as.character(ylevels)

      # Fill in states that were present
      for (i in seq_along(states_present)) {
        original_state <- states_present[i]
        sop_tx_full[, as.character(original_state)] <- sop_result$sop_tx[, i]
        sop_ctrl_full[, as.character(original_state)] <- sop_result$sop_ctrl[,
          i
        ]
      }

      return(list(
        sop_tx = sop_tx_full,
        sop_ctrl = sop_ctrl_full
      ))
    } else {
      # All states present - return as is
      return(sop_result)
    }
  }

  # Apply analysis function to bootstrap samples with JIT materialization
  sop_results <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = data,
    id_var = varnames$id,
    parallel = parallel,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "dplyr", "stats"),
    globals = c(
      "model",
      "times",
      "ylevels",
      "absorb",
      "varnames",
      "n_times",
      "n_states"
    )
  )

  # Extract and reshape results into long format
  # Format: boot_id | time | tx | state_1 | state_2 | ... | state_k
  result_list <- list()

  for (i in seq_along(sop_results)) {
    sop_i <- sop_results[[i]]

    if (!is.null(sop_i)) {
      # Convert treatment SOPs to tibble with time and state columns
      tx_df <- tibble::as_tibble(sop_i$sop_tx)
      colnames(tx_df) <- paste0("state_", ylevels)
      tx_df$time <- times
      tx_df$tx <- 1
      tx_df$boot_id <- i

      # Convert control SOPs to tibble with time and state columns
      ctrl_df <- tibble::as_tibble(sop_i$sop_ctrl)
      colnames(ctrl_df) <- paste0("state_", ylevels)
      ctrl_df$time <- times
      ctrl_df$tx <- 0
      ctrl_df$boot_id <- i

      # Combine treatment and control
      result_list[[i]] <- dplyr::bind_rows(tx_df, ctrl_df)
    }
  }

  # Combine all bootstrap iterations
  result <- dplyr::bind_rows(result_list) |>
    dplyr::select(boot_id, time, tx, dplyr::everything())

  return(result)
}


#' Use bootstrap to compute confidence intervals for the time in target state(s)
#'
#' Performs group bootstrap resampling and computes confidence intervals for the
#' time spent in specified target state(s). Uses parallel computation via
#' future.callr for efficiency.
#'
#' @param model A fitted model object from \code{rms::orm} or \code{VGAM::vglm}.
#' See details for specifications of these models.
#' @param data A data frame containing the patient data
#' @param n_boot Number of bootstrap samples
#' @param workers Number of workers used for parallelization. Default is
#'   parallel::detectCores() - 1
#' @param parallel Whether parallelization should be used (default TRUE)
#' @param ylevels States in the data
#' @param absorb Absorbing state
#' @param times Time points in the data. Default is \code{1:max(data[["time"]])}
#' @param target_states Integer vector of target state(s) to calculate time in.
#'   Default is 1 (typically "home" or baseline state). Can specify multiple
#'   states, e.g., c(1, 2) to sum time across states 1 and 2.
#' @param varnames List of variable names in the data
#'
#' @return A vector containing estimates of time in the target state
#'
#' @details
#' Uses group bootstrap resampling (resampling by patient ID) to preserve the
#' within-patient correlation structure. For each bootstrap iteration:
#' \enumerate{
#'   \item Extracts the bootstrap sample and creates unique IDs
#'   \item Checks which states are present and relevels factors to consecutive integers
#'   \item Adjusts \code{ylevels} and \code{absorb} parameters if states are missing
#'   \item Updates the datadist in the global environment (safe with future.callr)
#'   \item Refits the model using \code{update()} with the bootstrap sample
#'   \item Passes the refitted model to \code{\link{taooh}} to calculate the estimand
#' }
#'
#' \strong{Handling missing states:} When certain states are absent from a
#' bootstrap sample (e.g., if state 3 is missing), the function automatically:
#' \itemize{
#'   \item Relevels \code{y} and \code{yprev} to consecutive integers (1, 2, 4, 5, 6 becomes 1, 2, 3, 4, 5)
#'   \item Updates \code{ylevels} to the new consecutive sequence
#'   \item Adjusts the \code{absorb} parameter to the new position of the absorbing state
#' }
#' This prevents model fitting failures while maintaining the ordinal structure.
#'
#' Parallelization is handled via future.callr::callr strategy, which provides
#' isolated R processes for each worker, making it safe to modify the global
#' environment (for datadist) without conflicts.
#'
#' \strong{Important:} The \code{yprev} variable should be a factor before
#' fitting the model for proper handling of state levels.
#'
#' orm should be fitted with \code{x = TRUE, y = TRUE} to enable model updating with
#' bootstrap samples. vglm should use factor in the model formula (such as
#' \code{y ~ time + tx+ factor(yprev)}) and \code{model = TRUE} (see VGAM::vglm documentation)
#' to be compatible with \code{soprobMarkovOrdm()}. Otherwise soprobMarkovOrdm will
#' complain about "yprev" not being a factor.
#'
#' @keywords bootstrap time alive and out of hospital
#'
#' @importFrom future.callr callr
#' @importFrom future plan
#' @importFrom furrr future_map_dbl
#' @importFrom rms datadist
#' @importFrom rms rcs
#' @importFrom stats ave update
#'
#' @examples
#' \dontrun{
#' # Fit initial model
#' d <- rms::datadist(my_data)
#' options(datadist = "d")
#' m <- orm(y ~ tx + yprev + time, data = my_data, x = TRUE, y = TRUE)
#'
#' # Bootstrap for time at home (state 1)
#' bs_results <- taooh_bootstrap2(
#'   model = m,
#'   data = my_data,
#'   n_boot = 1000,
#'   target_states = 1
#' )
#'
#' # Bootstrap for time in states 1 or 2
#' bs_results <- taooh_bootstrap2(
#'   model = m,
#'   data = my_data,
#'   n_boot = 1000,
#'   target_states = c(1, 2)
#' )
#' }
#' @export

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Expand to support vgam models in addition to orm

taooh_bootstrap2 <- function(
  model,
  data,
  n_boot,
  workers = NULL,
  parallel = TRUE,
  ylevels = 1:6,
  absorb = 6,
  times = NULL,
  target_states = 1,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx"),
  quantiles = 0.95
) {
  # Check model class
  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  # Create matrix containing sampled IDs
  n_data <- length(unique(data[[varnames$id]]))

  sample_IDs <- matrix(
    sample(unique(data[[varnames$id]]), size = n_data * n_boot, replace = TRUE),
    ncol = n_boot
  )

  # Also split the original data by id to access each block later
  df_split <- split(data, data[[varnames$id]])

  # Loop arguments
  if (is.null(workers)) {
    workers <- parallel::detectCores() - 1
  }

  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  # Apply the taooh() function to the bootstrap samples in parallel.
  if (parallel) {
    plan(future.callr::callr, workers = workers)
  } else {
    plan("sequential")
  }

  # Loop through the resamples
  result <- future_map_dbl(as.data.frame(sample_IDs), function(x) {
    # Create resample based on selected IDs
    boot_data <- do.call(
      rbind,
      lapply(x, function(id) {
        df_split[[as.character(id)]]
      })
    )

    # Assign new unique IDs
    boot_data$block_count <- ave(
      boot_data[[varnames$id]],
      boot_data[[varnames$id]],
      boot_data[[varnames$tvarname]],
      FUN = seq_along
    )

    boot_data$new_id <- factor(paste(
      boot_data[[varnames$id]],
      boot_data$block_count,
      sep = "_"
    ))

    # Handle missing states in bootstrap sample by releveling
    # Get actual states present in bootstrap sample
    states_in_y <- unique(boot_data$y)
    states_in_yprev <- unique(boot_data[[varnames$pvarname]])
    states_present <- sort(unique(c(states_in_y, states_in_yprev)))

    # Check if all original states are present
    original_states <- if (is.factor(data$y)) {
      as.numeric(levels(data$y))
    } else {
      ylevels
    }

    missing_states <- setdiff(original_states, states_present)

    if (length(missing_states) > 0) {
      # Some states are missing - need to relevel to consecutive integers
      # Create mapping from old to new state numbers
      state_mapping <- setNames(seq_along(states_present), states_present)

      # Relevel y and yprev to consecutive integers
      boot_data$y <- state_mapping[as.character(boot_data$y)]
      boot_data[[varnames$pvarname]] <- state_mapping[
        as.character(boot_data[[varnames$pvarname]])
      ]

      # Update ylevels to new consecutive sequence
      boot_ylevels <- seq_along(states_present)

      # Update absorbing state if it was in the original states
      if (absorb %in% states_present) {
        boot_absorb <- state_mapping[as.character(absorb)]
      } else {
        # If absorbing state is missing, use the highest state
        boot_absorb <- max(boot_ylevels)
        warning(
          "Absorbing state ",
          absorb,
          " not present in bootstrap sample. ",
          "Using state ",
          boot_absorb,
          " as absorbing state."
        )
      }
    } else {
      # All states present - use original values
      boot_ylevels <- ylevels
      boot_absorb <- absorb
    }

    # Now convert to factors with the appropriate levels
    boot_data$y <- factor(boot_data$y, levels = boot_ylevels)
    boot_data[[varnames$pvarname]] <- factor(
      boot_data[[varnames$pvarname]],
      levels = boot_ylevels
    )
    boot_data$yprev <- factor(boot_data$yprev)

    # Update datadist and assign to global environment
    # This works because we use future.callr which has isolated processes
    if (inherits(model, "orm")) {
      dd <- rms::datadist(boot_data)
      assign("dd", dd, envir = .GlobalEnv)
      options(datadist = "dd")
    }

    # Refit model with bootstrap data
    m_boot <- tryCatch(
      update(model, data = boot_data),
      error = function(e) {
        warning("Bootstrap iteration failed: ", e$message)
        return(NULL)
      }
    )

    # Calculate taooh with refitted model
    if (!is.null(m_boot)) {
      tryCatch(
        taooh(
          model = m_boot,
          data = boot_data,
          times = times,
          ylevels = boot_ylevels,
          absorb = boot_absorb,
          target_states = target_states,
          varnames = list(
            tvarname = varnames$tvarname,
            pvarname = varnames$pvarname,
            id = "new_id", # Use the new_id for bootstrap samples
            tx = varnames$tx
          )
        )[1],
        error = function(e) {
          warning("taooh calculation failed: ", e$message)
          return(NA_real_)
        }
      )
    } else {
      return(NA_real_)
    }
  })

  # Free memory
  plan("sequential")

  # Get confidence intervals for time in target state

  return(list(
    t_est = mean(result),
    ci = quantile(result, c(1 - quantiles, quantiles), na.rm = TRUE),
    errors = sum(is.na(result)),
    all_results = result
  ))
}
