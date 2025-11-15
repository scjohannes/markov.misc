#' Compute time in target state(s) for treatment groups
#'
#' Calculates the expected time spent in specified target state(s) for treatment
#' and control groups using a proportional odds Markov model. Uses an existing
#' fitted model and state occupancy probabilities to estimate the total time in
#' the target state(s).
#'
#' @param model A fitted model object from \code{rms::orm}. The model is used
#'   as-is without refitting.
#' @param data A data frame containing patient trajectory data. If NULL, attempts
#'   to extract data from the model object (requires model fitted with
#'   \code{x = TRUE, y = TRUE}).
#' @param times Time points in the data. Default is \code{1:max(data[["time"]])}
#' @param ylevels States in the data (e.g., 1:6)
#' @param absorb Absorbing state (e.g., 6 for death)
#' @param target_states Integer vector of target state(s) to calculate time in.
#'   Default is 1 (typically "home" or baseline state). Can specify multiple
#'   states, e.g., c(1, 2) to sum time across states 1 and 2.
#' @param varnames List of variable names in the data with components:
#'   \itemize{
#'     \item tvarname: time variable name (default "time")
#'     \item pvarname: previous state variable name (default "yprev")
#'     \item id: patient identifier (default "id")
#'     \item tx: treatment indicator (default "tx")
#'   }
#'
#' @return A named numeric vector with three elements:
#'   \itemize{
#'     \item delta_taooh: difference in time in target state(s) (treatment - control)
#'     \item SOP_tx: expected time in target state(s) for treatment group
#'     \item SOP_ctrl: expected time in target state(s) for control group
#'   }
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Takes a pre-fitted model (no refitting occurs in this function)
#'   \item Computes state occupancy probabilities for each individual at each time point
#'   \item Sums probabilities across target states and time points
#'   \item Averages within treatment groups to get expected time in target state(s)
#' }
#'
#' This function assumes the model has already been fitted on the appropriate
#' data. For bootstrap inference, see \code{\link{taooh_bootstrap}}, which handles
#' model refitting for each bootstrap sample.
#'
#' @keywords time alive out of hospital state occupancy
#'
#' @importFrom Hmisc soprobMarkovOrdm
#' @importFrom stats aggregate
#'
#' @examples
#' \dontrun{
#' # Fit initial model
#' d <- rms::datadist(my_data)
#' options(datadist = "d")
#' m <- orm(y ~ tx + yprev + time, data = my_data, x = TRUE, y = TRUE)
#'
#' # Calculate time at home (state 1) for a dataset
#' result <- taooh(
#'   model = m,
#'   data = my_data,
#'   ylevels = 1:6,
#'   absorb = 6,
#'   target_states = 1
#' )
#'
#' # Calculate time in states 1 or 2 combined
#' result <- taooh(
#'   model = m,
#'   data = my_data,
#'   target_states = c(1, 2)
#' )
#' }
#' @export

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Handle missing states in the data (change ylevels and absorb depending on the data)
# - Expand to support vgam models in addition to orm

taooh <- function(
  model,
  data = NULL,
  times = NULL,
  ylevels = 1:6,
  absorb = 6,
  target_states = 1,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Check model class
  if (!inherits(model, "orm")) {
    stop(
      "model must be an orm object from rms package. vgam support coming soon."
    )
  }

  # If no data provided, extract from model
  if (is.null(data)) {
    data <- model$x
    if (is.null(data)) {
      stop("No data provided and model was not fitted with x=TRUE")
    }
    # Add response variable
    data$y <- model$y
  }

  # Check that yprev is a factor
  if (!is.factor(data[[varnames$pvarname]])) {
    warning(
      "Variable '",
      varnames$pvarname,
      "' should be a factor but is not. "
    )
  }

  # Determine id variable
  id_var <- varnames$id

  # Set times if not provided
  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  # Generate covariate data.frame to predict on
  X <- data[!duplicated(data[[id_var]]), ] # first row of each individual

  # Run soprobMarkovOrdm to get state probability predictions for each individual
  n_ids <- length(unique(data[[id_var]]))
  n_times <- max(times)
  n_states <- length(target_states)

  # Initialize array to store SOPs for target states
  sop_array <- array(
    dim = c(n_times, n_ids, n_states),
    dimnames = list(NULL, unique(data[[id_var]]), target_states)
  )

  for (i in 1:n_ids) {
    sop_full <- soprobMarkovOrdm(
      object = model,
      data = X[i, ],
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = varnames$tvarname,
      pvarname = varnames$pvarname,
      gap = 1
    )

    # Extract SOPs for target states
    for (s_idx in seq_along(target_states)) {
      state <- target_states[s_idx]
      sop_array[, i, s_idx] <- sop_full[, state]
    }
  }

  # Sum across target states to get total time in any target state
  sop_mat <- apply(sop_array, c(1, 2), sum)
  colnames(sop_mat) <- unique(data[[id_var]])

  # Indicate treatment and control ids
  id_tx <- aggregate(
    data[[varnames$tx]],
    by = list(data[[id_var]]),
    FUN = unique,
    simplify = TRUE
  )

  # Compute time in target state(s) by treatment group
  SOP_tx <- sum(rowMeans(sop_mat[,
    colnames(sop_mat) %in% id_tx[id_tx[, 2] == 1, 1],
    drop = FALSE
  ]))
  SOP_ctrl <- sum(rowMeans(sop_mat[,
    colnames(sop_mat) %in% id_tx[id_tx[, 2] == 0, 1],
    drop = FALSE
  ]))

  return(c(
    delta_taooh = SOP_tx - SOP_ctrl,
    SOP_tx = SOP_tx,
    SOP_ctrl = SOP_ctrl
  ))
}


#' Use bootstrap to compute confidence intervals for the time in target state(s)
#'
#' Performs group bootstrap resampling and computes confidence intervals for the
#' time spent in specified target state(s). Uses parallel computation via
#' future.callr for efficiency.
#'
#' @param model A fitted model object from \code{rms::orm}. Should be fitted
#'   with \code{x = TRUE, y = TRUE} to enable model updating with bootstrap samples.
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
#' @return An rsample bootstrap object with an additional column containing the
#'   bootstrap estimates of time in target state(s)
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
#' @keywords bootstrap time alive and out of hospital
#'
#' @importFrom rsample group_bootstraps analysis
#' @importFrom future.callr callr
#' @importFrom future plan
#' @importFrom furrr future_map
#' @importFrom rms datadist
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
  parallel = TRUE,
  ylevels = 1:6,
  absorb = 6,
  times = NULL,
  target_states = 1,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Check model class
  if (!inherits(model, "orm")) {
    stop(
      "model must be an orm object from rms package. vgam support coming soon."
    )
  }

  # Bootstrap samples
  resample <- group_bootstraps(
    data,
    group = id,
    times = n_boot,
    apparent = FALSE
  )

  # Arguments
  if (is.null(workers)) {
    workers <- parallel::detectCores() - 1
  }

  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  # Apply the taooh() function to the bootstrap samples in parallel.
  if (parallel) {
    plan(callr, workers = workers)
  } else {
    plan("sequential")
  }

  bs_SOP <- resample |>
    dplyr::mutate(
      models = future_map(splits, \(.x) {
        # Extract bootstrap sample
        boot_data <- analysis(.x)

        # Define new id variable (unique to each bootstrap draw)
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
            # If absorbing state is missing, absorbing state becomes NULL
            boot_absorb <- NULL
            warning(
              "Absorbing state ",
              absorb,
              " not present in bootstrap sample. ",
              "Set to NULL."
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

        # Update datadist and assign to global environment
        # This works because we use future.callr which has isolated processes
        dd <- rms::datadist(boot_data)
        assign("dd", dd, envir = .GlobalEnv)
        options(datadist = "dd")

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
    )

  plan("sequential")

  return(bs_SOP)
}


#' Bootstrap confidence intervals for model coefficients
#'
#' Performs group bootstrap resampling to compute confidence intervals for all
#' model parameters (intercepts and coefficients). This is a general-purpose
#' bootstrap function that works with any \code{rms::orm} model.
#'
#' @param model A fitted model object from \code{rms::orm}. Should be fitted
#'   with \code{x = TRUE, y = TRUE} to enable model updating with bootstrap samples.
#' @param data A data frame containing the patient data. This is required (cannot be NULL)
#'   because the ID variable is needed for group bootstrap but is typically not included
#'   in the model formula.
#' @param n_boot Number of bootstrap samples
#' @param workers Number of workers used for parallelization. Default is
#'   parallel::detectCores() - 1
#' @param parallel Whether parallelization should be used (default TRUE)
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
#' Uses group bootstrap resampling (resampling by patient ID) to preserve the
#' within-patient correlation structure. For each bootstrap iteration:
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
#' @importFrom rsample group_bootstraps analysis
#' @importFrom future.callr callr
#' @importFrom future plan
#' @importFrom furrr future_map
#' @importFrom rms datadist
#' @importFrom stats ave update coef
#' @importFrom dplyr select row_number mutate
#' @importFrom tidyr unnest_wider
#' @importFrom rlang sym
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
  parallel = TRUE,
  id_var = "id"
) {
  # Check model class
  if (!inherits(model, "orm")) {
    stop(
      "model must be an orm object from rms package. vgam support coming soon."
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

  # Bootstrap samples
  resample <- group_bootstraps(
    data,
    group = !!rlang::sym(id_var),
    times = n_boot,
    apparent = FALSE
  )

  # Arguments
  if (is.null(workers)) {
    workers <- parallel::detectCores() - 1
  }

  # Apply bootstrap in parallel
  if (parallel) {
    plan(callr, workers = workers)
  } else {
    plan("sequential")
  }

  bs_coefs <- resample |>
    dplyr::mutate(
      boot_id = dplyr::row_number(),
      coefs = future_map(splits, \(.x) {
        # Extract bootstrap sample
        boot_data <- analysis(.x)

        # Define new id variable (unique to each bootstrap draw)
        boot_data$block_count <- ave(
          boot_data[[id_var]],
          boot_data[[id_var]],
          FUN = seq_along
        )
        boot_data$new_id <- factor(paste(
          boot_data[[id_var]],
          boot_data$block_count,
          sep = "_"
        ))

        # Handle missing factor levels by releveling to consecutive values
        # Identify factor columns in original data
        factor_cols <- names(data)[sapply(data, is.factor)]

        for (col in factor_cols) {
          if (col %in% names(boot_data)) {
            # Get levels present in bootstrap sample
            levels_present <- sort(unique(boot_data[[col]]))
            original_levels <- levels(data[[col]])
            missing_levels <- setdiff(original_levels, levels_present)

            if (
              length(missing_levels) > 0 &&
                is.numeric(type.convert(original_levels, as.is = TRUE))
            ) {
              # Numeric factor - relevel to consecutive integers
              # This is important for ordered factors like y and yprev
              level_mapping <- setNames(
                seq_along(levels_present),
                levels_present
              )
              boot_data[[col]] <- level_mapping[as.character(boot_data[[col]])]
              boot_data[[col]] <- factor(
                boot_data[[col]],
                levels = seq_along(levels_present)
              )
            } else {
              # Non-numeric factor - just use present levels
              boot_data[[col]] <- factor(
                boot_data[[col]],
                levels = levels_present
              )
            }
          }
        }

        # Update datadist and assign to global environment
        # This works because we use future.callr which has isolated processes
        dd <- rms::datadist(boot_data)
        assign("dd", dd, envir = .GlobalEnv)
        options(datadist = "dd")

        # Refit model with bootstrap data
        m_boot <- tryCatch(
          update(model, data = boot_data),
          error = function(e) {
            warning("Bootstrap iteration failed: ", e$message)
            return(NULL)
          }
        )

        # Extract coefficients
        if (!is.null(m_boot)) {
          coefs <- coef(m_boot)
          return(as.list(coefs))
        } else {
          return(NULL)
        }
      })
    )

  plan("sequential")

  # Unnest the coefficients into columns
  result <- bs_coefs |>
    dplyr::select(boot_id, coefs) |>
    tidyr::unnest_wider(coefs)

  return(result)
}
