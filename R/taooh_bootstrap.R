#' Compute time in target state(s) for treatment groups
#'
#' Calculates the expected time spent in specified target state(s) for treatment
#' and control groups using a proportional odds Markov model. Uses an existing
#' fitted model and state occupancy probabilities to estimate the total time in
#' the target state(s).
#'
#' @param model A fitted model object from \code{rms::orm} or \code{VGAM::vglm}.
#'   The model is used as-is without refitting.
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


#' Compute standardized state occupancy probabilities under treatment and control
#'
#' Calculates marginalized (standardized) state occupancy probabilities using
#' g-computation. For each individual, predictions are made under both treatment
#' (tx=1) and control (tx=0), then averaged across the population to obtain the
#' marginal state occupancy probabilities for each state at each time point.
#'
#' @param model A fitted model object from \code{rms::orm}. The model is used
#'   as-is without refitting.
#' @param data A data frame containing patient trajectory data. If NULL, attempts
#'   to extract data from the model object (requires model fitted with
#'   \code{x = TRUE, y = TRUE}).
#' @param times Time points in the data. Default is \code{1:max(data[["time"]])}
#' @param ylevels States in the data (e.g., 1:6)
#' @param absorb Absorbing state (e.g., 6 for death)
#' @param varnames List of variable names in the data with components:
#'   \itemize{
#'     \item tvarname: time variable name (default "time")
#'     \item pvarname: previous state variable name (default "yprev")
#'     \item id: patient identifier (default "id")
#'     \item tx: treatment indicator (default "tx")
#'   }
#'
#' @return A list with two components:
#'   \itemize{
#'     \item sop_tx: Matrix of state occupancy probabilities under treatment
#'       (rows = time points, columns = states)
#'     \item sop_ctrl: Matrix of state occupancy probabilities under control
#'       (rows = time points, columns = states)
#'   }
#'
#' @details
#' This function implements g-computation (standardization) to estimate the
#' marginal causal effect of treatment:
#' \enumerate{
#'   \item For each individual, create two counterfactual datasets: one with tx=1
#'     (treatment) and one with tx=0 (control)
#'   \item For each counterfactual dataset, predict state occupancy probabilities
#'     at each time point using \code{soprobMarkovOrdm}
#'   \item Average predictions across all individuals to obtain marginal
#'     (population-averaged) state occupancy probabilities
#' }
#'
#' Unlike \code{\link{taooh}}, which averages within observed treatment groups,
#' this function estimates what would happen if the entire population received
#' treatment vs. if the entire population received control, thus providing
#' a causal estimate under the assumption of no unmeasured confounding.
#'
#' The returned matrices can be used for:
#' \itemize{
#'   \item Plotting standardized state occupancy curves
#'   \item Computing treatment effects for specific states and time periods
#'   \item Calculating summary measures like total time in target states
#' }
#'
#' @keywords standardization g-computation state occupancy causal inference
#'
#' @importFrom Hmisc soprobMarkovOrdm
#'
#' @examples
#' \dontrun{
#' # Fit initial model
#' d <- rms::datadist(my_data)
#' options(datadist = "d")
#' m <- orm(y ~ tx + yprev + time, data = my_data, x = TRUE, y = TRUE)
#'
#' # Get standardized state occupancy probabilities
#' sop_std <- standardize_sops(
#'   model = m,
#'   data = my_data,
#'   ylevels = 1:6,
#'   absorb = 6
#' )
#'
#' # Plot results
#' library(ggplot2)
#' library(tidyr)
#' sop_df <- bind_rows(
#'   as.data.frame(sop_std$sop_tx) |> mutate(time = row_number(), tx = "Treatment"),
#'   as.data.frame(sop_std$sop_ctrl) |> mutate(time = row_number(), tx = "Control")
#' ) |>
#'   pivot_longer(cols = -c(time, tx), names_to = "state", values_to = "probability")
#'
#' ggplot(sop_df, aes(x = time, y = probability, color = state)) +
#'   geom_line() +
#'   facet_wrap(~tx)
#'
#' # Calculate treatment effect for time in state 1 (home)
#' sum(sop_std$sop_tx[, "1"]) - sum(sop_std$sop_ctrl[, "1"])
#' }
#' @export

standardize_sops <- function(
  model,
  data = NULL,
  times = NULL,
  ylevels = 1:6,
  absorb = 6,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Check model class
  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
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
  tx_var <- varnames$tx

  # Set times if not provided
  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  # Generate covariate data.frame to predict on (one row per individual)
  X <- data[!duplicated(data[[id_var]]), ]

  # Number of individuals, time points, and states
  n_ids <- nrow(X)
  n_times <- max(times)
  n_states <- length(ylevels)

  # Initialize arrays to store SOPs under treatment and control
  # Dimensions: [time, individual, state]
  sop_array_tx <- array(
    dim = c(n_times, n_ids, n_states),
    dimnames = list(NULL, unique(data[[id_var]]), as.character(ylevels))
  )

  sop_array_ctrl <- array(
    dim = c(n_times, n_ids, n_states),
    dimnames = list(NULL, unique(data[[id_var]]), as.character(ylevels))
  )

  # Loop through each individual
  for (i in 1:n_ids) {
    # Create counterfactual dataset with tx = 1 (treatment)
    X_tx <- X[i, , drop = FALSE]
    X_tx[[tx_var]] <- 1

    # Create counterfactual dataset with tx = 0 (control)
    X_ctrl <- X[i, , drop = FALSE]
    X_ctrl[[tx_var]] <- 0

    # Predict SOPs under treatment
    sop_full_tx <- soprobMarkovOrdm(
      object = model,
      data = X_tx,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = varnames$tvarname,
      pvarname = varnames$pvarname,
      gap = 1
    )

    # Predict SOPs under control
    sop_full_ctrl <- soprobMarkovOrdm(
      object = model,
      data = X_ctrl,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = varnames$tvarname,
      pvarname = varnames$pvarname,
      gap = 1
    )

    # Store predictions for all states
    sop_array_tx[, i, ] <- sop_full_tx
    sop_array_ctrl[, i, ] <- sop_full_ctrl
  }

  # Average across individuals to get marginal SOPs
  # Result: matrix with rows = time points, columns = states
  sop_tx <- apply(sop_array_tx, c(1, 3), mean)
  sop_ctrl <- apply(sop_array_ctrl, c(1, 3), mean)

  # Ensure column names are character versions of ylevels
  colnames(sop_tx) <- as.character(ylevels)
  colnames(sop_ctrl) <- as.character(ylevels)

  return(list(
    sop_tx = sop_tx,
    sop_ctrl = sop_ctrl
  ))
}


#' Use bootstrap to compute confidence intervals for the time in target state(s)
#'
#' Performs group bootstrap resampling and computes confidence intervals for the
#' time spent in specified target state(s). Uses parallel computation via
#' future.callr for efficiency.
#'
#' @param model A fitted model object from \code{rms::orm} or \code{VGAM::vglm}.
#'   For \code{orm}, should be fitted with \code{x = TRUE, y = TRUE} to enable
#'   model updating with bootstrap samples.
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
#' bootstrap function that works with any \code{rms::orm} or \code{VGAM::vglm} model.
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


#' Bootstrap confidence intervals for standardized state occupancy probabilities
#'
#' Performs group bootstrap resampling to compute confidence intervals for
#' standardized state occupancy probabilities (SOPs) across all states and time
#' points. Uses g-computation to estimate marginal (population-averaged) SOPs
#' under treatment and control.
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
#' \strong{Handling missing states:} This is the key feature of this function.
#' When a state (e.g., state 3) is missing from a bootstrap sample:
#' \itemize{
#'   \item The model is fit on releveled data (states 1,2,4,5,6 become 1,2,3,4,5)
#'   \item SOPs are predicted for the releveled states
#'   \item Results are mapped back to original state numbering (1,2,3,4,5,6)
#'   \item Missing state 3 is assigned probability 0 at all time points
#' }
#'
#' This ensures that all bootstrap iterations return SOPs for all original states,
#' which is required for computing valid confidence intervals. The zero padding
#' is mathematically correct: if a state never appears in the bootstrap sample,
#' the probability of occupying that state is truly zero.
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
#' @importFrom rsample group_bootstraps analysis
#' @importFrom future.callr callr
#' @importFrom future plan
#' @importFrom furrr future_map
#' @importFrom rms datadist
#' @importFrom stats ave update
#' @importFrom dplyr mutate select row_number bind_rows filter
#' @importFrom tidyr pivot_longer starts_with
#' @importFrom tibble as_tibble
#' @importFrom rlang sym
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
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Check model class
  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  # Bootstrap samples
  resample <- group_bootstraps(
    data,
    group = !!rlang::sym(varnames$id),
    times = n_boot,
    apparent = FALSE
  )

  # Set up workers
  if (is.null(workers)) {
    workers <- parallel::detectCores() - 1
  }

  # Set times if not provided
  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  n_times <- length(times)
  n_states <- length(ylevels)

  # Apply bootstrap in parallel
  if (parallel) {
    plan(callr, workers = workers)
  } else {
    plan("sequential")
  }

  bs_results <- resample |>
    dplyr::mutate(
      boot_id = dplyr::row_number(),
      sops = future_map(splits, \(.x) {
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

        # Relevel factors to handle missing states
        releveled <- relevel_factors_consecutive(
          data = boot_data,
          factor_cols = c("y", varnames$pvarname),
          original_data = data,
          ylevels = ylevels,
          absorb = absorb
        )

        boot_data <- releveled$data
        boot_ylevels <- releveled$ylevels
        boot_absorb <- releveled$absorb
        missing_states <- releveled$missing_levels

        # Get states present in original numbering (for mapping back later)
        states_present <- setdiff(ylevels, as.numeric(missing_states))

        # Update datadist
        dd <- rms::datadist(boot_data)
        assign("dd", dd, envir = .GlobalEnv)
        options(datadist = "dd")

        # Refit model
        m_boot <- tryCatch(
          update(model, data = boot_data),
          error = function(e) {
            warning("Bootstrap iteration failed: ", e$message)
            return(NULL)
          }
        )

        if (is.null(m_boot)) {
          return(NULL)
        }

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
            sop_tx_full[, as.character(original_state)] <- sop_result$sop_tx[,
              i
            ]
            sop_ctrl_full[, as.character(
              original_state
            )] <- sop_result$sop_ctrl[, i]
          }

          return(list(
            sop_tx = sop_tx_full,
            sop_ctrl = sop_ctrl_full
          ))
        } else {
          # All states present - return as is
          return(sop_result)
        }
      })
    )

  plan("sequential")

  # Extract and reshape results into long format
  # Format: boot_id | time | tx | state_1 | state_2 | ... | state_n
  result_list <- list()

  for (i in seq_len(nrow(bs_results))) {
    sop_i <- bs_results$sops[[i]]

    if (!is.null(sop_i)) {
      # Convert treatment SOPs to tibble with time and state columns
      tx_df <- tibble::as_tibble(sop_i$sop_tx)
      colnames(tx_df) <- paste0("state_", ylevels)
      tx_df$time <- times
      tx_df$tx <- 1
      tx_df$boot_id <- bs_results$boot_id[i]

      # Convert control SOPs to tibble with time and state columns
      ctrl_df <- tibble::as_tibble(sop_i$sop_ctrl)
      colnames(ctrl_df) <- paste0("state_", ylevels)
      ctrl_df$time <- times
      ctrl_df$tx <- 0
      ctrl_df$boot_id <- bs_results$boot_id[i]

      # Combine treatment and control
      result_list[[i]] <- dplyr::bind_rows(tx_df, ctrl_df)
    }
  }

  # Combine all bootstrap iterations
  result <- dplyr::bind_rows(result_list) |>
    dplyr::select(boot_id, time, tx, dplyr::everything())

  return(result)
}
