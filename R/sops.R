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
  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  # Check data
  if (is.null(data)) {
    stop(
      "Provide data used for model fitting."
    )
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
