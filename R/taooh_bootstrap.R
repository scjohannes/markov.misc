#' Compute time in target state(s) for treatment groups
#'
#' Calculates the expected time spent in specified target state(s) for treatment
#' and control groups using a proportional odds Markov model. The function fits
#' an ordinal regression model and uses state occupancy probabilities to estimate
#' the total time in the target state(s).
#'
#' @param data A data frame containing patient trajectory data, or an rsample
#'   splits object from bootstrap resampling
#' @param formula Model formula for orm (ordinal regression model)
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
#'   \item Fits a proportional odds model to the trajectory data
#'   \item Computes state occupancy probabilities for each individual at each time point
#'   \item Sums probabilities across target states and time points
#'   \item Averages within treatment groups to get expected time in target state(s)
#' }
#'
#' When data is an rsample splits object (from bootstrap), it extracts the
#' analysis set and creates unique IDs for each resample. When data is a regular
#' data frame, it uses the IDs as provided.
#'
#' @keywords time alive out of hospital state occupancy
#'
#' @importFrom rms orm
#' @importFrom Hmisc soprobMarkovOrdm
#' @importFrom rsample analysis
#' @importFrom stats ave
#' @importFrom stats aggregate
#'
#' @examples
#' \dontrun{
#' # Calculate time at home (state 1) for a dataset
#' result <- taooh(
#'   data = my_data,
#'   formula = y ~ tx + yprev + time,
#'   ylevels = 1:6,
#'   absorb = 6,
#'   target_states = 1
#' )
#'
#' # Calculate time in states 1 or 2 combined
#' result <- taooh(
#'   data = my_data,
#'   formula = y ~ tx + yprev + time,
#'   target_states = c(1, 2)
#' )
#' }
#' @export

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Handle missing states in the data (change ylevels and absorb depanding on the data)
# - Use different bootstrapping procedure.
#     - Up to now, some resamples contain less unique IDs than the original data.
#       Unclear, why this can happen.

taooh <- function(
  data,
  formula,
  times = NULL,
  ylevels = 1:6,
  absorb = 6,
  target_states = 1,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Check if data is a splits object (bootstrap) or regular data frame
  if (inherits(data, "rsplit")) {
    data <- analysis(data)

    # Define new id variable (unique to each bootstrap draw)
    data$block_count <- ave(
      data[[varnames$id]],
      data[[varnames$id]],
      data[[varnames$tvarname]],
      FUN = seq_along
    )
    data$new_id <- factor(paste(
      data[[varnames$id]],
      data$block_count,
      sep = "_"
    ))
    id_var <- "new_id"
  } else {
    # Use original ID for regular data frame
    id_var <- varnames$id
  }

  # Set times if not provided
  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  # 1. Fit model
  m <- orm(formula, data = data)

  # 2. Generate covariate data.frame to predict on
  X <- data[!duplicated(data[[id_var]]), ] # first row of each individual

  # 3. Run soprobMarkovOrdm to get state probability predictions for each individual
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
      object = m,
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

  # 4. Indicate treatment and control ids
  id_tx <- aggregate(
    data[[varnames$tx]],
    by = list(data[[id_var]]),
    FUN = unique,
    simplify = TRUE
  )

  # 5. Compute time in target state(s) by treatment group
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
#' @param data A data frame containing the patient data
#' @param n_boot Number of bootstrap samples
#' @param formula Model formula for orm
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
#' within-patient correlation structure. Parallelization is handled via
#' future.callr::callr strategy, which provides better isolation and stability
#' than multisession.
#'
#' @keywords bootstrap time alive and out of hospital
#'
#' @importFrom rsample group_bootstraps
#' @importFrom future.callr callr
#' @importFrom future plan
#' @importFrom furrr future_map
#' @importFrom stats as.formula
#'
#' @examples
#' \dontrun{
#' # Bootstrap for time at home (state 1)
#' bs_results <- taooh_bootstrap(
#'   data = my_data,
#'   n_boot = 1000,
#'   formula = y ~ tx + yprev + time,
#'   target_states = 1
#' )
#'
#' # Bootstrap for time in states 1 or 2
#' bs_results <- taooh_bootstrap(
#'   data = my_data,
#'   n_boot = 1000,
#'   formula = y ~ tx + yprev + time,
#'   target_states = c(1, 2)
#' )
#' }
#' @export

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Handle missing states in the data (change ylevels and absorb depending on the data)
# - Use different bootstrapping procedure.
#     - Up to now, some resamples contain less unique IDs than the original data.
#       Unclear, why this can happen.

taooh_bootstrap <- function(
  data,
  n_boot,
  formula,
  workers = NULL,
  parallel = TRUE,
  ylevels = 1:6,
  absorb = 6,
  times = NULL,
  target_states = 1,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  # Bootstrap samples
  resample <- group_bootstraps(
    data,
    group = id,
    times = n_boot,
    apparent = FALSE
  )

  # Arguments
  formula <- as.formula(formula)
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
    plan(sequential)
  }

  bs_SOP <- resample |>
    dplyr::mutate(
      models = future_map(splits, \(.x) {
        taooh(
          .x,
          formula = formula,
          times = times,
          ylevels = ylevels,
          absorb = absorb,
          target_states = target_states,
          varnames = varnames
        )[1]
      })
    )

  plan(sequential)

  return(bs_SOP)
}
