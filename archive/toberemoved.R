#' Calculate State Occupation Probabilities for First-Order Markov Models
#'
#' Estimates state occupation probabilities over time by iterating a transition
#' matrix derived from a fitted model object (e.g., `vglm`, `rms`, or `rmsb`).
#' This function supports linear time iteration as well as complex, non-linear
#' time specifications (e.g., splines) via a covariate lookup table.
#'
#' @param object A fitted model object. Supported classes include:
#'   \code{"lrm"}, \code{"orm"} (from package `rms`),
#'   \code{"blrm"} (from package `rmsb`), and
#'   \code{"vglm"}, \code{"vgam"} (from package `VGAM`).
#' @param data A data frame containing the baseline covariates for the prediction.
#'   Should generally contain a single row representing the subject profile.
#'   Must contain the variable specified in \code{tvarname} (usually initialized to 0 or start time).
#' @param times A numeric vector of time points to iterate over.
#'   The function calculates probabilities for each time point in this sequence.
#' @param ylevels A character vector defining the names of the outcome levels (states).
#'   These must match the levels used in the fitted `object`.
#' @param absorb (Optional) A character vector of absorbing states (states from which
#'   transitions out are impossible). Defaults to \code{NULL}.
#' @param tvarname A character string specifying the name of the time variable in the
#'   model formula. Defaults to \code{"time"}.
#' @param pvarname A character string specifying the name of the previous state variable
#'   used in the model formula. Defaults to \code{"yprev"}.
#' @param gap (Optional) A character string specifying the name of the variable representing
#'   the time gap (delta time) between observations, if used in the model. Defaults to \code{NULL}.
#' @param t_covs (Optional) A data frame used for non-linear time handling (e.g., splines
#'   or basis functions).
#'   \itemize{
#'     \item **Structure:** The number of rows in \code{t_covs} must exactly match the length of \code{times}.
#'     \item **Columns:** Column names must match the specific basis variables used in the model formula (e.g., \code{t1}, \code{t2}).
#'     \item **Usage:** At step \code{i}, the values from the \code{i}-th row of \code{t_covs} are injected into the prediction data.
#'   }
#'   This allows for `vglm` models where basis functions are supplied as separate columns rather than calculated on the fly.
#'
#' @details
#' The function iterates through the vector \code{times}. At each time point \eqn{t}, it:
#' 1. Updates the time variable (and `t_covs` if supplied).
#' 2. Predicts the transition probability matrix \eqn{T(t)} based on the covariates and the previous state.
#' 3. Updates the state occupation probability vector \eqn{P(t)} using the matrix multiplication:
#'    \eqn{P(t) = P(t-1) \times T(t)}.
#'
#' If the model is Bayesian (\code{rmsb}), the output preserves the posterior draws.
#'
#' @return
#' An array of state probabilities.
#' \itemize{
#'   \item **Frequentist fit:** A matrix of dimension \code{[length(times) x length(ylevels)]}.
#'   \item **Bayesian fit:** An array of dimension \code{[n_draws x length(times) x length(ylevels)]}.
#' }
#'
#' @noRd
#' @keywords internal
soprob_markov_old <- function(
  object,
  data,
  times,
  ylevels,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  gap = NULL,
  t_covs = NULL
) {
  # --- Initial Checks ---
  cl <- class(object)[1]
  ftypes <- c(
    lrm = "rms",
    orm = "rms",
    blrm = "rmsb",
    vglm = "vgam",
    vgam = "vgam"
  )
  ftype <- ftypes[cl]
  if (is.na(ftype)) {
    stop(paste(
      "object must be a fit from one of:",
      paste(ftypes, collapse = " ")
    ))
  }

  # --- Define Prediction Function ---
  prd <- switch(
    ftype,
    rms = function(obj, data) predict(obj, data, type = "fitted.ind"),
    vgam = function(obj, data) VGAM::predict(obj, data, type = "response"),
    rmsb = function(obj, data) {
      predict(obj, data, type = "fitted.ind", posterior.summary = "all")
    }
  )

  if (pvarname %nin% names(data)) {
    stop(paste(pvarname, "is not in data"))
  }
  if (length(absorb) && (pvarname %in% absorb)) {
    stop("initial state cannot be an absorbing state")
  }

  # --- Check Time Covariates (New Logic) ---
  if (!is.null(t_covs)) {
    if (!is.data.frame(t_covs)) {
      stop("'t_covs' must be a data frame")
    }
    if (nrow(t_covs) != length(times)) {
      stop("The number of rows in 't_covs' must match the length of 'times'")
    }
  }

  # --- Setup Output Arrays ---
  nd <- if (ftype == "rmsb" && length(object$draws)) nrow(object$draws) else 0
  if ((nd == 0) != (ftype != "rmsb")) {
    stop("model fit inconsistent with having posterior draws")
  }

  k <- length(ylevels)
  s <- length(times)

  P <- if (nd == 0) {
    array(
      NA,
      c(s, k),
      dimnames = list(as.character(times), as.character(ylevels))
    )
  } else {
    array(
      NA,
      c(nd, s, k),
      dimnames = list(
        paste("draw", 1:nd),
        as.character(times),
        as.character(ylevels)
      )
    )
  }

  # --- Initialize Time 1 ---
  data[[tvarname]] <- times[1]
  if (length(gap)) {
    data[[gap]] <- times[1]
  }

  # Inject basis functions for Time 1
  if (!is.null(t_covs)) {
    for (nm in names(t_covs)) {
      data[[nm]] <- t_covs[1, nm]
    }
  }

  data <- as.data.frame(data)

  # --- Prediction at Time 1 ---
  p <- prd(object, data)
  if (nd == 0) {
    P[1, ] <- p
  } else {
    P[, 1, ] <- p
  }

  # --- Setup Transition Matrices ---
  rnameprev <- paste("t-1", ylevels)
  rnameprevna <- paste("t-1", setdiff(ylevels, absorb))
  if (length(absorb)) {
    rnamepreva <- paste("t-1", absorb)
    cnamea <- paste("t", absorb)
  }

  cp <- matrix(
    0,
    nrow = k,
    ncol = k,
    dimnames = list(rnameprev, paste("t", ylevels))
  )
  if (length(absorb)) {
    cp[cbind(rnamepreva, cnamea)] <- 1
  }

  # --- Prepare Expansion Data ---
  data <- as.list(data)
  yna <- setdiff(ylevels, absorb)
  data[[pvarname]] <- yna
  edata <- expand.grid(data) # creates rows for every previous state

  # --- Iterate Through Time ---
  for (it in 2:s) {
    # Update standard time variable
    edata[[tvarname]] <- times[it]

    # Update gap if necessary
    if (length(gap)) {
      edata[[gap]] <- times[it] - times[it - 1]
    }

    # NEW: Update basis function variables for current time 'it'
    if (!is.null(t_covs)) {
      for (nm in names(t_covs)) {
        # We assign the scalar value from the lookup table to the whole column in edata
        edata[[nm]] <- t_covs[it, nm]
      }
    }

    # Calculate transition probabilities
    pp <- prd(object, edata)

    # Matrix Multiplication (Markov Step)
    if (nd == 0) {
      cp[rnameprevna, ] <- pp
      P[it, ] <- t(cp) %*% P[it - 1, ] # State t = TransMatrix * State t-1
    } else {
      for (idraw in 1:nd) {
        cp[rnameprevna, ] <- pp[idraw, , ]
        P[idraw, it, ] <- t(cp) %*% P[idraw, it - 1, ]
      }
    }
  }

  P
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
#' @noRd
#' @keywords internal
standardize_sops_old <- function(
  model,
  data = NULL,
  times = NULL,
  ylevels = factor(1:6),
  absorb = 6,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx"),
  t_covs = NULL
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
  }
  # Check that yprev is a factor
  if (!is.factor(data[[varnames$pvarname]])) {
    warning(
      "Variable '",
      varnames$pvarname,
      "' should be a factor but is not. Converting to factor."
    )
    data[[varnames$pvarname]] <- factor(data[[varnames$pvarname]])
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
    sop_full_tx <- soprob_markov_ord_m(
      object = model,
      data = X_tx,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = varnames$tvarname,
      pvarname = varnames$pvarname,
      gap = 1,
      t_covs = t_covs
    )

    # Predict SOPs under control
    sop_full_ctrl <- soprob_markov_ord_m(
      object = model,
      data = X_ctrl,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = varnames$tvarname,
      pvarname = varnames$pvarname,
      gap = 1,
      t_covs = t_covs
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


#' Compute time in target state(s) for treatment groups
#'
#' Calculates the expected time spent in specified target state(s) for treatment
#' and control groups using a proportional odds Markov model. Uses an existing
#' fitted model and state occupancy probabilities to estimate the total time in
#' the target state(s).
#'
#' @param model A fitted model object from \code{rms::orm} or \code{VGAM::vglm}.
#'   The model is used as-is without refitting.
#' @param data A data frame containing patient baseline_data. If NULL, attempts
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
#' data. For bootstrap inference, see \code{taooh_bootstrap2}, which handles
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

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Handle missing states in the data (change ylevels and absorb depending on the data)
# - Expand to support vgam models in addition to orm

#' @keywords internal
#' @export

taooh <- function(
  model,
  data = NULL,
  times = NULL,
  ylevels = as.character(1:6),
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
      "' should be a factor but is not. Converting to factor."
    )
    data[[varnames$pvarname]] <- factor(data[[varnames$pvarname]])
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

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)
# - Handle missing states in the data (change ylevels and absorb depending on the data)
# - Expand to support vgam models in addition to orm

#' @noRd
#' @keywords internal
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

# Planned extensions:
# - If beta == TRUE, stop the function and return beta treatment (for beta bootstrapping)

#' @keywords internal
#' @export

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
  result <- future_map_dbl(
    as.data.frame(sample_IDs),
    function(x) {
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
    },
    .options = furrr_options(packages = c("rms", "markov.misc"))
  )

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
