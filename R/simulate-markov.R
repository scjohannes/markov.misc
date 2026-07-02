#' Simulate Individual Patient Trajectories Using Markov Model
#'
#' This function generates individual patient trajectories over time based on
#' a proportional odds or partial proportional odds model with customizable
#' intercepts, linear predictor function, and baseline patient characteristics.
#'
#' @param baseline_data A data frame containing baseline patient characteristics.
#'   Must include columns: `id`, `yprev` (initial state), and any covariates
#'   needed by `lp_function`. Typically includes `tx`, `age`, `sofa`, etc.
#'   Default is `violet_baseline`, a dataset of 10,000 patients derived from
#'   the VIOLET trial (see `?violet_baseline` for details).
#' @param follow_up_time Integer. Number of time periods to simulate (default: 60).
#' @param intercepts Numeric vector of intercepts for the proportional odds model.
#'   Should have length = (number of states - 1). Default values are from VIOLET
#'   study with 6-state expansion:
#'   c(-9.353, -4.294, -1.389, -0.556, 3.127)
#' @param lp_function A function that calculates the linear predictor for each
#'   patient at each time point. Should accept parameters:
#'   - yprev: previous state
#'   - t: current time
#'   - tx: treatment indicator
#'   - Additional named arguments matching columns in `baseline_data`
#'   - parameter: treatment effect (log odds ratio)
#'   - extra_params: named vector of coefficients
#'   Must return either one numeric value, used for every threshold in a
#'   proportional odds model, or a numeric vector with length equal to
#'   `length(intercepts)`, used as threshold-specific linear predictors in a
#'   partial proportional odds model. Default is `lp_violet`, which implements
#'   the VIOLET study model (see `?lp_violet` for details).
#' @param extra_params Named numeric vector of model coefficients used by
#'   `lp_function`. Default values from VIOLET study include time, spline,
#'   age, sofa, and previous state effects with interactions.
#' @param parameter Numeric. Treatment effect on log odds ratio scale (default: 0,
#'   meaning odds ratio = 1).
#' @param absorbing_states Integer vector of states that are absorbing (once
#'   entered, patients cannot leave). Default is 6 (death).
#' @param seed Integer. Random seed for reproducibility (default: NULL, no seed set).
#' @param covariate_names Character vector of covariate names from `baseline_data`
#'   to pass to `lp_function`. Default is c("age", "sofa", "tx"). The function
#'   will automatically detect and pass these to `lp_function`.
#'
#' @return A data frame with columns:
#'   - id: patient identifier
#'   - time: time point (1 to follow_up_time)
#'   - y: observed state at this time
#'   - yprev: state at previous time point
#'   - all columns from baseline_data
#'
#' @details
#' The function implements a discrete-time Markov model where transition
#' probabilities are determined by a proportional odds or partial proportional
#' odds model. At each time step:
#' 1. Calculate linear predictor for each active patient
#' 2. Apply intercepts to get cumulative probabilities
#' 3. Convert to state probabilities
#' 4. Sample next state for each patient
#' 5. Patients in absorbing states remain there
#'
#' States are assumed to be ordered integers (e.g., 1=Home, 2=Hospital mild,
#' 3=Hospital oxygen, 4=Hospital NIV, 5=Ventilator, 6=Death). By default the
#' model uses a proportional odds structure where higher intercepts represent
#' more severe health states. If `lp_function` returns threshold-specific
#' linear predictors, the induced cumulative probabilities must remain
#' nondecreasing across thresholds; otherwise the function stops because the
#' implied state probabilities would be negative.
#'
#' @examples
#' \dontrun{
#' # Using all defaults (violet_baseline data, lp_violet function, VIOLET params)
#' trajectories <- sim_trajectories_markov(
#'   parameter = log(0.8),  # OR = 0.8 for treatment effect
#'   seed = 12345
#' )
#'
#' # Using default function with custom baseline data
#' custom_baseline <- data.frame(
#'   id = 1:100,
#'   yprev = sample(2:5, 100, replace = TRUE),
#'   tx = rbinom(100, 1, 0.5),
#'   age = rnorm(100, 60, 15),
#'   sofa = rpois(100, 5)
#' )
#'
#' trajectories_custom <- sim_trajectories_markov(
#'   baseline_data = custom_baseline,
#'   follow_up_time = 30,
#'   seed = 12345
#' )
#'
#' # Define a custom linear predictor function
#' my_custom_lp <- function(yprev, t, age, sofa, tx, parameter = 0, extra_params) {
#'   # Custom implementation
#'   tx_effect <- parameter * tx
#'   time_effect <- extra_params["time"] * t
#'   # ... your custom logic ...
#'   return(tx_effect + time_effect)
#' }
#'
#' trajectories_custom_lp <- sim_trajectories_markov(
#'   lp_function = my_custom_lp,
#'   parameter = log(0.75),
#'   seed = 12345
#' )
#'
#' # Partial proportional odds: treatment affects each threshold differently
#' my_ppo_lp <- function(yprev, t, age, sofa, tx, parameter = 0, extra_params) {
#'   base_lp <- lp_violet(
#'     yprev = yprev,
#'     t = t,
#'     age = age,
#'     sofa = sofa,
#'     tx = 0,
#'     parameter = parameter,
#'     extra_params = extra_params
#'   )
#'   base_lp + tx * parameter * c(0.25, 0.5, 1, 1, 1)
#' }
#'
#' trajectories_ppo <- sim_trajectories_markov(
#'   lp_function = my_ppo_lp,
#'   parameter = log(0.75),
#'   seed = 12345
#' )
#' }
#'
#' @seealso
#' \code{\link{lp_violet}} for the default linear predictor function
#'
#' \code{\link{violet_baseline}} for the default baseline dataset
#'
#' @importFrom stats plogis
#'
#' @export
sim_trajectories_markov <- function(
  baseline_data = violet_baseline,
  follow_up_time = 60,
  intercepts = c(-9.353417, -4.294121, -1.389221, -0.555688, 3.127056),
  lp_function = lp_violet,
  extra_params = c(
    "time" = -0.738194,
    "time'" = 0.7464006,
    "age" = 0.010321,
    "sofa" = 0.046901,
    "yprev=1" = -8.518344,
    "yprev=3" = 0,
    "yprev=4" = 1.315332,
    "yprev=5" = 6.576662,
    "yprev=1 * time" = 0,
    "yprev=3 * time" = 0,
    "yprev=4 * time" = 0,
    "yprev=5 * time" = 0
  ),
  parameter = 0,
  absorbing_states = 6,
  seed = NULL,
  covariate_names = c("age", "sofa", "tx")
) {
  # Input validation
  if (!is.data.frame(baseline_data)) {
    stop("baseline_data must be a data frame")
  }

  required_cols <- c("id", "yprev")
  missing_cols <- setdiff(required_cols, names(baseline_data))
  if (length(missing_cols) > 0) {
    stop(
      "baseline_data must contain columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.numeric(intercepts) || length(intercepts) < 1) {
    stop("intercepts must be a numeric vector with at least one element")
  }

  if (!is.function(lp_function)) {
    stop("lp_function must be a function")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Setup
  N <- nrow(baseline_data)
  times <- 1:follow_up_time
  n_states <- length(intercepts) + 1

  # Initialize state matrix
  # Rows = patients, Columns = time points (0 to follow_up_time)
  state_matrix <- matrix(
    as.numeric(baseline_data$yprev),
    nrow = N,
    ncol = length(times) + 1,
    dimnames = list(baseline_data$id, 0:max(times))
  )

  # Detect which covariates from baseline_data need to be passed to lp_function
  available_covariates <- intersect(covariate_names, names(baseline_data))
  if (length(available_covariates) == 0) {
    warning("None of the specified covariate_names found in baseline_data")
  }

  # Main simulation loop
  for (t in times) {
    # Current states
    yprev <- state_matrix[, t]

    # Identify patients not in absorbing states
    active_idx <- which(!yprev %in% absorbing_states)

    # If all patients are in absorbing states, fill remaining time and exit
    if (length(active_idx) == 0) {
      if (t <= max(times)) {
        for (i in (t + 1):ncol(state_matrix)) {
          state_matrix[, i] <- state_matrix[, t]
        }
      }
      break
    }

    # Extract data for active patients
    yprev_active <- yprev[active_idx]
    X_active <- baseline_data[active_idx, , drop = FALSE]

    # Build argument list for lp_function
    # This allows dynamic passing of covariates
    # We pass yprev as a factor to ensure compatibility with lp functions that expect it
    lp_args <- list(
      yprev = factor(yprev_active, levels = 1:n_states)
    )

    # MoreArgs for scalar values that don't vary across patients
    more_args <- list(
      t = t,
      parameter = parameter,
      extra_params = extra_params
    )

    # Add covariates dynamically
    for (cov_name in available_covariates) {
      lp_args[[cov_name]] <- X_active[[cov_name]]
    }

    # Special handling for treatment variable (might be tx_val in function)
    if ("tx" %in% names(X_active) && !"tx" %in% names(lp_args)) {
      lp_args[["tx"]] <- X_active[["tx"]]
    }

    # Calculate linear predictor for each active patient. Each call may return
    # one shared PO value or one value per threshold for PPO simulation.
    lp_list <- do.call(
      mapply,
      c(
        list(FUN = lp_function, SIMPLIFY = FALSE, USE.NAMES = FALSE),
        lp_args,
        list(MoreArgs = more_args)
      )
    )

    n_thresholds <- length(intercepts)
    lp_values <- vapply(
      seq_along(lp_list),
      function(i) {
        normalize_markov_lp(lp_list[[i]], n_thresholds, active_idx[i], t)
      },
      numeric(n_thresholds)
    )
    lp_matrix <- matrix(
      as.numeric(lp_values),
      nrow = length(active_idx),
      ncol = n_thresholds,
      byrow = TRUE
    )

    # Calculate transition probabilities using scalar PO or threshold-specific
    # PPO linear predictors.
    cum_probs <- plogis(
      matrix(
        intercepts,
        nrow = length(active_idx),
        ncol = n_thresholds,
        byrow = TRUE
      ) +
        lp_matrix
    )

    if (n_thresholds > 1 && any(t(apply(cum_probs, 1, diff)) < -1e-12)) {
      stop(
        "lp_function produced threshold-specific linear predictors that ",
        "create decreasing cumulative probabilities. This implies negative ",
        "state probabilities; use threshold-specific effects that preserve ",
        "nondecreasing cumulative probabilities across intercepts."
      )
    }

    # Convert cumulative probabilities to individual state probabilities
    prob_matrix <- cbind(cum_probs, 1) - cbind(0, cum_probs)

    # States in descending order (highest to lowest)
    states_vec <- n_states:1

    # Sample next state for each active patient
    y_new_active <- apply(prob_matrix, 1, function(p) {
      sample(states_vec, size = 1, prob = p)
    })

    # Update state matrix
    # First, carry forward previous state for all patients
    state_matrix[, t + 1] <- state_matrix[, t]
    # Then update only active patients with new states
    state_matrix[active_idx, t + 1] <- y_new_active
  }

  # Convert matrix to long format data frame.
  result <- matrix_to_long(state_matrix, value_name = "y")
  id_key <- as.character(baseline_data$id)
  id_match <- match(result$id, id_key)
  result$id <- baseline_data$id[id_match]
  result$time <- as.integer(result$time)
  result <- left_join_preserve_order(result, baseline_data, by = "id")
  id_order <- match(as.character(result$id), id_key)
  result <- result[order(id_order, result$time), , drop = FALSE]
  result$yprev <- ave(result$y, result$id, FUN = function(x) {
    c(NA, utils::head(x, -1))
  })
  result <- result[result$time > 0, , drop = FALSE] # Remove time 0.
  rownames(result) <- NULL

  return(result)
}

normalize_markov_lp <- function(lp, n_thresholds, active_index, time) {
  if (!is.numeric(lp) || is.matrix(lp) || is.array(lp)) {
    stop("lp_function must return a numeric scalar or numeric vector")
  }

  lp <- as.numeric(lp)
  if (!(length(lp) %in% c(1, n_thresholds))) {
    stop(
      "lp_function must return either length 1 or length(intercepts) (",
      n_thresholds,
      ") for each patient; got length ",
      length(lp),
      " for active patient index ",
      active_index,
      " at time ",
      time
    )
  }

  if (any(!is.finite(lp))) {
    stop("lp_function must return finite numeric values")
  }

  if (length(lp) == 1) {
    rep(lp, n_thresholds)
  } else {
    lp
  }
}
