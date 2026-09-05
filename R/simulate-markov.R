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

    n_thresholds <- length(intercepts)
    if (identical(lp_function, lp_violet)) {
      lp <- do.call(lp_function, c(lp_args, more_args))
      lp_matrix <- matrix(
        lp,
        nrow = length(active_idx),
        ncol = n_thresholds
      )
    } else {
      # Arbitrary user functions retain scalar invocation and validation.
      lp_list <- do.call(
        mapply,
        c(
          list(FUN = lp_function, SIMPLIFY = FALSE, USE.NAMES = FALSE),
          lp_args,
          list(MoreArgs = more_args)
        )
      )
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
    }

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

    if (
      n_thresholds > 1 &&
        any(
          cum_probs[, -1L, drop = FALSE] -
            cum_probs[, -n_thresholds, drop = FALSE] <
            -1e-12
        )
    ) {
      stop(
        "lp_function produced threshold-specific linear predictors that ",
        "create decreasing cumulative probabilities. This implies negative ",
        "state probabilities; use threshold-specific effects that preserve ",
        "nondecreasing cumulative probabilities across intercepts."
      )
    }

    # Convert cumulative probabilities to individual state probabilities
    prob_matrix <- cbind(cum_probs, 1) - cbind(0, cum_probs)

    # Probability columns run from the highest to the lowest state.
    y_new_active <- n_states + 1L - sample_categorical_rows(prob_matrix)

    # Update state matrix
    # First, carry forward previous state for all patients
    state_matrix[, t + 1] <- state_matrix[, t]
    # Then update only active patients with new states
    state_matrix[active_idx, t + 1] <- y_new_active
  }

  # Construct the patient-major long result directly.
  patient_rows <- rep(seq_len(N), each = follow_up_time)
  result <- data.frame(
    id = baseline_data$id[patient_rows],
    time = rep(seq_len(follow_up_time), times = N),
    y = as.vector(t(state_matrix[, -1L, drop = FALSE])),
    check.names = FALSE
  )
  baseline_columns <- setdiff(names(baseline_data), "id")
  result <- cbind(
    result,
    baseline_data[patient_rows, baseline_columns, drop = FALSE]
  )
  result$yprev <- as.vector(t(
    state_matrix[, -(follow_up_time + 1L), drop = FALSE]
  ))
  rownames(result) <- NULL

  return(result)
}

#' Simulate ACTT-2 Markov Trajectories
#'
#' Generates ACTT-2-like eight-state ordinal trajectories from a fitted
#' proportional-odds transition model. The wrapper samples baseline states and
#' treatment assignments, evaluates the model's restricted cubic spline in day,
#' and delegates trajectory generation to [sim_trajectories_markov()].
#'
#' @param n_patients Integer. Number of patients to simulate (default: 1000).
#' @param treatment_prob Numeric. Probability of assignment to Baricitinib plus
#'   Remdesivir, represented by `tx = 1` (default: 0.5). Placebo plus Remdesivir
#'   is represented by `tx = 0`.
#' @param treatment_effect Numeric. Day-1 log-odds coefficient for Baricitinib
#'   plus Remdesivir (`tx = 1`) relative to Placebo plus Remdesivir (`tx = 0`).
#'   The default -0.0657 is the active-arm equivalent of the reported reverse
#'   contrast.
#' @param treatment_effect_decay Nonnegative exponential decay rate per day.
#'   The effect on day `t` is `treatment_effect *
#'   exp(-treatment_effect_decay * (t - 1))`. The default 0 gives a constant
#'   effect.
#' @param seed Integer. Random seed for reproducibility (default: `NULL`).
#' @param ... Additional arguments passed to [sim_trajectories_markov()]. These
#'   can override ACTT-2 defaults such as `follow_up_time`. Overrides are applied
#'   after the wrapper generates its baseline data; supplying `baseline_data`
#'   replaces those generated baseline rows.
#'
#' @return A data frame with the same columns returned by
#'   [sim_trajectories_markov()]. The default output contains 28 rows per
#'   patient and states numbered 1 through 8.
#'
#' @details
#' Initial previous states are sampled from states 4 through 7 with
#' probabilities 0.1374637, 0.5459826, 0.2090997, and 0.1074540. They represent
#' the states immediately before the first generated day-1 outcome. State 8 is
#' absorbing by default.
#'
#' The transition model uses the reported `rms::orm()` thresholds for `y >= 2`
#' through `y >= 8`, reversed to match the cumulative-probability ordering used
#' by [sim_trajectories_markov()]. Shared previous-state and spline coefficients
#' retain their reported order and signs. The treatment effect is supplied
#' separately through `treatment_effect`. The reported coefficients are rounded,
#' so this function reproduces the supplied model summary rather than an
#' unavailable full-precision fit.
#'
#' @examples
#' actt2 <- sim_actt2_markov(n_patients = 20, seed = 123)
#' head(actt2)
#' decaying <- sim_actt2_markov(
#'   n_patients = 20,
#'   treatment_effect = log(0.8),
#'   treatment_effect_decay = -log(0.01) / 14,
#'   seed = 123
#' )
#'
#' @export
sim_actt2_markov <- function(
  n_patients = 1000,
  treatment_prob = 0.5,
  treatment_effect = -0.0657,
  treatment_effect_decay = 0,
  seed = NULL,
  ...
) {
  validate_markov_wrapper_inputs(n_patients, treatment_prob)

  validate_markov_treatment_effect(treatment_effect, treatment_effect_decay)

  n_patients <- as.integer(n_patients)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  actt2 <- actt2_markov_parameters()
  baseline_data <- data.frame(
    id = seq_len(n_patients),
    yprev = sample(
      actt2$baseline_states,
      size = n_patients,
      replace = TRUE,
      prob = actt2$baseline_probabilities
    ),
    tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
  )

  args <- utils::modifyList(
    list(
      baseline_data = baseline_data,
      follow_up_time = 28L,
      intercepts = rev(actt2$orm_intercepts),
      lp_function = actt2_markov_lp,
      extra_params = c(
        actt2$extra_params,
        treatment_effect_decay = treatment_effect_decay
      ),
      parameter = treatment_effect,
      absorbing_states = 8L,
      seed = NULL,
      covariate_names = "tx"
    ),
    list(...)
  )

  do.call(sim_trajectories_markov, args)
}


#' Simulate 60-day ACTT-2 calibrated Markov trajectories
#'
#' Generates 60-day eight-state trajectories from a full proportional-odds
#' Markov model calibrated to 100,000 patients generated by
#' [sim_trajectories_brownian_gap()]. The calibration model was
#' `y ~ tx + yprev + rms::rcs(time, 7)` and used 5,565,851 non-post-absorption
#' transitions.
#'
#' @param n_patients Number of patients to simulate. Must be a positive integer.
#' @param treatment_prob Probability of assignment to the active arm (`tx = 1`).
#'   Must be a numeric scalar between 0 and 1.
#' @param treatment_effect Shared proportional-odds coefficient applied to the
#'   active arm. The default is the fitted null-calibration estimate
#'   -0.0039178706.
#' @param seed Optional random-number seed. The seed is set once before baseline
#'   states and treatment assignments are sampled, and trajectory generation
#'   continues from the same random-number stream.
#' @param ... Advanced overrides passed to [sim_trajectories_markov()]. These can
#'   override defaults such as `follow_up_time`. Supplying `baseline_data`
#'   replaces the generated baseline rows.
#'
#' @return A data frame with the same columns returned by
#'   [sim_trajectories_markov()]. The default output contains 60 rows per
#'   patient and states numbered 1 through 8.
#'
#' @details
#' Baseline states 4 through 7 are sampled with probabilities 0.04797, 0.46254,
#' 0.35816, and 0.13133, estimated from the 100,000-patient Brownian-gap
#' calibration cohort. State 8 is absorbing.
#'
#' The fitted restricted cubic spline uses knots at days 2, 11, 20, 30, 39, 49,
#' and 59. Printed `rms::orm()` thresholds are reversed only at the boundary to
#' [sim_trajectories_markov()]; all fitted slopes retain their reported signs.
#' Treatment is coded as `tx = 0` for placebo/control and `tx = 1` for active
#' treatment. Because higher states are worse, a negative treatment coefficient
#' is protective.
#'
#' This full proportional-odds transition model is a compact approximation to
#' the Brownian-gap source process, not an occupancy-equivalent replacement. In
#' the calibration check, source and wrapper day-60 mortality were 9.58 and
#' 15.55 percent, respectively, and the maximum absolute state-by-day difference
#' was 25.18 percentage points. Use [sim_actt2_brownian()] when the source
#' Brownian dynamics are required. Follow-up beyond day 60 is extrapolation
#' outside the fitted time support.
#'
#' @examples
#' trajectories <- sim_actt2_markov_60day(n_patients = 20, seed = 123)
#' head(trajectories)
#'
#' @export
sim_actt2_markov_60day <- function(
  n_patients = 1000,
  treatment_prob = 0.5,
  treatment_effect = -0.00344081587038234,
  seed = NULL,
  ...
) {
  validate_markov_wrapper_inputs(n_patients, treatment_prob)
  validate_markov_treatment_effect(treatment_effect, 0)

  n_patients <- as.integer(n_patients)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  calibration <- actt2_markov_60day_parameters()
  baseline_data <- data.frame(
    id = seq_len(n_patients),
    yprev = sample(
      calibration$baseline_states,
      size = n_patients,
      replace = TRUE,
      prob = calibration$baseline_probabilities
    ),
    tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
  )

  args <- utils::modifyList(
    list(
      baseline_data = baseline_data,
      follow_up_time = 60L,
      intercepts = rev(calibration$orm_intercepts),
      lp_function = actt2_markov_60day_lp,
      extra_params = calibration$extra_params,
      parameter = treatment_effect,
      absorbing_states = 8L,
      seed = NULL,
      covariate_names = "tx"
    ),
    list(...)
  )

  do.call(sim_trajectories_markov, args)
}


#' Simulate ACTT-1 ordinal patient trajectories
#'
#' @description
#' Simulates individual 8-state ACTT-1 trajectories from a simplified Bayesian
#' Markov proportional-odds model. The wrapper uses posterior mean coefficients,
#' the reported baseline-state counts, and the ACTT-2 restricted cubic spline
#' knots because the original ACTT-1 knot locations are unavailable.
#'
#' @param n_patients Number of patients to simulate. Must be a positive integer.
#' @param treatment_prob Probability of assignment to Remdesivir. Must be a
#'   numeric scalar between 0 and 1.
#' @param treatment_effect Day-1 proportional-odds coefficient applied to
#'   Remdesivir (`tx = 1`) relative to placebo (`tx = 0`). The default -0.0657
#'   matches [sim_actt2_markov()].
#' @param treatment_effect_decay Nonnegative exponential decay rate per day.
#'   The effect on day `t` is `treatment_effect *
#'   exp(-treatment_effect_decay * (t - 1))`. The default 0 gives a constant
#'   effect.
#' @param seed Optional random-number seed. The seed is set once before baseline
#'   states and treatment assignments are sampled, and trajectory generation
#'   continues from the same random-number stream.
#' @param ... Advanced overrides passed to [sim_trajectories_markov()]. These
#'   can override defaults such as `follow_up_time`. Overrides are applied after
#'   the wrapper generates its baseline data; supplying `baseline_data` replaces
#'   those generated baseline rows.
#'
#' @return A data frame with the same columns returned by
#'   [sim_trajectories_markov()]. The default output contains 28 rows per
#'   patient and states numbered 1 through 8.
#'
#' @details
#' Initial previous states are sampled from states 4 through 7 according to the
#' reported counts 138, 435, 193, and 285. These counts are normalized over the
#' 1,051 patients with a reported state and represent states immediately before
#' the first generated day-1 outcome. State 8 is absorbing by default.
#'
#' Treatment is coded as `tx = 0` for placebo and `tx = 1` for Remdesivir. The
#' user-specified effect replaces the fitted ACTT-1 treatment main effect and
#' treatment-by-day spline interactions. Age and sex are omitted because the
#' required joint baseline distribution is unavailable or the fitted effect was
#' judged negligible.
#'
#' The fitted constrained partial proportional-odds time deviation is retained.
#' For the cumulative cutoff `y >= k`, it adds `-0.0188 * day * k` to the linear
#' predictor. This threshold-specific term is essential for reproducing ACTT-1
#' mortality over follow-up even though its printed coefficient is small.
#'
#' The reported thresholds for `y >= 2` through `y >= 8` are reversed only when
#' passed to [sim_trajectories_markov()]; the shared previous-state and time
#' coefficients retain their reported signs. Coefficients are rounded, so this
#' function reproduces the supplied model summary rather than an unavailable
#' full-precision fit.
#'
#' @examples
#' actt1 <- sim_actt1_markov(n_patients = 20, seed = 123)
#' head(actt1)
#' decaying <- sim_actt1_markov(
#'   n_patients = 20,
#'   treatment_effect = log(0.8),
#'   treatment_effect_decay = -log(0.01) / 14,
#'   seed = 123
#' )
#'
#' @export
sim_actt1_markov <- function(
  n_patients = 1000,
  treatment_prob = 0.5,
  treatment_effect = -0.0657,
  treatment_effect_decay = 0,
  seed = NULL,
  ...
) {
  validate_markov_wrapper_inputs(n_patients, treatment_prob)
  validate_markov_treatment_effect(treatment_effect, treatment_effect_decay)

  n_patients <- as.integer(n_patients)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  actt1 <- actt1_markov_parameters()
  baseline_data <- data.frame(
    id = seq_len(n_patients),
    yprev = sample(
      actt1$baseline_states,
      size = n_patients,
      replace = TRUE,
      prob = actt1$baseline_probabilities
    ),
    tx = stats::rbinom(n_patients, size = 1L, prob = treatment_prob)
  )

  args <- utils::modifyList(
    list(
      baseline_data = baseline_data,
      follow_up_time = 28L,
      intercepts = rev(actt1$orm_intercepts),
      lp_function = actt1_markov_lp,
      extra_params = c(
        actt1$extra_params,
        treatment_effect_decay = treatment_effect_decay
      ),
      parameter = treatment_effect,
      absorbing_states = 8L,
      seed = NULL,
      covariate_names = "tx"
    ),
    list(...)
  )

  do.call(sim_trajectories_markov, args)
}


actt1_markov_parameters <- function() {
  baseline_counts <- c(138, 435, 193, 285)

  list(
    baseline_states = 4:7,
    baseline_probabilities = baseline_counts / sum(baseline_counts),
    orm_intercepts = c(
      -3.3535,
      -3.5896,
      -3.8395,
      -6.0379,
      -9.7465,
      -12.3361,
      -18.8334
    ),
    extra_params = c(
      "day" = -0.1110,
      "day'" = 1.6444,
      "day''" = -3.1537,
      "day'''" = 1.2387,
      "day''''" = 1.1910,
      "yprev=1" = 0,
      "yprev=2" = -0.0402,
      "yprev=3" = 3.9042,
      "yprev=4" = 5.0754,
      "yprev=5" = 8.0171,
      "yprev=6" = 11.5073,
      "yprev=7" = 16.1963,
      "day x f(y)" = -0.0188
    )
  )
}


actt1_markov_lp <- function(yprev, t, tx, parameter, extra_params) {
  basis <- actt_markov_rcs_basis(t)
  day_names <- c("day", "day'", "day''", "day'''", "day''''")
  time_effect <- drop(basis %*% extra_params[day_names])
  previous_state_effect <- extra_params[
    paste0("yprev=", as.character(yprev))
  ]
  treatment_effect <- decayed_markov_treatment_effect(
    parameter,
    t,
    extra_params
  )
  shared_effect <- time_effect + previous_state_effect + treatment_effect * tx
  cutoff_scores <- rev(2:8)
  partial_po_effect <- extra_params[["day x f(y)"]] * t * cutoff_scores

  as.numeric(shared_effect + partial_po_effect)
}


actt2_markov_parameters <- function() {
  list(
    baseline_states = 4:7,
    baseline_probabilities = c(0.1374637, 0.5459826, 0.2090997, 0.1074540),
    orm_intercepts = c(
      -1.7503,
      -2.4174,
      -2.6911,
      -4.8406,
      -8.9716,
      -12.6833,
      -19.5542
    ),
    extra_params = c(
      "day" = -0.3376,
      "day'" = 2.6884,
      "day''" = -4.9025,
      "day'''" = 1.4712,
      "day''''" = 2.1778,
      "yprev=1" = 0,
      "yprev=2" = -0.4585,
      "yprev=3" = 3.5573,
      "yprev=4" = 4.5199,
      "yprev=5" = 7.4148,
      "yprev=6" = 11.6247,
      "yprev=7" = 16.8083
    )
  )
}


actt2_markov_60day_parameters <- function() {
  time_knots <- c(2, 11, 20, 30, 39, 49, 59)
  extra_params <- structure(
    c(
      "time" = -0.0551480786641631,
      "time'" = 0.705650641016978,
      "time''" = -1.89732599383498,
      "time'''" = 1.72332047124566,
      "time''''" = -0.72372165657455,
      "time'''''" = 0.324725461513862,
      "yprev=1" = 0,
      "yprev=2" = 7.78190841234478,
      "yprev=3" = 9.7553019571941,
      "yprev=4" = 10.9217897817441,
      "yprev=5" = 12.789303156482,
      "yprev=6" = 15.4186523641647,
      "yprev=7" = 19.3967392254069
    ),
    time_knots = time_knots
  )

  list(
    baseline_states = 4:7,
    baseline_probabilities = c(0.04797, 0.46254, 0.35816, 0.13133),
    time_knots = time_knots,
    treatment_effect = -0.00344081587038234,
    drift_start_lambda = 6.026099,
    n_fit_transitions = 5590720L,
    orm_intercepts = c(
      -6.48590352305606,
      -8.98334851035646,
      -10.2826703630608,
      -11.0388174995576,
      -14.1517663647389,
      -16.9418587523501,
      -22.8326631030386
    ),
    extra_params = extra_params,
    source_occupancy_anchors = data.frame(
      time = c(1, 7, 14, 28, 42, 60),
      state1 = c(0.00026, 0.10199, 0.49575, 0.74997, 0.82855, 0.89781),
      state8 = c(0, 0.00131, 0.03383, 0.09147, 0.09583, 0.09583)
    )
  )
}


actt2_markov_60day_lp <- function(yprev, t, tx, parameter, extra_params) {
  basis <- markov_rcs_basis(t, attr(extra_params, "time_knots"))
  time_names <- c(
    "time",
    "time'",
    "time''",
    "time'''",
    "time''''",
    "time'''''"
  )
  time_effect <- drop(basis %*% extra_params[time_names])
  previous_state_effect <- extra_params[
    paste0("yprev=", as.character(yprev))
  ]

  as.numeric(time_effect + previous_state_effect + parameter * tx)
}


actt2_markov_lp <- function(yprev, t, tx, parameter, extra_params) {
  basis <- actt_markov_rcs_basis(t)
  day_names <- c("day", "day'", "day''", "day'''", "day''''")
  time_effect <- drop(basis %*% extra_params[day_names])
  previous_state_effect <- extra_params[
    paste0("yprev=", as.character(yprev))
  ]

  treatment_effect <- decayed_markov_treatment_effect(
    parameter,
    t,
    extra_params
  )

  as.numeric(time_effect + previous_state_effect + treatment_effect * tx)
}


actt_markov_rcs_basis <- function(day) {
  markov_rcs_basis(day, c(2, 6, 11, 16, 21, 27))
}


markov_rcs_basis <- function(day, knots) {
  n_knots <- length(knots)
  last_knot <- knots[n_knots]
  penultimate_knot <- knots[n_knots - 1L]
  scale <- (last_knot - knots[1L])^(2 / 3)
  nonlinear <- matrix(0, nrow = length(day), ncol = n_knots - 2L)

  for (column in seq_len(n_knots - 2L)) {
    nonlinear[, column] <-
      pmax((day - knots[column]) / scale, 0)^3 +
      ((penultimate_knot - knots[column]) *
        pmax((day - last_knot) / scale, 0)^3 -
        (last_knot - knots[column]) *
          pmax((day - penultimate_knot) / scale, 0)^3) /
        (last_knot - penultimate_knot)
  }

  cbind(day, nonlinear)
}


validate_markov_wrapper_inputs <- function(n_patients, treatment_prob) {
  if (
    length(n_patients) != 1L ||
      !is.numeric(n_patients) ||
      !is.finite(n_patients) ||
      n_patients <= 0 ||
      n_patients > .Machine$integer.max ||
      n_patients != as.integer(n_patients)
  ) {
    stop(simpleError(
      "n_patients must be a positive integer",
      call = sys.call(-1)
    ))
  }

  if (
    length(treatment_prob) != 1L ||
      !is.numeric(treatment_prob) ||
      !is.finite(treatment_prob) ||
      treatment_prob < 0 ||
      treatment_prob > 1
  ) {
    stop(simpleError(
      "treatment_prob must be a numeric scalar between 0 and 1",
      call = sys.call(-1)
    ))
  }
}


validate_markov_treatment_effect <- function(
  treatment_effect,
  treatment_effect_decay
) {
  if (
    length(treatment_effect) != 1L ||
      !is.numeric(treatment_effect) ||
      !is.finite(treatment_effect)
  ) {
    stop(simpleError(
      "treatment_effect must be a finite numeric scalar",
      call = sys.call(-1)
    ))
  }

  if (
    length(treatment_effect_decay) != 1L ||
      !is.numeric(treatment_effect_decay) ||
      !is.finite(treatment_effect_decay) ||
      treatment_effect_decay < 0
  ) {
    stop(simpleError(
      "treatment_effect_decay must be a nonnegative finite numeric scalar",
      call = sys.call(-1)
    ))
  }
}


decayed_markov_treatment_effect <- function(parameter, time, extra_params) {
  decay <- if ("treatment_effect_decay" %in% names(extra_params)) {
    extra_params[["treatment_effect_decay"]]
  } else {
    0
  }

  parameter * exp(-decay * (time - 1))
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
