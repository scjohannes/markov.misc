#' Default Linear Predictor Function for VIOLET-style Simulations
#'
#' A linear predictor function compatible with the VIOLET study parameters and
#' the `violet_baseline` dataset. This function calculates the linear predictor
#' for a proportional odds Markov transition model with effects for treatment,
#' covariates (age, SOFA), time (with spline), and previous state.
#'
#' @param yprev Numeric vector. Previous health state (1-6). State 2 is the
#'   reference level for previous state effects.
#' @param t Numeric scalar. Current time point.
#' @param age Numeric vector. Patient age in years.
#' @param sofa Numeric vector. Sequential Organ Failure Assessment (SOFA) score.
#' @param tx Numeric vector. Treatment indicator (0 = control, 1 = treatment).
#' @param parameter Numeric scalar. Treatment effect on log odds ratio scale
#'   (default: 0, meaning OR = 1).
#' @param extra_params Named numeric vector of model coefficients. Should include:
#'   - "time": Linear time effect
#'   - "time'": Spline knot effect for time > 2
#'   - "age": Age effect per year
#'   - "sofa": SOFA score effect per point
#'   - "yprev=1", "yprev=3", "yprev=4", "yprev=5": Previous state effects
#'     (yprev=2 is reference, so no coefficient needed)
#'   - "yprev=1 * time", "yprev=3 * time", etc.: Time interactions
#'
#' @return Numeric vector of linear predictor values, one per patient.
#'
#' @details
#' The linear predictor is calculated as:
#'
#' LP = β_tx × tx + β_age × age + β_sofa × sofa +
#'      β_time × t + β_time' × max(t-2, 0) +
#'      β_yprev + β_yprev×time × t
#'
#' where:
#' - Treatment effect (β_tx) comes from the `parameter` argument
#' - Covariate effects (β_age, β_sofa) come from `extra_params`
#' - Time has a linear effect plus a spline knot at t=2
#' - Previous state effects are relative to state 2 (reference)
#' - Previous state can have time-varying effects
#'
#' **State 2 as reference**: The model uses state 2 (Hospital - mild) as the
#' reference category for previous state. This means:
#' - yprev=1 (Home): Typically large negative effect (easier transitions)
#' - yprev=2 (Hospital mild): No effect (reference = 0)
#' - yprev=3-5 (Sicker states): Positive effects (harder transitions)
#'
#' **Time spline**: The function uses a linear spline with a knot at t=2,
#' implemented as: β_time × t + β_time' × max(t-2, 0). This allows the time
#' effect to change slope after day 2.
#'
#' @examples
#' # Example with default VIOLET parameters
#' lp_violet(
#'   yprev = c(2, 3, 5),
#'   t = 10,
#'   age = c(65, 70, 55),
#'   sofa = c(5, 8, 6),
#'   tx = c(0, 1, 1),
#'   parameter = log(0.8),  # OR = 0.8 for treatment
#'   extra_params = c(
#'     "time" = -0.738194,
#'     "time'" = 0.7464006,
#'     "age" = 0.010321,
#'     "sofa" = 0.046901,
#'     "yprev=1" = -8.518344,
#'     "yprev=3" = 0,
#'     "yprev=4" = 1.315332,
#'     "yprev=5" = 6.576662,
#'     "yprev=1 * time" = 0,
#'     "yprev=3 * time" = 0,
#'     "yprev=4 * time" = 0,
#'     "yprev=5 * time" = 0
#'   )
#' )
#'
#' \dontrun{
#' # Use with sim_trajectories_markov()
#' trajectories <- sim_trajectories_markov(
#'   baseline_data = violet_baseline,
#'   lp_function = lp_violet,
#'   parameter = log(0.75),  # OR = 0.75
#'   seed = 12345
#' )
#'
#' # Create a custom version with different time interactions
#' my_lp <- function(yprev, t, age, sofa, tx, parameter = 0, extra_params) {
#'   # Start with standard VIOLET LP
#'   base_lp <- lp_violet(yprev, t, age, sofa, tx, parameter, extra_params)
#'
#'   # Add custom modification
#'   # ... your custom code ...
#'
#'   return(base_lp)
#' }
#' }
#'
#' @seealso
#' \code{\link{sim_trajectories_markov}} for using this function in simulations
#'
#' \code{\link{violet_baseline}} for the default baseline dataset
#'
#' @export
lp_violet <- function(yprev, t, age, sofa, tx, parameter = 0, extra_params) {
  # Enforce yprev is a factor
  if (!is.factor(yprev)) {
    yprev <- factor(yprev)
  }

  # 1. Treatment Effect
  tx_effect <- parameter * tx

  # 2. Covariate Effects
  sofa_effect <- extra_params["sofa"] * sofa
  age_effect <- extra_params["age"] * age

  # 3. Time Effects (linear + spline knot at t=2)
  time_effect <- extra_params["time"] * t
  time_spline_effect <- extra_params["time'"] * pmax(t - 2, 0)

  # 4. Previous State Effect (relative to yprev=2, which is the reference)
  # Initialize to zero (reference level)
  yprev_effect <- rep(0, length(yprev))
  yprev_time_effect <- rep(0, length(yprev))

  # For each non-reference state, add its effect
  for (state in unique(yprev)) {
    if (state != 2) {
      # Skip reference level
      state_idx <- which(yprev == state)

      # Main effect of previous state
      yprev_coef_name <- paste0('yprev=', state)
      if (yprev_coef_name %in% names(extra_params)) {
        yprev_effect[state_idx] <- extra_params[[yprev_coef_name]]
      }

      # Time interaction for this previous state
      yprev_time_name <- paste0('yprev=', state, ' * time')
      if (yprev_time_name %in% names(extra_params)) {
        yprev_time_effect[state_idx] <- extra_params[[yprev_time_name]] * t
      }
    }
  }

  # 5. Sum all effects to get the linear predictor
  lp <- as.numeric(
    tx_effect +
      sofa_effect +
      age_effect +
      time_effect +
      time_spline_effect +
      yprev_effect +
      yprev_time_effect
  )

  return(lp)
}
