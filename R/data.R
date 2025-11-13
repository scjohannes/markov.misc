#' Default Baseline Dataset for Markov Simulations
#'
#' A dataset containing baseline patient characteristics derived from the VIOLET
#' trial data. This can be used as the default `baseline_data` argument in
#' `sim_trajectories_markov()` for simulating patient trajectories.
#'
#' @format A data frame with 250,000 rows and 5 variables:
#' \describe{
#'   \item{id}{Patient identifier (1 to 250,000)}
#'   \item{yprev}{Initial health state (2-5):
#'     \itemize{
#'       \item 2 = Hospital (mild)
#'       \item 3 = Hospital with oxygen
#'       \item 4 = Hospital with non-invasive ventilation (NIV)
#'       \item 5 = Mechanical ventilation
#'     }
#'   }
#'   \item{age}{Age in years}
#'   \item{sofa}{Sequential Organ Failure Assessment (SOFA) score at baseline}
#'   \item{tx}{Treatment assignment (0 = control, 1 = treatment)}
#' }
#'
#' @details
#' This dataset is derived from the `simlongord` dataset in the Hmisc package,
#' which contains simulated data based on the VIOLET trial (a randomized trial
#' of Vitamin C, Thiamine, and Steroids in patients with sepsis and acute
#' respiratory failure).
#'
#' **Data Processing:**
#' - Original 4-state model expanded to 6-state model by splitting the "Hospital"
#'   state into three severity levels (mild, nasal xygen, NIV)
#' - Proportion of mechanically ventilated patients reduced from ~30% to ~5%
#'   through resampling to better represent typical respiratory illness severity
#' - Baseline covariates (age, SOFA score) retained from original VIOLET data
#' - Treatment assignment randomly generated with equal probability
#'
#' **State Definitions:**
#' The model uses 6 ordered health states:
#' - State 1: Home (discharged)
#' - State 2: Hospital - mild
#' - State 3: Hospital - requiring nasal oxygen
#' - State 4: Hospital - requiring non-invasive ventilation (NIV)
#' - State 5: Hospital - requiring mechanical ventilation
#' - State 6: Death (absorbing state)
#'
#' @source
#' Based on `simlongord` dataset from the Hmisc package, created by Frank Harrell.
#' See `Hmisc::getHdata(simlongord)` for the original data.
#'
#' @references
#' Harrell, F. E. Jr. (2001) *Regression Modeling Strategies*.
#' Springer-Verlag, New York.
#'
#' @examples
#' # View structure of baseline data
#' head(violet_baseline)
#'
#' # Distribution of initial states
#' table(violet_baseline$yprev)
#'
#' \dontrun{
#' # Use as default baseline in simulation
#' library(markov.misc)
#'
#' # Define linear predictor function (simplified example)
#' my_lp <- function(yprev, t, age, sofa, tx, parameter = 0, extra_params) {
#'   # ... function body ...
#' }
#'
#' trajectories <- sim_trajectories_markov(
#'   baseline_data = violet_baseline,
#'   follow_up_time = 60,
#'   lp_function = my_lp
#' )
#' }
#'
"violet_baseline"
