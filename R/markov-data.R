#' Prepare Markov Data for Proportional Odds Modeling
#'
#' Prepares longitudinal Markov trajectory data for proportional odds regression
#' by removing absorbing states and converting variables to appropriate types.
#' Works with both VGAM::vglm and rms::orm models.
#'
#' @param data A data frame containing Markov trajectory data with columns
#'   `y` (current state) and `yprev` (previous state).
#' @param absorbing_state Integer or integer vector. The absorbing state(s) to
#'   filter out (default: 6 for death). Observations where `yprev` equals any of
#'   these values are removed. Can be a single value (e.g., 6) or a vector of
#'   multiple absorbing states (e.g., c(1, 6) for both home discharge and death).
#' @param ordered_response Logical. Should `y` be converted to ordered factor?
#'   (default: TRUE). Required for clm and vglm, not required for orm.
#' @param factor_previous Logical. Should `yprev` be converted to factor?
#'   Defaults to TRUE for backward compatibility and categorical previous-state
#'   effects. Set to FALSE when modeling the previous state linearly or with
#'   nonlinear numeric terms such as `rms::rcs(yprev, 6)`.
#'
#' @return A data frame with modified `y` and `yprev` columns and absorbing
#'   states removed.
#'
#' @details
#' This is a convenience function that performs standard data preparation for
#' proportional odds models on Markov data:
#' - Removes rows where patients were in an absorbing state at the previous time
#' - Converts current state to ordered factor (for cumulative models)
#' - Converts previous state to factor by default (for categorical transition
#'   effects)
#' - Can preserve previous state as numeric for linear or spline transition
#'   effects
#'
#' **Package compatibility**: This function prepares data for:
#' - **VGAM::vglm**: Use `cumulative(parallel = TRUE, reverse = TRUE)` default
#'   datasets.
#' - **rms::orm**: Automatically treats integers as ordered factors.
#' - **ordinal::clm**: Only accepts ordered factors.
#'
#' **Multiple absorbing states**: When modeling scenarios with multiple absorbing
#' states (e.g., discharge home and death), specify them as a vector. For example,
#' `absorbing_state = c(1, 6)` will remove all transitions from either state 1
#' (home) or state 6 (death).
#'
#' @examples
#' \dontrun{
#' # Prepare data with single absorbing state (death only)
#' model_data <- prepare_markov_data(trajectories, absorbing_state = 6)
#'
#' # Prepare data with multiple absorbing states (home and death)
#' model_data <- prepare_markov_data(trajectories, absorbing_state = c(1, 6))
#'
#' # Fit with VGAM (no robust SEs; use bootstrap for inference)
#' library(VGAM)
#' fit <- vglm(y ~ tx + rcs(time, 4) + yprev,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = model_data)
#'
#' # Keep yprev numeric to use one linear or nonlinear previous-state effect
#' model_data <- prepare_markov_data(
#'   trajectories,
#'   absorbing_state = 6,
#'   factor_previous = FALSE
#' )
#' fit <- vglm_markov(
#'   ordered(y) ~ tx + rcs(time, 4) + rcs(yprev, 6),
#'   family = cumulative(parallel = TRUE, reverse = TRUE),
#'   data = model_data
#' )
#'
#' # Fit with rms (supports robust SEs via robcov)
#' library(rms)
#' dd <- datadist(model_data)
#' options(datadist = "dd")
#' fit <- orm(y ~ tx + rcs(time, 4) + yprev, data = model_data, x = TRUE, y = TRUE)
#' fit_robust <- robcov(fit, cluster = model_data$id)
#'
#' # Fit with ordinal
#' library(ordinal)
#' fit <- clm(y ~ tx + time + yprev, data = model_data)
#' }
#'
#' @importFrom stats setNames aggregate model.matrix predict
#' @importFrom utils modifyList
#' @importFrom methods slot
#'
#' @export
prepare_markov_data <- function(
  data,
  absorbing_state = 6,
  ordered_response = TRUE,
  factor_previous = TRUE
) {
  if (!is.null(absorbing_state)) {
    data <- data[!data$yprev %in% absorbing_state, , drop = FALSE]
  }

  if (factor_previous) {
    data$yprev <- factor(data$yprev)
  }

  if (ordered_response) {
    data$y <- ordered(data$y)
  }

  data
}


#' Relevel Factors to Consecutive Integers
#'
#' Handles missing state levels in data by releveling state columns to
#' consecutive integers. This is useful when bootstrap samples or subsets are
#' missing certain levels, which can cause model fitting failures.
#'
#' @param data A data frame containing the data to relevel.
#' @param factor_cols Character vector of column names to relevel. Despite the
#'   historical name, these may be factor or numeric state columns. Numeric
#'   columns stay numeric after releveling.
#' @param original_data Optional data frame containing the original data before
#'   subsetting. Used to determine which levels are present in the full dataset.
#' @param ylevels Optional integer vector of original state levels (e.g., 1:6).
#'   If provided along with absorb, these will be updated to match the new levels.
#' @param absorb Optional integer specifying the absorbing state in the original
#'   levels. Will be updated to the new position if provided.
#'
#' @return A list with components:
#'   - data: The data frame with releveled factors
#'   - ylevels: Updated ylevels (if provided as input), or NULL
#'   - absorb: Updated absorb (if provided as input), or NULL
#'   - missing_levels: Character vector of levels that were missing
#'
#' @details
#' When certain state levels are absent from a dataset (e.g., in bootstrap
#' samples), this function:
#' 1. Identifies which levels are present
#' 2. Creates a mapping from old to new consecutive integers
#' 3. Relevels the specified state columns while preserving factor vs numeric
#'    type
#' 4. Updates ylevels and absorb parameters if provided
#'
#' This is particularly useful for ordered factors representing health states
#' in Markov models, where missing states would otherwise cause model fitting
#' failures.
#'
#' @examples
#' \dontrun{
#' # Relevel y and yprev in bootstrap sample
#' releveled <- relevel_factors_consecutive(
#'   data = boot_data,
#'   factor_cols = c("y", "yprev"),
#'   original_data = full_data,
#'   ylevels = 1:6,
#'   absorb = 6
#' )
#'
#' boot_data <- releveled$data
#' boot_ylevels <- releveled$ylevels
#' boot_absorb <- releveled$absorb
#' }
#'
#' @export
relevel_factors_consecutive <- function(
  data,
  factor_cols = c("y", "yprev"),
  original_data = NULL,
  ylevels = NULL,
  absorb = NULL
) {
  # Determine original levels
  if (!is.null(original_data)) {
    original_levels <- sort(unique(unlist(lapply(
      factor_cols,
      function(col) {
        if (col %in% names(original_data)) {
          if (is.factor(original_data[[col]])) {
            as.numeric(levels(original_data[[col]]))
          } else {
            unique(original_data[[col]])
          }
        }
      }
    ))))
  } else if (!is.null(ylevels)) {
    original_levels <- ylevels
  } else {
    # Infer from data
    original_levels <- sort(unique(unlist(lapply(
      factor_cols,
      function(col) {
        if (col %in% names(data)) {
          unique(data[[col]])
        }
      }
    ))))
  }

  # Get levels present in current data
  states_present <- sort(unique(unlist(lapply(
    factor_cols,
    function(col) {
      if (col %in% names(data)) {
        if (is.factor(data[[col]])) {
          as.numeric(as.character(data[[col]]))
        } else {
          unique(data[[col]])
        }
      }
    }
  ))))

  missing_levels <- setdiff(original_levels, states_present)

  # If no levels are missing, return data as-is
  if (length(missing_levels) == 0) {
    return(list(
      data = data,
      ylevels = ylevels,
      absorb = absorb,
      missing_levels = character(0)
    ))
  }

  # Create mapping from old to new state numbers
  state_mapping <- stats::setNames(seq_along(states_present), states_present)

  # Relevel each state column while preserving the column type. This matters for
  # previous-state predictors that are intentionally numeric, e.g. rcs(yprev, 6).
  for (col in factor_cols) {
    if (col %in% names(data)) {
      prototype <- if (!is.null(original_data) && col %in% names(original_data)) {
        original_data[[col]]
      } else {
        data[[col]]
      }

      mapped <- unname(state_mapping[as.character(data[[col]])])

      if (is.factor(prototype) || is.factor(data[[col]])) {
        data[[col]] <- factor(
          mapped,
          levels = seq_along(states_present),
          ordered = is.ordered(prototype) || is.ordered(data[[col]])
        )
      } else if (is.integer(prototype) || is.integer(data[[col]])) {
        data[[col]] <- as.integer(mapped)
      } else if (is.numeric(prototype) || is.numeric(data[[col]])) {
        data[[col]] <- as.numeric(mapped)
      } else {
        data[[col]] <- as.character(mapped)
      }
    }
  }

  # Update ylevels if provided
  new_ylevels <- if (!is.null(ylevels)) {
    seq_along(states_present)
  } else {
    NULL
  }

  # Update absorbing state if provided and present
  new_absorb <- if (!is.null(absorb)) {
    absorb_present <- absorb[absorb %in% states_present]
    absorb_missing <- setdiff(absorb, states_present)
    if (length(absorb_present) > 0) {
      if (length(absorb_missing) > 0) {
        warning(
          "Absorbing state(s) ",
          paste(absorb_missing, collapse = ", "),
          " not present in data and were dropped."
        )
      }
      unname(state_mapping[as.character(absorb_present)])
    } else {
      warning(
        "Absorbing state(s) ",
        paste(absorb, collapse = ", "),
        " not present in data. Set to NULL."
      )
      NULL
    }
  } else {
    NULL
  }

  return(list(
    data = data,
    ylevels = new_ylevels,
    absorb = new_absorb,
    missing_levels = as.character(missing_levels)
  ))
}
