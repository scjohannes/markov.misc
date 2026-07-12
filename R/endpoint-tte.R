#' Convert State Trajectory Data to Time-to-event data
#'
#' This function calculates time-to-event for usage in Cox or recurrent event
#' models.
#'
#' @param data A data frame containing trajectory data, typically the output from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`. Must contain
#'   columns: `id`, `time`, `y` (state), and `tx` (treatment).
#' @param covariates Character vector of additional covariate names to include in
#'   the output (default: c("age", "sofa")). The first value of each covariate
#'   per patient will be retained.
#'
#' @return A data frame with columns:
#'   - id: patient identifier
#'   - tx: treatment assignment
#'   - start: start of the interval, first value per participant is zero
#'   - stop: end of the interval
#'   - y: state at the end of the interval
#'   - any additional covariates specified
#'
#' @details
#' Each event change within a participant concludes a time interval indicated by
#' start and stop and the state at the end of the interval. Columns given in
#' covariates will be kept in the output.
#'
#' @export
states_to_tte <- function(
  data,
  covariates = c("age", "sofa")
) {
  # Input validation
  required_cols <- c("id", "time", "y", "tx", "yprev")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate Andersen-Gill count data format ----------------------------------

  # Indicate rows with state changes
  data$state_change <- data[["y"]] - data[["yprev"]]
  available_covs <- if (is.null(covariates)) {
    character(0)
  } else {
    intersect(covariates, colnames(data))
  }

  # Keep state changes and ...
  keep <- is.na(data[["yprev"]]) |
    (!is.na(data$state_change) & data$state_change != 0) |
    # ... the last entry of each individual, unless they died
    (data[["time"]] == ave(data[["time"]], data[["id"]], FUN = max) &
      as.character(data[["y"]]) != "6")
  data <- data[
    keep,
    ,
    drop = FALSE
  ]

  data[["start"]] <- ave(data$time, data$id, FUN = function(t) {
    c(0, t[-length(t)])
  })
  data[["stop"]] <- data[["time"]]

  # Add covariates if specified ------------------------------------------------
  keep_cols <- c("id", "tx", "yprev", "start", "stop", "y", available_covs)
  data <- data[, unique(keep_cols), drop = FALSE]

  return(result = data)
}


#' Convert Markov Data to Time-to-Event (Start-Stop) Format (deprecated version)
#'
#' Converts a longitudinal dataset of daily states into a collapsed
#' start-stop format suitable for survival analysis
#'
#' @param data A data frame containing trajectory data, typically the output from
#'   `sim_trajectories_markov()` or `sim_trajectories_brownian()`. Must contain
#'   columns: `id`, `time`, `y` (state), and `tx` (treatment).
#' @param covariates Character vector of additional covariate names to include in
#'   the output (default: c("age", "sofa")). The first value of each covariate
#'   per patient will be retained.
#' @param absorbing_state Integer. The absorbing state code (e.g., 6 for death).
#'   Rows where \code{yprev} equals the absorbing state are filtered out, as
#'   transitions from the absorbing state are not meaningful. Default is 6.
#'
#' @return A tibble in start-stop format with columns:
#'   \code{id}, \code{tx}, \code{start}, \code{stop}, \code{y}, \code{yprev},
#'   and specified covariates.
#'
#' @details
#' Each event change within a participant concludes a time interval indicated by
#' start and stop and the state at the end of the interval. Columns given in
#' covariates will be kept in the output.
#'
#' @examples
#' \dontrun{
#' # After simulating trajectories
#' trajectories <- sim_trajectories_markov(baseline_data, follow_up_time = 60,
#'                                         lp_function = my_lp)
#' tte_data <- states_to_tte_v2(trajectories)
#'
#' # With different covariates
#' tte_data <- states_to_tte(trajectories, covariates = c("age", "sofa", "baseline_severity"))
#' }
#' @export
states_to_tte_v2 <- function(
  data,
  covariates = c("age", "sofa"),
  absorbing_state = 6
) {
  # 1. Input Validation
  required_cols <- c("id", "time", "y", "yprev", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  if (
    is.factor(data$id) ||
      !is.numeric(data$time) ||
      anyNA(data[, c("id", "time", "y")])
  ) {
    return(states_to_tte_v2_reference(data, covariates, absorbing_state))
  }

  keep <- is.na(data$yprev) |
    as.character(data$yprev) != as.character(absorbing_state)
  data <- data[keep, , drop = FALSE]
  data <- data[order(data$id, data$time, method = "radix"), , drop = FALSE]
  available_covs <- intersect(covariates, names(data))

  if (nrow(data) == 0L) {
    return(states_to_tte_v2_reference(data, covariates, absorbing_state))
  }

  n <- nrow(data)
  new_id <- c(TRUE, data$id[-1L] != data$id[-n])
  state_change <- c(TRUE, as.character(data$y[-1L]) != as.character(data$y[-n]))
  run_start <- new_id | state_change
  start_rows <- which(run_start)
  end_rows <- c(start_rows[-1L] - 1L, n)
  raw_start <- numeric(length(start_rows))
  prior_rows <- start_rows - 1L
  has_prior <- !new_id[start_rows]
  raw_start[has_prior] <- data$time[prior_rows[has_prior]]

  result <- data.frame(
    id = data$id[start_rows],
    tx = data$tx[start_rows],
    start = raw_start,
    stop = data$time[end_rows],
    y = data$y[start_rows],
    yprev = data$yprev[start_rows],
    check.names = FALSE
  )
  for (cov in available_covs) {
    result[[cov]] <- data[[cov]][start_rows]
  }
  rownames(result) <- NULL
  result
}

states_to_tte_v2_reference <- function(data, covariates, absorbing_state) {
  keep <- is.na(data$yprev) |
    as.character(data$yprev) != as.character(absorbing_state)
  data <- data[keep, , drop = FALSE]
  data <- data[order(data$id, data$time), , drop = FALSE]
  available_covs <- intersect(covariates, names(data))

  bind_rows_fill(lapply(unique(data$id), function(id) {
    group_data <- data[data$id == id, , drop = FALSE]
    raw_start <- c(0, utils::head(group_data$time, -1))
    raw_stop <- group_data$time
    state_change <- c(TRUE, group_data$y[-1] != utils::head(group_data$y, -1))
    run_id <- cumsum(state_change)

    bind_rows_fill(lapply(unique(run_id), function(run) {
      rows <- run_id == run
      out <- data.frame(
        id = id,
        tx = group_data$tx[rows][1],
        start = min(raw_start[rows]),
        stop = max(raw_stop[rows]),
        y = group_data$y[rows][1],
        yprev = group_data$yprev[rows][1],
        check.names = FALSE
      )
      for (cov in available_covs) {
        out[[cov]] <- group_data[[cov]][rows][1]
      }
      out
    }))
  }))
}


#' Calculate Difference in Time Alive and Out Of Hospital
#'
#' Computes the true treatment effect as the difference in time spent in a
#' target state (typically "home" or "alive and out of hospital") between
#' treatment groups in simulated Markov trajectory data.
#'
#' @param data A data frame containing trajectory data with columns: `id`,
#'   `y` (state), a treatment column, and a time column.
#' @param target_state Optional vector of state values to summarize. If `NULL`
#'   (default), all observed states in `y` are included.
#' @param time_var Character string giving the name of the time column.
#'   Default is `"time"`.
#' @param txvarname Character string giving the name of the treatment
#'   column. Default is `"tx"`.
#' @param reference_level Optional reference treatment level. Differences are
#'   calculated for every other treatment level relative to this value. If
#'   `NULL` (default), the first factor level is used for factor treatments;
#'   otherwise the smallest observed treatment value is used.
#'
#' @return A data frame with the following columns:
#'   - `state`: State being summarized
#'   - `reference_level`: Reference treatment level
#'   - `treatment_level`: Non-reference treatment level being compared to the
#'     reference
#'   - `true_effect`: Difference in state occupation probability summed across
#'     all time points
#'   - `reference_mean_time`: Mean time spent in target state for the reference
#'     treatment
#'   - `reference_sd_time`: Standard deviation of time spent in target state
#'     for the reference treatment
#'   - `treatment_mean_time`: Mean time spent in target state for the
#'     comparison treatment
#'   - `treatment_sd_time`: Standard deviation of time spent in target state
#'     for the comparison treatment
#'
#' @details
#' This function calculates the true treatment effect in two ways:
#'
#' 1. **State Occupation Probability Method**: Calculates the proportion of
#'    patients in each selected state at each time point for each treatment
#'    group, then sums the difference across all time points relative to the
#'    reference treatment level.
#'
#' 2. **Individual Time in State**: For each patient, counts the total number
#'    of time periods spent in each selected state, then computes group-level
#'    means and standard deviations for each treatment level.
#'
#' The `true_effect` represents the area between the state occupation curves
#' and can be interpreted as the expected difference in total days spent in
#' each state over the follow-up period.
#'
#' @examples
#' \dontrun{
#' # Calculate effects for all observed states
#' trajectories <- sim_trajectories_markov(baseline_data, lp_function = my_lp)
#' effect_summary <- calc_time_in_state_diff(trajectories)
#'
#' # Use custom time and treatment columns
#' effect_summary <- calc_time_in_state_diff(
#'   trajectories,
#'   target_state = 1,
#'   time_var = "visit_day",
#'   txvarname = "arm",
#'   reference_level = "control"
#' )
#'
#' # View results
#' print(effect_summary)
#' }
#'
#' @importFrom stats aggregate sd
#'
#' @export
calc_time_in_state_diff <- function(
  data,
  target_state = NULL,
  time_var = "time",
  txvarname = "tx",
  reference_level = NULL
) {
  # Input validation
  required_cols <- c("id", "y", txvarname, time_var)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  data_standardized <- data.frame(
    id = data[["id"]],
    y = if (is.factor(data[["y"]])) {
      as.character(data[["y"]])
    } else {
      data[["y"]]
    },
    .time = data[[time_var]],
    .treatment = if (is.factor(data[[txvarname]])) {
      as.character(data[[txvarname]])
    } else {
      data[[txvarname]]
    },
    stringsAsFactors = FALSE
  )

  observed_states <- sort(unique(data_standardized[["y"]]))

  if (is.null(target_state)) {
    target_state <- observed_states
  }

  target_state <- if (is.factor(target_state)) {
    as.character(target_state)
  } else {
    target_state
  }

  missing_states <- setdiff(target_state, observed_states)
  if (length(missing_states) > 0) {
    stop(
      "target_state must be drawn from observed y values: ",
      paste(observed_states, collapse = ", ")
    )
  }

  treatment_levels <- sort(unique(data_standardized[[".treatment"]]))

  if (length(treatment_levels) < 2) {
    stop("txvarname must contain at least 2 observed treatment levels")
  }

  if (is.null(reference_level)) {
    reference_level <- treatment_levels[[1]]
  }

  if (is.factor(reference_level)) {
    reference_level <- as.character(reference_level)
  }

  if (!reference_level %in% treatment_levels) {
    stop(
      "reference_level must be one of: ",
      paste(treatment_levels, collapse = ", ")
    )
  }

  state_results <- lapply(target_state, function(state_value) {
    state_indicator <- data_standardized[["y"]] %in% state_value
    unique_times <- sort(unique(data_standardized[[".time"]]))

    sops <- stats::aggregate(
      state_indicator,
      by = list(
        .treatment = data_standardized[[".treatment"]],
        .time = data_standardized[[".time"]]
      ),
      FUN = mean
    )
    names(sops)[names(sops) == "x"] <- "sop"

    sops_complete <- merge(
      expand.grid(
        .treatment = treatment_levels,
        .time = unique_times,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      ),
      sops,
      by = c(".treatment", ".time"),
      all.x = TRUE,
      sort = FALSE
    )
    sops_complete[["sop"]][is.na(sops_complete[["sop"]])] <- 0

    reference_sops <- sops_complete[
      sops_complete[[".treatment"]] == reference_level,
    ]
    reference_sops <- reference_sops[, c(".time", "sop"), drop = FALSE]
    names(reference_sops)[names(reference_sops) == "sop"] <- "reference_sop"

    comparison_sops <- sops_complete[
      sops_complete[[".treatment"]] != reference_level,
    ]
    comparison_sops <- comparison_sops[,
      c(".treatment", ".time", "sop"),
      drop = FALSE
    ]
    names(comparison_sops) <- c("treatment_level", ".time", "treatment_sop")

    estimand_data <- merge(
      comparison_sops,
      reference_sops,
      by = ".time",
      all.x = TRUE,
      sort = FALSE
    )
    estimand_data[["reference_sop"]][is.na(estimand_data[[
      "reference_sop"
    ]])] <- 0
    estimand_data[["effect"]] <- estimand_data[["treatment_sop"]] -
      estimand_data[["reference_sop"]]

    estimand <- stats::aggregate(
      effect ~ treatment_level,
      data = estimand_data,
      FUN = sum
    )
    names(estimand)[names(estimand) == "effect"] <- "true_effect"

    patient_time <- stats::aggregate(
      state_indicator,
      by = list(
        id = data_standardized[["id"]],
        .treatment = data_standardized[[".treatment"]]
      ),
      FUN = sum
    )
    names(patient_time)[names(patient_time) == "x"] <- "time_in_state"

    mean_summary <- stats::aggregate(
      time_in_state ~ .treatment,
      data = patient_time,
      FUN = mean
    )
    sd_summary <- stats::aggregate(
      time_in_state ~ .treatment,
      data = patient_time,
      FUN = stats::sd
    )
    time_summary <- merge(
      mean_summary,
      sd_summary,
      by = ".treatment",
      suffixes = c("_mean", "_sd"),
      sort = FALSE
    )
    names(time_summary) <- c(".treatment", "mean_time", "sd_time")

    reference_summary <- time_summary[
      time_summary[[".treatment"]] == reference_level,
    ]
    comparison_summary <- time_summary[
      time_summary[[".treatment"]] != reference_level,
    ]

    result <- merge(
      comparison_summary,
      estimand,
      by.x = ".treatment",
      by.y = "treatment_level",
      all.x = TRUE,
      sort = FALSE
    )

    data.frame(
      state = rep(state_value, nrow(result)),
      reference_level = rep(reference_level, nrow(result)),
      treatment_level = result[[".treatment"]],
      true_effect = result[["true_effect"]],
      reference_mean_time = rep(reference_summary[["mean_time"]], nrow(result)),
      reference_sd_time = rep(reference_summary[["sd_time"]], nrow(result)),
      treatment_mean_time = result[["mean_time"]],
      treatment_sd_time = result[["sd_time"]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, state_results)
}
