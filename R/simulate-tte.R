#' Simulate Patient Trajectories Using Latent Time-to-Event Model
#'
#' Generates individual patient trajectories using a latent time-to-event mechanism
#' with recurrent events. Each patient transitions through ordered states based on
#' state-specific event rates, with support for treatment effects via hazard ratios.
#'
#' @param baseline_data A data frame containing baseline patient characteristics.
#'   Must include columns: `id`, `tx` (treatment indicator), and `state` (initial state).
#'   Default is `NULL` (user must provide). Each row represents one patient.
#' @param baseline_states Integer vector. States that patients can start in at time 0.
#'   Default is NULL, which allows starting in any state. If provided, only these
#'   states will be assigned at baseline.
#' @param prob Probabilities for each baseline_state.
#' @param n Integer. Number of patients to simulate when `baseline_data` is
#'   `NULL`. Ignored when `baseline_data` is supplied.
#' @param states Integer vector. Ordered states in the model (default: 1:6).
#'   States should be numbered consecutively.
#' @param absorbing_states Integer vector. States that are absorbing (once entered,
#'   patients cannot leave). Default is 6 (typically death).
#' @param follow_up_time Integer. Number of days to simulate per patient (default: 60).
#' @param param Numeric vector. Baseline event rates (hazards) for each state.
#'   Length must equal `length(states)`. These represent the rate parameter
#'   lambda for
#'   exponential waiting times in the control group.
#' @param hazard_ratios List of numeric vectors. Treatment effects as multiplicative
#'   factors on the baseline rates for each positive treatment arm. `tx = 0`
#'   uses the control rates in `param`, `tx = 1` uses `hazard_ratios[[1]]`,
#'   `tx = 2` uses `hazard_ratios[[2]]`, and so on. Length of each vector
#'   must equal `length(states)`.
#'   Default is `list(c(1, 1, 1, 1, 1, 1))` (no treatment effect).
#' @param b Numeric. Autoregressive coefficient for recurrent events within a state
#'   (default: 0). The rate for event j in a given state is:
#'   rate_j = param + b * (j - 1). Positive b means events become more likely over time
#'   (Poisson process acceleration); b = 0 means independent events. If a vector is given,
#'   it must be the same length as the states vector (event-type-specific acceleration).
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#'
#' @return A data frame (tibble) with columns:
#'   - id: patient identifier
#'   - tx: treatment assignment (0 = control, 1 = treatment)
#'   - time: time/day (1 to follow_up_time)
#'   - y: observed state at time
#'   - yprev: state at previous time point
#'
#' @details
#' This function implements a latent time-to-event data generating mechanism:
#' 1. For each patient, state, and treatment arm, waiting times to recurrent events
#'    are generated from exponential distributions with state-specific rates.
#' 2. Rates are adjusted by the corresponding hazard ratio (treatment effect).
#' 3. Event times are cumulated to determine when state transitions occur.
#' 4. Events occurring after follow-up are censored (administrative censoring).
#' 5. Each patient's trajectory is expanded to include all days 1 to follow_up_time,
#'    with states carried forward between event times (last observation carried forward, LOCF).
#' 6. Once an absorbing state is entered, no further transitions occur.
#'
#' Patients in absorbing states at baseline are handled by advancing their first
#' absorbing event to day 1.
#'
#' @examples
#' \dontrun{
#' # Basic example with control parameters
#' baseline <- data.frame(
#'   id = 1:100,
#'   tx = rbinom(100, 1, 0.5),
#'   state = sample(2:5, 100, replace = TRUE)
#' )
#'
#' trajectories <- sim_trajectories_tte(
#'   baseline_data = baseline,
#'   param = c(0.05, 0.0035, 0.0025, 0.002, 0.002, 0.005),
#'   hazard_ratios = list(c(0.8, 1, 1, 1, 1, 1)),
#'   follow_up_time = 30,
#'   seed = 12345
#' )
#'
#' # Generate baseline_data within the function
#' test_traj <- sim_trajectories_tte(
#'   baseline_data = NULL,
#' states = 1:6,
#'   baseline_states = c(2:5),
#'   prob = c(0.55, 0.2, 0.15, 0.1),
#'   n = 10000,
#'   absorbing_states = 6,
#'   follow_up_time = 60,
#'   param = c(0.05, 0.003, 0.001, 0.001, 0.001, 0.0015),
#'   hazard_ratios = list(c(1.145, 1, 1, 1, 1, 1)),
#'   b = 0,
#'   seed = NULL
#' )
#' plot_sops(test_traj)
#' calc_time_in_state_diff(test_traj)
#'
#' # With multiple treatment arms
#' trajectories_multi <- sim_trajectories_tte(
#'   baseline_data = baseline,
#'   param = c(0.05, 0.0035, 0.0025, 0.002, 0.002, 0.005),
#'   hazard_ratios = list(
#'     c(0.8, 1, 1, 1, 1, 1),  # Treatment arm 1
#'     c(0.9, 1, 1, 1, 1, 1)   # Treatment arm 2
#'   ),
#'   follow_up_time = 30,
#'   seed = 12345
#' )
#' }
#'
#' @seealso
#' \code{\link{recurr_event}} for the underlying recurrent event generation
#'
#' @importFrom stats rnorm
#'
#' @export

sim_trajectories_tte <- function(
  baseline_data,
  baseline_states = NULL,
  prob = c(0.55, 0.2, 0.15, 0.1),
  n = 1000,
  states = 1:6,
  absorbing_states = 6,
  follow_up_time = 60,
  param = c(0.05, 0.0035, 0.0025, 0.002, 0.002, 0.005),
  hazard_ratios = list(c(1, 1, 1, 1, 1, 1)),
  b = 0,
  seed = NULL
) {
  # Input validation
  if (!is.null(baseline_data) && !is.data.frame(baseline_data)) {
    stop("baseline_data must be a data frame or NULL")
  }

  if (is.null(baseline_data)) {
    baseline_data <- data.frame(
      id = 1:n,
      event_time = 0,
      tx = sample(c(0, 1), n, replace = TRUE),
      state = sample(
        if (!is.null(baseline_states)) baseline_states else states,
        n,
        replace = TRUE,
        prob = prob
      )
    )
  }

  required_cols <- c("id", "tx", "state", "event_time")
  missing_cols <- setdiff(required_cols, names(baseline_data))
  if (length(missing_cols) > 0) {
    stop(
      "baseline_data must contain columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!is.numeric(param) || length(param) != length(states)) {
    stop("param length must equal length(states)")
  }

  if (!is.list(hazard_ratios)) {
    stop("hazard_ratios must be a list of numeric vectors")
  }

  if (
    !all(vapply(
      hazard_ratios,
      function(x) {
        is.numeric(x) && length(x) == length(states)
      },
      logical(1)
    ))
  ) {
    stop(
      "Each element of hazard_ratios must have length equal to length(states)"
    )
  }

  if (any(baseline_data$state %in% absorbing_states)) {
    warning(
      "Some individuals are in absorbing state at baseline. ",
      "Consider changing the baseline distribution."
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  tx_levels <- sort(unique(baseline_data$tx))
  if (
    !is.numeric(tx_levels) ||
      anyNA(tx_levels) ||
      any(tx_levels < 0) ||
      any(tx_levels != floor(tx_levels))
  ) {
    stop("baseline_data$tx must contain non-negative integer treatment arm codes")
  }

  missing_hr <- tx_levels[tx_levels > length(hazard_ratios)]
  if (length(missing_hr) > 0) {
    stop(
      "hazard_ratios must provide one vector for each positive treatment arm. ",
      "Missing arm(s): ",
      paste(missing_hr, collapse = ", ")
    )
  }

  params_for_tx <- function(tx) {
    if (tx == 0) {
      return(param)
    }
    param * hazard_ratios[[tx]]
  }

  # Expand b vector for each event type
  if (length(b) == 1) {
    # Fixed acceleration
    b <- rep(b, times = length(states))
  } else {
    # Event-type-specific acceleration: check input
    if (length(b) != length(states)) {
      stop(
        "Autoregressive parameter b must have length one or equal to length(states)"
      )
    }
  }

  # Generate recurrent events for each treatment arm by state combination
  state_changes <- list()
  m <- 1

  for (tx_val in tx_levels) {
    ids_tx <- baseline_data$id[baseline_data$tx == tx_val]
    tx_param <- params_for_tx(tx_val)

    for (state_pos in seq_along(states)) {
      i <- states[state_pos]
      state_change <- recurr_event(
        id = ids_tx,
        param = tx_param[state_pos],
        b = b[state_pos],
        follow_up = follow_up_time,
        max_events = NULL
      )
      if (nrow(state_change) == 0L) {
        state_changes[[m]] <- data.frame(
          id = numeric(0),
          event_time = numeric(0),
          state = numeric(0),
          tx = numeric(0)
        )
      } else {
        state_changes[[m]] <- cbind(state_change, state = i, tx = tx_val)
      }
      m <- m + 1
    }
  }

  # Combine event times, add baseline, and order
  state_changes_long <- rbind(
    do.call(rbind, state_changes),
    baseline_data
  )

  state_changes_long <- state_changes_long[
    order(state_changes_long[, "id"], state_changes_long[, "event_time"]),
  ]
  state_changes_long[, "time"] <- floor(state_changes_long[, "event_time"])

  # Split by patient ID
  split_state_changes_long <- split(
    as.data.frame(state_changes_long),
    state_changes_long[, "id"]
  )

  # LOCF helper function
  locf <- function(v) {
    ind <- which(!is.na(v))
    if (length(ind) == 0L) {
      return(v)
    }
    fi <- findInterval(seq_along(v), ind)
    res <- v[ind][pmax(fi, 1L)]
    res[fi == 0L] <- NA
    res
  }

  # Process each patient: remove post-absorbing events, handle duplicates, expand days, and LOCF
  state_changes_list <- lapply(split_state_changes_long, function(x) {
    this_id <- unique(x$id[!is.na(x$id)])
    if (length(this_id) != 1L) {
      this_id <- this_id[1L]
    }

    this_tx <- unique(x$tx[!is.na(x$tx)])
    if (length(this_tx) != 1L) {
      this_tx <- this_tx[1L]
    }

    # Absorbing state at baseline: set to day 1 (observed next day)
    if (any(x$time == 0 & x$state %in% absorbing_states)) {
      x$time[x$time == 0 & x$state %in% absorbing_states] <- 1
    }

    # Remove events after absorbing state
    absorbing_selector <- ave(
      x[, "state"] %in% absorbing_states,
      x[, "id"],
      FUN = cumsum
    )
    absorbing_selector <- c(0, absorbing_selector[-length(absorbing_selector)])
    x <- x[absorbing_selector == 0, , drop = FALSE]

    # Keep last event of each day
    x <- x[!duplicated(x$time, fromLast = TRUE), , drop = FALSE]

    # Expand to all days (0:follow_up_time)
    complete_times <- 0:follow_up_time
    x <- merge(
      data.frame(time = complete_times),
      x,
      by = "time",
      all.x = TRUE
    )

    # Restore id and tx
    x$id <- this_id
    x$tx <- this_tx

    # Convert state to numeric and apply LOCF
    if (is.factor(x$state)) {
      x$state <- as.integer(as.character(x$state))
    }
    x$y <- locf(x$state)

    # Generate yprev
    x$yprev <- c(NA, x$y[-nrow(x)])

    # Drop day 0 (keep only days 1 to follow_up_time)
    x <- x[x$time > 0, , drop = FALSE]

    # Select relevant columns
    x[, c("id", "tx", "time", "y", "yprev"), drop = FALSE]
  })

  state_changes_long <- do.call(rbind, state_changes_list)
  rownames(state_changes_long) <- NULL

  state_changes_long
}
