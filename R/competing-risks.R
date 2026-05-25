#' Format Time-to-Event Output for Competing Risks
#'
#' Prepares the output from `states_to_tte()` for competing-risks analysis
#' (e.g., Fine-Gray or Cause-Specific Cox) by identifying transitions from
#' a risk state (e.g., hospital) to an absorbing state (e.g., discharge or death).
#'
#' @description
#' This function transforms "state occupancy" data into "time-to-event" data.
#' It looks ahead to the subsequent interval to detect when a subject transitions
#' from a generic state to an event state. The interval ending at that transition
#' is assigned the corresponding event status. Rows representing the absorbing
#' states themselves are excluded, as the subject is no longer at risk.
#'
#' @param data A data frame produced by `states_to_tte()` containing at least
#'   columns `id`, `start`, `stop`, `y` (state), and `tx`.
#' @param event_status Integer code identifying the event of interest (e.g., discharge)
#'   in the `y` column (default: 1).
#' @param death_status Integer code identifying the competing event (e.g., death)
#'   in the `y` column (default: 6).
#' @param version_old Logical. If TRUE (default), uses the original implementation.
#'   If FALSE, uses an updated implementation with different logic for status assignment.
#'
#' @return A data frame where each row represents an interval at risk. The output
#'   includes a `status` column:
#'   \itemize{
#'     \item \code{0}: Censored (no event observed or lost to follow-up)
#'     \item \code{1}: Event of interest occurred at the end of this interval
#'     \item \code{2}: Competing event occurred at the end of this interval
#'   }
#'   Rows where `y` matches `event_status` or `death_status` are removed.
#'
#' @details
#' The function performs the following steps for each subject (`id`):
#' \enumerate{
#'   \item **Sorts** the data by time (`start`) to ensure chronological order.
#'   \item **Looks ahead** to the state (`y`) of the *next* interval.
#'   \item **Assigns status**:
#'     \itemize{
#'       \item If the *next* state is `event_status`, the current interval gets `status = 1`.
#'       \item If the *next* state is `death_status`, the current interval gets `status = 2`.
#'       \item Otherwise, `status = 0`.
#'     }
#'   \item **Filters**: Removes rows where the subject is already in an absorbing
#'     state (`event_status` or `death_status`), as they are no longer at risk.
#'   \item **Truncates**: Ensures only the trajectory up to the first observed event
#'     is retained per subject.
#' }
#'
#' @examples
#' \dontrun{
#' tte <- states_to_tte(trajectories)
#' # Prepare for Fine-Gray model (1=Discharge vs 2=Death)
#' fg_ready <- format_competing_risks(tte, event_status = 1, death_status = 6)
#' }
#'
#' @export
format_competing_risks <- function(
  data,
  event_status = 1,
  death_status = 6,
  version_old = TRUE
) {
  # Basic validation
  required_cols <- c("id", "start", "stop", "y", "tx")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!version_old) {
    data <- data[order(data$id, data$start), , drop = FALSE]
    data <- bind_rows_fill(lapply(unique(data$id), function(id) {
      group_data <- data[data$id == id, , drop = FALSE]
      next_y <- c(group_data$y[-1], NA)
      status <- ifelse(
        next_y == event_status,
        1L,
        ifelse(next_y == death_status, 2L, 0L)
      )
      status[is.na(status)] <- 0L
      cum_events <- cumsum(status > 0)
      prev_events <- c(0, utils::head(cum_events, -1))
      group_data$status <- status
      group_data <- group_data[prev_events == 0, , drop = FALSE]
      group_data
    }))
    data <- reorder_columns(data, c("id", "tx", "start", "stop", "status", "y"))
    return(data)
  }

  if (version_old) {
    # Split data into id groups
    data_split <- split(data, data[["id"]])

    # For each individual ...
    data_split <- lapply(data_split, FUN = function(d) {
      # ...if any event occurred in this group, select the according rows...
      event_or_death <- d[["y"]] == event_status | d[["y"]] == death_status
      if (any(event_or_death)) {
        d <- d[event_or_death, ]
        # ...else select the last row
      } else {
        d <- d[nrow(d), ]
      }
      # introduce status variable, which indicates event (1) or death (2)
      d[["status"]] <- factor(ifelse(
        d[["y"]] == event_status,
        1,
        ifelse(d[["y"]] == death_status, 2, 0)
      ))
      d <- d[1, ]
      d[["etime"]] <- as.numeric(d[["stop"]])
      return(d)
    })
    return(do.call(rbind, data_split))
  }
}
