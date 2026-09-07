#' Generate recurrent event times for individuals (helper for latent TTE DGM)
#'
#' Simulate recurrent event times in long format for a set of individuals.
#' For each individual the waiting times between events are drawn from
#' exponential distributions with state-specific rates and then cumulated to
#' event times. Event times that exceed the administrative follow-up are
#' censored (removed).
#'
#' @param id Optional numeric vector of subject identifiers. If supplied, these
#'   identifiers are repeated `max_events` times internally (one block of
#'   rates per subject). If missing, identifiers 1:n are generated where
#'   `n` is provided or inferred.
#' @param dist Character. Name of the waiting-time distribution. Currently only
#'   "Exponential" is supported (default).
#' @param param Numeric. Baseline rate parameter for the waiting-time distribution.
#'   Must have length 1, or length equal to the number of individuals in `id`.
#'   When a vector is supplied, each individual receives their own baseline rate.
#' @param b Numeric scalar. Relative factor applied to the rate for every subsequent
#'   event (autoregressive Poisson / accelerating rates). Rate for event j is
#'   lambda_i * b^(j - 1). Default `1` (constant rates).
#' @param follow_up Numeric scalar. Administrative follow-up time. Event times
#'   \code{>= follow_up} are censored and not returned. Default `60`.
#' @param max_events Optional integer. Maximum number of events to generate per
#'   individual. If NULL (default) the function chooses a suitable `max_events`
#'   between 3 and 20 so that the probability of observing all events within
#'   `follow_up` is very small (threshold 1e-4). When provided, the function
#'   uses the given value.
#'
#' @return A two-column numeric matrix / data.frame with columns
#'   - \code{id}: subject identifier (numeric)
#'   - \code{event_time}: event time (numeric) for events strictly less than \code{follow_up}
#'
#' @details
#' Algorithm summary:
#' 1. Construct a vector of per-event rates of length \code{max_events} for each
#'    subject: \eqn{rate_{ij} = lambda_i * b^(j - 1)}.
#' 2. If \code{max_events} is not supplied, choose the smallest value in
#'    3:20 for which the probability of experiencing all \code{max_events}
#'    within \code{follow_up} is < 1e-4. When rates are constant the gamma
#'    distribution (pgamma) is used; otherwise \code{sdprisk::phypoexp} is used.
#' 3. Simulate waiting times with \code{rexp(n * max_events, rate = rates_rep)}.
#'    Rates are arranged in subject-major order so that each subject receives
#'    their own block of per-event rates.
#' 4. For each subject compute the cumulative sums of waiting times to obtain
#'    event times, and retain only those event times < \code{follow_up}.
#'
#' Notes and caveats:
#' - Event times exactly equal to \code{follow_up} are excluded (strict <).
#'
#' @examples
#' # small example
#' recurr_event(id = 1:3, param = 0.1, b = 1, follow_up = 10, max_events = 5)
#'
#' # let function choose max_events automatically
#' recurr_event(id = 1:10, param = 0.05, b = 0.99, follow_up = 30)
#'
#' @seealso \code{\link[sdprisk]{phypoexp}} for hypoexponential CDF used internally
#' @keywords recurr_event
#' @importFrom sdprisk phypoexp
#' @export

recurr_event <- function(
  id,
  dist = "Exponential",
  param,
  b = 1,
  follow_up = 60,
  max_events = NULL
) {
  if (missing(id)) {
    id <- 1L
  }
  n <- length(id)

  if (!identical(dist, "Exponential")) {
    stop("Only dist = 'Exponential' is supported.")
  }

  lambda_i <- param
  if (length(lambda_i) == 1) {
    lambda_i <- rep(lambda_i, n)
  } else if (length(lambda_i) != n) {
    stop("Length of param must be one or equal to the number of participants.")
  }

  #______________________________________________________________________________#
  #____Find max_events and corresponding rates___________________________________#
  #______________________________________________________________________________#
  # Find max_events so that the probability of experiencing all events within
  # the study interval is very small (here set to 0.0001, can be changed below)
  if (is.null(max_events)) {
    lambda_i_min <- min(lambda_i)
    for (candidate_max_events in 3:20) {
      candidate_rates <- lambda_i_min * b^(0:(candidate_max_events - 1))
      if (any(candidate_rates <= 0)) {
        stop("Rates must be strictly positive. Check param and b.")
      }

      if (length(unique(candidate_rates)) == 1) {
        p <- stats::pgamma(
          follow_up,
          shape = candidate_max_events,
          rate = unique(candidate_rates)
        )
      } else {
        p <- phypoexp(follow_up, rate = candidate_rates)
      }

      if (p < 0.0001) {
        max_events <- candidate_max_events
        break()
      }

      if (candidate_max_events == 20) {
        warning(paste0(
          "max_events = 20, but probability of experiencing all events within time ",
          follow_up,
          " is still ",
          p,
          ". To lower this probability, try different follow-up-time or lambda_i or override max_events."
        ))
        max_events <- candidate_max_events
      }
    }
  }

  rate_mat <- outer(lambda_i, 0:(max_events - 1), function(lam, j) lam * b^j)
  if (any(rate_mat <= 0)) {
    stop("Rates must be strictly positive. Check param and b.")
  }
  rates_rep <- as.vector(t(rate_mat))

  #______________________________________________________________________________#
  #____ Model waiting times between events and store as vector   ________________#
  #______________________________________________________________________________#
  id_rep <- rep(id, each = max_events)

  waiting_times_orig <- rexp(n * max_events, rate = rates_rep)

  event_times <- ave(waiting_times_orig, id_rep, FUN = cumsum)
  names(event_times) <- id_rep

  censored_event_times <- event_times[event_times < follow_up]

  return(cbind(
    id = as.numeric(names(censored_event_times)),
    event_time = censored_event_times
  ))
}
