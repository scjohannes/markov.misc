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
#' @param param Numeric. Parameter vector for the waiting-time distribution.
#'   The first element is taken as the baseline rate lambda used to
#'   construct the per-event rates. Only `param[1]` is used by the current
#'   implementation.
#' @param b Numeric scalar. Increment added to the rate for every subsequent
#'   event (autoregressive Poisson / accelerating rates). Rate for event j is
#'   lambda + b * (j - 1). Default `0` (constant rates).
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
#' 1. Construct a vector of per-event rates of length \code{max_events}:
#'    \eqn{rate_j = lambda + b * (j - 1)}.
#' 2. If \code{max_events} is not supplied, choose the smallest value in
#'    3:20 for which the probability of experiencing all \code{max_events}
#'    within \code{follow_up} is < 1e-4. When rates are constant the gamma
#'    distribution (pgamma) is used; otherwise \code{sdprisk::phypoexp} is used.
#' 3. Simulate waiting times with \code{rexp(n * max_events, rate = rates)}.
#'    Rates are recycled to length \code{n * max_events} so that each subject
#'    receives the block of per-event rates in order.
#' 4. For each subject compute the cumulative sums of waiting times to obtain
#'    event times, and retain only those event times < \code{follow_up}.
#'
#' Notes and caveats:
#' - Event times exactly equal to \code{follow_up} are excluded (strict <).
#'
#' @examples
#' # small example
#' recurr_event(id = 1:3, param = 0.1, b = 0, follow_up = 10, max_events = 5)
#'
#' # let function choose max_events automatically
#' recurr_event(id = 1:10, param = 0.05, b = 0.01, follow_up = 30)
#'
#' @seealso \code{\link[sdprisk]{phypoexp}} for hypoexponential CDF used internally
#' @keywords recurr_event
#' @importFrom sdprisk phypoexp
#' @export

recurr_event <- function(
  id,
  dist = "Exponential",
  param,
  b = 0,
  follow_up = 60,
  max_events = NULL
) {
  # prepare variables
  if (missing(id)) {
    id <- 1L
  }
  n <- length(id)

  lambda_i <- param

  #______________________________________________________________________________#
  #____Find max_events and corresponding rates___________________________________#
  #______________________________________________________________________________#
  # Find max_events so that the probability of experiencing all events within
  # the study interval is very small (here set to 0.0001, can be changed below)
  if (is.null(max_events)) {
    for (max_events in 3:20) {
      rates <- numeric(max_events)

      if (length(param) == 1) {
        for (j in 1:max_events) {
          rates[j] <- lambda_i + b * (j - 1)
        }
      } else {
        if (length(param) == max_events) {
          rates <- param
        } else {
          stop("Length of param must be one or equal to max_events.")
        }
      }

      if (length(unique(rates)) == 1) {
        p <- stats::pgamma(
          follow_up,
          shape = length(rates),
          rate = unique(rates)
        )
      } else {
        p <- phypoexp(follow_up, rate = rates)
      }

      if (p < 0.0001) {
        break()
      }

      if (max_events == 20) {
        warning(paste0(
          "max_events = 20, but propability of experiencing all events within time ",
          follow_up,
          " is still ",
          p,
          ". To lower this probability, try different follow-up-time or lambda_i or override max_events."
        ))
      }
    }
  } else {
    # define rates for waiting time prior to event i
    rates <- numeric(max_events)
    for (i in 1:max_events) {
      rates[i] <- lambda_i + b * (i - 1)
    }
  }

  #______________________________________________________________________________#
  #____ Model waiting times between events and store as vector   ________________#
  #______________________________________________________________________________#
  # Assign patient ID for later dataframe
  id <- rep(id, each = max_events)

  # generate the waiting times between events for each individual
  rates_rep <- rep(rates, times = n)
  waiting_times_orig <- rexp(n * max_events, rate = rates_rep)

  # vector of cumulative sums of individuals, i.e., time points of state change
  event_times <- ave(waiting_times_orig, id, FUN = cumsum)
  names(event_times) <- id

  #   - censor intervals which exceed the follow-up (administrative censoring)
  censored_event_times <- event_times[event_times < follow_up]

  return(cbind(
    id = as.numeric(names(censored_event_times)),
    event_time = censored_event_times
  ))
}
