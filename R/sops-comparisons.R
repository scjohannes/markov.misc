# Average SOP comparisons.

#' Average Comparisons From Markov State Occupancy Probabilities
#'
#' Computes average contrasts between counterfactual state occupancy summaries.
#' The first value in `variables[[1]]` is the reference level, and every later
#' value is compared to it. Reverse the order of values to reverse the contrast.
#'
#' @param model A fitted Markov transition model supported by [avg_sops()].
#' @param newdata Optional data frame of standardization profiles. See
#'   [avg_sops()].
#' @param variables A named list with one counterfactual variable. The variable
#'   must contain at least two values, for example `list(tx = c(0, 1))`.
#' @param metric Character scalar. One of `"sop"`, `"time_in_state"`, or
#'   `"time_benefit"`.
#' @param states State selection for `"sop"` and `"time_in_state"`. `NULL`
#'   returns one result per state. An atomic vector is treated as one lumped
#'   state set. A named list returns one result per named state set.
#' @param comparison Character scalar. `"difference"` computes comparison level
#'   minus reference level. `"ratio"` computes comparison level divided by
#'   reference level. `"time_benefit"` only supports `"difference"`.
#' @param by Optional character vector of variables to stratify by after
#'   standardization.
#' @param times Required visit-scale time points. See [avg_sops()].
#' @param ylevels A vector of state levels. See [avg_sops()].
#' @param absorb The absorbing state. See [avg_sops()].
#' @param time_map Optional named numeric vector or data frame mapping visit
#'   labels to real elapsed times. Used by `"time_in_state"` and
#'   `"time_benefit"`.
#' @param origin_time Optional real time for an empirical baseline anchor. See
#'   [interpolate_sops()].
#' @param xout Optional numeric real-time grid when `time_map` is supplied.
#' @param origin Origin handling when `origin_time` is supplied. See
#'   [interpolate_sops()].
#' @param time_unit Optional label stored in output.
#' @param refit_data Optional full longitudinal data used only by refit-bootstrap
#'   inference. It is not used for point estimates. See [avg_sops()].
#' @param id_var Name of the patient ID variable. See [avg_sops()].
#' @param tvarname Name of the time variable in the model. See [avg_sops()].
#' @param pvarname Name of the previous state variable in the model. See
#'   [avg_sops()].
#' @param p2varname Optional second previous-state variable. See [avg_sops()].
#' @param gap Name of the time gap variable, if used. See [avg_sops()].
#' @param t_covs Optional time-varying covariate lookup table for explicit
#'   precomputed time-basis columns. See [avg_sops()].
#' @param include_re Logical. For `rmsb::blrm()` fits with `cluster()`, include
#'   fitted random-effect draws in posterior predictions for known IDs.
#' @param n_draws Integer number of posterior draws to sample for `blrm`, or
#'   `NULL` to use all stored draws.
#' @param seed Optional random seed for reproducible posterior draw sampling.
#' @param posterior_summary Summary used for `blrm` posterior comparison draws:
#'   `"mean"` or `"median"`.
#' @param conf_level Confidence level for `blrm` posterior intervals.
#' @param return_draws Logical. For `blrm`, store posterior comparison draws.
#'   Frequentist uncertainty draws are stored by [inferences()].
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A data frame of class `markov_avg_comparisons`.
#'
#' @details
#' `metric = "sop"` and `metric = "time_in_state"` are computed from marginal
#' SOPs. `metric = "time_benefit"` is computed before marginalization from
#' paired patient/profile-level counterfactual state distributions, because it
#' is nonlinear in the two state distributions.
#'
#' Frequentist uncertainty follows the existing package workflow:
#' `avg_comparisons(...) |> inferences(...)`. Bayesian `blrm` models use their
#' posterior SOP draws directly and therefore return uncertainty from
#' `avg_comparisons()`.
#'
#' @seealso [avg_sops()], [sops()], [inferences()]
#'
#' @examples
#' \dontrun{
#' avg_comparisons(
#'   fit,
#'   variables = list(tx = c(0, 1)),
#'   metric = "time_in_state",
#'   states = "1",
#'   times = 1:30,
#'   ylevels = 1:6,
#'   absorb = 6
#' ) |>
#'   inferences(method = "simulation", n_sim = 500)
#' }
#'
#' @export
avg_comparisons <- function(
  model,
  newdata = NULL,
  variables,
  metric = c("sop", "time_in_state", "time_benefit"),
  states = NULL,
  comparison = c("difference", "ratio"),
  by = NULL,
  times,
  ylevels = NULL,
  absorb = NULL,
  time_map = NULL,
  origin_time = NULL,
  xout = NULL,
  origin = c("empirical_baseline", "none"),
  time_unit = NULL,
  refit_data = NULL,
  id_var = NULL,
  tvarname = "time",
  pvarname = "yprev",
  p2varname = NULL,
  gap = NULL,
  t_covs = NULL,
  include_re = FALSE,
  n_draws = 100L,
  seed = NULL,
  posterior_summary = c("mean", "median"),
  conf_level = 0.95,
  return_draws = FALSE,
  ...
) {
  metric <- match.arg(metric)
  comparison <- match.arg(comparison)
  origin <- match.arg(origin)
  posterior_summary <- match.arg(posterior_summary)

  if (missing(times) || is.null(times)) {
    stop("`times` must be supplied to `avg_comparisons()`.")
  }
  conf_level <- validate_conf_level(conf_level)
  variables <- validate_avg_comparison_variables(variables)
  validate_avg_comparison_metric(metric, states, comparison)

  if (metric == "time_benefit") {
    setup <- avg_comparison_setup(
      model = model,
      newdata = newdata,
      refit_data = refit_data,
      variables = variables,
      by = by,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      id_var = id_var,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      ...
    )
    tb <- avg_comparison_time_benefit_point(
      model = model,
      setup = setup,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = return_draws,
      comparison = comparison,
      time_map = time_map,
      origin_time = origin_time,
      xout = xout,
      origin = origin,
      time_unit = time_unit,
      ...
    )
    result <- tb$result
    setup <- tb$setup
  } else {
    avg <- avg_comparison_replay_avg_sops(
      model = model,
      newdata = newdata,
      refit_data = refit_data,
      variables = variables,
      by = by,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = tvarname,
      pvarname = pvarname,
      id_var = id_var,
      p2varname = p2varname,
      gap = gap,
      t_covs = t_covs,
      include_re = include_re,
      n_draws = n_draws,
      seed = seed,
      posterior_summary = posterior_summary,
      conf_level = conf_level,
      return_draws = inherits(model, "blrm") || isTRUE(return_draws),
      extra_args = list(...)
    )
    state_sets <- normalize_comparison_state_sets(states, attr(avg, "ylevels"))
    result <- avg_comparison_from_avg_sops(
      avg,
      metric = metric,
      state_sets = state_sets,
      comparison = comparison,
      time_map = time_map,
      origin_time = origin_time,
      xout = xout,
      origin = origin,
      time_unit = time_unit,
      return_draws = return_draws
    )
    setup <- avg_comparison_setup_from_sops(avg)
  }

  result <- set_avg_comparison_attrs(
    result = result,
    model = model,
    setup = setup,
    metric = metric,
    states = states,
    comparison = comparison,
    include_re = include_re,
    n_draws = n_draws,
    seed = seed,
    posterior_summary = posterior_summary,
    conf_level = conf_level,
    time_map = time_map,
    origin_time = origin_time,
    xout = xout,
    origin = origin,
    time_unit = time_unit,
    extra_args = list(...)
  )

  result
}
