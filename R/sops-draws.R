# SOP draw summaries and extraction.

#' Compute Confidence Intervals from Draws
#'
#' Computes confidence intervals from simulation or bootstrap draws.
#'
#' @param draws_df Data frame with draws (must have `estimate` and `draw_id` columns).
#' @param group_cols Character vector of grouping columns.
#' @param conf_level Confidence level (default 0.95).
#' @param conf_type Type of CI: "perc" (percentile) or "wald".
#'
#' @return Data frame with summary statistics.
#'
#' @keywords internal
compute_ci_from_draws <- function(
  draws_df,
  group_cols,
  conf_level = 0.95,
  conf_type = "perc"
) {
  alpha <- 1 - conf_level

  # Aggregate to get summary statistics
  agg_formula <- stats::as.formula(
    paste("estimate ~", paste(group_cols, collapse = " + "))
  )

  if (conf_type == "perc") {
    # Percentile confidence intervals
    summary_stats <- stats::aggregate(
      agg_formula,
      data = draws_df,
      FUN = function(x) {
        c(
          conf.low = stats::quantile(x, alpha / 2, na.rm = TRUE),
          conf.high = stats::quantile(x, 1 - alpha / 2, na.rm = TRUE),
          std.error = stats::sd(x, na.rm = TRUE)
        )
      }
    )
  } else if (conf_type == "wald") {
    # Wald confidence intervals
    summary_stats <- stats::aggregate(
      agg_formula,
      data = draws_df,
      FUN = function(x) {
        se <- stats::sd(x, na.rm = TRUE)
        critical <- abs(stats::qnorm(alpha / 2))
        mean_est <- mean(x, na.rm = TRUE)
        c(
          conf.low = mean_est - critical * se,
          conf.high = mean_est + critical * se,
          std.error = se
        )
      }
    )
  } else {
    stop("conf_type must be 'perc' or 'wald'")
  }

  # Fix aggregate's matrix column output
  mat <- summary_stats$estimate
  summary_stats$estimate <- NULL
  summary_stats$conf.low <- mat[, 1]
  summary_stats$conf.high <- mat[, 2]
  summary_stats$std.error <- mat[, 3]

  summary_stats
}


#' Extract Individual Draws from Inference Objects
#'
#' Extracts the individual draws (bootstrap samples, simulated values from MVN, etc.)
#' from an object returned by `inferences()` with `return_draws = TRUE`.
#' This function joins the draws back to the original point estimate object,
#' preserving all covariates, grouping variables, and summary statistics.
#'
#' @param object A SOP or average-comparison object returned with stored draws.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item draw_id: Simulation or bootstrap iteration number
#'     \item time: Time point
#'     \item state: State number
#'     \item draw: Draw-specific estimate of state occupation probability
#'     \item estimate: The original point estimate from the model
#'     \item conf.low, conf.high, std.error: Summary statistics from the point estimate object
#'     \item fwb_weight or score_weight: Optional draw-specific patient weight
#'       for eligible ungrouped stored-data `sops()` draws
#'     \item Additional columns from the original object (covariates, etc.)
#'   }
#'   Each row represents one draw for a specific time-state combination.
#'
#' @details
#' This function retrieves the draws from objects created by `inferences()`.
#' It performs a join between the draws and the point estimate object to ensure
#' that all metadata is preserved. To avoid conflict, the `estimate` column from
#' the draws is renamed to `draw`.
#'
#' Draw-specific patient weights are included only when they are meaningful for
#' manual downstream averaging: ungrouped `sops()` output from FWB or score
#' bootstrap on the stored empirical prediction cohort. They are not included for
#' `avg_sops()` or grouped `sops()` because those draws are already averaged, and
#' they are not included when the original SOP object used user-supplied
#' prediction profiles. The column names `fwb_weight` and `score_weight` are
#' reserved for these draw weights when they are attached.
#'
#' @examples
#' \dontrun{
#' # Create object with bootstrap draws
#' result <- avg_sops(
#'   model = fit,
#'   newdata = data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:30,
#'   ylevels = 1:6,
#'   absorb = 6,
#'   id_var = "id"
#' ) |>
#'   inferences(n_boot = 1000, return_draws = TRUE)
#'
#' # Extract draws
#' draws <- get_draws(result)
#'
#' # Plot distribution for state 1 at time 10 under treatment
#' library(ggplot2)
#' draws |>
#'   filter(time == 10, state == 1, tx == 1) |>
#'   ggplot(aes(x = draw)) +
#'   geom_histogram(bins = 30) +
#'   labs(title = "Distribution: P(State 1 | Time 10, Treatment)")
#'
#' # Compute mean time in state 1 with bootstrap CI
#' library(dplyr)
#' time_in_state_boot <- draws |>
#'   filter(state == 1, tx == 1) |>
#'   group_by(draw_id) |>
#'   summarise(total_time = sum(draw))
#'
#' quantile(time_in_state_boot$total_time, c(0.025, 0.5, 0.975))
#'
#' # Compare treatment effect on time in state
#' treatment_effect <- draws |>
#'   filter(state == 1) |>
#'   group_by(draw_id, tx) |>
#'   summarise(total_time = sum(draw), .groups = "drop") |>
#'   pivot_wider(names_from = tx, values_from = total_time,
#'               names_prefix = "tx") |>
#'   mutate(effect = tx1 - tx0)
#'
#' quantile(treatment_effect$effect, c(0.025, 0.5, 0.975))
#' }
#'
#' @seealso [inferences()], [avg_sops()]
#'
#' @export
get_draws <- function(object) {
  if (
    !inherits(
      object,
      c("markov_avg_sops", "markov_sops", "markov_avg_comparisons")
    )
  ) {
    stop(
      "get_draws() requires an object from inferences(). ",
      "Got: ",
      paste(class(object), collapse = ", ")
    )
  }

  # 1. Extract draws attribute
  draws <- attr(object, "bootstrap_draws")
  if (is.null(draws)) {
    draws <- attr(object, "simulation_draws")
  }
  if (is.null(draws)) {
    draws <- attr(object, "draws")
  }

  if (is.null(draws)) {
    method <- attr(object, "method")
    msg <- "No draws found. Run inferences() with return_draws = TRUE."
    if (!is.null(method)) {
      msg <- paste0(msg, " (Method used: '", method, "')")
    }
    stop(msg)
  }

  # 2. Prepare metadata from the original object
  meta <- as.data.frame(object)

  # Rename 'estimate' in draws to 'draw' to avoid conflict with point estimates
  if ("estimate" %in% names(draws)) {
    names(draws)[names(draws) == "estimate"] <- "draw"
  }

  # 3. Identify join keys
  # Keys are columns present in both draws and meta (excluding estimate/draw)
  common_cols <- intersect(names(draws), names(meta))
  keys <- setdiff(common_cols, c("draw", "estimate"))

  if (length(keys) == 0) {
    return(draws)
  }

  # 4. Join metadata back to draws
  draws <- merge(draws, meta, by = keys, all.x = TRUE, sort = FALSE)

  return(draws)
}
