#' Tidy Bootstrap Coefficient Estimates
#'
#' Summarizes bootstrap coefficient estimates by computing confidence intervals
#' and point estimates. Typically used with output from
#' \code{\link{bootstrap_model_coefs}}.
#'
#' @param boot_coefs A data frame or tibble containing bootstrap coefficient
#'   estimates, typically the output from \code{\link{bootstrap_model_coefs}}.
#'   Should have a column for bootstrap iteration ID (default: "boot_id") and
#'   one column per model coefficient.
#' @param id_col Name of the column containing bootstrap iteration IDs
#'   (default: "boot_id"). This column will be excluded from the summary.
#' @param probs Numeric vector of length 2 specifying the lower and upper
#'   quantiles for confidence intervals. Default is c(0.025, 0.975) for 95% CI.
#'   Can be changed to other levels, e.g., c(0.05, 0.95) for 90% CI.
#' @param estimate Character string specifying the point estimate to use:
#'   "median" (default), "mean", or NULL to omit point estimate. This is
#'   computed separately from the quantiles in \code{probs}.
#' @param na.rm Logical indicating whether to remove NA values before computing
#'   statistics. Default is FALSE, which will result in NA output if any
#'   bootstrap iteration failed. Set to TRUE to compute statistics on
#'   non-missing values only.
#'
#' @return A tibble with one row per coefficient containing:
#'   \itemize{
#'     \item term: Coefficient name
#'     \item estimate: Point estimate (median or mean, if requested)
#'     \item lower: Lower confidence limit (first value in \code{probs})
#'     \item upper: Upper confidence limit (second value in \code{probs})
#'   }
#'
#' @details
#' This function computes quantile-based confidence intervals from bootstrap
#' coefficient estimates. The default settings provide:
#' \itemize{
#'   \item Median as point estimate
#'   \item 2.5th percentile as lower bound
#'   \item 97.5th percentile as upper bound
#' }
#'
#' This corresponds to the percentile bootstrap confidence interval method,
#' which is appropriate when the bootstrap distribution is roughly symmetric
#' and unbiased.
#'
#' For highly skewed or biased bootstrap distributions, consider using
#' bias-corrected and accelerated (BCa) intervals instead (not currently
#' implemented in this function).
#'
#' **Usage with assess_operating_characteristics()**: When using this function
#' within fitting functions for \code{\link{assess_operating_characteristics}},
#' you'll need to add additional columns and rename \code{lower}/\code{upper}
#' to \code{conf_low}/\code{conf_high}. See examples below.
#'
#' @examples
#' \dontrun{
#' # After running bootstrap
#' boot_coefs <- bootstrap_model_coefs(
#'   model = m,
#'   data = my_data,
#'   n_boot = 1000
#' )
#'
#' # Get 95% CI with median
#' tidy_boot <- tidy_bootstrap_coefs(boot_coefs)
#'
#' # Get 90% CI with mean
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   probs = c(0.05, 0.95),
#'   estimate = "mean"
#' )
#'
#' # Get 99% CI with median
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   probs = c(0.005, 0.995)
#' )
#'
#' # Remove NA values if some bootstrap iterations failed
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   na.rm = TRUE
#' )
#'
#' # Format for assess_operating_characteristics()
#' tidy_boot <- tidy_bootstrap_coefs(
#'   boot_coefs,
#'   probs = c(0.025, 0.975),
#'   estimate = "median"
#' ) |>
#'   mutate(
#'     iter = 1,
#'     analysis = "markov",
#'     se_type = "boot",
#'     std_error = NULL,
#'     statistic = NULL,
#'     p_value = NULL,
#'     conf_low = lower,
#'     conf_high = upper,
#'     .before = 1
#'   ) |>
#'   select(iter, analysis, se_type, term, estimate, std_error,
#'          statistic, p_value, conf_low, conf_high)
#' }
#'
#' @importFrom stats quantile median
#'
#' @export
tidy_bootstrap_coefs <- function(
  boot_coefs,
  id_col = "boot_id",
  probs = c(0.025, 0.975),
  estimate = "median",
  na.rm = TRUE
) {
  # Input validation
  if (!is.data.frame(boot_coefs)) {
    stop("boot_coefs must be a data frame or tibble")
  }

  if (!id_col %in% names(boot_coefs)) {
    stop("id_col '", id_col, "' not found in boot_coefs")
  }

  if (
    !is.numeric(probs) || length(probs) != 2 || any(probs < 0) || any(probs > 1)
  ) {
    stop(
      "probs must be a numeric vector of length 2 with values between 0 and 1"
    )
  }

  if (!is.null(estimate) && !estimate %in% c("median", "mean")) {
    stop("estimate must be 'median', 'mean', or NULL")
  }

  if (!is.logical(na.rm) || length(na.rm) != 1) {
    stop("na.rm must be a single logical value (TRUE or FALSE)")
  }

  # Ensure probs are sorted
  probs <- sort(probs)

  coef_cols <- setdiff(names(boot_coefs), id_col)
  result_long <- bind_rows_fill(lapply(coef_cols, function(term) {
    x <- boot_coefs[[term]]
    out <- data.frame(term = term, check.names = FALSE)
    if (!is.null(estimate)) {
      out$estimate <- if (estimate == "median") {
        stats::median(x, na.rm = na.rm)
      } else {
        mean(x, na.rm = na.rm)
      }
    }
    out$lower <- as.numeric(stats::quantile(x, probs[1], na.rm = na.rm))
    out$upper <- as.numeric(stats::quantile(x, probs[2], na.rm = na.rm))
    out$n_iter <- length(x)
    out
  }))

  return(result_long)
}
