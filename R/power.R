#' Sample Patients from Arrow Dataset
#'
#' Draws a random sample of patients from an Arrow dataset with specified
#' treatment allocation ratio.
#'
#' @param data_path Character. Path to the parquet file or Arrow dataset.
#' @param sample_size Integer. Total number of patients to sample.
#' @param allocation_ratio Numeric. Proportion of patients assigned to treatment
#'   (default: 0.5 for 1:1 allocation).
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#'
#' @return A tibble containing the sampled data with columns from the original
#'   dataset.
#'
#' @details
#' This function efficiently samples from large Arrow datasets without loading
#' the entire dataset into memory. It:
#' 1. Identifies unique patient IDs in each treatment group
#' 2. Randomly samples the specified number from each group
#' 3. Filters and collects only the selected patients
#'
#' The function assumes the dataset has an `id` column for patient identifiers
#' and a `tx` column for treatment assignment (0 = control, 1 = treatment).
#'
#' @examples
#' \dontrun{
#' # Sample 200 patients with 1:1 allocation
#' sample_data <- sample_from_arrow(
#'   data_path = "sim_data/trajectories.parquet",
#'   sample_size = 200,
#'   seed = 123
#' )
#'
#' # Sample with 2:1 treatment:control allocation
#' sample_data <- sample_from_arrow(
#'   data_path = "sim_data/trajectories.parquet",
#'   sample_size = 300,
#'   allocation_ratio = 0.67,
#'   seed = 456
#' )
#' }
#'
#' @importFrom arrow open_dataset
#' @importFrom dplyr filter select distinct collect pull
#'
#' @export
sample_from_arrow <- function(
  data_path,
  sample_size,
  allocation_ratio = 0.5,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ds <- arrow::open_dataset(data_path)

  # Get unique IDs by treatment group
  tx_ids <- ds |>
    dplyr::filter(tx == 1) |>
    dplyr::select(id) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::pull(id)

  soc_ids <- ds |>
    dplyr::filter(tx == 0) |>
    dplyr::select(id) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::pull(id)

  # Sample IDs based on allocation ratio
  n_tx <- round(allocation_ratio * sample_size)
  n_soc <- sample_size - n_tx

  tx_sample <- sample(tx_ids, n_tx, replace = FALSE)
  soc_sample <- sample(soc_ids, n_soc, replace = FALSE)
  selected_ids <- c(tx_sample, soc_sample)

  # Load only selected patients
  result <- ds |>
    dplyr::filter(id %in% selected_ids) |>
    dplyr::collect()

  return(result)
}


#' Prepare Markov Data for Proportional Odds Modeling
#'
#' Prepares longitudinal Markov trajectory data for proportional odds regression
#' by removing absorbing states and converting variables to appropriate types.
#'
#' @param data A data frame containing Markov trajectory data with columns
#'   `y` (current state) and `yprev` (previous state).
#' @param absorbing_state Integer. The absorbing state to filter out (default: 6
#'   for death). Observations where `yprev` equals this value are removed.
#' @param ordered_response Logical. Should `y` be converted to ordered factor?
#'   (default: TRUE). Required for most proportional odds models.
#' @param factor_previous Logical. Should `yprev` be converted to factor?
#'   (default: TRUE).
#'
#' @return A data frame with modified `y` and `yprev` columns and absorbing
#'   states removed.
#'
#' @details
#' This is a convenience function that performs standard data preparation for
#' proportional odds models on Markov data:
#' - Removes rows where patients were in an absorbing state at the previous time
#' - Converts current state to ordered factor (for cumulative models)
#' - Converts previous state to factor (for including as predictor)
#'
#' Users can then fit models using their preferred package (VGAM, rms, ordinal).
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' model_data <- prepare_markov_data(trajectories)
#'
#' # Fit with VGAM
#' library(VGAM)
#' fit <- vglm(y ~ tx + rcs(time, 4) + yprev,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = model_data)
#'
#' # Fit with ordinal
#' library(ordinal)
#' fit <- clm(y ~ tx + time + yprev, data = model_data)
#' }
#'
#' @importFrom dplyr filter mutate
#'
#' @export
prepare_markov_data <- function(
  data,
  absorbing_state = 6,
  ordered_response = TRUE,
  factor_previous = TRUE
) {
  result <- data |>
    dplyr::filter(yprev != absorbing_state)

  if (factor_previous) {
    result <- result |>
      dplyr::mutate(yprev = factor(yprev))
  }

  if (ordered_response) {
    result <- result |>
      dplyr::mutate(y = ordered(y))
  }

  result
}


#' Extract Tidy Coefficients from Proportional Odds Models
#'
#' Extracts coefficients with standard errors, p-values, and confidence intervals
#' from various proportional odds model objects, with standardized output format.
#'
#' @param fit A fitted proportional odds model object. Supported classes:
#'   "vglm" (VGAM), "orm" (rms), "clm" (ordinal).
#' @param se_type Character. Type of standard errors: "naive" or "robust"
#'   (default: "naive"). Robust SEs only available for rms::orm with robcov.
#' @param conf_level Numeric. Confidence level for intervals (default: 0.95).
#'
#' @return A tibble with columns:
#'   - term: Coefficient name (standardized: "tx", "yprev2", etc.)
#'   - estimate: Point estimate (log odds ratio scale)
#'   - std_error: Standard error
#'   - statistic: Test statistic (z-score)
#'   - p_value: Two-sided p-value
#'   - conf_low: Lower confidence limit
#'   - conf_high: Upper confidence limit
#'
#' @details
#' This function provides a standardized interface for extracting coefficients
#' across different proportional odds packages:
#' - **VGAM::vglm**: Extracts from vglm objects (no robust SEs)
#' - **rms::orm**: Extracts from orm objects (robust SEs if robcov was applied)
#' - **ordinal::clm**: Extracts from clm objects (no robust SEs)
#'
#' Intercepts are automatically excluded. Previous state coefficients are
#' standardized (e.g., "yprev2" for all packages).
#'
#' @examples
#' \dontrun{
#' # VGAM
#' library(VGAM)
#' fit <- vglm(y ~ tx + yprev,
#'             family = cumulative(parallel = TRUE),
#'             data = model_data)
#' tidy_po(fit)
#'
#' # rms with robust SEs
#' library(rms)
#' fit <- orm(y ~ tx + yprev, data = model_data)
#' fit_robust <- robcov(fit, cluster = model_data$id)
#' tidy_po(fit_robust, se_type = "robust")
#'
#' # ordinal
#' library(ordinal)
#' fit <- clm(y ~ tx + yprev, data = model_data)
#' tidy_po(fit)
#' }
#'
#' @importFrom dplyr case_when mutate filter
#' @importFrom tibble tibble as_tibble
#' @importFrom stats coef vcov pnorm confint
#'
#' @export
tidy_po <- function(fit, se_type = c("naive", "robust"), conf_level = 0.95) {
  se_type <- match.arg(se_type)
  alpha <- 1 - conf_level
  z_crit <- stats::qnorm(1 - alpha / 2)

  # VGAM method
  if (inherits(fit, "vglm")) {
    if (se_type == "robust") {
      warning("VGAM does not support robust standard errors. Using naive SEs.")
      se_type <- "naive"
    }

    coefs <- stats::coef(fit)
    ses <- sqrt(diag(stats::vcov(fit)))

    # Filter to non-intercept terms
    non_intercepts <- !grepl("\\(Intercept\\)", names(coefs))

    result <- tibble::tibble(
      term = names(coefs)[non_intercepts],
      estimate = coefs[non_intercepts],
      std_error = ses[non_intercepts],
      statistic = estimate / std_error,
      p_value = 2 * stats::pnorm(-abs(statistic)),
      conf_low = estimate - z_crit * std_error,
      conf_high = estimate + z_crit * std_error
    )
  } else if (inherits(fit, "orm")) {
    # rms::orm method
    coefs <- stats::coef(fit)
    ses <- sqrt(diag(stats::vcov(fit)))

    # Filter out intercepts (they start with 'y>=')
    non_intercepts <- !grepl("^y>=", names(coefs))

    result <- tibble::tibble(
      term = names(coefs)[non_intercepts],
      estimate = coefs[non_intercepts],
      std_error = ses[non_intercepts],
      statistic = estimate / std_error,
      p_value = 2 * stats::pnorm(-abs(statistic)),
      conf_low = estimate - z_crit * std_error,
      conf_high = estimate + z_crit * std_error
    )
  } else if (inherits(fit, "clm")) {
    # ordinal::clm method
    if (se_type == "robust") {
      warning(
        "ordinal::clm does not support robust standard errors. Using naive SEs."
      )
      se_type <- "naive"
    }

    summ <- summary(fit)
    coef_table <- summ$coefficients

    result <- tibble::as_tibble(coef_table, rownames = "term") |>
      dplyr::filter(!grepl("\\|", term)) |> # Remove threshold parameters
      dplyr::mutate(
        estimate = Estimate,
        std_error = `Std. Error`,
        statistic = `z value`,
        p_value = `Pr(>|z|)`,
        conf_low = estimate - z_crit * std_error,
        conf_high = estimate + z_crit * std_error
      ) |>
      dplyr::select(
        term,
        estimate,
        std_error,
        statistic,
        p_value,
        conf_low,
        conf_high
      )
  } else {
    stop("fit must be a vglm, orm, or clm object")
  }

  result
}


#' Run One Power Simulation Iteration
#'
#' Executes a complete simulation iteration: sampling patients, fitting models,
#' and extracting treatment effects across multiple analysis methods.
#'
#' @param iter_num Integer. Iteration number for tracking.
#' @param data_paths Named list of paths to datasets. Must include "markov" and
#'   optionally "t_test" and "drs". Example:
#'   `list(markov = "ha_or_0.8.parquet", t_test = "t_data_or_0.8.parquet")`.
#' @param sample_size Integer. Total number of patients to sample (default: 250).
#' @param allocation_ratio Numeric. Proportion assigned to treatment (default: 0.5).
#' @param seed Integer. Base random seed (default: 123). Actual seed will be
#'   `seed + iter_num`.
#' @param fit_functions Named list of functions that fit models and return tidy
#'   results. Each function should accept `data` and `iter` arguments and return
#'   a tibble with columns: `iter`, `analysis`, `term`, `estimate`, `std_error`,
#'   `p_value`, `conf_low`, `conf_high`.
#'
#' @return A tibble combining results from all analysis functions with columns:
#'   - iter: Iteration number
#'   - analysis: Analysis type
#'   - term: Coefficient name
#'   - estimate: Point estimate
#'   - std_error: Standard error
#'   - statistic: Test statistic
#'   - p_value: P-value
#'   - conf_low, conf_high: Confidence interval
#'
#' @details
#' This function provides a flexible framework for power simulations by accepting
#' custom fitting functions. This allows complete control over model specification,
#' package choice, and fitting options.
#'
#' Each fitting function should:
#' 1. Accept `data` (the sampled dataset) and `iter` (iteration number)
#' 2. Fit the model
#' 3. Return a tibble with results in standardized format
#'
#' @examples
#' \dontrun{
#' # Define custom fitting functions
#' fit_markov <- function(data, iter) {
#'   model_data <- prepare_markov_data(data)
#'   fit <- vglm(y ~ tx + rcs(time, 4) + yprev,
#'               family = cumulative(parallel = TRUE, reverse = TRUE),
#'               data = model_data)
#'   tidy_po(fit) |>
#'     filter(term == "tx") |>
#'     mutate(iter = iter, analysis = "markov", .before = 1)
#' }
#'
#' fit_ttest <- function(data, iter) {
#'   res <- t.test(y ~ tx, data = data)
#'   tibble(
#'     iter = iter,
#'     analysis = "t_test",
#'     term = "tx",
#'     estimate = -diff(res$estimate),
#'     std_error = res$stderr,
#'     p_value = res$p.value,
#'     conf_low = -res$conf.int[2],
#'     conf_high = -res$conf.int[1]
#'   )
#' }
#'
#' # Run iteration
#' result <- run_power_iteration(
#'   iter_num = 1,
#'   data_paths = list(
#'     markov = "sim_data/ha_or_0.8.parquet",
#'     t_test = "sim_data/t_data_or_0.8.parquet"
#'   ),
#'   sample_size = 200,
#'   fit_functions = list(markov = fit_markov, t_test = fit_ttest)
#' )
#' }
#'
#' @importFrom dplyr bind_rows
#'
#' @export
run_power_iteration <- function(
  iter_num,
  data_paths,
  sample_size = 250,
  allocation_ratio = 0.5,
  seed = 123,
  fit_functions
) {
  set.seed(seed + iter_num)

  results_list <- list()

  # Iterate over each analysis
  for (analysis_name in names(fit_functions)) {
    if (!analysis_name %in% names(data_paths)) {
      warning(
        "No data path provided for analysis: ",
        analysis_name,
        ". Skipping."
      )
      next
    }

    # Sample data
    sampled_data <- sample_from_arrow(
      data_path = data_paths[[analysis_name]],
      sample_size = sample_size,
      allocation_ratio = allocation_ratio
    )

    # Fit model and extract results
    fit_fn <- fit_functions[[analysis_name]]
    results_list[[analysis_name]] <- fit_fn(
      data = sampled_data,
      iter = iter_num
    )
  }

  dplyr::bind_rows(results_list)
}

#' Summarize Power Simulation Results
#'
#' Aggregates results from multiple simulation iterations to compute power,
#' bias, coverage, mean standard errors, and Monte Carlo standard errors.
#'
#' @param results A data frame from `run_power_iteration()` containing columns:
#'   `iter`, `analysis`, `term`, `estimate`, `std_error`, `p_value`, `conf_low`,
#'   `conf_high`.
#' @param true_value Numeric. True parameter value on the scale of estimation
#'   (e.g., log OR for treatment effect). Used to compute bias and coverage
#'   (default: NULL).
#' @param alpha Numeric. Significance level for power calculation (default: 0.05).
#'
#' @return A tibble with columns:
#'   - analysis: Analysis method
#'   - term: Coefficient name
#'   - mean_estimate: Mean point estimate across iterations
#'   - bias: Mean estimate minus true value (if true_value provided)
#'   - mean_se: Mean standard error
#'   - sd_estimate: Empirical standard deviation of estimates
#'   - power: Proportion of iterations with p < alpha
#'   - power_mcse: Monte Carlo standard error of power estimate
#'   - coverage: Proportion of CIs containing true value (if true_value provided)
#'   - n_iters: Number of iterations
#'
#' @details
#' This function computes key performance metrics for simulation studies:
#' - **Power**: Proportion rejecting null at given alpha level
#' - **Bias**: Difference between mean estimate and true value
#' - **Coverage**: Proportion of confidence intervals containing true value
#' - **Monte Carlo SE**: Uncertainty in estimates due to finite simulations
#'
#' Results are grouped by analysis method and coefficient term.
#'
#' @examples
#' \dontrun{
#' # Run multiple iterations
#' results <- purrr::map_dfr(
#'   1:1000,
#'   ~run_power_iteration(
#'     iter_num = .,
#'     data_paths = paths,
#'     fit_functions = fit_fns
#'   )
#' )
#'
#' # Summarize with true log(OR) = log(0.8)
#' summary <- summarize_power_results(
#'   results,
#'   true_value = log(0.8)
#' )
#' }
#'
#' @importFrom dplyr group_by summarise ungroup mutate
#'
#' @export
summarize_power_results <- function(results, true_value = NULL, alpha = 0.05) {
  if (!is.null(true_value)) {
    # Calculate coverage before summarizing
    results <- results |>
      dplyr::mutate(
        covers_truth = conf_low <= true_value & conf_high >= true_value
      )
  }

  summary <- results |>
    dplyr::group_by(analysis, term) |>
    dplyr::summarise(
      mean_estimate = mean(estimate),
      mean_se = mean(std_error),
      sd_estimate = sd(estimate),
      power = mean(p_value < alpha),
      power_mcse = jackknife_mcse(p_value < alpha),
      coverage = if (!is.null(true_value)) mean(covers_truth) else NA_real_,
      n_iters = n(),
      .groups = "drop"
    )

  if (!is.null(true_value)) {
    summary <- summary |>
      dplyr::mutate(
        bias = mean_estimate - true_value,
        .after = mean_estimate
      )
  }

  summary
}
