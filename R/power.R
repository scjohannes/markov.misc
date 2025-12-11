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


#' Extract Tidy Coefficients from Proportional Odds, Cox, and GLM Models
#'
#' Extracts coefficients with standard errors, p-values, and confidence intervals
#' from various model objects, with standardized output format.
#' Supports VGAM::vglm, rms::orm, ordinal::clm, survival::coxph, and stats::glm models.
#'
#' @param fit A fitted model object. Supported classes:
#'   "vglm" (VGAM), "orm" (rms), "clm" (ordinal), "coxph" (survival), "glm" (stats).
#'   For robust standard errors with orm models, apply `rms::robcov()` to the model
#'   before passing to this function.
#' @param alternative Character. Specifies the alternative hypothesis for p-value
#'   calculation: "two.sided" (default), "greater" (H₁: β > 0), or "less" (H₁: β < 0).
#'   Affects both p-values and confidence interval bounds.
#' @param conf_level Numeric. Confidence level for intervals (default: 0.95).
#'
#' @return A tibble with columns:
#'   - term: Coefficient name (standardized: "tx", "yprev2", etc.)
#'   - estimate: Point estimate (log odds ratio scale)
#'   - std_error: Standard error (robust if `robcov()` was applied to orm model)
#'   - statistic: Test statistic (z-score)
#'   - p_value: P-value (one-sided or two-sided based on `alternative`)
#'   - conf_low: Lower confidence limit (may be -Inf for one-sided tests)
#'   - conf_high: Upper confidence limit (may be Inf for one-sided tests)
#'
#' @details
#' This function provides a standardized interface for extracting coefficients
#' across different packages. Intercepts are automatically excluded from the output.
#'
#' Supported models:
#' - **VGAM::vglm**: Extracts naive SEs using `vcov()`.
#' - **rms::orm**: Automatically detects if `robcov()` was applied.
#' - **ordinal::clm**: Extracts naive SEs from model summary.
#' - **survival::coxph**: Extracts SEs (robust if cluster/weights used).
#' - **stats::glm**: Extracts SEs using `vcov()`. Useful for binary logistic regression.
#'
#' @examples
#' \dontrun{
#' # VGAM (naive SEs only)
#' library(VGAM)
#' fit <- vglm(y ~ tx + yprev,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = model_data)
#' tidy_po(fit)  # Two-sided test (default)
#' tidy_po(fit, alternative = "greater")  # One-sided: H₁: β > 0
#'
#' # rms with naive SEs
#' library(rms)
#' fit <- orm(y ~ tx + yprev, data = model_data, x = TRUE, y = TRUE)
#' tidy_po(fit)  # Naive SEs
#'
#' # rms with robust SEs
#' fit_robust <- robcov(fit, cluster = model_data$id)
#' tidy_po(fit_robust)  # Automatically uses robust SEs
#'
#' # ordinal (naive SEs only)
#' library(ordinal)
#' fit <- clm(y ~ tx + yprev, data = model_data)
#' tidy_po(fit, alternative = "less")  # One-sided: H₁: β < 0
#' }
#'
#' @importFrom dplyr case_when mutate filter
#' @importFrom tibble tibble as_tibble
#' @importFrom stats coef vcov pnorm confint
#'
#' @export
tidy_po <- function(
  fit,
  alternative = c("two.sided", "greater", "less"),
  conf_level = 0.95
) {
  alternative <- match.arg(alternative)
  alpha <- 1 - conf_level

  # Calculate critical value based on alternative
  if (alternative == "two.sided") {
    z_crit <- stats::qnorm(1 - alpha / 2)
  } else {
    z_crit <- stats::qnorm(1 - alpha)
  }

  # --- VGAM method ---
  if (inherits(fit, "vglm")) {
    coefs <- stats::coef(fit)
    ses <- sqrt(diag(stats::vcov(fit)))

    # Filter to non-intercept terms
    non_intercepts <- !grepl("\\(Intercept\\)", names(coefs))

    result <- tibble::tibble(
      term = names(coefs)[non_intercepts],
      estimate = coefs[non_intercepts],
      std_error = ses[non_intercepts],
      statistic = estimate / std_error
    )

    # --- rms::orm method ---
  } else if (inherits(fit, "orm")) {
    # rms::orm method
    # Note: fit$var is only populated after robcov() is applied
    # For models without robcov(), fit$var is NULL and we use vcov() instead
    coefs <- fit$coefficients

    if (!is.null(fit$var)) {
      # Robust SEs available
      vcov_matrix <- fit$var
      ses <- sqrt(diag(vcov_matrix))
      non_intercepts <- !grepl("^y>=", names(coefs))

      result <- tibble::tibble(
        term = names(coefs)[non_intercepts],
        estimate = coefs[non_intercepts],
        std_error = ses[non_intercepts],
        statistic = estimate / std_error
      )
    } else {
      # Naive SEs
      vcov_matrix <- stats::vcov(fit)
      vcov_names <- rownames(vcov_matrix)

      # Extract only the coefficients that are in the vcov matrix
      # and exclude intercepts
      non_intercept_names <- vcov_names[!grepl("^y>=", vcov_names)]

      result <- tibble::tibble(
        term = non_intercept_names,
        estimate = coefs[non_intercept_names],
        std_error = sqrt(diag(vcov_matrix))[non_intercept_names],
        statistic = estimate / std_error
      )
    }

    # --- ordinal::clm method ---
  } else if (inherits(fit, "clm")) {
    summ <- summary(fit)
    coef_table <- summ$coefficients

    result <- tibble::as_tibble(coef_table, rownames = "term") |>
      dplyr::filter(!grepl("\\|", term)) |>
      dplyr::mutate(
        estimate = Estimate,
        std_error = `Std. Error`,
        statistic = `z value`
      ) |>
      dplyr::select(term, estimate, std_error, statistic)

    # --- survival::coxph method (Supports Fine-Gray) ---
  } else if (inherits(fit, "coxph")) {
    # vcov(fit) automatically returns the robust variance if weights/cluster were used

    coefs <- stats::coef(fit)
    # Handle cases where coefficients might be NA (rare but possible in rank deficient models)
    if (any(is.na(coefs))) {
      warning("Some coefficients in coxph model are NA")
    }

    vcov_mat <- stats::vcov(fit)
    ses <- sqrt(diag(vcov_mat))

    result <- tibble::tibble(
      term = names(coefs),
      estimate = coefs,
      std_error = ses,
      statistic = estimate / std_error
    )
    # --- stats::glm method ---
  } else if (inherits(fit, "glm")) {
    coefs <- stats::coef(fit)
    ses <- sqrt(diag(stats::vcov(fit)))

    # Filter out Intercepts
    non_intercepts <- !grepl("\\(Intercept\\)", names(coefs))

    result <- tibble::tibble(
      term = names(coefs)[non_intercepts],
      estimate = coefs[non_intercepts],
      std_error = ses[non_intercepts],
      statistic = estimate / std_error
    )
  } else {
    stop("fit must be a vglm, orm, clm, coxph, or glm object")
  }

  # --- Common Calculation for P-values and CIs ---
  # Apply to the 'result' tibble regardless of which method generated it
  result <- result |>
    dplyr::mutate(
      p_value = dplyr::case_when(
        alternative == "two.sided" ~ 2 * stats::pnorm(-abs(statistic)),
        alternative == "greater" ~ stats::pnorm(statistic, lower.tail = FALSE),
        alternative == "less" ~ stats::pnorm(statistic)
      ),
      conf_low = dplyr::case_when(
        alternative == "two.sided" ~ estimate - z_crit * std_error,
        alternative == "greater" ~ estimate - z_crit * std_error,
        alternative == "less" ~ -Inf
      ),
      conf_high = dplyr::case_when(
        alternative == "two.sided" ~ estimate + z_crit * std_error,
        alternative == "greater" ~ Inf,
        alternative == "less" ~ estimate + z_crit * std_error
      )
    )

  return(result)
}

#' Assess Operating Characteristics for One Iteration
#'
#' Executes a complete simulation iteration by sampling patients once, then
#' applying multiple fitting functions to the same sample. Bootstrap inference
#' is handled within individual fitting functions, not by this orchestrator.
#'
#' @param iter_num Integer. Iteration number for tracking.
#' @param data_paths Named list of paths to datasets. Names must match those in
#'   fit_functions. Example:
#'   `list(markov = "ha_or_0.8.parquet", t_test = "t_data_or_0.8.parquet")`.
#' @param sample_size Integer. Total number of patients to sample (default: 250).
#' @param allocation_ratio Numeric. Proportion assigned to treatment (default: 0.5).
#' @param seed Integer. Base random seed (default: 123). Actual seed will be
#'   `seed + iter_num`.
#' @param output_path Character. Path to directory where additional results will
#'   be saved (default: NULL). If NULL and fit functions return additional results,
#'   an error will be raised. A subdirectory `details/` will be created under this
#'   path, with further subdirectories for each analysis type.
#' @param fit_functions Named list of functions that fit models and return results.
#'   Each function receives `data` (sampled data) and `iter` (iteration number)
#'   arguments. Must return a **list** where:
#'   - First element: A one-row tibble/data frame with summary statistics containing
#'     columns: `iter`, `analysis`, `se_type`, `term`, `estimate`, `std_error`,
#'     `statistic`, `p_value`, `conf_low`, `conf_high`.
#'   - Remaining elements (optional): Additional results to save to disk (e.g.,
#'     bootstrap distributions, SOPs, detailed model output). These will be saved
#'     as RDS files in `output_path/details/analysis_name/analysis_name_iter_N.rds`.
#'
#' @return A tibble combining results from all analysis functions with columns:
#'   - iter: Iteration number
#'   - analysis: Analysis type
#'   - se_type: Type of standard error ("naive", "boot", etc.)
#'   - term: Coefficient name
#'   - estimate: Point estimate
#'   - std_error: Standard error
#'   - statistic: Test statistic
#'   - p_value: P-value
#'   - conf_low, conf_high: Confidence interval
#'
#' @details
#' This function implements the correct simulation structure for assessing
#' operating characteristics (power, Type I error, bias, coverage):
#'
#' 1. **Iteration-level sampling**: Draws a fixed sample of patients from the
#'    superpopulation for this iteration (using `seed + iter_num`)
#' 2. **Multiple analyses**: Applies all fitting functions to the *same* sample
#' 3. **Within-analysis inference**: Individual fitting functions handle their
#'    own inference strategy (naive SEs, robust SEs, bootstrap, etc.)
#' 4. **Save additional results**: If fitting functions return more than just
#'    summary statistics, additional results are saved to disk as RDS files
#'
#' **Key principle**: Each iteration gets a different sample from the superpopulation,
#' but all analyses within an iteration see the same sample. This allows fair
#' comparison of different methods on identical data.
#'
#' **Fitting function requirements**:
#' Each function should:
#' 1. Accept `data` (the sampled dataset) and `iter` (iteration number)
#' 2. Perform any necessary data preparation
#' 3. Fit the model (with or without bootstrap, as appropriate)
#' 4. Return a **list** where:
#'    - First element: One-row tibble with summary statistics
#'    - Remaining elements (optional): Additional data to save (bootstrap samples,
#'      SOPs, detailed output, etc.)
#'
#' **File saving**:
#' When fitting functions return lists with more than one element, all elements
#' after the first are saved to disk at:
#' `output_path/details/analysis_name/analysis_name_iter_N.rds`
#'
#' The saved object is a list containing `results[-1]` (all elements except the
#' first). Directories are created recursively as needed.
#'
#' **Bootstrap handling**: Bootstrap inference (if needed) should be implemented
#' within individual fitting functions, not by this orchestrator. This design:
#' - Keeps parallelization simple (across iterations only)
#' - Allows complete control over inference within each analysis
#' - Avoids nested parallelization issues
#'
#' **Parallelization strategy**: The typical workflow is to parallelize across
#' simulation iterations:
#' ```r
#' library(furrr)
#' plan(multisession, workers = 8)
#' results <- future_map_dfr(
#'   1:1000,
#'   ~assess_operating_characteristics(iter_num = ., ...),
#'   .options = furrr_options(seed = TRUE)
#' )
#' ```
#'
#' @examples
#' \dontrun{
#' # Define fitting functions that receive sampled data
#'
#' # Markov analysis with bootstrap - returns list with summary + bootstrap details
#' fit_markov_boot <- function(data, iter) {
#'   model_data <- prepare_markov_data(data, absorbing_state = 6)
#'
#'   # Fit model
#'   fit <- vglm(
#'     y ~ tx + rcs(time, 4) + yprev,
#'     family = cumulative(parallel = TRUE, reverse = TRUE),
#'     data = model_data
#'   )
#'
#'   # Bootstrap within this function
#'   boot_coef <- bootstrap_model_coefs(
#'     fit,
#'     data = model_data,
#'     n_boot = 100,
#'     id_var = "id"
#'   )
#'
#'   # Tidy bootstrap results and add required columns
#'   tidy_bootstrap_coefs(boot_coef, probs = c(0.025, 0.975), estimate = "median") |>
#'     filter(term == "tx") |>
#'     mutate(
#'       iter = iter,
#'       analysis = "markov",
#'       se_type = "boot",
#'       # Calculate bootstrap SE from confidence limits (approximate)
#'       std_error = NA_real_,
#'       statistic = NA_real_,
#'       p_value = NA_real_,
#'       conf_low = lower,
#'       conf_high = upper,
#'       .before = 1
#'     ) |>
#'     select(iter, analysis, se_type, term, estimate, std_error,
#'            statistic, p_value, conf_low, conf_high)
#'
#'   # Return list: summary + additional results to save
#'   list(
#'     summary_results,
#'     list(
#'       boot_coefs = boot_coef,
#'       model_fit = fit
#'     )
#'   )
#' }
#'
#' # T-test analysis (no additional results) - returns list with just summary
#' fit_ttest <- function(data, iter) {
#'   t_result <- t.test(y ~ tx, data = data)
#'
#'   summary_results <- tibble(
#'     iter = iter,
#'     analysis = "t_test",
#'     se_type = "naive",
#'     term = "tx",
#'     estimate = -diff(t_result$estimate),
#'     std_error = t_result$stderr,
#'     statistic = t_result$statistic,
#'     p_value = t_result$p.value,
#'     conf_low = -t_result$conf.int[2],
#'     conf_high = -t_result$conf.int[1]
#'   )
#'
#'   # Return list with just summary (no additional results to save)
#'   list(summary_results)
#' }
#'
#' # Parallelize across iterations - each gets different sample
#' library(furrr)
#' plan(multisession, workers = 8)
#' results <- future_map_dfr(
#'   1:1000,
#'   ~assess_operating_characteristics(
#'     iter_num = .,
#'     data_paths = list(
#'       markov = "ha_or_0.8.parquet",
#'       t_test = "t_data_or_0.8.parquet"
#'     ),
#'     output_path = "power_sim_results",
#'     sample_size = 250,
#'     fit_functions = list(
#'       markov = fit_markov_boot,
#'       t_test = fit_ttest
#'     ),
#'     seed = 123
#'   ),
#'   .options = furrr_options(seed = TRUE, packages = c("rms", "VGAM", "dplyr"))
#' )
#' }
#'
#' @importFrom dplyr bind_rows
#'
#' @export
assess_operating_characteristics <- function(
  iter_num,
  data_paths,
  output_path,
  fit_functions,
  sample_size = 250,
  allocation_ratio = 0.5,
  seed = 123,
  rerandomize = FALSE,
  id_var = "id",
  tx_var = "tx"
) {
  if (!is.null(seed)) {
    set.seed(seed + iter_num)
  }

  if (!is.null(output_path)) {
    out_dir <- paste0(output_path, "/details")
  }

  results_list <- list()

  # Sample data once per analysis type and apply fitting function
  for (analysis_name in names(fit_functions)) {
    if (!analysis_name %in% names(data_paths)) {
      warning(
        "No data path provided for analysis: ",
        analysis_name,
        ". Skipping."
      )
      next
    }

    # Sample data once for this iteration
    sampled_data <- sample_from_arrow(
      data_path = data_paths[[analysis_name]],
      sample_size = sample_size,
      allocation_ratio = allocation_ratio
      # Note: No seed argument - uses set.seed() from above
    )

    if (rerandomize) {
      ids <- unique(sampled_data[[id_var]])
      tx_ids <- sample(
        ids,
        size = round(allocation_ratio * length(ids)),
        replace = FALSE
      )
      sampled_data[[tx_var]] <- ifelse(
        sampled_data[[id_var]] %in% tx_ids,
        1,
        0
      )
    }

    # Apply fitting function to this specific sample
    # All of the fit functions must return a list
    # The first element should always be a one row tibble (or df) of summary
    # results
    # The second element contains everything else (may be multiplpe data frames)
    # The second element should be saved to disk
    results <- fit_functions[[analysis_name]](
      data = sampled_data,
      iter = iter_num
    )

    if (is.data.frame(results)) {
      stop(
        "Fitting function must return a list with the first element as a data frame of results. You likely used an old fit function."
      )
    }

    # Save additional results to disk if present
    if (length(results) > 1) {
      if (is.null(output_path)) {
        stop(
          "output_path is NULL, but fitting function returned additional results to save. ",
          analysis_name
        )
      } else {
        # Create directory recursively (including parent directories)
        analysis_dir <- file.path(out_dir, analysis_name)
        if (!dir.exists(analysis_dir)) {
          dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
        }

        # Save all additional results (everything except first element)
        saveRDS(
          results[-1],
          file = file.path(
            analysis_dir,
            paste0(analysis_name, "_iter_", iter_num, ".rds")
          )
        )
      }
    }

    # Store summary results (first element)
    results_list[[analysis_name]] <- results[[1]]
  }

  dplyr::bind_rows(results_list)
}

#' @rdname assess_operating_characteristics
#' @export
run_power_iteration <- assess_operating_characteristics

#' Summarize Operating Characteristics from Simulation Results
#'
#' Aggregates results from multiple simulation iterations to compute power,
#' Type I error, bias, coverage, mean standard errors, and Monte Carlo standard errors.
#'
#' @param results A data frame from `assess_operating_characteristics()` containing
#'   columns: `iter`, `analysis`, `term`, `estimate`, `std_error`, `p_value`,
#'   `conf_low`, `conf_high`.
#' @param true_value Numeric. True parameter value on the scale of estimation
#'   (e.g., log OR for treatment effect). Used to compute bias and coverage
#'   (default: NULL). If true_value = 0, the "power" metric actually represents
#'   Type I error rate.
#' @param alpha Numeric. Significance level (default: 0.05).
#'
#' @return A tibble with columns:
#'   - analysis: Analysis method
#'   - term: Coefficient name
#'   - mean_estimate: Mean point estimate across iterations
#'   - bias: Mean estimate minus true value (if true_value provided)
#'   - mean_se: Mean standard error
#'   - sd_estimate: Empirical standard deviation of estimates
#'   - power: Proportion of iterations with p < alpha (Type I error if true_value = 0)
#'   - power_mcse: Monte Carlo standard error of power estimate
#'   - coverage: Proportion of CIs containing true value (if true_value provided)
#'   - n_iters: Number of iterations
#'
#' @details
#' This function computes key operating characteristics for simulation studies:
#' - **Power/Type I error**: Proportion rejecting null at given alpha level
#'   - When true_value ≠ 0: measures power (ability to detect true effects)
#'   - When true_value = 0: measures Type I error (false positive rate)
#' - **Bias**: Difference between mean estimate and true value
#' - **Coverage**: Proportion of confidence intervals containing true value
#' - **Monte Carlo SE**: Uncertainty in estimates due to finite simulations
#'
#' Results are grouped by analysis method and coefficient term.
#'
#' @examples
#' \dontrun{
#' # Run multiple iterations for power analysis
#' results <- purrr::map_dfr(
#'   1:1000,
#'   ~assess_operating_characteristics(
#'     iter_num = .,
#'     data_paths = paths,
#'     fit_functions = fit_fns
#'   )
#' )
#'
#' # Summarize with true log(OR) = log(0.8)
#' summary <- summarize_oc_results(
#'   results,
#'   true_value = log(0.8)
#' )
#'
#' # For Type I error assessment
#' summary_t1e <- summarize_oc_results(
#'   results,
#'   true_value = 0  # Null hypothesis
#' )
#' }
#'
#' @importFrom dplyr group_by summarise ungroup mutate
#'
#' @export
summarize_oc_results <- function(results, true_value = NULL, alpha = 0.05) {
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

#' @rdname summarize_oc_results
#' @export
summarize_power_results <- summarize_oc_results
