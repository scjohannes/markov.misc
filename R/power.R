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
#' Works with both VGAM::vglm and rms::orm models.
#'
#' @param data A data frame containing Markov trajectory data with columns
#'   `y` (current state) and `yprev` (previous state).
#' @param absorbing_state Integer. The absorbing state to filter out (default: 6
#'   for death). Observations where `yprev` equals this value are removed.
#' @param ordered_response Logical. Should `y` be converted to ordered factor?
#'   (default: TRUE). Required for clm and vglm, not requierd for orm.
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
#' **Package compatibility**: This function prepares data for:
#' - **VGAM::vglm**: Use `cumulative(parallel = TRUE, reverse = TRUE)` default
#'   datasets.
#' - **rms::orm**: Automatically treats integers as ordered factors.
#' - **ordinal::clm**: Only accepts ordered factors.
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' model_data <- prepare_markov_data(trajectories)
#'
#' # Fit with VGAM (no robust SEs; use bootstrap for inference)
#' library(VGAM)
#' fit <- vglm(y ~ tx + rcs(time, 4) + yprev,
#'             family = cumulative(parallel = TRUE, reverse = TRUE),
#'             data = model_data)
#'
#' # Fit with rms (supports robust SEs via robcov)
#' library(rms)
#' dd <- datadist(model_data)
#' options(datadist = "dd")
#' fit <- orm(y ~ tx + rcs(time, 4) + yprev, data = model_data, x = TRUE, y = TRUE)
#' fit_robust <- robcov(fit, cluster = model_data$id)
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
#' Supports both VGAM::vglm, rms::orm, and ordinal::clm models.
#'
#' @param fit A fitted proportional odds model object. Supported classes:
#'   "vglm" (VGAM), "orm" (rms), "clm" (ordinal). For robust standard errors with
#'   orm models, apply `rms::robcov()` to the model before passing to this function.
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
#' across different proportional odds packages:
#'
#' - **VGAM::vglm**: Extracts naive SEs using `vcov()`. For robust inference
#'   accounting for clustering, use bootstrap via `assess_operating_characteristics()`.
#'
#' - **rms::orm**: Automatically detects whether `robcov()` was applied. If
#'   `fit$var` exists (populated by `robcov()`), uses robust SEs. Otherwise,
#'   uses `vcov(fit)` for naive SEs. To get robust SEs:
#'   `fit_robust <- robcov(fit, cluster = data$id); tidy_po(fit_robust)`
#'
#' - **ordinal::clm**: Extracts naive SEs from model summary. For robust inference
#'   accounting for clustering, use bootstrap via `assess_operating_characteristics()`.
#'
#' Intercepts are automatically excluded. Previous state coefficients are
#' standardized (e.g., "yprev2" for all packages).
#'
#' **Note**: The function automatically determines the SE type for orm models.
#' For VGAM and ordinal models, bootstrap-based inference is recommended when
#' accounting for within-patient correlation is critical.
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

  # Calculate critical value and confidence interval bounds based on alternative
  if (alternative == "two.sided") {
    z_crit <- stats::qnorm(1 - alpha / 2)
  } else {
    z_crit <- stats::qnorm(1 - alpha)
  }

  # VGAM method
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

    # Calculate p-values based on alternative hypothesis
    result <- result |>
      dplyr::mutate(
        p_value = dplyr::case_when(
          alternative == "two.sided" ~ 2 * stats::pnorm(-abs(statistic)),
          alternative == "greater" ~ stats::pnorm(
            statistic,
            lower.tail = FALSE
          ),
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
  } else if (inherits(fit, "orm")) {
    # rms::orm method
    # Note: fit$var is only populated after robcov() is applied
    # For models without robcov(), fit$var is NULL and we use vcov() instead
    coefs <- fit$coefficients

    if (!is.null(fit$var)) {
      # Robust SEs available (robcov was applied)
      # fit$var contains the full covariance matrix for all coefficients
      vcov_matrix <- fit$var
      ses <- sqrt(diag(vcov_matrix))

      # Filter out intercepts (they start with 'y>=')
      non_intercepts <- !grepl("^y>=", names(coefs))

      result <- tibble::tibble(
        term = names(coefs)[non_intercepts],
        estimate = coefs[non_intercepts],
        std_error = ses[non_intercepts],
        statistic = estimate / std_error
      )
    } else {
      # Naive SEs (use vcov)
      # vcov() for orm returns a reduced matrix that may not include all intercepts
      # We need to match coefficients to the vcov matrix by names
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

    # Calculate p-values based on alternative hypothesis
    result <- result |>
      dplyr::mutate(
        p_value = dplyr::case_when(
          alternative == "two.sided" ~ 2 * stats::pnorm(-abs(statistic)),
          alternative == "greater" ~ stats::pnorm(
            statistic,
            lower.tail = FALSE
          ),
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
  } else if (inherits(fit, "clm")) {
    # ordinal::clm method
    summ <- summary(fit)
    coef_table <- summ$coefficients

    result <- tibble::as_tibble(coef_table, rownames = "term") |>
      dplyr::filter(!grepl("\\|", term)) |> # Remove threshold parameters
      dplyr::mutate(
        estimate = Estimate,
        std_error = `Std. Error`,
        statistic = `z value`
      )

    # Calculate p-values based on alternative hypothesis (override clm's two-sided p-value)
    result <- result |>
      dplyr::mutate(
        p_value = dplyr::case_when(
          alternative == "two.sided" ~ 2 * stats::pnorm(-abs(statistic)),
          alternative == "greater" ~ stats::pnorm(
            statistic,
            lower.tail = FALSE
          ),
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


#' Assess Operating Characteristics for One Iteration
#'
#' Executes a complete simulation iteration: sampling patients, fitting models,
#' and extracting treatment effects across multiple analysis methods. This function
#' can be used to assess both power (when true effect ≠ 0) and Type I error
#' (when true effect = 0).
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
#' @param inference_type Character. Type of inference: "coef" for coefficient-based
#'   (default) or "state" for state-occupancy-based. State-based inference always
#'   uses bootstrap and requires additional parameters (see below).
#' @param use_bootstrap Logical. Whether to use bootstrap for coefficient-based
#'   inference (default: FALSE). Always TRUE for state-based inference.
#' @param n_bootstrap Integer. Number of bootstrap samples (required when
#'   use_bootstrap = TRUE or inference_type = "state").
#' @param id_var Character. Name of the ID variable for group bootstrap (default: "id").
#' @param workers Integer. Number of parallel workers for bootstrap (default: NULL).
#'   When NULL, bootstrap runs sequentially. Set to > 1 to parallelize bootstrap
#'   within a single iteration (useful when running a single analysis, not when
#'   parallelizing across iterations). **Note**: Typically you parallelize across
#'   calls to `assess_operating_characteristics()` rather than within each call.
#' @param ylevels Integer vector. State levels for state-based inference (default: 1:6).
#' @param absorb Integer. Absorbing state for state-based inference (default: 6).
#' @param target_states Integer vector. Target states for state-based inference (default: 1).
#' @param varnames List. Variable names for state-based inference (default:
#'   list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")).
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
#' This function provides a flexible framework for assessing operating characteristics
#' (power, Type I error, bias, coverage) by accepting custom fitting functions.
#' This allows complete control over model specification, package choice, and
#' fitting options.
#'
#' Each fitting function should:
#' 1. Accept `data` (the sampled dataset) and `iter` (iteration number)
#' 2. Fit the model
#' 3. Return a tibble with results in standardized format
#'
#' **Inference types**:
#' 1. **Coefficient-based** (inference_type = "coef", default):
#'    - Fits models using fit_functions and extracts coefficients
#'    - Can use bootstrap or naive SEs depending on use_bootstrap
#'    - Suitable for testing treatment effects via regression coefficients
#'
#' 2. **State-based** (inference_type = "state"):
#'    - Calculates treatment effect as difference in time spent in target state(s)
#'    - Always uses bootstrap (required for CIs on state occupancy differences)
#'    - Requires a fitted model in fit_functions that returns a model object
#'    - Uses state occupancy probabilities via soprobMarkovOrdm (for orm) or similar
#'
#' **Bootstrap inference**: When `use_bootstrap = TRUE` or `inference_type = "state"`,
#' the function performs group bootstrap resampling. This is particularly useful for:
#' - vglm models (VGAM package) which don't support robust standard errors
#' - State-based inference (only way to get CIs on state occupancy differences)
#' - Accounting for within-patient correlation in coefficient estimates
#'
#' **Parallelization strategy**: The typical workflow is to parallelize across
#' simulation iterations, not within them:
#' ```r
#' # Recommended: Parallelize across iterations
#' library(furrr)
#' plan(multisession, workers = 8)
#' results <- future_map_dfr(1:1000, ~assess_operating_characteristics(
#'   iter_num = ., ...
#' ))
#' ```
#' Only set `workers > 1` when running a single analysis with many bootstrap samples
#' and NOT parallelizing across iterations.
#'
#' @examples
#' \dontrun{
#' # Coefficient-based inference (default)
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
#' result_coef <- assess_operating_characteristics(
#'   iter_num = 1,
#'   data_paths = list(markov = "sim_data/ha_or_0.8.parquet"),
#'   sample_size = 200,
#'   fit_functions = list(markov = fit_markov),
#'   inference_type = "coef",
#'   use_bootstrap = TRUE,
#'   n_bootstrap = 500
#' )
#'
#' # State-based inference (always uses bootstrap)
#' fit_markov_model <- function(data, iter) {
#'   model_data <- prepare_markov_data(data)
#'   d <- datadist(model_data)
#'   options(datadist = "d")
#'   orm(y ~ tx + rcs(time, 4) + yprev, data = model_data, x = TRUE, y = TRUE)
#' }
#'
#' result_state <- assess_operating_characteristics(
#'   iter_num = 1,
#'   data_paths = list(markov = "sim_data/ha_or_0.8.parquet"),
#'   sample_size = 200,
#'   fit_functions = list(markov = fit_markov_model),
#'   inference_type = "state",
#'   n_bootstrap = 500,
#'   target_states = 1
#' )
#' }
#'
#' @importFrom dplyr bind_rows group_by summarise mutate
#' @importFrom rsample group_bootstraps analysis
#' @importFrom future.callr callr
#' @importFrom future plan
#' @importFrom furrr future_map
#' @importFrom stats ave pnorm quantile
#'
#' @export
assess_operating_characteristics <- function(
  iter_num,
  data_paths,
  sample_size = 250,
  allocation_ratio = 0.5,
  seed = 123,
  fit_functions,
  inference_type = c("coef", "state"),
  use_bootstrap = FALSE,
  n_bootstrap = NULL,
  id_var = "id",
  workers = NULL,
  ylevels = 1:6,
  absorb = 6,
  target_states = 1,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx")
) {
  inference_type <- match.arg(inference_type)
  set.seed(seed + iter_num)

  # State-based inference always requires bootstrap
  if (inference_type == "state") {
    use_bootstrap <- TRUE
  }

  # Validate bootstrap parameters
  if (use_bootstrap && is.null(n_bootstrap)) {
    stop(
      "n_bootstrap must be specified when use_bootstrap = TRUE or inference_type = 'state'"
    )
  }

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

    if (inference_type == "state") {
      # State-based inference using bootstrap
      # Fit model on original sample
      fitted_model <- fit_functions[[analysis_name]](
        data = sampled_data,
        iter = iter_num
      )

      # Create bootstrap samples
      bs_samples <- rsample::group_bootstraps(
        sampled_data,
        group = !!rlang::sym(id_var),
        times = n_bootstrap,
        apparent = FALSE
      )

      # Set up parallel processing (only if workers specified)
      if (!is.null(workers) && workers > 1) {
        future::plan(future.callr::callr, workers = workers)
      } else {
        future::plan("sequential")
      }

      # Bootstrap state occupancy calculations
      bs_results <- bs_samples |>
        dplyr::mutate(
          state_diff = furrr::future_map(splits, \(.x) {
            boot_data <- rsample::analysis(.x)

            # Create unique IDs
            boot_data$block_count <- stats::ave(
              boot_data[[id_var]],
              boot_data[[id_var]],
              FUN = seq_along
            )
            boot_data$new_id <- factor(paste(
              boot_data[[id_var]],
              boot_data$block_count,
              sep = "_"
            ))

            # Relevel factors
            releveled <- relevel_factors_consecutive(
              data = boot_data,
              factor_cols = c("y", varnames$pvarname),
              original_data = sampled_data,
              ylevels = ylevels,
              absorb = absorb
            )

            boot_data <- releveled$data
            boot_ylevels <- releveled$ylevels
            boot_absorb <- releveled$absorb

            # Update datadist for orm models
            if (inherits(fitted_model, "orm")) {
              dd <- rms::datadist(boot_data)
              assign("dd", dd, envir = .GlobalEnv)
              options(datadist = "dd")
            }

            # Refit model
            m_boot <- tryCatch(
              stats::update(fitted_model, data = boot_data),
              error = function(e) {
                warning("Bootstrap model fit failed: ", e$message)
                return(NULL)
              }
            )

            # Calculate state occupancy difference
            if (!is.null(m_boot)) {
              tryCatch(
                taooh(
                  model = m_boot,
                  data = boot_data,
                  times = 1:max(boot_data[[varnames$tvarname]]),
                  ylevels = boot_ylevels,
                  absorb = boot_absorb,
                  target_states = target_states,
                  varnames = list(
                    tvarname = varnames$tvarname,
                    pvarname = varnames$pvarname,
                    id = "new_id",
                    tx = varnames$tx
                  )
                )[1], # Extract delta_taooh
                error = function(e) {
                  warning("State occupancy calculation failed: ", e$message)
                  return(NA_real_)
                }
              )
            } else {
              return(NA_real_)
            }
          })
        )

      future::plan("sequential")

      # Summarize bootstrap results
      bs_vector <- unlist(bs_results$state_diff)
      bs_vector <- bs_vector[!is.na(bs_vector)]

      if (length(bs_vector) > 0) {
        results_list[[analysis_name]] <- tibble::tibble(
          iter = iter_num,
          analysis = analysis_name,
          term = "state_diff",
          estimate = mean(bs_vector),
          std_error = sd(bs_vector),
          statistic = mean(bs_vector) / sd(bs_vector),
          p_value = 2 * stats::pnorm(-abs(mean(bs_vector) / sd(bs_vector))),
          conf_low = quantile(bs_vector, 0.025),
          conf_high = quantile(bs_vector, 0.975)
        )
      }
    } else if (use_bootstrap) {
      # Coefficient-based inference with bootstrap
      # Create bootstrap samples
      bs_samples <- rsample::group_bootstraps(
        sampled_data,
        group = !!rlang::sym(id_var),
        times = n_bootstrap,
        apparent = FALSE
      )

      # Apply fitting function to each bootstrap sample
      bs_results <- list()
      for (b in seq_len(nrow(bs_samples))) {
        boot_data <- rsample::analysis(bs_samples$splits[[b]])

        # Create unique IDs for bootstrap sample
        boot_data$block_count <- stats::ave(
          boot_data[[id_var]],
          boot_data[[id_var]],
          FUN = seq_along
        )
        boot_data$new_id <- factor(paste(
          boot_data[[id_var]],
          boot_data$block_count,
          sep = "_"
        ))

        # Relevel factors if needed
        if ("y" %in% names(boot_data) && "yprev" %in% names(boot_data)) {
          releveled <- relevel_factors_consecutive(
            data = boot_data,
            factor_cols = c("y", "yprev"),
            original_data = sampled_data
          )
          boot_data <- releveled$data
        }

        # Fit model on bootstrap sample
        bs_results[[b]] <- tryCatch(
          fit_functions[[analysis_name]](
            data = boot_data,
            iter = iter_num
          ),
          error = function(e) {
            warning("Bootstrap iteration ", b, " failed: ", e$message)
            NULL
          }
        )
      }

      # Combine bootstrap results
      bs_combined <- dplyr::bind_rows(bs_results)

      # Calculate bootstrap statistics
      if (nrow(bs_combined) > 0) {
        bs_summary <- bs_combined |>
          dplyr::group_by(analysis, term) |>
          dplyr::summarise(
            estimate = mean(estimate, na.rm = TRUE),
            std_error = sd(estimate, na.rm = TRUE),
            conf_low = stats::quantile(estimate, 0.025, na.rm = TRUE),
            conf_high = stats::quantile(estimate, 0.975, na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::mutate(
            statistic = estimate / std_error,
            p_value = 2 * stats::pnorm(-abs(statistic)),
            iter = iter_num,
            .before = 1
          )

        results_list[[analysis_name]] <- bs_summary
      }
    } else {
      # Standard approach without bootstrap (coefficient-based only)
      fit_fn <- fit_functions[[analysis_name]]
      results_list[[analysis_name]] <- fit_fn(
        data = sampled_data,
        iter = iter_num
      )
    }
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
