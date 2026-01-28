#' Calculate State Occupation Probabilities for First-Order Markov Models
#'
#' Estimates state occupation probabilities over time by iterating a transition
#' matrix derived from a fitted model object (e.g., `vglm`, `rms`, or `rmsb`).
#' This function supports linear time iteration as well as complex, non-linear
#' time specifications (e.g., splines) via a covariate lookup table.
#'
#' @param object A fitted model object. Supported classes include:
#'   \code{"lrm"}, \code{"orm"} (from package `rms`),
#'   \code{"blrm"} (from package `rmsb`), and
#'   \code{"vglm"}, \code{"vgam"} (from package `VGAM`).
#' @param data A data frame containing the baseline covariates for the prediction.
#'   Rows represent unique patients. Columns must contain baseline covariates and
#'   initial values for time-varying variables.
#' @param times A numeric vector of time points to iterate over.
#'   The function calculates probabilities for each time point in this sequence.
#' @param ylevels A character vector defining the names of the outcome levels (states).
#'   These must match the levels used in the fitted `object`.
#' @param absorb (Optional) A character vector of absorbing states (states from which
#'   transitions out are impossible). Defaults to \code{NULL}.
#' @param tvarname A character string specifying the name of the time variable in the
#'   model formula. Defaults to \code{"time"}.
#' @param pvarname A character string specifying the name of the previous state variable
#'   used in the model formula. Defaults to \code{"yprev"}.
#' @param gap (Optional) A character string specifying the name of the variable representing
#'   the time gap (delta time) between observations, if used in the model. Defaults to \code{NULL}.
#' @param t_covs (Optional) A data frame used for non-linear time handling (e.g., splines
#'   or basis functions).
#'   \itemize{
#'     \item **Structure:** The number of rows in \code{t_covs} must exactly match the length of \code{times}.
#'     \item **Columns:** Column names must match the specific basis variables used in the model formula (e.g., \code{t1}, \code{t2}).
#'     \item **Usage:** At step \code{i}, the values from the \code{i}-th row of \code{t_covs} are injected into the prediction data.
#'   }
#'   This allows for `vglm` models where basis functions are supplied as separate columns rather than calculated on the fly.
#'
#' @details
#'
#' \strong{1. Data Expansion:}
#' We construct a "long" expansion dataset at every time point. For \eqn{N} patients
#' and \eqn{K} non-absorbing states, the expansion dataset contains \eqn{N \times K} rows.
#' \itemize{
#'   \item Rows \eqn{1 \dots N}: All patients assuming \eqn{y_{prev} = \text{State } 1}
#'   \item Rows \eqn{(N+1) \dots 2N}: All patients assuming \eqn{y_{prev} = \text{State } 2}
#'   \item ... and so on.
#' }
#' This allows a single \code{predict()} call to generate transition probabilities for the entire
#' cohort for all possible previous states in one step.
#'
#' \strong{2. Element-wise weighted sum}
#' We use an element-wise weighted sum approach based on the
#' Law of Total Probability:
#' \deqn{P(S_t = k) = \sum_{j} P(S_{t-1} = j) \times P(S_t = k | S_{t-1} = j)}
#'
#' where:
#' \itemize{
#'   \item \eqn{S_t}: State occupied at time \eqn{t}.
#'   \item \eqn{S_{t-1}}: State occupied at time \eqn{t-1}.
#'   \item \eqn{k}: The target state at time \eqn{t}.
#'   \item \eqn{j}: The origin state at time \eqn{t-1}.
#' }
#'
#' Let \eqn{\mathbf{v}_{prev, j}} be a vector of length \eqn{N} containing the probability that each
#' patient was in state \eqn{j} at time \eqn{t-1}.
#' Let \eqn{\mathbf{M}_{trans, j \to k}} be a vector of length \eqn{N} containing the transition
#' probability from \eqn{j} to \eqn{k} for each patient (derived from the batched prediction).
#'
#' The probability of being in state \eqn{k} at time \eqn{t} is updated as:
#' \deqn{\mathbf{v}_{curr, k} = \sum_{j} (\mathbf{v}_{prev, j} \odot \mathbf{M}_{trans, j \to k})}
#' where \eqn{\odot} denotes element-wise multiplication.
#'
#' \strong{3. Absorbing States:}
#' If an absorbing state \eqn{a} is present, the update logic handles the accumulation of probability mass:
#' \deqn{P(S_t = a) = P(S_{t-1} = a) + \sum_{j \neq a} P(S_{t-1} = j) \times P(S_t = a | S_{t-1} = j)}
#'
#' @return
#' An array of state probabilities.
#' \itemize{
#'   \item **Frequentist fit:** An array of dimension \code{[n_patients x n_times x n_states]}.
#'   \item **Bayesian fit:** An array of dimension \code{[n_draws x n_patients x n_times x n_states]}.
#' }
#'
#' @export
soprob_markov <- function(
  object,
  data,
  times,
  ylevels,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  gap = NULL,
  t_covs = NULL
) {
  # --- 1. Initial Checks & Setup ---
  cl <- class(object)[1]
  ftypes <- c(
    lrm = "rms",
    orm = "rms",
    blrm = "rmsb",
    vglm = "vgam",
    vgam = "vgam",
    robcov_vglm = "robcov"
  )
  ftype <- ftypes[cl]

  if (is.na(ftype)) {
    stop("Object class not supported")
  }

  # Define prediction function
  prd <- switch(
    ftype,
    rms = function(obj, d) stats::predict(obj, d, type = "fitted.ind"),
    vgam = function(obj, d) VGAM::predict(obj, d, type = "response"),
    rmsb = function(obj, d) {
      stats::predict(obj, d, type = "fitted.ind", posterior.summary = "all")
    },
    robcov = function(obj, d) {
      if (is.null(obj$vglm_fit)) {
        stop(
          "robcov_vglm object does not contain the original vglm fit. ",
          "Please re-run robcov_vglm() with the latest version of the package."
        )
      }
      VGAM::predict(obj$vglm_fit, d, type = "response")
    }
  )

  # Prepare dimensions
  n_pat <- nrow(data)
  n_times <- length(times)
  n_states <- length(ylevels)
  yna <- setdiff(ylevels, absorb) # Non-absorbing states
  n_yna <- length(yna)

  # Check Bayesian draws
  nd <- if (ftype == "rmsb" && length(object$draws)) nrow(object$draws) else 0

  # Initialize Output Array
  # Structure: [Patients, Time, States]
  if (nd == 0) {
    P <- array(
      0,
      dim = c(n_pat, n_times, n_states),
      dimnames = list(rownames(data), times, ylevels)
    )
  } else {
    P <- array(
      0,
      dim = c(nd, n_pat, n_times, n_states),
      dimnames = list(1:nd, rownames(data), times, ylevels)
    )
  }

  # --- 2. Time 1 Initialization (Start) ---
  # Update time variables for T1
  data[[tvarname]] <- times[1]
  if (!is.null(gap)) {
    data[[gap]] <- times[1]
  }

  # Inject T1 basis functions
  if (!is.null(t_covs)) {
    for (nm in names(t_covs)) {
      data[[nm]] <- t_covs[1, nm]
    }
  }

  # Predict probabilities at T1
  p_t1 <- prd(object, data) # Returns [n_pat x n_states] or [nd x n_pat x n_states]

  if (nd == 0) {
    P[, 1, ] <- p_t1
  } else {
    P[,, 1, ] <- p_t1
  }

  # --- 3. Prepare Expansion Data for Transitions ---
  # We create a long dataframe where every patient is repeated for every possible PREVIOUS state.
  # Order: Patient 1 (State A), Patient 2 (State A)... Patient 1 (State B)...
  # This ordering is crucial for the vectorized update logic later.

  edata_base <- data[rep(1:n_pat, times = n_yna), , drop = FALSE]

  # Assign the previous state variable
  # We repeat each state n_pat times

  if (is.factor(data[[pvarname]])) {
    # If the original variable is a factor, ensure the new column is also a factor
    # with the same levels
    edata_base[[pvarname]] <- factor(
      rep(yna, each = n_pat),
      levels = levels(data[[pvarname]])
    )
  } else {
    edata_base[[pvarname]] <- rep(yna, each = n_pat)
  }

  # --- 4. Iterate Through Time (Vectorized over Patients) ---
  for (it in 2:n_times) {
    # Update time in the expanded dataset
    edata_base[[tvarname]] <- times[it]

    if (!is.null(gap)) {
      edata_base[[gap]] <- times[it] - times[it - 1]
    }

    if (!is.null(t_covs)) {
      for (nm in names(t_covs)) {
        edata_base[[nm]] <- t_covs[it, nm]
      }
    }

    # Get Transition Probabilities
    # This returns matrix: [ (n_pat * n_yna) x n_states ]
    # Rows 1:N are transitions from State 1, Rows N+1:2N from State 2, etc.
    trans_probs <- prd(object, edata_base)

    # --- 5. The Update Step (Markov) ---
    # Formula: P(S_t = k) = Sum_over_j [ P(S_t-1 = j) * P(S_t=k | S_t-1=j) ]

    if (nd == 0) {
      # Initialize current time probabilities with 0
      p_current <- matrix(0, nrow = n_pat, ncol = n_states)
      colnames(p_current) <- ylevels

      # Add contribution from Non-Absorbing Previous States
      for (i in seq_along(yna)) {
        prev_state_name <- yna[i]

        # Extract prob of being in this previous state for all patients
        # Vector of length n_pat
        prob_prev <- P[, it - 1, prev_state_name]

        # Extract transition probs GIVEN this previous state
        # Rows correspond to the block for this state in edata_base
        row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
        probs_transition <- trans_probs[row_indices, ]

        # Weighted sum: multiply column-wise by the probability of being in prev state
        p_current <- p_current + (probs_transition * prob_prev)
      }

      # Add contribution from Absorbing States (if any)
      # Absorbing states transition to themselves with Prob 1
      if (!is.null(absorb)) {
        for (a_state in absorb) {
          # If you were in absorb state at t-1, you are in absorb state at t
          p_current[, a_state] <- p_current[, a_state] + P[, it - 1, a_state]
        }
      }

      P[, it, ] <- p_current
    } else {
      # --- Bayesian Block (Complex due to extra dimension) ---
      # trans_probs is [nd x (n_pat*n_yna) x n_states]

      # Iterate draws (hard to vectorize draws + patients + states simultaneously without memory issues)
      for (d in 1:nd) {
        p_current_d <- matrix(0, nrow = n_pat, ncol = n_states)

        for (i in seq_along(yna)) {
          prev_state_name <- yna[i]
          prob_prev <- P[d, , it - 1, prev_state_name] # [n_pat]

          row_indices <- ((i - 1) * n_pat + 1):(i * n_pat)
          probs_transition <- trans_probs[d, row_indices, ] # [n_pat x n_states]

          p_current_d <- p_current_d + (probs_transition * prob_prev)
        }

        if (!is.null(absorb)) {
          for (a_state in absorb) {
            # total dead at t = new deaths + already dead
            p_current_d[, a_state] <- p_current_d[, a_state] +
              P[d, , it - 1, a_state]
          }
        }
        P[d, , it, ] <- p_current_d
      }
    }
  }

  return(P)
}

#' Compute standardized state occupancy probabilities (Vectorized)
#'
#' Updates the standardize_sops wrapper to use the optimized
#' soprob_markov_vectorized function.
#'
#' @param model A fitted model object (orm, vglm, rmsb).
#' @param data A data frame containing patient trajectory data.
#' @param times Time points to predict.
#' @param ylevels States in the data.
#' @param absorb Absorbing state name.
#' @param varnames List of variable names: tvarname, pvarname, id, tx, gap (optional).
#' @param t_covs Time-dependent covariate lookup table.
#'
#' @return A list with two components:
#'   \item{sop_tx}{Matrix (Time x State) or Array (Draws x Time x State) for Treatment}
#'   \item{sop_ctrl}{Matrix (Time x State) or Array (Draws x Time x State) for Control}
#'
#' @export
standardize_sops <- function(
  model,
  data = NULL,
  times = NULL,
  ylevels = factor(1:6),
  absorb = 6,
  varnames = list(
    tvarname = "time",
    pvarname = "yprev",
    id = "id",
    tx = "tx",
    gap = NULL
  ),
  t_covs = NULL
) {
  # --- 1. Setup & Data Extraction ---
  if (
    !inherits(model, "orm") &&
      !inherits(model, "vglm") &&
      !inherits(model, "rmsb") &&
      !inherits(model, "robcov_vglm")
  ) {
    stop("model must be an orm, rmsb, vglm, or robcov_vglm object.")
  }

  if (is.null(data)) {
    data <- model$x
    if (is.null(data)) {
      stop("No data provided and model was not fitted with x=TRUE")
    }
    data$y <- model$y
  }

  # Variable mapping
  id_var <- varnames$id
  tx_var <- varnames$tx
  tvar <- varnames$tvarname
  pvar <- varnames$pvarname
  gap_var <- varnames$gap %||% NULL # Helper if varnames$gap is missing

  if (is.null(times)) {
    times <- 1:max(data[[tvar]], na.rm = TRUE)
  }

  # Check factors
  # We need to make this more sophisticated, e.g., checking levels match model
  if (!is.factor(data[[pvar]])) {
    # Enforce factor conversion for robustness
    if (getOption("markov.misc.verbose", default = FALSE)) {
      warning(paste(pvar, "should be a factor. Converting automatically."))
    }
    data[[pvar]] <- factor(data[[pvar]])
  }

  # --- 2. Create Counterfactual Cohorts ---
  # Extract one row per patient (baseline)
  X_base <- data[!duplicated(data[[id_var]]), ]

  # Create Treated Cohort (Everyone tx=1)
  X_tx <- X_base
  X_tx[[tx_var]] <- 1

  # Create Control Cohort (Everyone tx=0)
  X_ctrl <- X_base
  X_ctrl[[tx_var]] <- 0

  # --- 3. Vectorized Prediction ---
  # These calls replace the loop.
  # Returns: [Patients x Time x States] (Freq) OR [Draws x Patients x Time x States] (Bayes)

  res_tx <- soprob_markov(
    object = model,
    data = X_tx,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvar,
    pvarname = pvar,
    gap = gap_var,
    t_covs = t_covs
  )

  res_ctrl <- soprob_markov(
    object = model,
    data = X_ctrl,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvar,
    pvarname = pvar,
    gap = gap_var,
    t_covs = t_covs
  )

  # --- 4. Marginalize (Average over Patients) ---

  # Helper to average over the patient dimension
  # Freq dims: [Pat, Time, State] -> Target: [Time, State] (Ave over dim 1)
  # Bayes dims: [Draw, Pat, Time, State] -> Target: [Draw, Time, State] (Ave over dim 2)

  calc_marginal <- function(arr) {
    ndims <- length(dim(arr))
    if (ndims == 3) {
      # Frequentist: Average over rows (Patients)
      # Result: [Time x States]
      return(apply(arr, c(2, 3), mean, na.rm = TRUE))
    } else if (ndims == 4) {
      # Bayesian: Average over dim 2 (Patients)
      # Result: [Draws x Time x States]
      return(apply(arr, c(1, 3, 4), mean, na.rm = TRUE))
    } else {
      stop("Unexpected array dimensions returned from vectorized function.")
    }
  }

  sop_tx_marg <- calc_marginal(res_tx)
  sop_ctrl_marg <- calc_marginal(res_ctrl)

  # Ensure column names
  if (length(dim(sop_tx_marg)) == 2) {
    colnames(sop_tx_marg) <- as.character(ylevels)
    colnames(sop_ctrl_marg) <- as.character(ylevels)
  } else {
    dimnames(sop_tx_marg)[[3]] <- as.character(ylevels)
    dimnames(sop_ctrl_marg)[[3]] <- as.character(ylevels)
  }

  return(list(
    sop_tx = sop_tx_marg,
    sop_ctrl = sop_ctrl_marg
  ))
}

#' Compute Total Time in Target State(s)
#'
#' Calculates the expected total time spent in specified target state(s).
#' This function integrates state occupancy probabilities over time assuming unit steps.
#' It automatically adapts to the input format:
#' \itemize{
#'   \item **Patient-Level:** If input is an array from \code{soprob_markov}, it returns time-in-state for each patient.
#'   \item **Bootstrap-Level:** If input is a data frame from \code{\link{bootstrap_standardized_sops}}, it returns mean time-in-state per treatment group and the difference for each bootstrap sample.
#' }
#'
#' @param sops Input object.
#'   \itemize{
#'     \item **Array:** `[Patients x Time x States]` (Frequentist) or `[Draws x Patients x Time x States]` (Bayes).
#'     \item **Data Frame:** Output from `bootstrap_standardized_sops` containing columns `boot_id`, `time`, `tx`, and `state_*`.
#'   }
#' @param target_states Vector of target state(s) to include in the time calculation.
#'   Can be integer indices or character names (e.g., \code{1} or \code{c("Home", "Rehab")}).
#'   For bootstrap outputs, these must match the suffix of the `state_` columns (e.g., if column is `state_1`, use `1`).
#'
#' @return
#' \itemize{
#'   \item **Input = Array (Frequentist):** A named numeric vector of length \code{n_patients}.
#'   \item **Input = Array (Bayesian):** A matrix of dimension \code{[n_draws x n_patients]}.
#'   \item **Input = Bootstrap DF:** A tibble with columns:
#'     \itemize{
#'       \item `boot_id`: Bootstrap iteration
#'       \item `SOP_tx`: Mean time in target state (Treatment)
#'       \item `SOP_ctrl`: Mean time in target state (Control)
#'       \item `delta`: Difference (Treatment - Control)
#'     }
#' }
#'
#' @examples
#' \dontrun{
#' # --- Scenario 1: Patient-Level Estimates ---
#' sops_arr <- soprob_markov(model, data, times = 1:30, ylevels = 1:6)
#' days_home <- time_in_state(sops_arr, target_states = 1)
#'
#' # --- Scenario 2: Bootstrap Inference ---
#' bs_res <- bootstrap_standardized_sops(model, data, n_boot=100)
#' # Get distribution of treatment effects
#' bs_effects <- time_in_state(bs_res, target_states = 1)
#' }
#'
#' @keywords time-in-state auc bootstrap
#' @export
time_in_state <- function(sops, target_states = 1) {
  # =========================================================================
  # BRANCH 1: Bootstrap Data Frame Input
  # =========================================================================
  if (inherits(sops, "data.frame")) {
    # 1. Validation
    if (!all(c("boot_id", "time", "tx") %in% colnames(sops))) {
      stop("Input data frame must contain 'boot_id', 'time', and 'tx' columns.")
    }

    # Identify state columns in the dataframe
    # The bootstrap function outputs columns like "state_1", "state_2", etc.
    state_cols_available <- grep("^state_", colnames(sops), value = TRUE)

    if (length(state_cols_available) == 0) {
      stop(
        "Input data frame does not contain any columns starting with 'state_'."
      )
    }

    # 2. Map target_states to column names
    target_suffixes <- as.character(target_states)
    target_cols <- paste0("state_", target_suffixes)

    # Check for missing columns
    missing_cols <- setdiff(target_cols, colnames(sops))
    if (length(missing_cols) > 0) {
      stop(paste(
        "The following target state columns were not found in the bootstrap output:",
        paste(missing_cols, collapse = ", ")
      ))
    }

    # 3. Sum probabilities across target states (Row-wise)
    # If multiple states are targeted (e.g. Home + Rehab), we sum their probabilities first
    if (length(target_cols) == 1) {
      prob_target <- sops[[target_cols]]
    } else {
      prob_target <- rowSums(sops[, target_cols, drop = FALSE])
    }

    # Add temporary column for aggregation
    # We use base R aggregation for dependency minimization, or dplyr if available
    # Using dplyr approach as it's cleaner and likely present with the bootstrap workflow

    # We construct a minimal DF to aggregate
    agg_df <- sops[, c("boot_id", "tx", "time")]
    agg_df$prob <- prob_target

    # 4. Sum over time (Area Under Curve)
    # Group by boot_id and tx
    # Result: One row per boot_id per tx group
    res_grouped <- stats::aggregate(
      prob ~ boot_id + tx,
      data = agg_df,
      FUN = sum
    )

    # 5. Pivot to Wide Format (Tx vs Ctrl)
    # Split into two dataframes
    res_tx <- res_grouped[res_grouped$tx == 1, c("boot_id", "prob")]
    res_ctrl <- res_grouped[res_grouped$tx == 0, c("boot_id", "prob")]

    # Merge
    res_wide <- merge(
      res_tx,
      res_ctrl,
      by = "boot_id",
      suffixes = c("_tx", "_ctrl")
    )

    # Calculate difference
    colnames(res_wide)[colnames(res_wide) == "prob_tx"] <- "SOP_tx"
    colnames(res_wide)[colnames(res_wide) == "prob_ctrl"] <- "SOP_ctrl"
    res_wide$delta <- res_wide$SOP_tx - res_wide$SOP_ctrl

    # Return tibble if tibble package is available, otherwise DF
    if (requireNamespace("tibble", quietly = TRUE)) {
      return(tibble::as_tibble(res_wide))
    } else {
      return(res_wide)
    }
  }

  # =========================================================================
  # BRANCH 2: Array Input (Original Logic)
  # =========================================================================

  # --- 1. Detect Dimensions ---
  dims <- dim(sops)
  if (is.null(dims)) {
    stop("Input 'sops' must be an array or data frame.")
  }

  ndim <- length(dims)
  dnames <- dimnames(sops)

  if (ndim == 3) {
    # Frequentist: [Patients, Time, States]
    state_dim <- 3
  } else if (ndim == 4) {
    # Bayesian: [Draws, Patients, Time, States]
    state_dim <- 4
  } else {
    stop("Input array must be 3D or 4D.")
  }

  # --- 2. Resolve Target States ---
  state_names <- dnames[[state_dim]]
  if (is.null(state_names)) {
    state_names <- as.character(1:dims[state_dim])
  }

  target_chars <- as.character(target_states)
  missing_states <- setdiff(target_chars, state_names)
  if (length(missing_states) > 0) {
    stop(paste(
      "Target states not found in input array:",
      paste(missing_states, collapse = ", ")
    ))
  }

  # --- 3. Filter & Aggregate States ---
  subset_args <- rep(list(TRUE), ndim)
  subset_args[[state_dim]] <- target_chars

  sops_slice <- do.call("[", c(list(sops), subset_args, drop = FALSE))
  prob_in_target <- apply(sops_slice, (1:ndim)[-state_dim], sum)

  # --- 4. Integrate Over Time ---
  if (ndim == 3) {
    # Freq Input -> [Pat, Time] -> Apply over rows
    res <- rowSums(prob_in_target)
  } else {
    # Bayes Input -> [Draws, Pat, Time] -> Apply over Draws & Pat
    res <- apply(prob_in_target, c(1, 2), sum)
  }

  return(res)
}
#' Calculate Individual State Occupation Probabilities
#'
#' Computes individual-level state occupation probabilities (SOPs) for each row
#' in a dataset.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`).
#' @param newdata Optional. A data frame of new data for prediction. If NULL,
#'   uses the data used to fit the model.
#' @param times A numeric vector of time points to estimate.
#' @param ylevels A vector of state levels. If NULL, attempts to infer from model.
#' @param absorb The absorbing state.
#' @param tvarname Name of the time variable in the model.
#' @param pvarname Name of the previous state variable in the model.
#' @param gap Name of the time gap variable (if used).
#' @param t_covs Optional time-varying covariate lookup table, e.g., spline
#' basis functions.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class `markov_sops` containing:
#'   \item{rowid}{Row identifier from newdata}
#'   \item{time}{Time point}
#'   \item{state}{State name}
#'   \item{estimate}{Probability of being in the state}
#'   Plus all columns from `newdata`.
#'
#' @details
#' This function wraps `soprob_markov()` and converts its array output to a tidy
#' data frame. The output contains one row per patient-time-state combination.
#'
#' For computing marginal/standardized SOPs (G-computation), use `avg_sops()`
#' instead, which creates counterfactual datasets and averages over individuals.
#'
#' @seealso [avg_sops()] for marginal SOPs, [soprob_markov()] for the underlying
#'   computation.
#'
#' @importFrom stats model.frame
#' @export
sops <- function(
  model,
  newdata = NULL,
  times = NULL,
  ylevels = NULL,
  absorb = NULL,
  tvarname = "time",
  pvarname = "yprev",
  gap = NULL,
  t_covs = NULL,
  ...
) {
  # --- 1. Setup & Defaults ---
  if (is.null(newdata)) {
    # newdata <- model$x
    # if (is.null(newdata)) {
    #   newdata <- tryCatch(model.frame(model), error = function(e) NULL)
    # }
    # if (is.null(newdata)) {
    stop("Please provide `newdata` (should be baseline data).")
    # }
  }

  # Ensure rowid exists for tracking
  if (!"rowid" %in% names(newdata)) {
    newdata$rowid <- seq_len(nrow(newdata))
  }

  if (is.null(times)) {
    if (is.null(tvarname) || !tvarname %in% names(newdata)) {
      stop("`times` must be specified if `tvarname` is not in data.")
    }
    times <- sort(unique(newdata[[tvarname]]))
  }

  if (is.null(ylevels)) {
    # Try to infer from model
    if (inherits(model, "vglm")) {
      # VGAM stores response levels
      ylevels <- model@extra$colnames.y
      if (is.null(ylevels)) stop("`ylevels` cannot be NULL")
    } else if (inherits(model, "robcov_vglm")) {
      # robcov_vglm stores extra slot in list
      ylevels <- model$extra$colnames.y
      if (is.null(ylevels)) stop("`ylevels` cannot be NULL")
    } else if (inherits(model, "orm")) {
      # rms orm stores levels
      ylevels <- model$yunique
      if (is.null(ylevels)) stop("`ylevels` cannot be NULL")
    } else {
      stop("`ylevels` cannot be NULL")
    }
  }

  # --- 2. Compute SOPs (Vectorized) ---
  sops_array <- soprob_markov(
    object = model,
    data = newdata,
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    gap = gap,
    t_covs = t_covs
  )

  # --- 3. Tidy the Output ---
  # Array dims: [n_pat, n_times, n_states]
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]

  # Flatten array (Reminder for myself: R is column-major, iterates dim1 fastest)
  # as.vector order: [1,1,1], [2,1,1]... [N,1,1], [1,2,1]... [N,T,S]
  probs_flat <- as.vector(sops_array)

  # Construct indices matching as.vector order
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)

  # Build result by repeating newdata rows

  result <- newdata[idx_pat, , drop = FALSE]
  result$time <- idx_time
  result$state <- idx_state
  result$estimate <- probs_flat
  rownames(result) <- NULL

  # Store attributes for downstream use
  attr(result, "model") <- model
  attr(result, "tvarname") <- tvarname
  attr(result, "pvarname") <- pvarname
  attr(result, "ylevels") <- ylevels
  attr(result, "absorb") <- absorb
  attr(result, "gap") <- gap
  attr(result, "t_covs") <- t_covs
  attr(result, "call_args") <- list(
    times = times,
    ylevels = ylevels,
    absorb = absorb,
    tvarname = tvarname,
    pvarname = pvarname,
    gap = gap,
    t_covs = t_covs
  )
  # Store original newdata for inferences() (needed for simulation-based inference)
  attr(result, "newdata_orig") <- newdata

  class(result) <- c("markov_sops", class(result))
  return(result)
}


#' Calculate Averaged State Occupation Probabilities (Marginal Effects)
#'
#' Computes standardized (marginal) state occupation probabilities using
#' G-computation. Creates counterfactual cohorts by setting all individuals
#' to each level of the treatment variable and averaging over the covariate
#' distribution.
#'
#' @param model A fitted model object (e.g., `vglm`, `orm`).
#' @param newdata Data frame for prediction. For **simulation inference**, pass
#'   baseline data only (one row per patient). For **bootstrap inference**, you
#'   must pass the full longitudinal dataset (all time points) since the model
#'   needs to be refit on bootstrap samples. If NULL, extracts from model.
#' @param variables A named list specifying the variable(s) to standardize over.
#'   E.g., `list(tx = c(0, 1))` creates counterfactual datasets for treatment
#'   and control.
#' @param by Optional character vector of additional variables to group by,
#' after standardization.
#' @param times Numeric vector of time points. If NULL, inferred from data.
#' @param id_var Name of the patient ID variable (default "id"). Required for
#'   bootstrap inference.
#' @param ... Additional arguments passed to `sops()` (e.g., `ylevels`, `absorb`,
#'   `tvarname`, `pvarname`, `t_covs`).
#'
#' @return A data frame of class `markov_avg_sops` with columns:
#'   \item{time}{Time point}
#'   \item{state}{State level}
#'   \item{(variables)}{Value of standardization variable (e.g., tx)}
#'   \item{estimate}{Average probability across individuals}
#'
#' @details
#' This function implements G-computation (standardization) for Markov SOPs:
#'
#' 1. **Counterfactual Creation**: For each value in `variables`, creates a
#'    copy of `newdata` with that variable set to the specified value.
#'
#' 2. **Individual Prediction**: Calls `sops()` to compute individual-level
#'    SOPs for each patient under each counterfactual scenario.
#'
#' 3. **Marginalization**: Averages individual SOPs across patients within
#'    each time-state-treatment combination, yielding population-average
#'    (marginal) probabilities.
#'
#' The result represents the expected SOP if the entire population received
#' treatment vs. control, averaged over the observed covariate distribution.
#'
#' @seealso [sops()] for individual-level SOPs, [inferences()] for bootstrap
#'   uncertainty, [standardize_sops()] for the underlying implementation.
#'
#' @examples
#' \dontrun{
#' # For simulation inference: use baseline data (one row per patient)
#' baseline_data <- data |> filter(time == 1)
#' result <- avg_sops(
#'   model = fit,
#'   newdata = baseline_data,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6
#' ) |> inferences(method = "simulation", n_sim = 500)
#'
#' # For bootstrap inference: must use full data (all time points)
#' result_boot <- avg_sops(
#'   model = fit,
#'   newdata = data,  # Full longitudinal data
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = 1:6,
#'   absorb = 6
#' ) |> inferences(method = "bootstrap", n_sim = 500)
#' }
#'
#' @export
avg_sops <- function(
  model,
  newdata = NULL,
  variables = NULL,
  by = NULL,
  times = NULL,
  id_var = "id",
  ...
) {
  # --- 1. Input Validation ---
  if (is.null(variables)) {
    stop(
      "`variables` is required for G-computation. ",
      "Specify the treatment variable, e.g., `variables = list(tx = c(0, 1))`."
    )
  }

  if (is.null(newdata)) {
    # this creates too many issues
    # newdata <- model$x
    # if (is.null(newdata)) {
    #   newdata <- tryCatch(model.frame(model), error = function(e) NULL)
    # }
    # if (is.null(newdata)) {
    stop("Provide newdata or ensure model stores data (x = TRUE).")
    # }
  }

  # Validate id_var exists
  if (!id_var %in% names(newdata)) {
    stop("ID variable '", id_var, "' not found in data.")
  }

  # Validate variables exist in data
  missing_vars <- setdiff(names(variables), names(newdata))
  if (length(missing_vars) > 0) {
    stop("Variables not in data: ", paste(missing_vars, collapse = ", "))
  }

  # --- 2. Extract Baseline Data (One Row Per Patient) ---
  # For standardization, we need unique patient baseline covariates
  # This matches standardize_sops() behavior
  baseline_data <- newdata[!duplicated(newdata[[id_var]]), ]

  # Could consider not checking ID, but taking tvarname as an argument and then
  # filtering for min.
  # newdata[newdata[, tvarname] == 1, ]

  # Ensure rowid exists
  if (!"rowid" %in% names(baseline_data)) {
    baseline_data$rowid <- seq_len(nrow(baseline_data))
  }

  # --- 3. Create Counterfactual Datasets ---
  # For each combination in variables, create a copy of baseline_data with
  # the variable(s) set to that value

  if (!is.list(variables)) {
    var_list <- list()
    for (i in seq_along(variables)) {
      var_list[[variables[i]]] <- unique(baseline_data[[variables[i]]])
    }
  } else {
    var_list <- variables
  }

  grid <- do.call(expand.grid, var_list)

  # expanded_data_list <- vector("list", nrow(grid))
  # for (i in seq_len(nrow(grid))) {
  #   dt_copy <- baseline_data
  #   for (v in names(grid)) {s
  #     dt_copy[[v]] <- grid[i, v]
  #   }
  #   expanded_data_list[[i]] <- dt_copy
  # }

  # newdata_expanded <- do.call(rbind, expanded_data_list)

  newdata_expanded <- create_counterfactual_data(baseline_data, grid, var_list)

  # --- 4. Compute Individual SOPs ---
  sops_ind <- sops(model, newdata = newdata_expanded, times = times, ...)

  # --- 5. Aggregate (Marginalize) ---
  # Group by time, state, and the variables used for standardization
  group_cols <- c("time", "state", names(var_list))

  # Add optional 'by' variables
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }

  # Validate grouping columns exist
  missing_groups <- setdiff(group_cols, names(sops_ind))
  if (length(missing_groups) > 0) {
    stop("Grouping variables missing: ", paste(missing_groups, collapse = ", "))
  }

  # Aggregate
  agg_formula <- stats::as.formula(
    paste("estimate ~", paste(group_cols, collapse = " + "))
  )
  result <- stats::aggregate(
    agg_formula,
    data = sops_ind,
    FUN = mean,
    na.rm = TRUE
  )

  # --- 5. Store Attributes ---
  # Copy attributes from sops_ind
  attr(result, "model") <- attr(sops_ind, "model")
  attr(result, "call_args") <- attr(sops_ind, "call_args")
  attr(result, "tvarname") <- attr(sops_ind, "tvarname")
  attr(result, "pvarname") <- attr(sops_ind, "pvarname")
  attr(result, "ylevels") <- attr(sops_ind, "ylevels")
  attr(result, "absorb") <- attr(sops_ind, "absorb")
  attr(result, "gap") <- attr(sops_ind, "gap")
  attr(result, "t_covs") <- attr(sops_ind, "t_covs")

  # Specific attributes for avg_sops/inferences
  attr(result, "avg_args") <- list(
    variables = var_list,
    by = by,
    times = times,
    id_var = id_var
  )
  # Store ORIGINAL newdata for bootstrap (not the expanded counterfactual)
  attr(result, "newdata_orig") <- newdata

  class(result) <- c("markov_avg_sops", class(result))
  return(result)
}


#' Inference for State Occupation Probabilities
#'
#' Adds confidence intervals to SOP objects using simulation-based (MVN) or
#' bootstrap methods. The default method is simulation, which draws coefficient
#' vectors from a multivariate normal distribution centered at the point
#' estimates.
#'
#' @param object A `markov_avg_sops` object from `avg_sops()` or a
#'   `markov_sops` object from `sops()`.
#' @param method Character. Inference method:
#'   \itemize{
#'     \item `"simulation"` (default): Draws from MVN(beta_hat, Sigma). Fast,
#'       does not refit models. Works for both individual and averaged SOPs.
#'     \item `"bootstrap"`: Resamples patients with replacement, refits model.
#'  Only works for `markov_avg_sops` objects.
#'   }
#' @param n_sim Number of simulation draws (for simulation) or bootstrap
#'   iterations (for bootstrap). Default is 1000.
#' @param vcov Optional custom variance-covariance matrix. If provided,
#'   overrides the vcov extracted from the model.
#' @param parallel Logical. Use parallel processing? Default is FALSE.
#' @param workers Number of parallel workers. If NULL, uses detectCores() - 1.
#' @param conf_level Confidence level for intervals (default 0.95).
#' @param conf_type Type of confidence interval (simulation method only):
#'   \itemize{
#'     \item `"perc"` (default): Percentile-based intervals from the simulation
#'       distribution.
#'     \item `"wald"`: Uses simulation standard errors with normal quantiles.
#'   }
#' @param return_draws Logical. If TRUE, stores all individual simulation/bootstrap
#'   draws as an attribute. Extract with `get_draws()`. Default is TRUE
#' @param update_datadist Logical. Whether to update datadist for rms models
#'   during bootstrap. Default is TRUE.
#' @param use_coefstart Logical. Use original coefficients as starting values
#'   for bootstrap refitting. Default is FALSE.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object with added columns:
#'   \item{conf.low}{Lower confidence bound}
#'   \item{conf.high}{Upper confidence bound}
#'   \item{std.error}{Standard error from simulation/bootstrap}
#'
#'   If `return_draws = TRUE`, the object also has a `"simulation_draws"` or
#'   `"bootstrap_draws"` attribute containing all individual draws. Extract
#'   with `get_draws()`.
#'
#' @details
#' ## Simulation Method
#'
#' The simulation method works as follows:
#' 1. Extract coefficient vector beta_hat and (robust) variance-covariance Sigma
#' 2. Draw n_sim coefficient vectors from MVN(beta_hat, Sigma)
#' 3. For each draw, replace model coefficients and compute SOPs
#' 4. Compute confidence intervals from the empirical distribution
#'
#' - Works for both individual-level (`sops()`) and averaged (`avg_sops()`) SOPs
#'
#' ## Bootstrap Method
#'
#' The bootstrap method resamples patients and refits the model:
#' 1. Sample patient IDs with replacement
#' 2. Refit model on bootstrap sample (handles missing states through releveling)
#' 3. Compute SOPs using G-computation
#' 4. Compute percentile-based confidence intervals
#' Bootstrap requires the full longitudinal dataset (all time
#' points) in the original `newdata` passed to `avg_sops()`, not just baseline
#' data. This is because the model must be refit on each bootstrap sample.
#'
#'
#' This design ensures consistency: the same vcov is used for both point
#' estimates and inference, regardless of how `inferences()` is called.
#'
#' @seealso [avg_sops()], [sops()], [get_draws()], [robcov_vglm()],
#'   [set_coef()]
#'
#' @examples
#' \dontrun{
#' # Step 1: Fit model and wrap with cluster-robust vcov
#' fit <- vglm(y ~ time + tx + yprev, family = cumulative(reverse = TRUE,
#' parallel = TRUE), data = data)
#' fit_robust <- robcov_vglm(fit, cluster = data$id)
#'
#' # Step 2a: Simulation inference (use baseline data only)
#' baseline_data <- data |> filter(time == 1)
#' result_sim <- avg_sops(
#'   model = fit_robust,
#'   newdata = baseline_data,  # Baseline only for simulation
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6",
#'   id_var = "id"
#' ) |>
#'   inferences(method = "simulation", n_sim = 1000)
#'
#' # Step 2b: Bootstrap inference (must use full data)
#' result_boot <- avg_sops(
#'   model = fit_robust,
#'   newdata = data,  # Full longitudinal data for bootstrap
#'   variables = list(tx = c(0, 1)),
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6",
#'   id_var = "id"
#' ) |>
#'   inferences(method = "bootstrap", n_sim = 500)
#'
#' # Individual-level SOPs with simulation inference
#' ind_result <- sops(
#'   model = fit_robust,
#'   newdata = baseline_data,
#'   times = 1:60,
#'   ylevels = factor(1:6),
#'   absorb = "6") |>
#'   inferences(n_sim = 500)
#'
#' # Extract draws for custom analyses
#' get_draws(ind_result)
#'
#' # Compute treatment effect on time in state
#' library(dplyr)
#' library(tidyr)
#' treatment_effect <- draws |>
#'   filter(state == 1) |>
#'   group_by(draw_id, tx) |>
#'   summarise(total_time = sum(draw), .groups = "drop") |>
#'   pivot_wider(names_from = tx, values_from = total_time) |>
#'   mutate(effect = `1` - `0`)
#' quantile(treatment_effect$effect, c(0.025, 0.5, 0.975))
#' }
#'
#' @seealso [get_draws()] to extract individual draws from any inference method
#'
#' @export
inferences <- function(
  object,
  method = "simulation",
  n_sim = 1000,
  vcov = NULL,
  parallel = FALSE,
  workers = NULL,
  conf_level = 0.95,
  conf_type = "perc",
  return_draws = TRUE,
  update_datadist = TRUE,
  use_coefstart = FALSE,
  ...
) {
  # --- Input Validation ---
  if (!inherits(object, c("markov_avg_sops", "markov_sops"))) {
    stop(
      "inferences() requires a 'markov_avg_sops' or 'markov_sops' object. ",
      "Got: ",
      paste(class(object), collapse = ", ")
    )
  }

  method <- match.arg(method, choices = c("simulation", "bootstrap"))

  # # Handle backward compatibility: n_boot overrides n_sim for bootstrap
  # if (!is.null(n_boot) && method == "bootstrap") {
  #   n_sim <- n_boot
  # }

  # --- Dispatch to Method-Specific Implementation ---
  if (method == "simulation") {
    inferences_simulation(
      object = object,
      n_sim = n_sim,
      vcov = vcov,
      conf_level = conf_level,
      conf_type = conf_type,
      parallel = parallel,
      workers = workers,
      return_draws = return_draws
    )
  } else if (method == "bootstrap") {
    inferences_bootstrap(
      object = object,
      n_boot = n_sim,
      parallel = parallel,
      workers = workers,
      conf_level = conf_level,
      return_draws = return_draws,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )
  }
}


# =============================================================================
# SIMULATION-BASED INFERENCE (MVN)
# =============================================================================

#' Simulation-Based Inference Using Multivariate Normal
#'
#' Internal function that implements MVN simulation-based inference for SOPs.
#' Draws coefficient vectors from a multivariate normal distribution centered
#' at the point estimates with the (robust) variance-covariance matrix.
#'
#' @param object A `markov_avg_sops` or `markov_sops` object.
#' @param n_sim Number of simulation draws.
#' @param vcov Custom variance-covariance matrix (optional, overrides model vcov).
#' @param conf_level Confidence level.
#' @param conf_type Type of confidence interval ("perc" or "wald").
#' @param parallel Use parallel processing?
#' @param workers Number of parallel workers.
#' @param return_draws Store individual draws?
#'
#' @return The input object with added confidence intervals.
#'
#' @keywords internal
inferences_simulation <- function(
  object,
  n_sim,
  vcov,
  conf_level,
  conf_type,
  parallel,
  workers,
  return_draws
) {
  # Check for mvtnorm package
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop(
      "Package 'mvtnorm' is required for simulation-based inference.\n",
      "Install with: install.packages('mvtnorm')"
    )
  }

  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  call_args <- attr(object, "call_args")
  avg_args <- attr(object, "avg_args")

  tvarname <- attr(object, "tvarname")
  pvarname <- attr(object, "pvarname")
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
  t_covs <- attr(object, "t_covs")

  # For avg_sops objects
  is_avg <- inherits(object, "markov_avg_sops")
  if (is_avg) {
    variables <- avg_args$variables
    by <- avg_args$by
    times <- avg_args$times
    id_var <- avg_args$id_var
  } else {
    # For individual sops objects
    times <- call_args$times
    variables <- NULL
    by <- NULL
    id_var <- NULL
  }

  if (is.null(model)) {
    stop("Model not stored in object. Cannot perform simulation inference.")
  }

  # --- 2. Get Coefficients and VCov from Model ---
  # The vcov is extracted directly from the model object:
  # - For robcov_vglm: uses the stored cluster-robust vcov
  # - For orm with rms::robcov(): uses the stored robust vcov
  # - For plain vglm/orm: uses model-based vcov
  # Users who want cluster-robust SEs should wrap their model with
  # robcov_vglm() or rms::robcov() BEFORE passing to sops()/avg_sops().
  beta_hat <- get_coef(model)

  if (!is.null(vcov) && is.matrix(vcov)) {
    # User-provided custom vcov matrix (for sensitivity analyses)
    Sigma <- vcov
  } else {
    # Extract vcov from model (robust if model was wrapped with robcov_vglm)
    Sigma <- get_vcov_robust(model)
  }

  # Validate dimensions
  if (length(beta_hat) != nrow(Sigma) || length(beta_hat) != ncol(Sigma)) {
    stop(
      "Dimension mismatch: coefficients (",
      length(beta_hat),
      ") vs ",
      "vcov matrix (",
      nrow(Sigma),
      " x ",
      ncol(Sigma),
      ")"
    )
  }

  # --- 3. Draw n_sim Coefficient Vectors from MVN ---
  beta_draws <- mvtnorm::rmvnorm(n_sim, mean = beta_hat, sigma = Sigma)

  # --- 4. Prepare Prediction Data (COMPUTED ONCE!) ---
  if (is_avg) {
    # For avg_sops: create counterfactual datasets
    baseline_data <- newdata_orig[!duplicated(newdata_orig[[id_var]]), ]
    grid <- do.call(expand.grid, variables)
    newdata_pred <- create_counterfactual_data(baseline_data, grid, variables)
    n_cf <- nrow(grid)
    n_each <- nrow(baseline_data)
  } else {
    # For individual sops: use data directly
    newdata_pred <- newdata_orig
    grid <- NULL
    n_cf <- 1
    n_each <- nrow(newdata_pred)
  }

  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 5. Detect Fast Path Eligibility ---
  # The fast path uses pre-computed design matrix decompositions for VGLM models
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  use_fast_path <- inherits(model_chk, "vglm")

  if (use_fast_path) {
    # --- FAST PATH: Pre-build components once, then run efficient Markov loop ---
    components <- tryCatch(
      markov_msm_build(
        model = model_chk,
        data = newdata_pred,
        t_covs = t_covs,
        ylevels = ylevels,
        pvarname = pvarname
      ),
      error = function(e) {
        warning(
          "Fast path build failed, falling back to slow path: ",
          e$message
        )
        NULL
      }
    )

    if (!is.null(components)) {
      # Get constraint list for computing Gamma from beta draws
      C_list <- VGAM::constraints(model_chk)

      # Define fast analysis function
      analysis_fn <- function(i) {
        # Compute effective coefficients from beta draw
        Gamma_i <- compute_Gamma(beta_draws[i, ], C_list)

        # Run fast Markov simulation
        sops_array <- tryCatch(
          markov_msm_run(components, Gamma_i, times, absorb),
          error = function(e) {
            warning("markov_msm_run failed in draw ", i, ": ", e$message)
            return(NULL)
          }
        )

        if (is.null(sops_array)) {
          return(NULL)
        }

        # Marginalize if needed (for avg_sops)
        if (is_avg) {
          result <- marginalize_sops_array(
            sops_array = sops_array,
            grid = grid,
            times = times,
            ylevels = ylevels,
            variables = variables,
            n_cf = n_cf,
            n_each = n_each
          )
        } else {
          # Individual-level: convert array to data frame
          result <- array_to_df_individual(
            sops_array,
            times,
            ylevels,
            newdata_pred
          )
        }

        result
      }
    } else {
      # Fast path build failed, fall back to slow path
      use_fast_path <- FALSE
    }
  }

  if (!use_fast_path) {
    # --- SLOW PATH: Use soprob_markov with coefficient replacement ---
    analysis_fn <- function(i) {
      # Replace coefficients
      model_i <- set_coef(model, beta_draws[i, ])

      # Compute SOPs
      sops_array <- tryCatch(
        soprob_markov(
          object = model_i,
          data = newdata_pred,
          times = times,
          ylevels = ylevels,
          absorb = absorb,
          tvarname = tvarname,
          pvarname = pvarname,
          t_covs = t_covs
        ),
        error = function(e) {
          warning("soprob_markov failed in draw ", i, ": ", e$message)
          return(NULL)
        }
      )

      if (is.null(sops_array)) {
        return(NULL)
      }

      # Marginalize if needed (for avg_sops)
      if (is_avg) {
        result <- marginalize_sops_array(
          sops_array = sops_array,
          grid = grid,
          times = times,
          ylevels = ylevels,
          variables = variables,
          n_cf = n_cf,
          n_each = n_each
        )
      } else {
        # Individual-level: convert array to data frame
        result <- array_to_df_individual(
          sops_array,
          times,
          ylevels,
          newdata_pred
        )
      }

      result
    }
  }

  # --- 6. Apply Across All Draws ---
  if (parallel) {
    # Setup parallel
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("Package 'furrr' is required for parallel processing")
    }

    if (is.null(workers)) {
      workers <- max(1, parallel::detectCores() - 1)
    }

    future::plan(future.callr::callr, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)

    # Build globals list dynamically based on fast/slow path
    globals_list <- c(
      "model",
      "beta_draws",
      "newdata_pred",
      "times",
      "ylevels",
      "absorb",
      "tvarname",
      "pvarname",
      "t_covs",
      "is_avg",
      "grid",
      "variables",
      "n_cf",
      "n_each",
      "use_fast_path"
    )
    if (use_fast_path) {
      globals_list <- c(globals_list, "components", "C_list")
    }

    sim_results <- furrr::future_map(
      seq_len(n_sim),
      analysis_fn,
      .options = furrr::furrr_options(
        seed = TRUE,
        globals = globals_list,
        packages = c("rms", "VGAM", "dplyr", "stats", "markov.misc")
      ),
      .progress = FALSE
    )
  } else {
    sim_results <- lapply(seq_len(n_sim), analysis_fn)
  }

  # --- 7. Combine Results ---
  sim_results <- Filter(Negate(is.null), sim_results)

  if (length(sim_results) == 0) {
    stop("All simulation draws failed.")
  }

  # Add draw_id
  for (i in seq_along(sim_results)) {
    sim_results[[i]]$draw_id <- i
  }
  draws_df <- dplyr::bind_rows(sim_results)

  # --- 8. Compute Confidence Intervals ---
  if (is_avg) {
    group_cols <- c("time", "state", names(variables))
    if (!is.null(by)) {
      group_cols <- unique(c(group_cols, by))
    }
  } else {
    group_cols <- c("rowid", "time", "state")
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = draws_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )

  # --- 9. Merge with Original Object ---
  final_result <- merge(object, summary_stats, by = group_cols, all.x = TRUE)

  # Restore attributes
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(final_result, a) <- attr(object, a)
    }
  }
  class(final_result) <- class(object)

  # Add metadata
  attr(final_result, "n_sim") <- n_sim
  attr(final_result, "n_successful") <- length(sim_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "conf_type") <- conf_type
  attr(final_result, "method") <- "simulation"

  if (return_draws) {
    attr(final_result, "simulation_draws") <- draws_df
  }

  final_result
}


# =============================================================================
# BOOTSTRAP-BASED INFERENCE
# =============================================================================

#' Bootstrap-Based Inference for SOPs
#'
#' Internal function that implements bootstrap inference for SOPs. Resamples
#' patients with replacement, refits the model, and computes SOPs for each
#' bootstrap sample.
#'
#' @param object A `markov_avg_sops` object.
#' @param n_boot Number of bootstrap iterations.
#' @param parallel Use parallel processing?
#' @param workers Number of parallel workers.
#' @param conf_level Confidence level.
#' @param return_draws Store individual bootstrap draws?
#' @param update_datadist Update datadist for rms models?
#' @param use_coefstart Use original coefficients as starting values?
#'
#' @return The input object with added confidence intervals.
#'
#' @keywords internal
inferences_bootstrap <- function(
  object,
  n_boot,
  parallel,
  workers,
  conf_level,
  return_draws,
  update_datadist,
  use_coefstart
) {
  # Bootstrap only supports avg_sops for now
  if (!inherits(object, "markov_avg_sops")) {
    stop(
      "Bootstrap inference currently only supports 'markov_avg_sops' objects. ",
      "For individual-level SOPs, use method = 'simulation'."
    )
  }

  # --- 1. Extract Stored Attributes ---
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  call_args <- attr(object, "call_args")
  avg_args <- attr(object, "avg_args")

  tvarname <- attr(object, "tvarname")
  pvarname <- attr(object, "pvarname")
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
  t_covs <- attr(object, "t_covs")

  variables <- avg_args$variables
  by <- avg_args$by
  times <- avg_args$times
  id_var <- avg_args$id_var

  if (is.null(newdata_orig)) {
    stop("Original newdata not stored. Cannot perform bootstrap.")
  }

  if (!id_var %in% names(newdata_orig)) {
    stop("ID variable '", id_var, "' not found in data.")
  }

  # --- 2. Validate Data is Longitudinal (Not Just Baseline) ---
  # Check if tvarname exists and has multiple time points per patient
  if (tvarname %in% names(newdata_orig)) {
    rows_per_patient <- newdata_orig |>
      dplyr::group_by(!!rlang::sym(id_var)) |>
      dplyr::summarise(n_rows = dplyr::n(), .groups = "drop") |>
      dplyr::pull(n_rows)

    if (all(rows_per_patient == 1)) {
      stop(
        "Bootstrap inference requires full longitudinal data (all time points), ",
        "but the data passed to avg_sops() appears to be baseline only ",
        "(one row per patient).\\n\\n",
        "For bootstrap: avg_sops(model, newdata = data, ...)\\n",
        "For simulation: avg_sops(model, newdata = data |> filter(time == 1), ...)"
      )
    }
  }

  # Dimensions for result expansion
  n_times <- length(times)
  n_states <- length(ylevels)

  # --- 2. Identify Factor Columns ---
  factor_cols <- c("y", pvarname)
  factor_cols <- intersect(factor_cols, names(newdata_orig))

  # --- 3. Generate Bootstrap ID Samples ---
  boot_ids <- fast_group_bootstrap(
    data = newdata_orig,
    id_var = id_var,
    n_boot = n_boot
  )

  # --- 4. Define Analysis Function ---
  analysis_fn <- function(boot_data) {
    # A. Relevel factors and refit model
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = factor_cols,
      original_data = newdata_orig,
      ylevels = ylevels,
      absorb = absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )

    m_boot <- boot_result$model
    boot_data <- boot_result$data
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(NULL)
    }

    # Get states present in original numbering (for mapping back later)
    states_present <- setdiff(ylevels, as.numeric(missing_states))

    # B. Compute standardized SOPs using G-computation on bootstrap data
    # Extract baseline data (one row per patient)
    baseline_boot <- boot_data[!duplicated(boot_data[["new_id"]]), ]

    # Create counterfactual datasets for each variable value
    grid <- do.call(expand.grid, variables)
    cf_data_list <- vector("list", nrow(grid))
    for (i in seq_len(nrow(grid))) {
      dt_copy <- baseline_boot
      for (v in names(grid)) {
        dt_copy[[v]] <- grid[i, v]
      }
      cf_data_list[[i]] <- dt_copy
    }
    newdata_cf <- do.call(rbind, cf_data_list)

    # Compute individual SOPs for counterfactual data
    sops_array <- tryCatch(
      soprob_markov(
        object = m_boot,
        data = newdata_cf,
        times = times,
        ylevels = factor(boot_ylevels),
        absorb = boot_absorb,
        tvarname = tvarname,
        pvarname = pvarname,
        t_covs = t_covs
      ),
      error = function(e) {
        warning("soprob_markov failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(sops_array)) {
      return(NULL)
    }

    # C. Marginalize (average) across patients for each counterfactual
    # sops_array is [n_pat, n_times, n_boot_states]
    n_pat <- dim(sops_array)[1]
    n_boot_states <- dim(sops_array)[3]
    n_cf <- nrow(grid)
    n_each <- n_pat %/% n_cf

    # Average within each counterfactual group
    avg_sops_list <- vector("list", n_cf)
    for (cf_i in seq_len(n_cf)) {
      start_idx <- (cf_i - 1) * n_each + 1
      end_idx <- cf_i * n_each

      # Subset and average [n_each x n_times x n_boot_states]
      sops_cf <- sops_array[start_idx:end_idx, , , drop = FALSE]
      avg_sops_mat <- apply(sops_cf, c(2, 3), mean) # [n_times x n_boot_states]

      avg_sops_list[[cf_i]] <- avg_sops_mat
    }

    # D. Expand back to original state space with zero-padding
    if (length(missing_states) > 0) {
      for (cf_i in seq_len(n_cf)) {
        avg_sops_mat <- avg_sops_list[[cf_i]]
        full_mat <- matrix(0, nrow = n_times, ncol = n_states)

        # Fill in states that were present
        for (i in seq_along(states_present)) {
          original_state <- states_present[i]
          original_idx <- which(ylevels == original_state)
          full_mat[, original_idx] <- avg_sops_mat[, i]
        }

        avg_sops_list[[cf_i]] <- full_mat
      }
    }

    # E. Format as data frame
    result_list <- vector("list", n_cf)
    for (cf_i in seq_len(n_cf)) {
      avg_sops_mat <- avg_sops_list[[cf_i]]

      df <- expand.grid(time = times, state = ylevels)
      df$estimate <- as.vector(avg_sops_mat)

      # Add variable values for this counterfactual
      for (v in names(grid)) {
        df[[v]] <- grid[cf_i, v]
      }

      result_list[[cf_i]] <- df
    }

    dplyr::bind_rows(result_list)
  }

  # --- 5. Apply to Bootstrap Samples ---
  boot_results <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = newdata_orig,
    id_var = id_var,
    parallel = parallel,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "dplyr", "stats"),
    globals = c(
      "model",
      "variables",
      "times",
      "ylevels",
      "absorb",
      "pvarname",
      "tvarname",
      "t_covs",
      "n_times",
      "n_states",
      "update_datadist",
      "factor_cols",
      "use_coefstart"
    )
  )

  # --- 6. Combine and Compute Summary Statistics ---
  boot_results <- Filter(Negate(is.null), boot_results)

  if (length(boot_results) == 0) {
    stop("All bootstrap iterations failed.")
  }

  # Add draw_id and combine
  for (i in seq_along(boot_results)) {
    boot_results[[i]]$draw_id <- i
  }
  boot_df <- dplyr::bind_rows(boot_results)

  # Compute confidence intervals
  group_cols <- c("time", "state", names(variables))
  if (!is.null(by)) {
    group_cols <- unique(c(group_cols, by))
  }

  summary_stats <- compute_ci_from_draws(
    draws_df = boot_df,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = "perc" # Bootstrap always uses percentile
  )

  # Merge with original object
  final_result <- merge(object, summary_stats, by = group_cols, all.x = TRUE)

  # Restore attributes
  for (a in names(attributes(object))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(final_result, a) <- attr(object, a)
    }
  }
  class(final_result) <- class(object)

  # Add bootstrap metadata
  attr(final_result, "n_boot") <- n_boot
  attr(final_result, "n_successful") <- length(boot_results)
  attr(final_result, "conf_level") <- conf_level
  attr(final_result, "method") <- "bootstrap"

  # Store full bootstrap draws if requested
  if (return_draws) {
    attr(final_result, "bootstrap_draws") <- boot_df
  }

  final_result
}


# =============================================================================
# HELPER FUNCTIONS FOR INFERENCE
# =============================================================================

#' Create Counterfactual Datasets for G-Computation
#'
#' Creates copies of baseline data with treatment variable set to each level.
#'
#' @param baseline_data Data frame with one row per patient.
#' @param grid Data frame of variable combinations.
#' @param variables Named list of variable values.
#'
#' @return Data frame with counterfactual data stacked.
#'
#' @keywords internal
create_counterfactual_data <- function(baseline_data, grid, variables) {
  cf_data_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt_copy <- baseline_data
    for (v in names(grid)) {
      dt_copy[[v]] <- grid[i, v]
    }
    cf_data_list[[i]] <- dt_copy
  }
  do.call(rbind, cf_data_list)
}


#' Marginalize SOPs Array Over Patients
#'
#' Averages individual-level SOPs to get population-average (marginal) SOPs.
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param grid Data frame of variable combinations.
#' @param times Vector of time points.
#' @param ylevels Vector of state levels.
#' @param variables Named list of variables.
#' @param n_cf Number of counterfactual scenarios.
#' @param n_each Number of patients per scenario.
#'
#' @return Data frame with marginalized SOPs.
#'
#' @keywords internal
marginalize_sops_array <- function(
  sops_array,
  grid,
  times,
  ylevels,
  variables,
  n_cf,
  n_each
) {
  n_times <- length(times)
  n_states <- length(ylevels)

  # Average within each counterfactual group
  avg_sops_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    start_idx <- (cf_i - 1) * n_each + 1
    end_idx <- cf_i * n_each

    sops_cf <- sops_array[start_idx:end_idx, , , drop = FALSE]
    avg_sops_mat <- apply(sops_cf, c(2, 3), mean)

    avg_sops_list[[cf_i]] <- avg_sops_mat
  }

  # Format as data frame
  result_list <- vector("list", n_cf)
  for (cf_i in seq_len(n_cf)) {
    avg_sops_mat <- avg_sops_list[[cf_i]]

    df <- expand.grid(time = times, state = ylevels)
    df$estimate <- as.vector(avg_sops_mat)

    for (v in names(grid)) {
      df[[v]] <- grid[cf_i, v]
    }

    result_list[[cf_i]] <- df
  }

  dplyr::bind_rows(result_list)
}


#' Convert SOPs Array to Individual-Level Data Frame
#'
#' @param sops_array Array of dimensions (n_pat x n_times x n_states).
#' @param times Vector of time points.
#' @param ylevels Vector of state levels.
#' @param newdata Original data with rowid.
#'
#' @return Data frame with individual SOPs.
#'
#' @keywords internal
array_to_df_individual <- function(sops_array, times, ylevels, newdata) {
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]

  # Flatten array
  probs_flat <- as.vector(sops_array)

  # Construct indices
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)

  # Build result
  result <- data.frame(
    rowid = if ("rowid" %in% names(newdata)) {
      newdata$rowid[idx_pat]
    } else {
      idx_pat
    },
    time = idx_time,
    state = idx_state,
    estimate = probs_flat
  )

  result
}


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
#' @param object An object returned by `inferences()` with `return_draws = TRUE`.
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item draw_id: Simulation or bootstrap iteration number
#'     \item time: Time point
#'     \item state: State number
#'     \item draw: Draw-specific estimate of state occupation probability
#'     \item estimate: The original point estimate from the model
#'     \item conf.low, conf.high, std.error: Summary statistics from the point estimate object
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
  if (!inherits(object, c("markov_avg_sops", "markov_sops"))) {
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
  if (requireNamespace("dplyr", quietly = TRUE)) {
    draws <- dplyr::left_join(draws, meta, by = keys)
  } else {
    draws <- merge(draws, meta, by = keys, all.x = TRUE, sort = FALSE)
  }

  return(draws)
}


# =============================================================================
# FAST PATH HELPER FUNCTIONS
# =============================================================================
# These internal functions provide the fast path for simulation-based inference
# on VGLM models by pre-computing design matrix decompositions and bypassing
# VGAM's predict() function.

#' Build Fast Markov Simulation Components (Internal)
#'
#' Pre-calculates design matrices and interaction structures for fast Markov
#' simulation. This is the "heavy lifting" step that should be done ONCE before
#' running thousands of simulation draws.
#'
#' @param model Fitted vglm model
#' @param data Baseline data frame (one row per patient)
#' @param t_covs Time-dependent covariate lookup (optional)
#' @param ylevels State labels
#' @param pvarname Name of previous state variable
#' @param ... Ignored
#'
#' @return A list containing pre-calculated matrices and metadata.
#'
#' @keywords internal
markov_msm_build <- function(
  model,
  data,
  t_covs = NULL,
  ylevels = 1:6,
  pvarname = "yprev",
  ...
) {
  n_pat <- nrow(data)
  n_states <- length(ylevels)
  M <- n_states - 1

  # Need Gamma structure to know which columns to keep
  Gamma_template <- get_effective_coefs(model)
  common_cols <- colnames(Gamma_template)

  # Prepare data
  if (!pvarname %in% names(data)) {
    data[[pvarname]] <- factor(ylevels[1], levels = ylevels)
  }
  tvar <- "time"
  if (!tvar %in% names(data)) {
    data[[tvar]] <- 0
  }

  # Terms object
  tt <- stats::terms(model)
  tt <- stats::delete.response(tt)

  # Helper to get Design Matrix (X) only
  get_X <- function(d) {
    X <- stats::model.matrix(tt, data = d)
    common <- intersect(colnames(X), common_cols)
    if (length(common) == 0) {
      return(NULL)
    }
    X[, common, drop = FALSE]
  }

  # A. X_0 (Base: Time=0, Prev=Ref)
  d_0 <- data
  if (!is.null(t_covs)) {
    for (n in names(t_covs)) {
      d_0[[n]] <- 0
    }
  }

  # if we don't want to treat yprev as factor, this becomes a problem
  d_0[[pvarname]] <- factor(ylevels[1], levels = ylevels)
  X_0 <- get_X(d_0)

  # B. Time Slopes
  X_slopes <- list()
  time_vars <- if (!is.null(t_covs)) names(t_covs) else character(0)

  for (tv in time_vars) {
    d_tv <- d_0
    d_tv[[tv]] <- 1
    X_tv <- get_X(d_tv)

    if (!is.null(X_tv) && !is.null(X_0)) {
      X_diff <- X_tv - X_0
      X_slopes[[tv]] <- X_diff
    }
  }

  # C. Prev Main Effects
  X_prev <- list()
  for (k in 1:n_states) {
    d_k <- d_0
    d_k[[pvarname]] <- factor(ylevels[k], levels = ylevels)
    X_k <- get_X(d_k)

    if (!is.null(X_k) && !is.null(X_0)) {
      X_prev[[k]] <- X_k - X_0
    }
  }

  # D. Interactions (Prev x Time)
  X_interactions <- list()

  for (k in 1:n_states) {
    X_interactions[[k]] <- list()

    d_k_base <- d_0
    d_k_base[[pvarname]] <- factor(ylevels[k], levels = ylevels)
    X_k <- X_prev[[k]]

    for (tv in time_vars) {
      d_tv_k <- d_k_base
      d_tv_k[[tv]] <- 1
      X_tv_k <- get_X(d_tv_k)
      X_slope <- X_slopes[[tv]]

      if (
        !is.null(X_tv_k) && !is.null(X_0) && !is.null(X_slope) && !is.null(X_k)
      ) {
        X_int <- (X_tv_k - X_0) - X_slope - X_k
        if (max(abs(X_int)) > 1e-10) {
          X_interactions[[k]][[tv]] <- X_int
        } else {
          X_interactions[[k]][[tv]] <- NULL
        }
      }
    }
  }

  # Baseline Y indices (for T=1 setup)
  y_base_fac <- factor(data[[pvarname]], levels = ylevels)
  y_base_idx <- as.integer(y_base_fac)

  list(
    X_0 = X_0,
    X_slopes = X_slopes,
    X_prev = X_prev,
    X_interactions = X_interactions,
    y_base_idx = y_base_idx,
    t_covs = t_covs,
    n_pat = n_pat,
    n_states = n_states,
    M = M,
    ylevels = ylevels,
    col_names = common_cols
  )
}


#' Run Fast Markov Simulation (Internal)
#'
#' Runs the Markov loop using pre-calculated components and a specific
#' coefficient vector.
#'
#' @param components List returned by `markov_msm_build`
#' @param Gamma Effective coefficient matrix with dimensions M x P (M thresholds by P predictors)
#' @param times Vector of time points
#' @param absorb Absorbing state(s). Can be NULL for no absorbing states.
#'
#' @return Array of state probabilities with dimensions n_pat x n_times x n_states
#'
#' @keywords internal
markov_msm_run <- function(components, Gamma, times, absorb = NULL) {
  # Unpack
  X_0 <- components$X_0
  X_slopes <- components$X_slopes
  X_prev <- components$X_prev
  X_int <- components$X_interactions
  y_base_idx <- components$y_base_idx
  t_covs <- components$t_covs
  n_pat <- components$n_pat
  n_states <- components$n_states
  M <- components$M
  ylevels <- components$ylevels
  col_names <- components$col_names
  n_times <- length(times)

  # Align Gamma to components
  Gamma <- Gamma[, col_names, drop = FALSE]
  Gamma_t <- t(Gamma) # [P x M]

  # Helper: X * Gamma_t -> LP [N x M]
  calc_lp <- function(X) {
    if (is.null(X)) {
      return(matrix(0, n_pat, M))
    }
    X %*% Gamma_t
  }

  # Pre-calculate LPs
  LP_0 <- calc_lp(X_0)
  LP_slopes <- lapply(X_slopes, calc_lp)
  LP_prev <- lapply(X_prev, calc_lp)

  LP_int <- vector("list", n_states)
  for (k in 1:n_states) {
    LP_int[[k]] <- list()
    if (!is.null(t_covs)) {
      for (tv in names(t_covs)) {
        if (!is.null(X_int[[k]][[tv]])) {
          LP_int[[k]][[tv]] <- calc_lp(X_int[[k]][[tv]])
        }
      }
    }
  }

  # Simulation
  P_out <- array(0, dim = c(n_pat, n_times, n_states))
  dimnames(P_out)[[3]] <- ylevels

  absorb_idx <- if (!is.null(absorb)) {
    which(as.character(ylevels) %in% as.character(absorb))
  } else {
    integer(0)
  }
  non_absorb_idx <- setdiff(1:n_states, absorb_idx)
  time_vars <- if (!is.null(t_covs)) names(t_covs) else character(0)

  # T=1
  LP_t1 <- LP_0
  for (tv in time_vars) {
    val <- t_covs[1, tv]
    if (val != 0) LP_t1 <- LP_t1 + (val * LP_slopes[[tv]])
  }

  Effect_prev_base <- matrix(0, nrow = n_pat, ncol = M)
  for (k in 1:n_states) {
    mask <- (y_base_idx == k)
    if (any(mask)) Effect_prev_base[mask, ] <- LP_prev[[k]][mask, ]
  }
  LP_t1 <- LP_t1 + Effect_prev_base

  for (tv in time_vars) {
    val <- t_covs[1, tv]
    if (val != 0) {
      for (k in 1:n_states) {
        lp_i <- LP_int[[k]][[tv]]
        if (!is.null(lp_i)) {
          mask <- (y_base_idx == k)
          if (any(mask)) LP_t1[mask, ] <- LP_t1[mask, ] + (val * lp_i[mask, ])
        }
      }
    }
  }

  P_out[, 1, ] <- lp_to_probs(LP_t1, M)

  # Loop 2..T
  for (t_idx in 2:n_times) {
    LP_base_t <- LP_0
    for (tv in time_vars) {
      val <- t_covs[t_idx, tv]
      if (val != 0) LP_base_t <- LP_base_t + (val * LP_slopes[[tv]])
    }

    p_current <- matrix(0, nrow = n_pat, ncol = n_states)

    for (k in non_absorb_idx) {
      p_prev_k <- P_out[, t_idx - 1, k]
      if (max(p_prev_k) < 1e-12) {
        next
      }

      LP_k <- LP_base_t + LP_prev[[k]]

      for (tv in time_vars) {
        val <- t_covs[t_idx, tv]
        if (val != 0) {
          lp_i <- LP_int[[k]][[tv]]
          if (!is.null(lp_i)) LP_k <- LP_k + (val * lp_i)
        }
      }

      trans_k <- lp_to_probs(LP_k, M)
      p_current <- p_current + (trans_k * p_prev_k)
    }

    for (a in absorb_idx) {
      p_current[, a] <- p_current[, a] + P_out[, t_idx - 1, a]
    }
    P_out[, t_idx, ] <- p_current
  }

  P_out
}


#' Convert Cumulative Log-Odds to Probabilities (Internal)
#'
#' Converts a matrix of cumulative log-odds to category probabilities.
#'
#' @param eta Matrix of linear predictors with dimensions N x M (N observations by M thresholds)
#' @param M Number of thresholds (one less than the number of categories)
#' @return Matrix of probabilities with dimensions N x (M+1)
#'
#' @keywords internal
lp_to_probs <- function(eta, M) {
  cum_probs <- stats::plogis(eta) # [N x M]
  n_rows <- nrow(eta)
  probs <- matrix(0, nrow = n_rows, ncol = M + 1)

  probs[, 1] <- 1 - cum_probs[, 1]

  if (M > 1) {
    probs[, 2:M] <- cum_probs[, 1:(M - 1)] - cum_probs[, 2:M]
  }

  probs[, M + 1] <- cum_probs[, M]
  probs[probs < 0] <- 0

  probs
}


#' Compute Effective Coefficients from Beta Vector (Internal)
#'
#' Transforms a vector of VGLM coefficients into an effective coefficient matrix
#' by applying the constraint matrices.
#'
#' @param beta Coefficient vector from a VGLM fit
#' @param C_list Constraint matrices list from VGAM::constraints()
#' @return Effective coefficient matrix with dimensions M x P (M linear predictors by P terms)
#'
#' @keywords internal
compute_Gamma <- function(beta, C_list) {
  M <- nrow(C_list[[1]])
  P <- length(C_list)
  term_names <- names(C_list)

  G <- matrix(0, nrow = M, ncol = P)
  colnames(G) <- term_names
  rownames(G) <- paste0("eta", 1:M)

  curr <- 1
  for (j in seq_along(C_list)) {
    k <- ncol(C_list[[j]])
    chunk <- beta[curr:(curr + k - 1)]
    G[, j] <- C_list[[j]] %*% chunk
    curr <- curr + k
  }
  G
}
