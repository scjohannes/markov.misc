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
    vgam = "vgam"
  )
  ftype <- ftypes[cl]

  if (is.na(ftype)) {
    stop("Object class not supported")
  }

  # Define prediction function
  prd <- switch(
    ftype,
    rms = function(obj, d) predict(obj, d, type = "fitted.ind"),
    vgam = function(obj, d) VGAM::predict(obj, d, type = "response"),
    rmsb = function(obj, d) {
      predict(obj, d, type = "fitted.ind", posterior.summary = "all")
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
      !inherits(model, "rmsb")
  ) {
    stop("model must be an orm, rmsb, or vglm object.")
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

  # Note: ensure soprob_markov_vectorized is sourced/loaded
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
#'   \item **Patient-Level:** If input is an array from \code{\link{soprob_markov_vectorized}}, it returns time-in-state for each patient.
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
    res_grouped <- aggregate(prob ~ boot_id + tx, data = agg_df, FUN = sum)

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
