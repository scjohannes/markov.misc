#' Build Fast Markov Simulation Components
#'
#' Pre-calculates design matrices and interaction structures for fast Markov simulation.
#' This is the "heavy lifting" step that should be done ONCE before running
#' thousands of simulation draws.
#'
#' @param model Fitted vglm model
#' @param data Baseline data frame (one row per patient)
#' @param t_covs Time-dependent covariate lookup (optional)
#' @param ylevels State labels
#' @param pvarname Name of previous state variable
#' @param ... Ignored
#'
#' @return A list containing pre-calculated matrices and metadata.
#' @export
markov_msm_build <- function(
  model,
  data,
  t_covs = NULL,
  ylevels = 1:6,
  pvarname = "yprev",
  ...
) {
  # --- 1. Setup ---
  # If t_covs is NULL, we assume no time-varying covariates via lookup
  # Note: The model might still use linear time (which is updated manually if tvar is handled)
  # But here we focus on the parts driven by t_covs (splines etc)
  
  n_pat <- nrow(data)
  n_states <- length(ylevels)
  M <- n_states - 1 
  
  # Need Gamma structure to know which columns to keep, even if we don't use values yet
  Gamma_template <- get_effective_coefs(model)
  common_cols <- colnames(Gamma_template)
  
  # Prepare data
  if (!pvarname %in% names(data)) {
     data[[pvarname]] <- factor(ylevels[1], levels = ylevels)
  }
  tvar <- "time"
  if (!tvar %in% names(data)) data[[tvar]] <- 0
  
  # Terms object
  tt <- stats::terms(model)
  tt <- stats::delete.response(tt)
  
  # Helper to get Design Matrix (X) only
  get_X <- function(d) {
    X <- stats::model.matrix(tt, data = d)
    common <- intersect(colnames(X), common_cols)
    if (length(common) == 0) return(NULL)
    X[, common, drop = FALSE]
  }
  
  # --- 2. Decomposition ---
  
  # A. X_0 (Base: Time=0, Prev=Ref)
  d_0 <- data
  if (!is.null(t_covs)) {
    for (n in names(t_covs)) d_0[[n]] <- 0
  }
  d_0[[pvarname]] <- factor(ylevels[1], levels = ylevels)
  
  X_0 <- get_X(d_0)
  
  # B. Time Slopes
  X_slopes <- list()
  time_vars <- if (!is.null(t_covs)) names(t_covs) else character(0)
  
  for (tv in time_vars) {
    d_tv <- d_0
    d_tv[[tv]] <- 1
    X_tv <- get_X(d_tv)
    
    # Difference
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
    
    # Base for this state
    d_k_base <- d_0
    d_k_base[[pvarname]] <- factor(ylevels[k], levels = ylevels)
    X_k <- X_prev[[k]] # Already diffed from X_0
    
    for (tv in time_vars) {
      d_tv_k <- d_k_base
      d_tv_k[[tv]] <- 1
      X_tv_k <- get_X(d_tv_k)
      
      X_slope <- X_slopes[[tv]]
      
      # Int = Total - Base - Slope - PrevEffect
      # Int = (X_tv_k - X_0) - X_slope - X_k
      
      if (!is.null(X_tv_k) && !is.null(X_0) && !is.null(X_slope) && !is.null(X_k)) {
        X_int <- (X_tv_k - X_0) - X_slope - X_k
        
        # Optimization: Store NULL if all zero
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


#' Run Fast Markov Simulation
#'
#' Runs the Markov loop using pre-calculated components and a specific coefficient vector.
#'
#' @param components List returned by `markov_msm_build`
#' @param Gamma Effective coefficient matrix with dimensions M x P (M thresholds by P predictors)
#' @param times Vector of time points
#' @param absorb Absorbing state(s). Can be NULL for no absorbing states.
#'
#' @return Array of state probabilities with dimensions n_pat x n_times x n_states
#'
#' @export
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
    if (is.null(X)) return(matrix(0, n_pat, M))
    X %*% Gamma_t
  }
  
  # --- Pre-calculate LPs ---
  
  # LP_0
  LP_0 <- calc_lp(X_0)
  
  # LP Slopes
  LP_slopes <- lapply(X_slopes, calc_lp)
  
  # LP Prev Effects
  LP_prev <- lapply(X_prev, calc_lp)
  
  # LP Interactions
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
  
  # --- Simulation ---
  
  P_out <- array(0, dim = c(n_pat, n_times, n_states))
  dimnames(P_out)[[3]] <- ylevels
  
  absorb_idx <- if (!is.null(absorb)) {
    which(as.character(ylevels) %in% as.character(absorb))
  } else integer(0)
  non_absorb_idx <- setdiff(1:n_states, absorb_idx)
  time_vars <- if(!is.null(t_covs)) names(t_covs) else character(0)
  
  # T=1
  LP_t1 <- LP_0
  
  # Add Time Effects T=1
  for (tv in time_vars) {
    val <- t_covs[1, tv]
    if (val != 0) LP_t1 <- LP_t1 + (val * LP_slopes[[tv]])
  }
  
  # Add Prev Effects (Baseline)
  Effect_prev_base <- matrix(0, nrow=n_pat, ncol=M)
  for (k in 1:n_states) {
    mask <- (y_base_idx == k)
    if (any(mask)) Effect_prev_base[mask, ] <- LP_prev[[k]][mask, ]
  }
  LP_t1 <- LP_t1 + Effect_prev_base
  
  # Add Interactions (Baseline)
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
  for (t in 2:n_times) {
    # LP_base_t
    LP_base_t <- LP_0
    for (tv in time_vars) {
      val <- t_covs[t, tv]
      if (val != 0) LP_base_t <- LP_base_t + (val * LP_slopes[[tv]])
    }
    
    p_current <- matrix(0, nrow = n_pat, ncol = n_states)
    
    for (k in non_absorb_idx) {
      p_prev_k <- P_out[, t-1, k]
      if (max(p_prev_k) < 1e-12) next
      
      # LP_k
      LP_k <- LP_base_t + LP_prev[[k]]
      
      # Interactions
      for (tv in time_vars) {
        val <- t_covs[t, tv]
        if (val != 0) {
          lp_i <- LP_int[[k]][[tv]]
          if (!is.null(lp_i)) LP_k <- LP_k + (val * lp_i)
        }
      }
      
      trans_k <- lp_to_probs(LP_k, M)
      p_current <- p_current + (trans_k * p_prev_k)
    }
    
    for (a in absorb_idx) {
       p_current[, a] <- p_current[, a] + P_out[, t-1, a]
    }
    P_out[, t, ] <- p_current
  }
  
  P_out
}


#' Fast Individual State Occupation Probabilities
#' 
#' Computes individual-level SOPs using the fast Markov engine if possible, 
#' falling back to the standard implementation otherwise.
#'
#' @inheritParams sops
#' @export
sops_fast <- function(
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
  # --- 1. Setup & Defaults (Matched to sops) ---
  if (is.null(newdata)) {
    newdata <- model$x
    if (is.null(newdata)) {
      newdata <- tryCatch(stats::model.frame(model), error = function(e) NULL)
    }
    if (is.null(newdata)) stop("No newdata")
  }
  
  if (!"rowid" %in% names(newdata)) {
    newdata$rowid <- seq_len(nrow(newdata))
  }
  
  if (is.null(times)) {
    if (is.null(tvarname) || !tvarname %in% names(newdata)) {
      stop("`times` must be specified if `tvarname` is not in data.")
    }
    times <- sort(unique(newdata[[tvarname]]))
  }
  
  # Resolve ylevels
  if (is.null(ylevels)) {
    if (inherits(model, "vglm")) {
      ylevels <- model@extra$colnames.y
    } else if (inherits(model, "robcov_vglm")) {
      ylevels <- model$extra$colnames.y
    } else if (inherits(model, "orm")) {
      ylevels <- model$yunique
    }
    if (is.null(ylevels)) ylevels <- 1:6
  }

  # --- 2. Try Fast Path ---
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  is_vglm <- inherits(model_chk, "vglm")
  
  # We can use fast path if it's VGLM. 
  # t_covs is optional now, so we don't strict check it here.
  can_use_fast <- is_vglm
  
  if (can_use_fast) {
     # FAST PATH
     components <- tryCatch(
       markov_msm_build(model_chk, newdata, t_covs, ylevels, pvarname),
       error = function(e) NULL
     )
     
     if (!is.null(components)) {
       Gamma <- get_effective_coefs(model_chk)
       sops_array <- markov_msm_run(components, Gamma, times, absorb)
       
       # Convert array to DF (Tidy)
       result <- array_to_df_individual(sops_array, times, ylevels, newdata)
       
       # Attach Components for fast inference
       attr(result, "msm_components") <- components
       attr(result, "msm_baseline_n") <- nrow(newdata)
       attr(result, "msm_grid") <- NULL # No grid for individual sops
       
     } else {
       # Fallback if build failed
       can_use_fast <- FALSE
     }
  }
  
  if (!can_use_fast) {
    # SLOW PATH (Fallback to soprob_markov logic)
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
    
    result <- array_to_df_individual(sops_array, times, ylevels, newdata)
  }
  
  # Attributes
  attr(result, "model") <- model
  attr(result, "tvarname") <- tvarname
  attr(result, "pvarname") <- pvarname
  attr(result, "ylevels") <- ylevels
  attr(result, "absorb") <- absorb
  attr(result, "gap") <- gap
  attr(result, "t_covs") <- t_covs
  attr(result, "call_args") <- list(
    times = times, ylevels = ylevels, absorb = absorb, 
    tvarname = tvarname, pvarname = pvarname, gap = gap, t_covs = t_covs
  )
  attr(result, "newdata_orig") <- newdata
  
  class(result) <- c("markov_sops", class(result))
  return(result)
}


#' Fast Averaged SOPs (Vectorized & Pre-Calculated)
#'
#' A high-performance version of `avg_sops` that uses pre-calculated effective
#' coefficients to bypass repeated `predict()` calls.
#'
#' @inheritParams avg_sops
#' @param t_covs Optional data frame for non-linear time handling (e.g., splines).
#'   Rows must match the length of \code{times}, with columns matching the
#'   time-varying covariate names used in the model formula.
#' @param ylevels Character or integer vector of outcome state levels. If NULL,
#'   automatically extracted from the model object.
#'
#' @export
avg_sops_fast <- function(
  model,
  newdata = NULL,
  variables = NULL,
  times = NULL,
  id_var = "id",
  t_covs = NULL,
  ylevels = NULL,
  ...
) {
  # Input Validation
  if (is.null(variables)) stop("variables required")
  if (is.null(newdata)) newdata <- model$x
  if (is.null(newdata)) newdata <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (is.null(newdata)) stop("No newdata")
  
  # Save original newdata for inference (bootstrap requires full longitudinal data)
  newdata_orig <- newdata
  
  baseline_data <- newdata[!duplicated(newdata[[id_var]]), ]
  if (!"rowid" %in% names(baseline_data)) {
    baseline_data$rowid <- seq_len(nrow(baseline_data))
  }
  
  # Create Counterfactuals (Stacked)
  grid <- do.call(expand.grid, variables)
  full_data <- create_counterfactual_data(baseline_data, grid, variables)
  
  # Resolve ylevels
  if (is.null(ylevels)) {
     if (inherits(model, "vglm")) ylevels <- model@extra$colnames.y
     else if (inherits(model, "robcov_vglm")) ylevels <- model$extra$colnames.y
     else if (inherits(model, "orm")) ylevels <- model$yunique
     if (is.null(ylevels)) ylevels <- 1:6
  }
  
  # Check Fast Path
  model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
  can_use_fast <- inherits(model_chk, "vglm")
  
  # Extract args from ...
  args <- list(...)
  absorb <- args$absorb
  pvarname <- if (!is.null(args$pvarname)) args$pvarname else "yprev"
  
  res_array <- NULL
  components <- NULL
  
  if (can_use_fast) {
     components <- tryCatch(
       markov_msm_build(model_chk, full_data, t_covs, ylevels, pvarname),
       error = function(e) NULL
     )
     
     if (!is.null(components)) {
        Gamma <- get_effective_coefs(model_chk)
        res_array <- markov_msm_run(components, Gamma, times, absorb)
     } else {
        can_use_fast <- FALSE
     }
  }
  
  if (!can_use_fast) {
    # SLOW PATH: Use standard soprob_markov logic on expanded data
    res_array <- soprob_markov(
      object = model,
      data = full_data,
      times = times,
      ylevels = ylevels,
      absorb = absorb,
      tvarname = if(!is.null(args$tvarname)) args$tvarname else "time",
      pvarname = pvarname,
      gap = args$gap,
      t_covs = t_covs
    )
  }
  
  # Marginalize
  n_pat <- nrow(baseline_data)
  n_cf <- nrow(grid)
  
  res <- marginalize_sops_array(
    sops_array = res_array,
    grid = grid,
    times = times,
    ylevels = ylevels,
    variables = variables,
    n_cf = n_cf,
    n_each = n_pat
  )
  res <- as.data.frame(res)
  class(res) <- c("markov_avg_sops", class(res))
  
  # Attributes
  attr(res, "model") <- model
  attr(res, "ylevels") <- ylevels
  attr(res, "times") <- times
  attr(res, "newdata_orig") <- newdata_orig
  
  # Store components if available (only from fast path)
  if (!is.null(components)) {
    attr(res, "msm_components") <- components
    attr(res, "msm_grid") <- grid
    attr(res, "msm_baseline_n") <- n_pat
  }
  
  attr(res, "t_covs") <- t_covs
  attr(res, "absorb") <- absorb
  attr(res, "avg_args") <- list(
    variables = variables,
    id_var = id_var,
    times = times
  )
  # Standard args needed for reconstruction if fast inference used later on slow object
  attr(res, "call_args") <- args
  attr(res, "pvarname") <- pvarname
  
  return(res)
}

#' Fast Inference for State Occupation Probabilities
#'
#' Adds confidence intervals using either fast MVN simulation or bootstrap.
#'
#' @inheritParams inferences
#' @export
inferences_fast <- function(
  object, 
  method = "simulation",
  n_sim = 1000, 
  n_boot = NULL,
  vcov = NULL, 
  conf_level = 0.95,
  parallel = FALSE,
  workers = NULL,
  ...
) {
  # Argument handling
  if (!is.null(n_boot) && method == "bootstrap") {
    n_sim <- n_boot
  }
  
  method <- match.arg(method, choices = c("simulation", "bootstrap"))
  
  if (method == "simulation") {
    return(inferences_fast_simulation(
      object = object,
      n_sim = n_sim,
      vcov = vcov,
      conf_level = conf_level,
      parallel = parallel,
      workers = workers
    ))
  } else if (method == "bootstrap") {
    return(inferences_fast_bootstrap(
      object = object,
      n_boot = n_sim,
      conf_level = conf_level,
      parallel = parallel,
      workers = workers,
      ...
    ))
  }
}

#' Fast MVN Simulation Inference (Internal)
#' 
#' @keywords internal
inferences_fast_simulation <- function(
  object, 
  n_sim = 1000, 
  vcov = NULL, 
  conf_level = 0.95,
  parallel = FALSE,
  workers = NULL
) {
  # Support both types
  if (!inherits(object, c("markov_avg_sops", "markov_sops"))) {
    stop("Only markov_avg_sops or markov_sops objects supported")
  }
  
  is_avg <- inherits(object, "markov_avg_sops")
  
  # Try to retrieve components
  components <- attr(object, "msm_components")
  
  # If components missing, try to build them (Auto-Build)
  if (is.null(components)) {
    model <- attr(object, "model")
    model_chk <- if (inherits(model, "robcov_vglm")) model$vglm_fit else model
    
    # Can only build if VGLM
    if (!inherits(model_chk, "vglm")) {
       stop("Fast simulation inference requires a VGLM model (or robcov_vglm). For other models, use standard inferences().")
    }
    
    # Reconstruct Build Args
    newdata_orig <- attr(object, "newdata_orig")
    t_covs <- attr(object, "t_covs")
    ylevels <- attr(object, "ylevels")
    pvarname <- attr(object, "pvarname")
    
    # Prepare Data
    if (is_avg) {
      # Need to expand
      avg_args <- attr(object, "avg_args")
      id_var <- avg_args$id_var
      variables <- avg_args$variables
      
      baseline_data <- newdata_orig[!duplicated(newdata_orig[[id_var]]), ]
      grid <- do.call(expand.grid, variables)
      full_data <- create_counterfactual_data(baseline_data, grid, variables)
      
      # Set attributes needed for marginalization later
      attr(object, "msm_grid") <- grid
      attr(object, "msm_baseline_n") <- nrow(baseline_data)
      
    } else {
      full_data <- newdata_orig
    }
    
    # Build
    components <- markov_msm_build(
       model = model_chk, 
       data = full_data, 
       t_covs = t_covs, 
       ylevels = ylevels, 
       pvarname = pvarname
    )
    # Save back to object for future use? (Optional, but good)
    attr(object, "msm_components") <- components
  }
  
  # Setup for Simulation
  model <- attr(object, "model") # Original model object
  grid <- attr(object, "msm_grid")
  n_pat_base <- attr(object, "msm_baseline_n") # For avg
  n_pat_ind <- components$n_pat # For ind
  
  # Retrieve times robustly
  times <- attr(object, "times")
  if (is.null(times)) {
     if (is_avg) {
        times <- attr(object, "avg_args")$times
     } else {
        times <- attr(object, "call_args")$times
     }
  }
  if (is.null(times)) stop("Could not determine 'times' from object attributes.")
  
  # 1. Draw Coefficients
  beta_hat <- coef(model)
  if (is.null(vcov)) vcov <- get_vcov_robust(model)
  
  if (length(beta_hat) != nrow(vcov)) stop("Coef/Vcov mismatch")
  
  beta_draws <- mvtnorm::rmvnorm(n_sim, mean = beta_hat, sigma = vcov)
  
  # Helper to compute Gamma from beta vector (Defined locally)
  C_list <- VGAM::constraints(if (inherits(model, "robcov_vglm")) model$vglm_fit else model)
  compute_Gamma_local <- function(beta) {
      compute_Gamma(beta, C_list)
  }
  
  # 2. Simulation Loop
  run_draw <- function(i) {
    beta_i <- beta_draws[i, ]
    Gamma_i <- compute_Gamma_local(beta_i)
    
    # Run Simulation
    res_array <- markov_msm_run(components, Gamma_i, times)
    
    # Output Processing
    if (is_avg) {
       # Marginalize
       n_cf <- nrow(grid)
       draw_res <- vector("list", n_cf)
       for (cf in 1:n_cf) {
          idx1 <- (cf-1)*n_pat_base + 1
          idx2 <- cf*n_pat_base
          sub <- res_array[idx1:idx2, , , drop=FALSE]
          draw_res[[cf]] <- apply(sub, c(2,3), mean) # [Time x State]
       }
       return(draw_res)
    } else {
       # Individual - just flatten/return as vector for efficiency
       # Or return the array? Typically we want quantiles per row.
       # Returning the whole array might be huge.
       # Let's return vectors to match existing pattern.
       return(as.vector(res_array))
    }
  }
  
  if (parallel) {
    if (is.null(workers)) workers <- parallel::detectCores() - 1
    future::plan(future::multisession, workers = workers)
    sim_results <- furrr::future_map(1:n_sim, run_draw, .options = furrr::furrr_options(seed=TRUE))
  } else {
    sim_results <- lapply(1:n_sim, run_draw)
  }
  
  # 3. Aggregation for CIs
  flatten_draw_avg <- function(d_list) {
     unlist(lapply(d_list, as.vector))
  }
  
  if (is_avg) {
    draws_mat <- vapply(sim_results, flatten_draw_avg, numeric(nrow(object)))
  } else {
    draws_mat <- vapply(sim_results, identity, numeric(nrow(object)))
  }
  
  alpha <- (1 - conf_level)/2
  ci_lower <- apply(draws_mat, 1, quantile, probs = alpha)
  ci_upper <- apply(draws_mat, 1, quantile, probs = 1 - alpha)
  
  object$conf.low <- ci_lower
  object$conf.high <- ci_upper
  
  attr(object, "n_sim") <- n_sim
  
  object
}

#' Fast Bootstrap Inference (Internal)
#'
#' @keywords internal
inferences_fast_bootstrap <- function(
  object,
  n_boot = 1000,
  conf_level = 0.95,
  parallel = FALSE,
  workers = NULL,
  ...
) {
  if (!inherits(object, "markov_avg_sops")) stop("Only avg_sops objects supported for bootstrap")
  
  # Extract needed attributes
  model <- attr(object, "model")
  newdata_orig <- attr(object, "newdata_orig")
  ylevels <- attr(object, "ylevels")
  absorb <- attr(object, "absorb")
  times <- attr(object, "times")
  t_covs <- attr(object, "t_covs")
  avg_args <- attr(object, "avg_args")
  
  if (is.null(newdata_orig)) {
    stop("Original data not found. avg_sops_fast() must be run with newdata containing all time points.")
  }
  
  res <- bootstrap_standardized_sops(
    model = model,
    data = newdata_orig,
    n_boot = n_boot,
    workers = workers,
    parallel = parallel,
    ylevels = factor(ylevels),
    absorb = absorb,
    times = times,
    t_covs = t_covs,
    varnames = list(
      tvarname = "time",
      pvarname = "yprev",
      id = avg_args$id_var,
      tx = names(avg_args$variables)[1] 
    ),
    ...
  )
  
  boot_df <- res$sops
  
  if (requireNamespace("tidyr", quietly = TRUE) && requireNamespace("dplyr", quietly = TRUE)) {
    long_boot <- boot_df |>
      tidyr::pivot_longer(
        cols = dplyr::starts_with("state_"),
        names_to = "state",
        values_to = "estimate",
        names_prefix = "state_"
      )
    
    alpha <- (1 - conf_level) / 2
    
    ci_df <- long_boot |>
      dplyr::group_by(time, state, tx) |>
      dplyr::summarise(
        conf.low = quantile(estimate, alpha, na.rm = TRUE),
        conf.high = quantile(estimate, 1 - alpha, na.rm = TRUE),
        .groups = "drop"
      )
    
    ci_df$state <- as.factor(ci_df$state)
    
    tx_var <- names(avg_args$variables)[1]
    if (tx_var != "tx") {
      colnames(ci_df)[colnames(ci_df) == "tx"] <- tx_var
    }
    
    final <- merge(object, ci_df, by = c("time", "state", tx_var), all.x = TRUE)
    class(final) <- class(object)
    attr(final, "n_boot") <- n_boot
    return(final)
  } else {
    stop("tidyr and dplyr required for bootstrap inference")
  }
}

# Maintain compatibility wrapper
soprob_markov_fast <- function(model, data, times, t_covs, ylevels=1:6, absorb=NULL, ...) {
   components <- markov_msm_build(model, data, t_covs, ylevels, ...)
   Gamma <- get_effective_coefs(model)
   markov_msm_run(components, Gamma, times, absorb)
}


#' Helper: Convert Cumulative Log-Odds to Probabilities
#' 
#' Converts a matrix of cumulative log-odds (linear predictors) to category
#' probabilities for an ordinal model.
#' 
#' @param eta Matrix of linear predictors with dimensions N x M (N observations by M thresholds)
#' @param M Number of thresholds (one less than the number of categories)
#' @return Matrix of probabilities with dimensions N x (M+1) (N observations by M+1 categories)
#' 
#' @keywords internal
lp_to_probs <- function(eta, M) {
  cum_probs <- stats::plogis(eta) # [N x M]
  n_rows <- nrow(eta)
  probs <- matrix(0, nrow = n_rows, ncol = M + 1)
  
  probs[, 1] <- 1 - cum_probs[, 1]
  
  if (M > 1) {
    probs[, 2:M] <- cum_probs[, 1:(M-1)] - cum_probs[, 2:M]
  }
  
  probs[, M + 1] <- cum_probs[, M]
  probs[probs < 0] <- 0
  
  return(probs)
}

# --- Internal Helpers ---

create_counterfactual_data <- function(baseline_data, grid, variables) {
  cf_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    dt <- baseline_data
    for (v in names(grid)) dt[[v]] <- grid[i, v]
    cf_list[[i]] <- dt
  }
  do.call(rbind, cf_list)
}

marginalize_sops_array <- function(sops_array, grid, times, ylevels, variables, n_cf, n_each) {
  out_list <- vector("list", n_cf)
  n_states <- dim(sops_array)[3]
  use_ylevels <- if(is.null(ylevels)) 1:n_states else ylevels
  
  for (i in seq_len(n_cf)) {
    idx_start <- (i - 1) * n_each + 1
    idx_end <- i * n_each
    sub_arr <- sops_array[idx_start:idx_end, , , drop = FALSE]
    avg_mat <- apply(sub_arr, c(2, 3), mean)
    
    df <- expand.grid(time = times, state = use_ylevels)
    df$estimate <- as.vector(avg_mat)
    for (v in names(grid)) df[[v]] <- grid[i, v]
    out_list[[i]] <- df
  }
  do.call(rbind, out_list)
}

array_to_df_individual <- function(sops_array, times, ylevels, newdata) {
  n_pat <- dim(sops_array)[1]
  n_times <- dim(sops_array)[2]
  n_states <- dim(sops_array)[3]
  
  probs_flat <- as.vector(sops_array)
  
  idx_pat <- rep(seq_len(n_pat), times = n_times * n_states)
  idx_time <- rep(rep(times, each = n_pat), times = n_states)
  idx_state <- rep(ylevels, each = n_pat * n_times)
  
  result <- newdata[idx_pat, , drop = FALSE]
  result$time <- idx_time
  result$state <- idx_state
  result$estimate <- probs_flat
  rownames(result) <- NULL
  
  result
}

#' Compute Effective Coefficient Matrix (Gamma)
#'
#' Transforms a vector of VGLM coefficients into an effective coefficient matrix
#' by applying the constraint matrices. This is a key step in the fast path for
#' simulation-based inference.
#'
#' @param beta Numeric vector of model coefficients from a VGLM fit.
#' @param C_list List of constraint matrices from \code{VGAM::constraints(model)}.
#'   Each element corresponds to a model term and has dimensions M x k, where
#'   M is the number of linear predictors (thresholds) and k is the number of
#'   free parameters for that term.
#'
#' @return A matrix of effective coefficients with dimensions M x P, where M is
#'   the number of linear predictors and P is the number of model terms. Row
#'   names are "eta1", "eta2", etc. and column names are the term names.
#'
#' @details
#' In VGAM models with constraints (e.g., proportional odds with
#' \code{parallel = TRUE}), the relationship between raw coefficients and

#' linear predictors is mediated by constraint matrices. For each term j:
#' \deqn{G_{.,j} = C_j \times \beta_j}
#' where \eqn{C_j} is the constraint matrix and \eqn{\beta_j} is the subset
#' of coefficients for that term.
#'
#' @keywords internal
compute_Gamma <- function(beta, C_list) {
  M <- nrow(C_list[[1]])
  P <- length(C_list)
  term_names <- names(C_list)
  
  G <- matrix(0, nrow=M, ncol=P)
  colnames(G) <- term_names
  rownames(G) <- paste0("eta", 1:M)
  
  # Map beta to chunks
  curr <- 1
  for (j in seq_along(C_list)) {
     k <- ncol(C_list[[j]])
     chunk <- beta[curr:(curr+k-1)]
     G[, j] <- C_list[[j]] %*% chunk
     curr <- curr + k
  }
  G
}
