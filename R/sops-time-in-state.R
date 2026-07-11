# SOP time-in-state summaries.

trapezoid_auc <- function(time, value) {
  ok <- !is.na(time) & !is.na(value)
  time <- time[ok]
  value <- value[ok]
  if (length(time) < 2) {
    return(0)
  }
  agg <- stats::aggregate(value ~ time, FUN = sum)
  agg <- agg[order(agg$time), , drop = FALSE]
  if (nrow(agg) < 2) {
    return(0)
  }
  dt <- diff(agg$time)
  sum(dt * (utils::head(agg$value, -1L) + utils::tail(agg$value, -1L)) / 2)
}

time_in_state_tidy <- function(x, target_states, real_time = FALSE) {
  value_col <- if ("estimate" %in% names(x)) {
    "estimate"
  } else if ("draw" %in% names(x)) {
    "draw"
  } else {
    stop("Tidy SOP input must contain an `estimate` or `draw` column.")
  }

  keep <- as.character(x$state) %in% as.character(target_states)
  if (!any(keep)) {
    stop("Target states not found in SOP data frame.")
  }
  x <- x[keep, , drop = FALSE]
  measure_cols <- sop_measure_cols(x)
  group_cols <- setdiff(names(x), c("time", "state", measure_cols))

  agg_formula <- stats::as.formula(
    paste(value_col, "~", paste(c(group_cols, "time"), collapse = " + "))
  )
  by_time <- stats::aggregate(agg_formula, data = x, FUN = sum, na.rm = TRUE)

  if (!real_time) {
    if (length(group_cols) == 0) {
      out <- data.frame(total_time = sum(by_time[[value_col]], na.rm = TRUE))
      return(out)
    }
    total_formula <- stats::as.formula(
      paste(value_col, "~", paste(group_cols, collapse = " + "))
    )
    out <- stats::aggregate(
      total_formula,
      data = by_time,
      FUN = sum,
      na.rm = TRUE
    )
    names(out)[names(out) == value_col] <- "total_time"
    return(out)
  }

  groups <- split(
    seq_len(nrow(by_time)),
    split_key(by_time, group_cols),
    drop = TRUE
  )
  pieces <- lapply(groups, function(idx) {
    group <- by_time[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    meta$total_time <- trapezoid_auc(group$time, group[[value_col]])
    meta
  })
  bind_rows_fill(pieces)
}

time_in_state_bootstrap_df <- function(sops, target_states, real_time = FALSE) {
  if (!all(c("boot_id", "time", "tx") %in% colnames(sops))) {
    stop("Input data frame must contain 'boot_id', 'time', and 'tx' columns.")
  }

  state_cols_available <- grep("^state_", colnames(sops), value = TRUE)
  if (length(state_cols_available) == 0) {
    stop(
      "Input data frame does not contain any columns starting with 'state_'."
    )
  }

  target_suffixes <- as.character(target_states)
  target_cols <- paste0("state_", target_suffixes)
  missing_cols <- setdiff(target_cols, colnames(sops))
  if (length(missing_cols) > 0) {
    stop(paste(
      "The following target state columns were not found in the bootstrap output:",
      paste(missing_cols, collapse = ", ")
    ))
  }

  prob_target <- if (length(target_cols) == 1) {
    sops[[target_cols]]
  } else {
    rowSums(sops[, target_cols, drop = FALSE])
  }

  agg_df <- sops[, c("boot_id", "tx", "time")]
  agg_df$prob <- prob_target

  if (real_time) {
    groups <- split(
      seq_len(nrow(agg_df)),
      split_key(agg_df, c("boot_id", "tx"))
    )
    res_grouped <- bind_rows_fill(lapply(groups, function(idx) {
      group <- agg_df[idx, , drop = FALSE]
      data.frame(
        boot_id = group$boot_id[1],
        tx = group$tx[1],
        prob = trapezoid_auc(group$time, group$prob)
      )
    }))
  } else {
    res_grouped <- stats::aggregate(
      prob ~ boot_id + tx,
      data = agg_df,
      FUN = sum
    )
  }

  res_tx <- res_grouped[res_grouped$tx == 1, c("boot_id", "prob")]
  res_ctrl <- res_grouped[res_grouped$tx == 0, c("boot_id", "prob")]
  res_wide <- merge(
    res_tx,
    res_ctrl,
    by = "boot_id",
    suffixes = c("_tx", "_ctrl")
  )

  colnames(res_wide)[colnames(res_wide) == "prob_tx"] <- "SOP_tx"
  colnames(res_wide)[colnames(res_wide) == "prob_ctrl"] <- "SOP_ctrl"
  res_wide$delta <- res_wide$SOP_tx - res_wide$SOP_ctrl
  res_wide
}

#' Compute Total Time in Target State(s)
#'
#' Calculates the expected total time spent in specified target state(s).
#' By default this sums state occupancy probabilities over visit-scale unit
#' steps. When `time_map` or `origin_time` is supplied, it maps/interpolates to
#' real time and uses trapezoidal AUC. Already interpolated SOP output from
#' [interpolate_sops()] is integrated on its current real-time grid.
#' It automatically adapts to the input format:
#' \itemize{
#'   \item **Patient-Level:** If input is an array from \code{soprob_markov}, it returns time-in-state for each patient.
#'   \item **Legacy Wide Bootstrap-Level:** If input is a data frame with
#'     `boot_id`, `time`, `tx`, and `state_*` columns, it returns mean
#'     time-in-state per treatment group and the difference for each bootstrap
#'     sample.
#' }
#'
#' @param x Input object.
#'   \itemize{
#'     \item **Array:** `[Patients x Time x States]` (Frequentist) or `[Draws x Patients x Time x States]` (Bayes).
#'     \item **Data Frame:** Legacy wide bootstrap SOP data containing columns
#'       `boot_id`, `time`, `tx`, and `state_*`.
#'   }
#' @param target_states Vector of target state(s) to include in the time calculation.
#'   Can be integer indices or character names (e.g., \code{1} or \code{c("Home", "Rehab")}).
#'   For bootstrap outputs, these must match the suffix of the `state_` columns (e.g., if column is `state_1`, use `1`).
#' @param time_map Optional named numeric vector or data frame mapping visit
#'   labels to real elapsed times. Supplying this switches tidy SOP and array
#'   inputs to trapezoidal real-time AUC.
#' @param origin_time Optional real time for an empirical baseline anchor. For
#'   tidy outputs from [sops()] or [avg_sops()], this is passed to
#'   [interpolate_sops()].
#' @param target_times Optional numeric real-time grid for tidy SOP outputs when
#'   `time_map` is supplied. This controls the interpolation grid used for AUC,
#'   for example `target_times = 1:28` with `origin_time = 0` uses day 0 as an anchor
#'   but starts the AUC at day 1.
#' @param origin Origin handling for tidy SOP outputs when `origin_time` is
#'   supplied. See [interpolate_sops()].
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
#' sops_arr <- soprob_markov(model, data, times = 1:30, y_levels = 1:6)
#' days_home <- time_in_state(sops_arr, target_states = 1)
#'
#' # Real-time AUC after fitting on factor visit indices
#' real_days_home <- time_in_state(
#'   avg,
#'   target_states = 1,
#'   time_map = c("1" = 3, "2" = 7, "3" = 14, "4" = 28),
#'   origin_time = 0,
#'   target_times = 1:28
#' )
#'
#' # --- Scenario 2: Bootstrap Inference ---
#' avg_boot <- avg_sops(
#'   model,
#'   variables = list(tx = c(0, 1)),
#'   times = 1:30,
#'   y_levels = 1:6,
#'   id_var = "id"
#' ) |>
#'   inferences(method = "bootstrap", n_draws = 100)
#' boot_effects <- time_in_state(avg_boot, target_states = 1)
#' }
#'
#' @keywords time-in-state auc bootstrap
#' @export
time_in_state <- function(
  x,
  target_states = 1,
  time_map = NULL,
  origin_time = NULL,
  target_times = NULL,
  origin = c("empirical_baseline", "none")
) {
  sops <- x
  origin <- match.arg(origin)
  use_real_time <- !is.null(time_map) || !is.null(origin_time)

  if (inherits(sops, "markov_interpolated_sops")) {
    if (!is.null(time_map) || !is.null(origin_time)) {
      stop(
        "`markov_interpolated_sops` is already on the real-time scale; ",
        "omit `time_map` and `origin_time`."
      )
    }
    if (!is.null(target_times)) {
      target_times <- validate_sop_xout(
        target_times,
        min(sops$time),
        max(sops$time)
      )
      missing_xout <- setdiff(target_times, unique(sops$time))
      if (length(missing_xout) > 0) {
        stop(
          "`target_times` for already interpolated SOPs must be contained in ",
          "the existing `time` values."
        )
      }
      sops <- sops[sops$time %in% target_times, , drop = FALSE]
    }
    return(time_in_state_tidy(sops, target_states, real_time = TRUE))
  }

  if (inherits(sops, c("markov_sops", "markov_avg_sops"))) {
    if (!is.null(target_times) && !use_real_time) {
      stop(
        "`target_times` requires `time_map` or an already interpolated SOP object."
      )
    }
    if (use_real_time) {
      if (is.null(time_map)) {
        stop("`time_map` must be supplied for real-time AUC.")
      }
      sops <- interpolate_sops(
        sops,
        time_map = time_map,
        target_times = target_times,
        origin_time = origin_time,
        origin = origin
      )
      return(time_in_state_tidy(sops, target_states, real_time = TRUE))
    }
    return(time_in_state_tidy(sops, target_states, real_time = FALSE))
  }

  # =========================================================================
  # BRANCH 1: Bootstrap Data Frame Input
  # =========================================================================
  if (inherits(sops, "data.frame")) {
    if (!is.null(target_times)) {
      stop("`target_times` is only supported for tidy SOP outputs.")
    }
    if (use_real_time) {
      if (is.null(time_map)) {
        stop("`time_map` must be supplied for real-time AUC.")
      }
      time_map <- standardize_time_map(time_map)
      sops$time <- map_sop_time_values(sops$time, time_map)
    }
    return(time_in_state_bootstrap_df(
      sops,
      target_states = target_states,
      real_time = use_real_time
    ))
  }

  # =========================================================================
  # BRANCH 2: Array Input (Original Logic)
  # =========================================================================

  if (!is.null(target_times)) {
    stop("`target_times` is only supported for tidy SOP outputs.")
  }

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
  subset_args[[state_dim]] <- match(target_chars, state_names)

  sops_slice <- do.call("[", c(list(sops), subset_args, drop = FALSE))
  prob_in_target <- apply(sops_slice, (1:ndim)[-state_dim], sum)

  # --- 4. Integrate Over Time ---
  if (use_real_time) {
    if (is.null(time_map)) {
      stop("`time_map` must be supplied for real-time AUC.")
    }
    if (!is.null(origin_time) && origin == "empirical_baseline") {
      stop(
        "Empirical `origin_time` anchoring is only supported for tidy SOP outputs."
      )
    }
    time_names <- dnames[[state_dim - 1L]]
    if (is.null(time_names)) {
      time_names <- seq_len(dims[state_dim - 1L])
    }
    time_map <- standardize_time_map(time_map)
    real_time <- map_sop_time_values(time_names, time_map)

    if (ndim == 3) {
      return(apply(prob_in_target, 1, function(z) trapezoid_auc(real_time, z)))
    }
    return(apply(prob_in_target, c(1, 2), function(z) {
      trapezoid_auc(real_time, z)
    }))
  }

  if (ndim == 3) {
    # Freq Input -> [Pat, Time] -> Apply over rows
    res <- rowSums(prob_in_target)
  } else {
    # Bayes Input -> [Draws, Pat, Time] -> Apply over Draws & Pat
    res <- apply(prob_in_target, c(1, 2), sum)
  }

  return(res)
}
