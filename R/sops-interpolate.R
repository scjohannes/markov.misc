# SOP interpolation helpers.

standardize_time_map <- function(time_map) {
  if (is.null(time_map)) {
    stop("`time_map` must be supplied.")
  }

  if (is.data.frame(time_map)) {
    if (all(c("visit", "real_time") %in% names(time_map))) {
      visit <- time_map$visit
      real_time <- time_map$real_time
    } else if (all(c("visit", "time") %in% names(time_map))) {
      visit <- time_map$visit
      real_time <- time_map$time
    } else if (ncol(time_map) >= 2) {
      visit <- time_map[[1]]
      real_time <- time_map[[2]]
    } else {
      stop("`time_map` data frames must have at least two columns.")
    }
  } else if (is.atomic(time_map) && !is.null(names(time_map))) {
    visit <- names(time_map)
    real_time <- time_map
  } else {
    stop(
      "`time_map` must be a named numeric vector or a data frame with ",
      "visit and real-time columns."
    )
  }

  real_time <- suppressWarnings(as.numeric(real_time))
  visit <- as.character(visit)
  bad <- is.na(visit) |
    !nzchar(visit) |
    is.na(real_time) |
    !is.finite(real_time)
  if (any(bad)) {
    stop("`time_map` contains missing or non-finite visit/time values.")
  }
  if (anyDuplicated(visit)) {
    dup <- unique(visit[duplicated(visit)])
    stop(
      "`time_map` contains duplicate visit entries: ",
      paste(dup, collapse = ", ")
    )
  }

  data.frame(visit = visit, real_time = real_time, stringsAsFactors = FALSE)
}

map_sop_time_values <- function(times, time_map) {
  labels <- as.character(times)
  idx <- match(labels, time_map$visit)
  if (anyNA(idx)) {
    missing <- unique(labels[is.na(idx)])
    stop(
      "`time_map` is missing entries for visit value(s): ",
      paste(missing, collapse = ", ")
    )
  }
  time_map$real_time[idx]
}

validate_sop_xout <- function(target_times, lower, upper) {
  if (
    !is.numeric(target_times) ||
      anyNA(target_times) ||
      any(!is.finite(target_times))
  ) {
    stop("`target_times` must be a finite numeric vector.")
  }
  if (any(target_times < lower | target_times > upper)) {
    stop(
      "`target_times` must stay within the supported time range [",
      lower,
      ", ",
      upper,
      "]."
    )
  }
  sort(unique(target_times))
}

sop_measure_cols <- function(x) {
  intersect(
    c("estimate", "conf.low", "conf.high", "std.error", "draw"),
    names(x)
  )
}

split_key <- function(data, cols) {
  if (length(cols) == 0) {
    return(rep("all", nrow(data)))
  }
  do.call(interaction, c(data[, cols, drop = FALSE], drop = TRUE, sep = "\r"))
}

coerce_state_like <- function(labels, prototype) {
  if (is.factor(prototype)) {
    return(factor(
      labels,
      levels = levels(prototype),
      ordered = is.ordered(prototype)
    ))
  }
  if (is.integer(prototype)) {
    return(as.integer(labels))
  }
  if (is.numeric(prototype)) {
    return(as.numeric(labels))
  }
  labels
}

rows_match_values <- function(data, values) {
  keep <- rep(TRUE, nrow(data))
  for (nm in names(values)) {
    keep <- keep & as.character(data[[nm]]) == as.character(values[[nm]])
  }
  keep
}

baseline_rows_for_anchor <- function(x, id_var) {
  id_var <- id_var %||% attr(x, "id_var")
  newdata_orig <- attr(x, "newdata_orig")
  if (is.null(newdata_orig)) {
    stop(
      "`origin_time` requires SOP output with stored `newdata_orig` ",
      "attributes."
    )
  }

  if (isTRUE(attr(x, "newdata_supplied"))) {
    newdata_pred <- attr(x, "newdata_pred")
    if (inherits(x, "markov_sops") && !is.null(newdata_pred)) {
      return(newdata_pred)
    }
    return(newdata_orig)
  }

  if (!is.null(id_var) && id_var %in% names(newdata_orig)) {
    return(resolve_markov_prediction_data(
      newdata_orig,
      id_var = id_var,
      time_var = attr(x, "time_var") %||% "time"
    ))
  }

  newdata_pred <- attr(x, "newdata_pred")
  if (!is.null(newdata_pred)) {
    return(newdata_pred)
  }
  newdata_orig
}

state_distribution_anchor <- function(
  x,
  combos,
  group_cols,
  filter_cols,
  states,
  origin_time,
  baseline,
  p_var
) {
  n_out <- max(1L, nrow(combos)) * length(states)
  anchor <- x[rep(NA_integer_, n_out), , drop = FALSE]
  row <- 1L

  for (combo_i in seq_len(max(1L, nrow(combos)))) {
    combo <- if (length(group_cols) > 0) {
      combos[combo_i, group_cols, drop = FALSE]
    } else {
      data.frame()
    }

    baseline_i <- baseline
    use_filters <- intersect(filter_cols, names(combo))
    if (length(use_filters) > 0) {
      baseline_i <- baseline_i[
        rows_match_values(baseline_i, combo[use_filters]),
        ,
        drop = FALSE
      ]
    }
    if (nrow(baseline_i) == 0) {
      stop("No baseline rows are available for an empirical origin anchor.")
    }

    state_idx <- match(as.character(baseline_i[[p_var]]), states)
    probs <- tabulate(state_idx, nbins = length(states)) / nrow(baseline_i)
    rows <- row:(row + length(states) - 1L)

    for (nm in group_cols) {
      anchor[[nm]][rows] <- combo[[nm]][1]
    }
    anchor$state[rows] <- coerce_state_like(states, x$state)
    anchor$estimate[rows] <- probs
    if ("conf.low" %in% names(anchor)) {
      anchor$conf.low[rows] <- probs
    }
    if ("conf.high" %in% names(anchor)) {
      anchor$conf.high[rows] <- probs
    }
    if ("std.error" %in% names(anchor)) {
      anchor$std.error[rows] <- 0
    }
    if ("draw" %in% names(anchor)) {
      anchor$draw[rows] <- probs
    }
    anchor$.sop_real_time[rows] <- origin_time
    row <- row + length(states)
  }

  anchor
}

empirical_baseline_anchor <- function(x, origin_time) {
  p_var <- attr(x, "p_var") %||% "yprev"
  if (!"state" %in% names(x)) {
    stop("SOP output must contain a `state` column.")
  }
  if (!"estimate" %in% names(x)) {
    stop("SOP output must contain an `estimate` column.")
  }

  states <- as_state_labels(attr(x, "y_levels") %||% unique(x$state))

  if (inherits(x, "markov_avg_sops")) {
    avg_args <- attr(x, "avg_args")
    variables <- names(avg_args$variables %||% list())
    by <- avg_args$by %||% character()
    group_cols <- unique(c(variables, by))
    combos <- if (length(group_cols) > 0) {
      unique(x[, group_cols, drop = FALSE])
    } else {
      data.frame(.dummy = 1L)
    }
    if (!p_var %in% names(attr(x, "newdata_orig"))) {
      stop(
        "Previous-state variable `",
        p_var,
        "` not found in `newdata_orig`."
      )
    }
    baseline <- baseline_rows_for_anchor(x, attr(x, "avg_args")$id_var)
    return(state_distribution_anchor(
      x = x,
      combos = combos,
      group_cols = group_cols,
      filter_cols = by,
      states = states,
      origin_time = origin_time,
      baseline = baseline,
      p_var = p_var
    ))
  }

  if (p_var %in% names(x)) {
    value_cols <- sop_measure_cols(x)
    meta_cols <- setdiff(names(x), c("time", value_cols, ".sop_real_time"))
    anchor_meta <- unique(x[, meta_cols, drop = FALSE])
    anchor <- x[rep(NA_integer_, nrow(anchor_meta)), , drop = FALSE]
    for (nm in meta_cols) {
      anchor[[nm]] <- anchor_meta[[nm]]
    }
    anchor$estimate <- as.numeric(
      as.character(anchor$state) == as.character(anchor[[p_var]])
    )
    if ("conf.low" %in% names(anchor)) {
      anchor$conf.low <- anchor$estimate
    }
    if ("conf.high" %in% names(anchor)) {
      anchor$conf.high <- anchor$estimate
    }
    if ("std.error" %in% names(anchor)) {
      anchor$std.error <- 0
    }
    if ("draw" %in% names(anchor)) {
      anchor$draw <- anchor$estimate
    }
    anchor$.sop_real_time <- origin_time
    return(anchor)
  }

  by <- attr(x, "by") %||% character()
  group_cols <- by
  combos <- if (length(group_cols) > 0) {
    unique(x[, group_cols, drop = FALSE])
  } else {
    data.frame(.dummy = 1L)
  }
  baseline <- baseline_rows_for_anchor(x, NULL)
  if (!p_var %in% names(baseline)) {
    stop(
      "Previous-state variable `",
      p_var,
      "` not found in `newdata_orig`."
    )
  }
  state_distribution_anchor(
    x = x,
    combos = combos,
    group_cols = group_cols,
    filter_cols = by,
    states = states,
    origin_time = origin_time,
    baseline = baseline,
    p_var = p_var
  )
}

interpolate_numeric_column <- function(time, value, target_times) {
  ok <- !is.na(time) & !is.na(value)
  time <- time[ok]
  value <- value[ok]

  if (length(time) == 0) {
    return(rep(NA_real_, length(target_times)))
  }

  agg <- stats::aggregate(value ~ time, FUN = mean)
  agg <- agg[order(agg$time), , drop = FALSE]
  if (nrow(agg) == 1L) {
    out <- rep(NA_real_, length(target_times))
    out[target_times == agg$time] <- agg$value
    return(out)
  }

  stats::approx(
    x = agg$time,
    y = agg$value,
    xout = target_times,
    rule = 1,
    ties = "ordered"
  )$y
}

normalize_interpolated_sops <- function(x) {
  value_col <- if ("estimate" %in% names(x)) {
    "estimate"
  } else if ("draw" %in% names(x)) {
    "draw"
  } else {
    return(x)
  }

  key_cols <- setdiff(names(x), c("state", sop_measure_cols(x)))
  key <- split_key(x, key_cols)
  sums <- ave(x[[value_col]], key, FUN = function(z) sum(z, na.rm = TRUE))
  scale <- !is.na(sums) & sums > 0
  x[[value_col]][scale] <- x[[value_col]][scale] / sums[scale]
  x
}

sop_draw_attr_name <- function(x) {
  for (nm in c("draws", "draws", "draws")) {
    if (!is.null(attr(x, nm))) {
      return(nm)
    }
  }
  NULL
}

sop_has_uncertainty_cols <- function(x) {
  any(c("conf.low", "conf.high", "std.error") %in% names(x))
}

warn_missing_interpolation_draws <- function(x) {
  if (!sop_has_uncertainty_cols(x)) {
    return(invisible(NULL))
  }

  if (identical(attr(x, "method"), "posterior")) {
    warning(
      "`interpolate_sops()` found uncertainty columns but no stored posterior ",
      "draws. It will interpolate interval endpoints. For draw-level ",
      "interpolation, rerun `sops()` or `avg_sops()` with `return_draws = TRUE`.",
      call. = FALSE
    )
  } else {
    warning(
      "`interpolate_sops()` found uncertainty columns but no stored draws. ",
      "It will interpolate interval endpoints. For draw-level interpolation, ",
      "rerun `inferences(..., return_draws = TRUE)`.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

sop_draw_value_col <- function(draws) {
  if ("estimate" %in% names(draws)) {
    return("estimate")
  }
  if ("draw" %in% names(draws)) {
    return("draw")
  }
  NULL
}

repeat_anchor_for_draws <- function(anchor, draws, value_col) {
  if (is.null(anchor) || !"draw_id" %in% names(draws)) {
    return(NULL)
  }

  draw_ids <- sort(unique(draws$draw_id))
  anchor_value <- anchor$estimate
  anchor <- anchor[,
    intersect(names(anchor), c(names(draws), ".sop_real_time")),
    drop = FALSE
  ]
  anchor[[value_col]] <- anchor_value

  bind_rows_fill(lapply(draw_ids, function(draw_id) {
    anchor_i <- anchor
    anchor_i$draw_id <- draw_id
    anchor_i
  }))
}

interpolate_sop_draws <- function(
  draws,
  time_map,
  target_times,
  anchor,
  normalize
) {
  if (
    !is.data.frame(draws) ||
      !all(c("draw_id", "time", "state") %in% names(draws))
  ) {
    return(NULL)
  }
  value_col <- sop_draw_value_col(draws)
  if (is.null(value_col)) {
    return(NULL)
  }

  work <- as.data.frame(draws)
  work$.sop_real_time <- map_sop_time_values(work$time, time_map)
  draw_anchor <- repeat_anchor_for_draws(anchor, work, value_col)
  work <- bind_rows_fill(list(draw_anchor, work))

  group_cols <- setdiff(names(work), c("time", value_col, ".sop_real_time"))
  groups <- split(seq_len(nrow(work)), split_key(work, group_cols), drop = TRUE)
  pieces <- lapply(groups, function(idx) {
    group <- work[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    out <- meta[rep(1L, length(target_times)), , drop = FALSE]
    out$time <- target_times
    out[[value_col]] <- interpolate_numeric_column(
      time = group$.sop_real_time,
      value = group[[value_col]],
      target_times = target_times
    )
    out
  })

  result <- bind_rows_fill(pieces)
  result <- result[, intersect(names(draws), names(result)), drop = FALSE]
  if (isTRUE(normalize)) {
    result <- normalize_interpolated_sops(result)
  }
  result
}

interpolated_ci_from_draws <- function(draws, result, conf_level, conf_type) {
  value_col <- sop_draw_value_col(draws)
  if (is.null(value_col) || !"draw_id" %in% names(draws)) {
    return(result)
  }

  draws_for_ci <- draws
  if (value_col != "estimate") {
    draws_for_ci$estimate <- draws_for_ci[[value_col]]
  }
  group_cols <- setdiff(names(draws), c(value_col, "draw_id"))
  ci <- compute_ci_from_draws(
    draws_df = draws_for_ci,
    group_cols = group_cols,
    conf_level = conf_level,
    conf_type = conf_type
  )

  join_cols <- intersect(group_cols, names(result))
  ci <- ci[, c(join_cols, "conf.low", "conf.high", "std.error"), drop = FALSE]

  result$.sop_order <- seq_len(nrow(result))
  result_base <- result[,
    setdiff(names(result), c("conf.low", "conf.high", "std.error")),
    drop = FALSE
  ]
  out <- left_join_preserve_order(result_base, ci, by = join_cols)
  out <- out[order(out$.sop_order), , drop = FALSE]
  out$.sop_order <- NULL
  rownames(out) <- NULL
  out
}

#' Interpolate SOPs from Visit Time to Real Time
#'
#' Maps visit-scale SOP output from [sops()] or [avg_sops()] to a real elapsed
#' time scale and linearly interpolates probabilities within patient, state,
#' draw, treatment, and strata groups.
#'
#' @param x A `markov_sops` or `markov_avg_sops` object.
#' @param time_map A named numeric vector mapping visit labels to real times, or
#'   a data frame with visit and real-time columns. Data frames may use columns
#'   `visit` and `real_time`, `visit` and `time`, or their first two columns.
#' @param target_times Optional numeric real-time grid. If `NULL`, uses mapped visit
#'   times, plus `origin_time` when an empirical origin anchor is requested.
#' @param origin_time Optional real time for an empirical baseline anchor.
#' @param origin Origin handling. `"empirical_baseline"` adds an anchor from the
#'   stored `newdata_orig` and previous-state variable; `"none"` does not.
#' @param normalize Logical. If `TRUE`, normalize interpolated estimates so
#'   state probabilities sum to one within each time and group.
#'
#' @return An interpolated data frame with the same core columns as `x`, where
#'   `time` is on the real-time scale.
#'
#' @details
#' The recommended workflow for irregular assessment schedules is:
#' 1. Recode assessment days to visit indices, for example day 3, 7, 14, and 28
#'    to visit `1:4`.
#' 2. Fit the Markov model with the visit index as a factor.
#' 3. Compute visit-scale SOPs with [sops()] or [avg_sops()].
#' 4. Use `interpolate_sops()` or `time_in_state(..., time_map = ...)` for
#'    real-day summaries.
#'
#' With RCT standardization from [avg_sops()], the empirical baseline anchor is
#' shared across treatment counterfactual groups.
#'
#' If `x` contains stored simulation, bootstrap, or posterior draws, the draws
#' are interpolated too. Interval columns are then recomputed from the
#' interpolated draws, with the empirical origin anchor treated as fixed.
#' If interval columns are present but stored draws are unavailable,
#' `interpolate_sops()` warns and falls back to interpolating interval endpoints.
#'
#' @examples
#' \dontrun{
#' avg <- avg_sops(
#'   fit,
#'   newdata = baseline,
#'   variables = list(tx = c(0, 1)),
#'   times = NULL
#' )
#' interpolate_sops(
#'   avg,
#'   time_map = c("1" = 3, "2" = 7, "3" = 14, "4" = 28),
#'   target_times = 0:28,
#'   origin_time = 0
#' )
#' }
#'
#' @export
interpolate_sops <- function(
  x,
  time_map,
  target_times = NULL,
  origin_time = NULL,
  origin = c("empirical_baseline", "none"),
  normalize = TRUE
) {
  if (!inherits(x, c("markov_sops", "markov_avg_sops"))) {
    stop("`x` must be a `markov_sops` or `markov_avg_sops` object.")
  }
  if (!all(c("time", "state") %in% names(x))) {
    stop("`x` must contain `time` and `state` columns.")
  }

  origin <- match.arg(origin)
  time_map <- standardize_time_map(time_map)
  real_time <- map_sop_time_values(x$time, time_map)
  mapped_range <- range(time_map$real_time)

  use_origin <- !is.null(origin_time) && origin == "empirical_baseline"
  if (!is.null(origin_time)) {
    if (
      !is.numeric(origin_time) || length(origin_time) != 1 || is.na(origin_time)
    ) {
      stop("`origin_time` must be a single finite numeric value.")
    }
    if (!is.finite(origin_time)) {
      stop("`origin_time` must be a single finite numeric value.")
    }
  }

  lower <- if (use_origin) origin_time else mapped_range[1]
  upper <- mapped_range[2]
  if (lower > upper) {
    stop("`origin_time` must not be greater than the largest mapped time.")
  }

  if (is.null(target_times)) {
    target_times <- sort(unique(c(
      if (use_origin) origin_time else NULL,
      real_time
    )))
  }
  target_times <- validate_sop_xout(target_times, lower, upper)

  work <- x
  work$.sop_real_time <- real_time

  anchor <- NULL
  if (use_origin) {
    anchor <- empirical_baseline_anchor(work, origin_time)
    work <- bind_rows_fill(list(anchor, work))
  } else {
    work <- as.data.frame(work)
  }

  value_cols <- sop_measure_cols(work)
  if (length(value_cols) == 0) {
    stop("`x` must contain an `estimate` or `draw` column to interpolate.")
  }
  group_cols <- setdiff(names(work), c("time", value_cols, ".sop_real_time"))
  groups <- split(seq_len(nrow(work)), split_key(work, group_cols), drop = TRUE)

  pieces <- lapply(groups, function(idx) {
    group <- work[idx, , drop = FALSE]
    meta <- group[1, group_cols, drop = FALSE]
    out <- meta[rep(1L, length(target_times)), , drop = FALSE]
    out$time <- target_times
    for (nm in value_cols) {
      out[[nm]] <- interpolate_numeric_column(
        time = group$.sop_real_time,
        value = group[[nm]],
        target_times = target_times
      )
    }
    out
  })

  result <- bind_rows_fill(pieces)
  result <- result[, intersect(names(x), names(result)), drop = FALSE]
  if (isTRUE(normalize)) {
    result <- normalize_interpolated_sops(result)
  }

  draw_attr <- sop_draw_attr_name(x)
  if (is.null(draw_attr)) {
    warn_missing_interpolation_draws(x)
  }
  interpolated_draws <- if (!is.null(draw_attr)) {
    interpolate_sop_draws(
      draws = attr(x, draw_attr),
      time_map = time_map,
      target_times = target_times,
      anchor = anchor,
      normalize = normalize
    )
  } else {
    NULL
  }
  if (!is.null(interpolated_draws)) {
    result <- interpolated_ci_from_draws(
      draws = interpolated_draws,
      result = result,
      conf_level = attr(x, "conf_level") %||% 0.95,
      conf_type = attr(x, "conf_type") %||% "perc"
    )
  }

  for (a in names(attributes(x))) {
    if (!a %in% c("names", "row.names", "class")) {
      attr(result, a) <- attr(x, a)
    }
  }
  if (!is.null(interpolated_draws)) {
    attr(result, draw_attr) <- interpolated_draws
  }
  attr(result, "time_map") <- time_map
  attr(result, "origin_time") <- if (use_origin) origin_time else NULL
  class(result) <- c("markov_interpolated_sops", class(x))
  result
}
