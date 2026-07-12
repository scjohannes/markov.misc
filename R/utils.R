# General package utilities.

# Helper for list access
`%||%` <- function(a, b) if (!is.null(a)) a else b

validate_conf_level <- function(conf_level, arg = "conf_level") {
  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    stop("`", arg, "` must be a single finite number between 0 and 1.")
  }

  conf_level
}

# Bind data frames row-wise while filling columns missing from individual inputs.
# This scoped helper is not a drop-in replacement for dplyr::bind_rows().
bind_rows_fill <- function(x) {
  x <- Filter(Negate(is.null), x)

  if (length(x) == 0) {
    return(data.frame())
  }

  x <- lapply(x, as.data.frame, stringsAsFactors = FALSE, optional = TRUE)
  all_cols <- unique(unlist(lapply(x, names), use.names = FALSE))

  x <- lapply(x, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    df[, all_cols, drop = FALSE]
  })

  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

# Left join while preserving left-hand order and repeated-key expansion.
# This scoped helper is not a drop-in replacement for dplyr::left_join().
left_join_preserve_order <- function(x, y, by) {
  y_cols <- setdiff(names(y), by)
  key_frame <- function(df) {
    do.call(paste, c(df[, by, drop = FALSE], sep = "\r"))
  }
  x_key <- key_frame(x)
  y_key <- key_frame(y)
  y_split <- split(seq_len(nrow(y)), y_key, drop = TRUE)
  y_rows <- lapply(x_key, function(key) y_split[[key]] %||% NA_integer_)
  counts <- lengths(y_rows)
  left_rows <- rep.int(seq_len(nrow(x)), counts)
  right_rows <- unlist(y_rows, use.names = FALSE)
  left <- x[left_rows, , drop = FALSE]
  right <- y[right_rows, y_cols, drop = FALSE]
  rownames(left) <- NULL
  rownames(right) <- NULL
  out <- cbind(left, right)
  rownames(out) <- NULL
  out
}

matrix_to_long <- function(
  mat,
  id_name = "id",
  time_name = "time",
  value_name = "value"
) {
  ids <- rownames(mat) %||% seq_len(nrow(mat))
  times <- colnames(mat) %||% seq_len(ncol(mat))

  out <- data.frame(
    rep(ids, each = ncol(mat)),
    rep(times, times = nrow(mat)),
    as.vector(t(mat)),
    check.names = FALSE
  )
  names(out) <- c(id_name, time_name, value_name)
  out
}

named_list_to_wide <- function(x, id = seq_along(x), id_name = "boot_id") {
  value_names <- unique(unlist(lapply(x, names), use.names = FALSE))
  out <- data.frame(id, check.names = FALSE)
  names(out) <- id_name

  for (name in value_names) {
    out[[name]] <- vapply(
      x,
      function(item) {
        if (is.null(item) || !name %in% names(item)) {
          return(NA_real_)
        }
        as.numeric(item[[name]][1])
      },
      numeric(1)
    )
  }

  out
}

pivot_state_columns_long <- function(
  data,
  values_to = "probability",
  names_to = "state",
  names_prefix = "state_"
) {
  state_cols <- grep(paste0("^", names_prefix), names(data), value = TRUE)
  id_cols <- setdiff(names(data), state_cols)

  if (length(state_cols) == 0) {
    stop("No state columns found.")
  }

  row_index <- rep(seq_len(nrow(data)), each = length(state_cols))
  out <- data[row_index, id_cols, drop = FALSE]
  rownames(out) <- NULL
  out[[names_to]] <- sub(
    paste0("^", names_prefix),
    "",
    rep(state_cols, times = nrow(data))
  )
  out[[values_to]] <- as.vector(t(as.matrix(data[, state_cols, drop = FALSE])))
  out
}

reorder_columns <- function(data, first) {
  first <- intersect(first, names(data))
  data[, c(first, setdiff(names(data), first)), drop = FALSE]
}

arrow_equal_expr <- function(field, value) {
  arrow::Expression$op(
    "==",
    arrow::Expression$field_ref(field),
    arrow::Expression$scalar(value)
  )
}

arrow_in_expr <- function(field, values) {
  values <- unique(values)

  if (length(values) == 0) {
    return(arrow_equal_expr(field, NA))
  }

  exprs <- lapply(values, function(value) arrow_equal_expr(field, value))
  Reduce(function(x, y) arrow::Expression$op("|", x, y), exprs)
}

arrow_collect_dataset <- function(dataset, filter = NULL, cols = NULL) {
  scanner <- dataset$NewScan()

  if (!is.null(filter)) {
    scanner$Filter(filter)
  }
  if (!is.null(cols)) {
    scanner$Project(cols)
  }

  as.data.frame(scanner$Finish()$ToTable())
}

terms_has_offset <- function(terms) {
  offset_terms <- attr(terms, "offset")
  !is.null(offset_terms) && length(offset_terms) > 0
}

model_uses_offset <- function(model) {
  model_chk <- if (inherits(model, "robcov_vglm")) {
    model$vglm_fit
  } else {
    model
  }

  terms_obj <- tryCatch(stats::terms(model_chk), error = function(e) NULL)
  if (!is.null(terms_obj) && terms_has_offset(terms_obj)) {
    return(TRUE)
  }

  offset_obj <- tryCatch(
    {
      if (inherits(model_chk, "vglm")) {
        methods::slot(model_chk, "offset")
      } else {
        model_chk$offset
      }
    },
    error = function(e) NULL
  )
  offset_is_zero <- FALSE
  if (!is.null(offset_obj) && length(offset_obj) > 0) {
    offset_num <- suppressWarnings(as.numeric(offset_obj))
    if (any(!is.na(offset_num) & offset_num != 0)) {
      return(TRUE)
    }
    offset_is_zero <- all(!is.na(offset_num) & offset_num == 0)
  }

  call_obj <- tryCatch(
    {
      if (inherits(model_chk, "vglm")) {
        methods::slot(model_chk, "call")
      } else {
        model_chk$call
      }
    },
    error = function(e) NULL
  )
  if (!is.null(call_obj)) {
    call_args <- as.list(call_obj)
    if (
      "offset" %in%
        names(call_args) &&
        !is.null(call_args[["offset"]]) &&
        !offset_is_zero
    ) {
      return(TRUE)
    }
  }

  FALSE
}

stop_unsupported_offset <- function() {
  stop(
    "Model offsets are not supported by markov.misc Markov SOP workflows. ",
    "Please remove offset() terms or the offset argument before fitting the model.",
    call. = FALSE
  )
}
