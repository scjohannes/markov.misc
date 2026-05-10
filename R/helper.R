# Helper for list access
`%||%` <- function(a, b) if (!is.null(a)) a else b

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

  pieces <- vector("list", nrow(x))
  for (i in seq_len(nrow(x))) {
    key <- x_key[i]
    y_rows <- y_split[[key]]

    if (is.null(y_rows)) {
      y_part <- y[NA_integer_, y_cols, drop = FALSE]
    } else {
      y_part <- y[y_rows, y_cols, drop = FALSE]
    }

    x_part <- x[rep(i, nrow(y_part)), , drop = FALSE]
    rownames(x_part) <- NULL
    rownames(y_part) <- NULL
    pieces[[i]] <- cbind(x_part, y_part)
  }

  bind_rows_fill(pieces)
}

matrix_to_long <- function(mat, id_name = "id", time_name = "time",
                           value_name = "value") {
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
    out[[name]] <- vapply(x, function(item) {
      if (is.null(item) || !name %in% names(item)) {
        return(NA_real_)
      }
      as.numeric(item[[name]][1])
    }, numeric(1))
  }

  out
}

pivot_state_columns_long <- function(data, values_to = "probability",
                                     names_to = "state",
                                     names_prefix = "state_") {
  state_cols <- grep(paste0("^", names_prefix), names(data), value = TRUE)
  id_cols <- setdiff(names(data), state_cols)

  if (length(state_cols) == 0) {
    stop("No state columns found.")
  }

  row_index <- rep(seq_len(nrow(data)), each = length(state_cols))
  out <- data[row_index, id_cols, drop = FALSE]
  rownames(out) <- NULL
  out[[names_to]] <- sub(paste0("^", names_prefix), "", rep(state_cols, times = nrow(data)))
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
