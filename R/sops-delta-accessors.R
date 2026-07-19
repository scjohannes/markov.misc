# Analytical SOP result accessors.

#' Analytical Delta-Method Memory Limit
#'
#' `markov.misc.delta_max_bytes` limits the size of guarded allocations used by
#' analytical delta-method SOP inference. The default is `256 * 1024^2` bytes
#' (256 MiB, or 268,435,456 bytes).
#'
#' The limit is checked for the SOP Jacobian and recursion workspace, retained
#' low-rank analytical state, standard-error workspaces, and covariance or
#' Jacobian blocks requested through [stats::vcov()] or [get_jacobian()]. Sizes
#' are measured in bytes. The checks estimate the principal numeric allocations
#' and are not a guarantee of total process memory use.
#'
#' Increase the option only after confirming that the reported allocation and
#' additional R overhead fit comfortably in available memory. Restore the prior
#' option after a temporary override.
#'
#' @name markov.misc.delta_max_bytes
#' @seealso [inferences()], [get_jacobian()], [stats::vcov()]
#' @examples
#' old_options <- options(
#'   markov.misc.delta_max_bytes = 512 * 1024^2 # 512 MiB
#' )
#' options(old_options)
NULL

delta_row_key_frame <- function(x) {
  key_cols <- setdiff(names(x), inference_columns())
  out <- as.data.frame(x)[, key_cols, drop = FALSE]
  rownames(out) <- NULL
  out
}

delta_row_key <- function(keys) {
  if (ncol(keys) == 0L) {
    return(as.character(seq_len(nrow(keys))))
  }
  do.call(
    paste,
    c(
      lapply(keys, function(value) {
        value <- as.character(value)
        value[is.na(value)] <- "<NA>"
        value
      }),
      sep = "\r"
    )
  )
}

delta_storage_bytes <- function(analytical) {
  component_bytes <- c(
    jacobian = if (is.null(analytical$jacobian)) {
      0
    } else {
      length(analytical$jacobian) * 8
    },
    average_jacobian = if (is.null(analytical$average_jacobian)) {
      0
    } else {
      length(analytical$average_jacobian) * 8
    },
    coefficient_vcov = if (is.null(analytical$coefficient_vcov)) {
      0
    } else {
      length(analytical$coefficient_vcov) * 8
    },
    influence = if (is.null(analytical$influence)) {
      0
    } else {
      length(analytical$influence) * 8
    }
  )
  list(
    components = component_bytes,
    stored = unname(sum(component_bytes)),
    limit = delta_max_bytes()
  )
}

delta_attach_byte_accounting <- function(analytical) {
  analytical$bytes <- delta_storage_bytes(analytical)
  delta_assert_bytes(
    analytical$bytes$stored,
    "Stored analytical inference state"
  )
  analytical
}

delta_analytical <- function(x) {
  if (
    !inherits(
      x,
      c("markov_sops", "markov_avg_sops", "markov_avg_comparisons")
    )
  ) {
    stop(
      "Analytical accessors require a SOP or average-comparison object.",
      call. = FALSE
    )
  }
  analytical <- attr(x, "analytical", exact = TRUE)
  if (is.null(analytical)) {
    stop(
      "No analytical inference state found. Run `inferences(method = ",
      "\"delta\")` first.",
      call. = FALSE
    )
  }
  analytical
}

delta_resolve_rows <- function(analytical, rows = NULL) {
  n <- length(analytical$row_key)
  if (is.null(rows)) {
    return(seq_len(n))
  }
  if (is.logical(rows)) {
    if (length(rows) != n || anyNA(rows)) {
      stop("Logical `rows` must be non-missing and match the result length.")
    }
    return(which(rows))
  }
  if (is.character(rows)) {
    index <- match(rows, analytical$row_key)
    if (anyNA(index)) {
      stop("Some requested analytical row keys were not found.")
    }
    return(index)
  }
  if (is.numeric(rows)) {
    if (
      anyNA(rows) ||
        any(!is.finite(rows)) ||
        any(rows != as.integer(rows)) ||
        any(rows < 1L) ||
        any(rows > n)
    ) {
      stop("Numeric `rows` must contain valid positive row indices.")
    }
    return(as.integer(rows))
  }
  stop("`rows` must be NULL, numeric, logical, or analytical row keys.")
}

delta_jacobian <- function(analytical) {
  jacobian <- analytical$jacobian %||% analytical$average_jacobian
  if (is.null(jacobian)) {
    stop("The analytical result does not contain an average Jacobian.")
  }
  jacobian
}

#' Extract an Analytical SOP Jacobian
#'
#' Extracts the raw-coefficient Jacobian retained by analytical delta-method
#' inference. Superpopulation-target results return the average Jacobian used in
#' the stacked influence function.
#'
#' @param x A SOP or average-comparison object returned by
#'   `inferences(method = "delta")`.
#' @param rows Optional numeric or logical row selection, or analytical row-key
#'   values. `NULL` returns every result row.
#'
#' @return A numeric matrix with one row per selected result cell and one column
#'   per fitted model coefficient.
#' @export
get_jacobian <- function(x, rows = NULL) {
  analytical <- delta_analytical(x)
  index <- delta_resolve_rows(analytical, rows)
  jacobian <- delta_jacobian(analytical)
  required <- as.double(length(index)) * as.double(ncol(jacobian)) * 8
  delta_assert_bytes(required, "The requested Jacobian block")
  out <- jacobian[index, , drop = FALSE]
  rownames(out) <- analytical$row_key[index]
  out
}

delta_vcov <- function(object, rows = NULL) {
  analytical <- delta_analytical(object)
  index <- delta_resolve_rows(analytical, rows)
  m <- length(index)

  if (identical(analytical$representation, "influence")) {
    n <- nrow(analytical$influence)
    if (n < 2L) {
      stop(
        "Superpopulation analytical covariance requires at least two patients."
      )
    }
    required <- (as.double(n) * m + as.double(m) * m) * 8
    delta_assert_bytes(required, "The requested analytical covariance block")
    influence <- analytical$influence[, index, drop = FALSE]
    out <- stats::cov(influence) / n
  } else {
    q <- ncol(analytical$jacobian)
    required <- (as.double(m) * q + as.double(m) * m) * 8
    delta_assert_bytes(required, "The requested analytical covariance block")
    jacobian <- analytical$jacobian[index, , drop = FALSE]
    coefficient_vcov <- analytical$coefficient_vcov
    out <- jacobian %*% coefficient_vcov %*% t(jacobian)
    out <- (out + t(out)) / 2
  }

  dimnames(out) <- list(
    analytical$row_key[index],
    analytical$row_key[index]
  )
  out
}

#' Analytical Covariance for Individual SOPs
#'
#' Materializes a selected covariance block from the low-rank analytical state
#' retained by `inferences(method = "delta")`.
#'
#' @param object A delta-inferred `markov_sops` object.
#' @param rows Optional result-row selection. See [get_jacobian()].
#' @param ... Reserved for the `stats::vcov()` generic.
#'
#' @return A symmetric covariance matrix for the selected result cells.
#' @export
vcov.markov_sops <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}

#' @rdname vcov.markov_sops
#' @export
vcov.markov_avg_sops <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}

#' @rdname vcov.markov_sops
#' @export
vcov.markov_avg_comparisons <- function(object, rows = NULL, ...) {
  delta_vcov(object, rows = rows)
}
