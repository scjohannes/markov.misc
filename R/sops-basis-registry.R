# Internal registry for fitted numeric basis transformations used by SOP plans.

rms_basis_registry <- function() {
  list(
    rcs = list(
      id = "rcs",
      helper = "rcs",
      assume = "rcspline",
      assume_code = 4L,
      pattern = "(^|:)\\s*(rms::)?rcs\\(",
      evaluate = evaluate_rcs_basis
    ),
    lsp = list(
      id = "lsp",
      helper = "lsp",
      assume = "lspline",
      assume_code = 3L,
      pattern = "(^|:)\\s*(rms::)?lsp\\(",
      evaluate = evaluate_lsp_basis
    )
  )
}

rms_basis_handler <- function(id) {
  handler <- rms_basis_registry()[[id]]
  if (is.null(handler)) {
    stop("Unknown RMS basis handler: ", id, call. = FALSE)
  }
  handler
}

rms_basis_handler_for_design <- function(assume, assume_code = NULL) {
  handlers <- rms_basis_registry()
  matched <- vapply(
    handlers,
    function(handler) {
      identical(as.character(assume), handler$assume) ||
        (!is.null(assume_code) &&
          identical(as.integer(assume_code), handler$assume_code))
    },
    logical(1)
  )
  if (!any(matched)) {
    return(NULL)
  }
  handlers[[which(matched)[1L]]]
}

rms_basis_handlers_for_term <- function(term) {
  handlers <- rms_basis_registry()
  names(handlers)[vapply(
    handlers,
    function(handler) grepl(handler$pattern, term, perl = TRUE),
    logical(1)
  )]
}

has_registered_rms_basis <- function(term) {
  length(rms_basis_handlers_for_term(term)) > 0L
}

new_compiled_rms_basis <- function(
  handler,
  variable,
  parameters,
  columns,
  nonlinear = NULL
) {
  handler <- rms_basis_handler(handler)
  parameters <- as.numeric(parameters)
  columns <- as.character(columns)
  if (!length(variable) || length(variable) != 1L || !nzchar(variable)) {
    stop("A compiled RMS basis requires one source variable.", call. = FALSE)
  }
  if (!length(parameters) || anyNA(parameters) || any(!is.finite(parameters))) {
    stop(
      "A compiled RMS basis requires finite fitted parameters.",
      call. = FALSE
    )
  }
  expected_columns <- length(parameters) + 1L
  if (handler$id == "rcs") {
    expected_columns <- length(parameters) - 1L
  }
  if (length(columns) != expected_columns) {
    stop(
      "Compiled ",
      handler$id,
      " basis expected ",
      expected_columns,
      " columns, got ",
      length(columns),
      ".",
      call. = FALSE
    )
  }
  if (is.null(nonlinear)) {
    nonlinear <- c(FALSE, rep(TRUE, expected_columns - 1L))
  }
  if (length(nonlinear) != expected_columns || anyNA(nonlinear)) {
    stop(
      "RMS basis nonlinear metadata is not aligned with its columns.",
      call. = FALSE
    )
  }
  structure(
    list(
      handler = handler$id,
      variable = variable,
      parameters = parameters,
      columns = columns,
      nonlinear = as.logical(nonlinear)
    ),
    class = "markov_rms_basis"
  )
}

evaluate_compiled_rms_basis <- function(basis, values) {
  if (!inherits(basis, "markov_rms_basis")) {
    stop("`basis` must be a compiled RMS basis.", call. = FALSE)
  }
  if (!is.numeric(values)) {
    stop("RMS numeric bases require numeric prediction values.", call. = FALSE)
  }
  handler <- rms_basis_handler(basis$handler)
  out <- handler$evaluate(values, basis$parameters)
  if (!is.matrix(out) || ncol(out) != length(basis$columns)) {
    stop("RMS basis handler returned an invalid matrix.", call. = FALSE)
  }
  colnames(out) <- basis$columns
  attr(out, "nonlinear") <- basis$nonlinear
  out
}

evaluate_rcs_basis <- function(values, knots) {
  evaluator <- utils::getFromNamespace("rcspline.eval", "Hmisc")
  as.matrix(evaluator(values, knots = knots, inclx = TRUE))
}

evaluate_lsp_basis <- function(values, knots) {
  out <- matrix(NA_real_, nrow = length(values), ncol = length(knots) + 1L)
  out[, 1L] <- values
  for (j in seq_along(knots)) {
    out[, j + 1L] <- pmax(values - knots[j], 0)
  }
  out
}
