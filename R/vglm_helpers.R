#' Fit a VGAM model with splines split for column-level constraints
#'
#' `VGAM::vglm()` groups multi-column inline terms such as `rms::rcs(time, 4)`
#' under one model term when processing constraints. This helper keeps the
#' user-facing formula unchanged, but splits the internal assignment metadata
#' for inline restricted cubic spline columns before calling `VGAM::vglm.fit()`.
#' That makes it possible to apply different constraint matrices to individual
#' spline basis columns while still fitting a model written with inline
#' `rcs()`.
#'
#' @param formula Model formula passed to `VGAM::vglm()`.
#' @param family A VGAM family object, e.g. `VGAM::cumulative()`.
#' @param data Model data.
#' @param constraints Optional VGAM constraints list. Names should match the
#'   column names in a full-PO fit returned by this helper. Factor terms may
#'   still need to be supplied at the formula-term level, as in VGAM itself.
#' @param control Control list from `VGAM::vglm.control()`.
#' @param model,x.arg,y.arg,qr.arg Passed through to `VGAM::vglm()` for the
#'   returned object.
#' @param ... Currently ignored.
#'
#' @return A fitted `vglm` object.
#' @export
vglm_split_rcs <- function(
  formula,
  family,
  data,
  constraints = NULL,
  control = VGAM::vglm.control(),
  model = FALSE,
  x.arg = TRUE,
  y.arg = TRUE,
  qr.arg = TRUE,
  ...
) {
  mf <- stats::model.frame(
    formula = formula,
    data = data,
    drop.unused.levels = TRUE
  )
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "any")
  x <- stats::model.matrix(mt, mf)

  split_assign <- split_rcs_assign(x, mt)
  attr(x, "assign") <- split_assign$assign
  attr(x, "orig.assign.lm") <- split_assign$orig_assign

  if (is.function(family)) {
    family <- family()
  }

  control$criterion <- control$criterion %||% "deviance"
  if (length(control$min.criterion) != 1) {
    control$min.criterion <- TRUE
  }

  fit <- VGAM::vglm.fit(
    x = x,
    y = y,
    w = stats::model.weights(mf) %||% rep_len(1, nrow(mf)),
    offset = stats::model.offset(mf) %||% 0,
    family = family,
    control = control,
    constraints = constraints,
    extra = list(),
    qr.arg = qr.arg,
    Terms = mt,
    function.name = "vglm"
  )

  ans <- VGAM::vglm(
    formula = formula,
    family = family,
    data = data,
    model = model,
    x.arg = x.arg,
    y.arg = y.arg,
    qr.arg = qr.arg
  )

  ans@assign <- attr(x, "assign")
  ans@call <- match.call()
  ans@coefficients <- fit$coefficients
  ans@constraints <- fit$constraints
  ans@criterion <- fit$crit.list
  ans@df.residual <- fit$df.residual
  ans@df.total <- fit$df.total
  ans@effects <- fit$effects
  ans@family <- fit$family
  ans@misc <- fit$misc
  ans@R <- fit$R
  ans@rank <- fit$rank
  ans@residuals <- as.matrix(fit$residuals)
  ans@x <- fit$x
  ans@y <- as.matrix(fit$y)
  ans@control <- fit$control
  ans@extra <- fit$extra
  ans@iter <- fit$iter
  ans@post <- fit$post
  if (length(fit$fitted.values)) {
    ans@fitted.values <- as.matrix(fit$fitted.values)
  }
  ans@predictors <- as.matrix(fit$predictors)
  if (length(fit$prior.weights)) {
    ans@prior.weights <- as.matrix(fit$prior.weights)
  }
  if (length(fit$weights)) {
    ans@weights <- as.matrix(fit$weights)
  }
  if (length(fit$qr)) {
    class(fit$qr) <- "list"
    ans@qr <- fit$qr
  }

  attr(ans, "markov_split_assign") <- TRUE
  ans
}

split_rcs_assign <- function(x, terms) {
  attrassigndefault <- getFromNamespace("attrassigndefault", "VGAM")
  orig_assign <- attr(x, "assign")
  assign <- attrassigndefault(x, terms)

  rcs_terms <- names(assign)[
    startsWith(names(assign), "rcs(") |
      startsWith(names(assign), "rms::rcs(")
  ]

  for (term in rcs_terms) {
    cols <- assign[[term]]
    pos <- match(term, names(assign))
    assign[term] <- NULL
    split_cols <- stats::setNames(as.list(cols), colnames(x)[cols])
    assign <- append(assign, split_cols, after = pos - 1)
  }

  list(assign = assign, orig_assign = orig_assign)
}
