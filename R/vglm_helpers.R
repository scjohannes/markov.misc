#' Fit a VGAM Markov Model
#'
#' `vglm.markov()` is the recommended VGAM fitting entrypoint for
#' `markov.misc` SOP workflows. It follows `VGAM::vglm()` closely and returns a
#' normal S4 `vglm` object, but marks the fit so the SOP prediction code can use
#' package-native Markov prediction. For inline restricted cubic spline terms
#' such as `rcs(time, 4)` or `rms::rcs(time, 4)`, the internal VGAM assignment
#' metadata is split by generated spline column. This makes it possible to give
#' individual spline basis columns separate constraint matrices for partial
#' proportional odds models.
#'
#' Non-spline terms, factor terms, and workflows that use explicit precomputed
#' basis columns keep ordinary VGAM assignment behavior. This function does not
#' compute robust covariance matrices; wrap the returned fit with
#' [robcov_vglm()] when robust covariance is needed.
#'
#' @inheritParams VGAM::vglm
#' @param family A VGAM family object, e.g. `VGAM::cumulative()`.
#' @param constraints Optional VGAM constraints list. For inline `rcs()` terms,
#'   names should match the column-level constraint names in a full proportional
#'   odds fit returned by `vglm.markov()`.
#'
#' @return A fitted S4 `vglm` object with an internal Markov marker attribute.
#'
#' @examples
#' \dontrun{
#' library(VGAM)
#' library(rms)
#'
#' fit_po <- vglm.markov(
#'   ordered(y) ~ rcs(time, 4) + tx + age + yprev,
#'   family = cumulative(reverse = TRUE, parallel = TRUE),
#'   data = data
#' )
#'
#' constraints <- fit_po@constraints
#' linear_time_term <- grep("^rcs\\(time, 4\\).*time$", names(constraints),
#'   value = TRUE
#' )[[1]]
#' n_thresholds <- nrow(constraints[[linear_time_term]])
#' constraints[[linear_time_term]] <- cbind(
#'   PO_effect = rep(1, n_thresholds),
#'   linear_deviation = seq_len(n_thresholds) - 1
#' )
#'
#' fit_ppo <- vglm.markov(
#'   ordered(y) ~ rcs(time, 4) + tx + age + yprev,
#'   family = cumulative(reverse = TRUE, parallel = FALSE),
#'   data = data,
#'   constraints = constraints
#' )
#' }
#'
#' @seealso [robcov_vglm()], [avg_sops()], [sops()]
#' @export
vglm.markov <- function(
  formula,
  family = stop("argument 'family' needs to be assigned"),
  data = list(),
  weights = NULL,
  subset = NULL,
  na.action,
  etastart = NULL,
  mustart = NULL,
  coefstart = NULL,
  control = VGAM::vglm.control(...),
  offset = NULL,
  model = FALSE,
  x.arg = TRUE,
  y.arg = TRUE,
  contrasts = NULL,
  constraints = NULL,
  extra = list(),
  form2 = NULL,
  qr.arg = TRUE,
  smart = TRUE,
  ...
) {
  dataname <- as.character(substitute(data))
  function.name <- "vglm"
  ocall <- match.call()

  setup_smart <- getFromNamespace("setup.smart", "VGAM")
  get_smart_prediction <- getFromNamespace("get.smart.prediction", "VGAM")
  wrapup_smart <- getFromNamespace("wrapup.smart", "VGAM")
  shadowvglm <- getFromNamespace("shadowvglm", "VGAM")
  eval_vcontrol <- make_vgam_vcontrol_eval()
  get_xlevels <- getFromNamespace(".getXlevels", "stats")

  if (smart) {
    setup_smart("write")
  }

  if (missing(data)) {
    data <- environment(formula)
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c(
      "formula", "data", "subset", "weights", "na.action",
      "etastart", "mustart", "offset"
    ),
    names(mf),
    0
  )
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  xlev <- get_xlevels(mt, mf)
  y <- stats::model.response(mf, "any")
  x <- if (!stats::is.empty.model(mt)) {
    stats::model.matrix(mt, mf, contrasts)
  } else {
    matrix(, NROW(y), 0)
  }

  split_assign <- split_rcs_assign(x, mt)
  attr(x, "assign") <- split_assign$assign
  attr(x, "orig.assign.lm") <- split_assign$orig_assign

  if (!is.null(form2)) {
    if (!is.null(subset)) {
      stop(
        "argument 'subset' cannot be used when argument 'form2' is used"
      )
    }
    retlist <- shadowvglm(
      formula = form2,
      family = family,
      data = data,
      na.action = na.action,
      control = VGAM::vglm.control(...),
      method = "vglm.fit",
      model = model,
      x.arg = x.arg,
      y.arg = y.arg,
      contrasts = contrasts,
      constraints = constraints,
      extra = extra,
      qr.arg = qr.arg
    )
    Ym2 <- retlist$Ym2
    Xm2 <- retlist$Xm2
    if (length(Ym2) && NROW(Ym2) != NROW(y)) {
      stop("number of rows of 'y' and 'Ym2' are unequal")
    }
    if (length(Xm2) && NROW(Xm2) != NROW(x)) {
      stop("number of rows of 'x' and 'Xm2' are unequal")
    }
  } else {
    Xm2 <- Ym2 <- NULL
  }

  offset <- stats::model.offset(mf)
  if (is.null(offset)) {
    offset <- 0
  }

  w <- stats::model.weights(mf)
  if (!length(w)) {
    w <- rep_len(1, nrow(mf))
  } else if (NCOL(w) == 1 && any(w < 0)) {
    stop("negative weights not allowed")
  }

  if (is.character(family)) {
    family <- get(family, envir = asNamespace("VGAM"), inherits = TRUE)
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "vglmff")) {
    stop("'family = ", family, "' is not a VGAM family function")
  }

  control <- eval_vcontrol(control, family, function.name, ...)
  if (length(methods::slot(family, "first"))) {
    eval(methods::slot(family, "first"))
  }

  fit <- VGAM::vglm.fit(
    x = x,
    y = y,
    w = w,
    offset = offset,
    Xm2 = Xm2,
    Ym2 = Ym2,
    etastart = etastart,
    mustart = mustart,
    coefstart = coefstart,
    family = family,
    control = control,
    constraints = constraints,
    extra = extra,
    qr.arg = qr.arg,
    Terms = mt,
    function.name = function.name,
    ...
  )

  fit$misc$dataname <- dataname
  if (smart) {
    fit$smart.prediction <- get_smart_prediction()
    wrapup_smart()
  }

  answer <- methods::new(
    Class = "vglm",
    assign = attr(x, "assign"),
    call = ocall,
    coefficients = fit$coefficients,
    constraints = fit$constraints,
    criterion = fit$crit.list,
    df.residual = fit$df.residual,
    df.total = fit$df.total,
    dispersion = 1,
    effects = fit$effects,
    family = fit$family,
    misc = fit$misc,
    model = if (model) mf else data.frame(),
    R = fit$R,
    rank = fit$rank,
    residuals = as.matrix(fit$residuals),
    ResSS = fit$ResSS,
    smart.prediction = as.list(fit$smart.prediction),
    terms = list(terms = mt)
  )

  if (!smart) {
    answer@smart.prediction <- list(smart.arg = FALSE)
  }
  if (qr.arg) {
    class(fit$qr) <- "list"
    methods::slot(answer, "qr") <- fit$qr
  }
  if (length(attr(x, "contrasts"))) {
    methods::slot(answer, "contrasts") <- attr(x, "contrasts")
  }
  if (length(fit$fitted.values)) {
    methods::slot(answer, "fitted.values") <- as.matrix(fit$fitted.values)
  }
  methods::slot(answer, "na.action") <- if (length(aaa <- attr(mf, "na.action"))) {
    list(aaa)
  } else {
    list()
  }
  if (length(offset)) {
    methods::slot(answer, "offset") <- as.matrix(offset)
  }
  if (length(fit$weights)) {
    methods::slot(answer, "weights") <- as.matrix(fit$weights)
  }
  if (x.arg) {
    methods::slot(answer, "x") <- fit$x
  }
  if (x.arg && length(Xm2)) {
    methods::slot(answer, "Xm2") <- Xm2
  }
  if (y.arg && length(Ym2)) {
    methods::slot(answer, "Ym2") <- as.matrix(Ym2)
  }
  if (!is.null(form2)) {
    methods::slot(answer, "callXm2") <- retlist$call
    answer@misc$Terms2 <- retlist$Terms2
  }
  answer@misc$formula <- formula
  answer@misc$form2 <- form2
  if (length(xlev)) {
    methods::slot(answer, "xlevels") <- xlev
  }
  if (y.arg) {
    methods::slot(answer, "y") <- as.matrix(fit$y)
  }
  methods::slot(answer, "control") <- fit$control
  methods::slot(answer, "extra") <- if (length(fit$extra)) {
    if (is.list(fit$extra)) {
      fit$extra
    } else {
      warning("'extra' is not a list, therefore placing 'extra' into a list")
      list(fit$extra)
    }
  } else {
    list()
  }
  methods::slot(answer, "iter") <- fit$iter
  methods::slot(answer, "post") <- fit$post
  fit$predictors <- as.matrix(fit$predictors)
  if (length(fit$misc$predictors.names) == ncol(fit$predictors)) {
    dimnames(fit$predictors) <- list(
      dimnames(fit$predictors)[[1]],
      fit$misc$predictors.names
    )
  }
  methods::slot(answer, "predictors") <- fit$predictors
  if (length(fit$prior.weights)) {
    methods::slot(answer, "prior.weights") <- as.matrix(fit$prior.weights)
  }

  attr(answer, "markov_vglm") <- TRUE
  attr(answer, "markov_split_assign") <- split_assign$has_rcs
  answer
}

split_rcs_assign <- function(x, terms) {
  attrassigndefault <- getFromNamespace("attrassigndefault", "VGAM")
  orig_assign <- attr(x, "assign")
  assign <- attrassigndefault(x, terms)

  rcs_terms <- names(assign)[has_inline_rcs(names(assign))]

  for (term in rcs_terms) {
    cols <- assign[[term]]
    pos <- match(term, names(assign))
    assign[term] <- NULL
    split_cols <- stats::setNames(as.list(cols), colnames(x)[cols])
    assign <- append(assign, split_cols, after = pos - 1)
  }

  list(
    assign = assign,
    orig_assign = orig_assign,
    has_rcs = length(rcs_terms) > 0
  )
}

has_inline_rcs <- function(term) {
  grepl("(^|:)\\s*(rms::)?rcs\\(", term)
}

make_vgam_vcontrol_eval <- function() {
  expr <- getFromNamespace("vcontrol.expression", "VGAM")
  eval(
    substitute(
      function(control, family, function.name, ...) {
        eval(EXPR)
        control
      },
      list(EXPR = expr)
    ),
    envir = asNamespace("VGAM")
  )
}
