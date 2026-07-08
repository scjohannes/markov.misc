#' Fit a VGAM Markov Model
#'
#' `vglm_markov()` is the recommended VGAM fitting entrypoint for
#' `markov.misc` SOP workflows. It follows `VGAM::vglm()` closely, stores the
#' row-aligned fitting data for downstream SOP inference, and marks the fit so
#' the SOP prediction code can use package-native Markov prediction. For inline
#' restricted cubic spline terms such as `rcs(time, 4)` or `rms::rcs(time, 4)`,
#' the internal VGAM assignment metadata is split by generated spline column.
#' This makes it possible to give individual spline basis columns separate
#' constraint matrices for partial proportional odds models. Inline splines may
#' be used for time or for numeric previous-state effects such as
#' `rcs(yprev, 6)`.
#'
#' `vglm_markov()` also makes the `rms` formula helpers `rcs()` and `%ia%`
#' available while building model frames and prediction design matrices. This
#' lets formulas use `time %ia% yprev` for a linear-time by previous-state
#' interaction without requiring `library(rms)` in the user's session.
#'
#' Non-spline terms, factor terms, and workflows that use explicit precomputed
#' basis columns keep ordinary VGAM assignment behavior. When `id_var` is
#' supplied, the returned object is wrapped with [robcov_vglm()] using that
#' patient/cluster ID. When `id_var` is omitted, a plain S4 `vglm` object is
#' returned; if `data` contains an `id` column, a warning reminds the user that
#' automatic cluster-robust covariance and FWB refits need `id_var`.
#'
#' Offsets are not supported in `markov.misc` SOP workflows, so `vglm_markov()`
#' rejects both `offset()` terms and the `offset` argument.
#'
#' @inheritParams VGAM::vglm
#' @param family A VGAM family object, e.g. `VGAM::cumulative()`.
#' @param id_var Optional character scalar naming the patient or cluster ID
#'   column in `data`. When supplied, [robcov_vglm()] is applied automatically.
#' @param constraints Optional VGAM constraints list. For inline `rcs()` terms,
#'   names should match the column-level constraint names in a full proportional
#'   odds fit returned by `vglm_markov()`.
#'
#' @return A fitted S4 `vglm` object with internal Markov marker attributes, or
#'   a `robcov_vglm` object when `id_var` is supplied.
#'
#' @examples
#' \dontrun{
#' library(VGAM)
#' library(rms)
#'
#' fit_po <- vglm_markov(
#'   ordered(y) ~ rcs(time, 4) + tx + age + yprev,
#'   family = cumulative(reverse = TRUE, parallel = TRUE),
#'   data = data,
#'   id_var = "id"
#' )
#'
#' fit_numeric_prev <- vglm_markov(
#'   ordered(y) ~ rcs(time, 4) + tx + age + rcs(yprev, 6),
#'   family = cumulative(reverse = TRUE, parallel = TRUE),
#'   data = prepare_markov_data(data, factor_previous = FALSE)
#' )
#'
#' fit_linear_time_by_prev <- vglm_markov(
#'   ordered(y) ~ rcs(time, 6) + tx + age + yprev + time %ia% yprev,
#'   family = cumulative(reverse = TRUE, parallel = TRUE),
#'   data = data
#' )
#'
#' # In contrast, rcs(time, 6) %ia% yprev follows rms semantics and interacts
#' # the spline basis columns with yprev when yprev is categorical.
#' fit_spline_basis_by_prev <- vglm_markov(
#'   ordered(y) ~ rcs(time, 6) + tx + age + yprev + rcs(time, 6) %ia% yprev,
#'   family = cumulative(reverse = TRUE, parallel = TRUE),
#'   data = data
#' )
#'
#' constraints <- fit_po$vglm_fit@constraints
#' linear_time_term <- grep("^rcs\\(time, 4\\).*time$", names(constraints),
#'   value = TRUE
#' )[[1]]
#' n_thresholds <- nrow(constraints[[linear_time_term]])
#' constraints[[linear_time_term]] <- cbind(
#'   PO_effect = rep(1, n_thresholds),
#'   linear_deviation = seq_len(n_thresholds) - 1
#' )
#'
#' fit_ppo <- vglm_markov(
#'   ordered(y) ~ rcs(time, 4) + tx + age + yprev,
#'   family = cumulative(reverse = TRUE, parallel = FALSE),
#'   data = data,
#'   constraints = constraints
#' )
#' }
#'
#' @seealso [robcov_vglm()], [avg_sops()], [sops()]
#' @export
vglm_markov <- function(
  formula,
  family = stop("argument 'family' needs to be assigned"),
  data = list(),
  id_var = NULL,
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
  data_was_supplied <- !missing(data)
  original_data <- if (data_was_supplied && is.data.frame(data)) data else NULL
  formula <- add_rms_formula_helpers(formula)
  if (!is.null(form2)) {
    form2 <- add_rms_formula_helpers(form2)
  }

  setup_smart <- utils::getFromNamespace("setup.smart", "VGAM")
  get_smart_prediction <- utils::getFromNamespace(
    "get.smart.prediction",
    "VGAM"
  )
  wrapup_smart <- utils::getFromNamespace("wrapup.smart", "VGAM")
  shadowvglm <- utils::getFromNamespace("shadowvglm", "VGAM")
  eval_vcontrol <- make_vgam_vcontrol_eval()
  get_xlevels <- utils::getFromNamespace(".getXlevels", "stats")

  smart_open <- FALSE
  if (smart) {
    setup_smart("write")
    smart_open <- TRUE
    on.exit(
      {
        if (smart_open) {
          wrapup_smart()
        }
      },
      add = TRUE
    )
  }

  if (missing(data)) {
    data <- environment(formula)
  }

  if (!is.null(id_var)) {
    if (is.null(original_data)) {
      stop("`id_var` requires `data` to be supplied as a data frame.")
    }
    id_var <- validate_markov_id_var(id_var, original_data, "data")
    warn_duplicate_markov_id_time(
      original_data,
      id_var,
      wrapper = "vglm_markov()"
    )
  } else if (!is.null(original_data) && "id" %in% names(original_data)) {
    warn_missing_markov_id_var("vglm_markov()")
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c(
      "formula",
      "data",
      "subset",
      "weights",
      "na.action",
      "etastart",
      "mustart",
      "offset"
    ),
    names(mf),
    0
  )
  mf <- mf[c(1, m)]
  mf$formula <- formula
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  fit_data <- markov_align_model_data(original_data, mf)

  mt <- attr(mf, "terms")
  if (terms_has_offset(mt) || !is.null(stats::model.offset(mf))) {
    stop_unsupported_offset()
  }

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
    smart_open <- FALSE
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
  methods::slot(answer, "na.action") <- if (
    length(aaa <- attr(mf, "na.action"))
  ) {
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
  answer <- markov_attach_model_data(answer, fit_data, id_var)

  if (!is.null(id_var)) {
    robust <- robcov_vglm(answer, cluster = fit_data[[id_var]])
    robust <- markov_attach_model_data(robust, fit_data, id_var)
    return(robust)
  }

  answer
}

add_rms_formula_helpers <- function(formula) {
  if (!inherits(formula, "formula")) {
    return(formula)
  }

  formula_env <- environment(formula)
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }

  helper_env <- new.env(parent = formula_env)

  if (!exists("rcs", envir = formula_env, inherits = TRUE)) {
    assign("rcs", utils::getFromNamespace("rcs", "rms"), envir = helper_env)
  }
  if (!exists("%ia%", envir = formula_env, inherits = TRUE)) {
    assign(
      "%ia%",
      utils::getFromNamespace("%ia%", "rms"),
      envir = helper_env
    )
  }

  environment(formula) <- helper_env
  formula
}

split_rcs_assign <- function(x, terms) {
  attrassigndefault <- utils::getFromNamespace("attrassigndefault", "VGAM")
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
  expr <- utils::getFromNamespace("vcontrol.expression", "VGAM")
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
