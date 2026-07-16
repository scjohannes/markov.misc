if (identical(Sys.getenv("MARKOV_MISC_BENCH_INSTALLED"), "true")) {
  library(markov.misc)
} else {
  devtools::load_all()
}

timings <- function(
  code,
  iterations = as.integer(Sys.getenv("MARKOV_MISC_BENCH_REPS", "7"))
) {
  code()
  elapsed <- numeric(iterations)
  for (iteration in seq_len(iterations)) {
    gc()
    elapsed[iteration] <- unname(system.time(code())[["elapsed"]])
  }
  c(
    median = stats::median(elapsed),
    p90 = unname(stats::quantile(elapsed, 0.9, names = FALSE))
  )
}

make_second_order_data <- function(
  patients = 300L,
  visits = 30L,
  states = 8L,
  seed = 741L
) {
  set.seed(seed)
  data <- expand.grid(
    id = seq_len(patients),
    time = seq_len(visits),
    KEEP.OUT.ATTRS = FALSE
  )
  data <- data[order(data$id, data$time), , drop = FALSE]
  data$tx <- rep(stats::rbinom(patients, 1L, 0.5), each = visits)
  data$y <- integer(nrow(data))
  for (id in seq_len(patients)) {
    rows <- (id - 1L) * visits + seq_len(visits)
    state <- sample.int(states - 1L, 1L)
    for (visit in seq_len(visits)) {
      state <- max(1L, min(states, state + sample(-1:1, 1L)))
      data$y[rows[visit]] <- state
    }
  }
  data$yprev <- ave(data$y, data$id, FUN = function(value) {
    c(value[1L], value[-length(value)])
  })
  data$ypprev <- ave(data$yprev, data$id, FUN = function(value) {
    c(value[1L], value[-length(value)])
  })
  data$y <- factor(data$y, levels = seq_len(states), ordered = TRUE)
  data$yprev <- factor(data$yprev, levels = seq_len(states))
  data$ypprev <- factor(data$ypprev, levels = seq_len(states))
  data
}

data <- make_second_order_data()
fit_data <- data[data$time > 2L, , drop = FALSE]
model <- rms::orm(
  y ~ tx + time + yprev + ypprev,
  data = fit_data,
  x = TRUE,
  y = TRUE
)
baseline <- data[data$time == 1L, , drop = FALSE]

namespace <- asNamespace("markov.misc")
set_internal <- function(name, value) {
  unlockBinding(name, namespace)
  assign(name, value, envir = namespace)
  lockBinding(name, namespace)
}

second_order_native <- get("markov_update_second_order_native", namespace)
second_order_reference <- function(
  previous,
  transition,
  older,
  current,
  absorb
) {
  dims <- dim(previous)
  out <- array(0, dims)
  observations <- dims[1L]
  states <- dims[2L]
  for (pair in seq_along(older)) {
    rows <- (pair - 1L) * observations + seq_len(observations)
    mass <- previous[, older[pair], current[pair]]
    for (target in seq_len(states)) {
      out[, current[pair], target] <- out[, current[pair], target] +
        transition[rows, target] * mass
    }
  }
  for (state in absorb) {
    out[, state, state] <- out[, state, state] +
      apply(previous[,, state, drop = FALSE], 1L, sum)
  }
  out
}

second_order <- timings(function() {
  soprob_markov(
    model,
    newdata = baseline,
    times = seq_len(30L),
    y_levels = seq_len(8L),
    p2_var = "ypprev"
  )
})
set_internal("markov_update_second_order_native", second_order_reference)
second_order_reference_time <- timings(function() {
  soprob_markov(
    model,
    newdata = baseline,
    times = seq_len(30L),
    y_levels = seq_len(8L),
    p2_var = "ypprev"
  )
})
set_internal("markov_update_second_order_native", second_order_native)

second_order_point <- avg_sops(
  model,
  newdata = baseline,
  variables = list(tx = c(0, 1)),
  times = seq_len(30L),
  y_levels = seq_len(8L),
  p2_var = "ypprev"
)
second_order_inference <- timings(function() {
  inferences(
    second_order_point,
    method = "mvn",
    n_draws = 5L,
    workers = 1L,
    seed = 744L,
    return_draws = FALSE
  )
})
compile_plan <- get("compile_sop_execution_plan", namespace)
set_internal(
  "compile_sop_execution_plan",
  function(...) stop("benchmark reference fallback")
)
second_order_inference_reference <- timings(function() {
  suppressWarnings(inferences(
    second_order_point,
    method = "mvn",
    n_draws = 5L,
    workers = 1L,
    seed = 744L,
    return_draws = FALSE
  ))
})
set_internal("compile_sop_execution_plan", compile_plan)

blrm <- model
class(blrm) <- unique(c("blrm", class(blrm)))
set.seed(742)
blrm$draws <- matrix(
  rep(stats::coef(model), each = 20L) +
    stats::rnorm(20L * length(stats::coef(model)), sd = 0.02),
  nrow = 20L,
  dimnames = list(NULL, names(stats::coef(model)))
)
blrm$pppo <- 0L
blrm$p <- length(stats::coef(model)) - model$non.slopes
blrm$tauInfo <- data.frame(name = character())
blrm$ylevels <- seq_len(8L)
blrm$yname <- "y"

# The fixture uses synthetic posterior draws attached to an ORM fit. Recreate
# the fitted numeric design directly so timing excludes MCMC fitting while
# retaining the same posterior propagation dimensions as a BLRM fit.
unlockBinding("blrm_design_matrix", namespace)
assign(
  "blrm_design_matrix",
  function(model, newdata, second = FALSE) {
    if (second) {
      return(matrix(numeric(), nrow = nrow(newdata), ncol = 0L))
    }
    out <- stats::model.matrix(
      ~ tx + time + yprev + ypprev,
      newdata
    )[, -1L, drop = FALSE]
    colnames(out) <- sub("^(yprev|ypprev)([0-9]+)$", "\\1=\\2", colnames(out))
    out
  },
  envir = namespace
)
lockBinding("blrm_design_matrix", namespace)

blrm_native <- get("blrm_probabilities_native", namespace)
blrm_reference <- function(
  base_eta,
  intercepts,
  threshold_eta = numeric(),
  threshold_scale = numeric()
) {
  thresholds <- ncol(intercepts)
  cumulative <- array(
    NA_real_,
    c(nrow(base_eta), ncol(base_eta), thresholds)
  )
  partial <- length(threshold_eta) > 0L
  for (threshold in seq_len(thresholds)) {
    eta <- sweep(base_eta, 1L, intercepts[, threshold], "+")
    if (partial) {
      eta <- eta + threshold_scale[threshold] * threshold_eta
    }
    cumulative[,, threshold] <- stats::plogis(eta)
  }
  out <- array(
    NA_real_,
    c(nrow(base_eta), ncol(base_eta), thresholds + 1L)
  )
  out[,, 1L] <- 1 - cumulative[,, 1L]
  out[,, 2L:thresholds] <-
    cumulative[,, seq_len(thresholds - 1L), drop = FALSE] -
    cumulative[,, 2L:thresholds, drop = FALSE]
  out[,, thresholds + 1L] <- cumulative[,, thresholds]
  out[out < 0] <- 0
  markov.misc:::normalize_probability_array(out)
}

blrm_avg <- timings(function() {
  avg_sops(
    blrm,
    newdata = baseline,
    variables = list(tx = c(0, 1)),
    times = seq_len(30L),
    y_levels = seq_len(8L),
    n_draws = 20L,
    seed = 743L
  )
})
set_internal("blrm_probabilities_native", blrm_reference)
blrm_reference_time <- timings(function() {
  avg_sops(
    blrm,
    newdata = baseline,
    variables = list(tx = c(0, 1)),
    times = seq_len(30L),
    y_levels = seq_len(8L),
    n_draws = 20L,
    seed = 743L
  )
})
set_internal("blrm_probabilities_native", blrm_native)

results <- data.frame(
  workflow = c(
    "second-order ORM",
    "second-order ORM + MVN",
    "BLRM posterior avg_sops"
  ),
  reference_median = c(
    second_order_reference_time[["median"]],
    second_order_inference_reference[["median"]],
    blrm_reference_time[["median"]]
  ),
  median = c(
    second_order[["median"]],
    second_order_inference[["median"]],
    blrm_avg[["median"]]
  ),
  speedup = c(
    second_order_reference_time[["median"]] / second_order[["median"]],
    second_order_inference_reference[["median"]] /
      second_order_inference[["median"]],
    blrm_reference_time[["median"]] / blrm_avg[["median"]]
  ),
  p90 = c(
    second_order[["p90"]],
    second_order_inference[["p90"]],
    blrm_avg[["p90"]]
  )
)
print(results)

required <- c(
  "second-order ORM + MVN" = 2,
  "BLRM posterior avg_sops" = 2.5
)
observed <- stats::setNames(results$speedup, results$workflow)
failed <- names(required)[observed[names(required)] < required]
if (length(failed) > 0L) {
  stop(
    "Missing-workflow benchmark gate failed: ",
    paste(failed, collapse = ", ")
  )
}
