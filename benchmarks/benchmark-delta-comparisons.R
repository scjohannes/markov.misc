# ACI-14: operator storage and propagation, excluding model fitting/recursion.
# Run from the repository root with Rscript benchmarks/benchmark-delta-comparisons.R.
devtools::load_all()

median_seconds <- function(code) {
  code()
  stats::median(replicate(3, unname(system.time(code())[["elapsed"]])))
}

benchmark_delta_comparisons <- function(visits, estimand) {
  avg <- expand.grid(
    time = seq_len(visits),
    state = as.character(1:4),
    tx = 0:2
  )
  avg$estimate <- rep(
    c(0.1, 0.2, 0.3, 0.4, 0.15, 0.15, 0.25, 0.45, 0.2, 0.1, 0.2, 0.5),
    each = visits
  )
  attr(avg, "y_levels") <- as.character(1:4)
  object <- expand.grid(
    time = if (estimand == "sop") seq_len(visits) else 1,
    state_set = as.character(1:4),
    comparison_level = 1:2
  )
  object$reference_level <- 0
  object$estimate <- rep(
    c(0.05, -0.05, -0.05, 0.05, 0.1, -0.1, -0.1, 0.1),
    each = if (estimand == "sop") visits else 1
  )
  args <- list(
    estimand = estimand,
    state_sets = as.list(stats::setNames(as.character(1:4), 1:4))
  )
  if (estimand == "time_in_state") {
    args$time_map <- stats::setNames(seq_len(visits), seq_len(visits))
    args$baseline_time <- 0
    args$target_times <- seq(0, visits, by = 0.5)
    object$estimate <- object$estimate * (visits - 0.5)
  }
  operator <- markov.misc:::delta_comparison_operator(
    object,
    avg,
    args,
    list(variables = list(tx = 0:2))
  )
  # Dense test reference; production must never expand the operator this way.
  dense <- matrix(0, nrow(object), nrow(avg))
  for (i in seq_along(operator)) {
    dense[i, operator[[i]]$index] <- operator[[i]]$weight
  }
  set.seed(14)
  jacobian <- matrix(stats::rnorm(nrow(avg) * 12), ncol = 12)
  influence <- matrix(stats::rnorm(nrow(avg) * 100), nrow = 100)
  compact_j <- function() {
    markov.misc:::delta_apply_comparison_operator(operator, jacobian)
  }
  compact_i <- function() {
    markov.misc:::delta_apply_comparison_operator(operator, influence, TRUE)
  }
  stopifnot(
    isTRUE(all.equal(compact_j(), dense %*% jacobian, tolerance = 1e-12)),
    isTRUE(all.equal(compact_i(), influence %*% t(dense), tolerance = 1e-12))
  )
  data.frame(
    visits = visits,
    estimand = estimand,
    results = nrow(object),
    sources = nrow(avg),
    dense_MiB = as.numeric(object.size(dense)) / 1024^2,
    compact_MiB = as.numeric(object.size(operator)) / 1024^2,
    dense_j_seconds = median_seconds(function() dense %*% jacobian),
    compact_j_seconds = median_seconds(compact_j),
    dense_i_seconds = median_seconds(function() influence %*% t(dense)),
    compact_i_seconds = median_seconds(compact_i)
  )
}

grid <- expand.grid(
  visits = c(20, 100, 250),
  estimand = c("sop", "time_in_state")
)
results <- lapply(seq_len(nrow(grid)), function(i) {
  benchmark_delta_comparisons(grid$visits[i], as.character(grid$estimand[i]))
})
print(do.call(rbind, results), row.names = FALSE, digits = 4)
