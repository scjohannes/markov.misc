if (identical(Sys.getenv("MARKOV_MISC_BENCH_INSTALLED"), "true")) {
  library(markov.misc)
} else {
  devtools::load_all()
}

timings <- function(code, iterations = 7L) {
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

benchmark_second_order_update <- function() {
  n <- 300L
  states <- 8L
  older <- rep(seq_len(states - 1L), times = states - 1L)
  current <- rep(seq_len(states - 1L), each = states - 1L)
  pairs <- length(older)
  previous <- array(stats::runif(n * states * states), c(n, states, states))
  transition <- matrix(stats::runif(n * pairs * states), n * pairs, states)

  reference <- function() {
    for (visit in seq_len(30L)) {
      out <- array(0, c(n, states, states))
      for (pair in seq_len(pairs)) {
        rows <- (pair - 1L) * n + seq_len(n)
        mass <- previous[, older[pair], current[pair]]
        for (target in seq_len(states)) {
          out[, current[pair], target] <- out[, current[pair], target] +
            transition[rows, target] * mass
        }
      }
      for (h in seq_len(states)) {
        out[, states, states] <- out[, states, states] + previous[, h, states]
      }
    }
    out
  }
  candidate <- function() {
    for (visit in seq_len(30L)) {
      out <- markov.misc:::markov_update_second_order_native(
        previous,
        transition,
        older,
        current,
        states
      )
    }
    out
  }
  stopifnot(isTRUE(all.equal(reference(), candidate(), tolerance = 1e-12)))
  list(reference = timings(reference), candidate = timings(candidate))
}

benchmark_blrm_conversion <- function() {
  draws <- 20L
  observations <- 300L * 7L
  thresholds <- 7L
  base_eta <- matrix(stats::rnorm(draws * observations), draws)
  intercepts <- matrix(stats::rnorm(draws * thresholds), draws)
  threshold_eta <- matrix(stats::rnorm(draws * observations), draws)
  threshold_scale <- seq(-0.4, 0.4, length.out = thresholds)

  reference <- function() {
    cumulative <- array(NA_real_, c(draws, observations, thresholds))
    for (threshold in seq_len(thresholds)) {
      cumulative[,, threshold] <- stats::plogis(
        sweep(base_eta, 1L, intercepts[, threshold], "+") +
          threshold_scale[threshold] * threshold_eta
      )
    }
    out <- array(NA_real_, c(draws, observations, thresholds + 1L))
    out[,, 1L] <- 1 - cumulative[,, 1L]
    out[,, 2L:thresholds] <-
      cumulative[,, seq_len(thresholds - 1L), drop = FALSE] -
      cumulative[,, 2L:thresholds, drop = FALSE]
    out[,, thresholds + 1L] <- cumulative[,, thresholds]
    out[out < 0] <- 0
    markov.misc:::normalize_probability_array(out)
  }
  candidate <- function() {
    markov.misc:::blrm_probabilities_native(
      base_eta,
      intercepts,
      threshold_eta,
      threshold_scale
    )
  }
  stopifnot(isTRUE(all.equal(reference(), candidate(), tolerance = 1e-12)))
  list(reference = timings(reference), candidate = timings(candidate))
}

set.seed(731)
components <- list(
  second_order_update = benchmark_second_order_update(),
  blrm_probability_conversion = benchmark_blrm_conversion()
)
results <- do.call(
  rbind,
  lapply(names(components), function(component) {
    value <- components[[component]]
    data.frame(
      component = component,
      reference_median = value$reference[["median"]],
      candidate_median = value$candidate[["median"]],
      speedup = value$reference[["median"]] / value$candidate[["median"]],
      candidate_p90 = value$candidate[["p90"]]
    )
  })
)
print(results)
