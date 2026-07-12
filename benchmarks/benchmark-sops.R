devtools::load_all()
source("tests/testthat/helper-simulation.R")

timings <- function(code, iterations = 7L) {
  values <- numeric(iterations)
  for (i in seq_len(iterations)) {
    gc()
    values[i] <- unname(system.time(code())["elapsed"])
  }
  c(
    median = stats::median(values),
    p90 = unname(stats::quantile(
      values,
      0.9,
      names = FALSE
    ))
  )
}

results <- list()
add_result <- function(
  workflow,
  candidate,
  reference = NULL,
  target = NA_real_
) {
  speedup <- if (is.null(reference)) {
    NA_real_
  } else {
    reference[["median"]] /
      candidate[["median"]]
  }
  results[[length(results) + 1L]] <<- data.frame(
    workflow = workflow,
    candidate_median = candidate[["median"]],
    candidate_p90 = candidate[["p90"]],
    reference_median = if (is.null(reference)) {
      NA_real_
    } else {
      reference[["median"]]
    },
    speedup = speedup,
    target = target
  )
}

data <- make_test_data(n_patients = 300, follow_up_time = 15, seed = 701)
baseline <- data[data$time == 1, , drop = FALSE]
time_covariates <- make_time_covariates(data, "time", "time_lin", "time_nlin_1")
model <- VGAM::vglm(
  ordered(y) ~ (time_lin + time_nlin_1) * tx + yprev,
  family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
  data = data
)
reference_sop <- timings(function() {
  soprob_markov_reference(
    model,
    baseline,
    times = 1:15,
    y_levels = 1:6,
    absorb = "6",
    time_covariates = time_covariates
  )
})
candidate_sop <- timings(function() {
  soprob_markov(
    model,
    baseline,
    times = 1:15,
    y_levels = 1:6,
    absorb = "6",
    time_covariates = time_covariates
  )
})
add_result("first-order PO SOP", candidate_sop, reference_sop, target = 5)

trajectory <- expand.grid(
  id = seq_len(500),
  time = seq_len(40),
  KEEP.OUT.ATTRS = FALSE
)
trajectory <- trajectory[order(trajectory$id, trajectory$time), ]
trajectory$y <- 1L + ((trajectory$id + trajectory$time) %% 5L)
trajectory$yprev <- ave(trajectory$y, trajectory$id, FUN = function(x) {
  c(NA_integer_, x[-length(x)])
})
trajectory$tx <- trajectory$id %% 2L
reference_tte <- timings(function() {
  states_to_tte_v2_reference(trajectory, NULL, 6)
})
candidate_tte <- timings(function() {
  states_to_tte_v2(trajectory, covariates = NULL, absorbing_state = 6)
})
add_result("states_to_tte_v2", candidate_tte, reference_tte, target = 10)

bootstrap_data <- expand.grid(
  id = seq_len(1000),
  time = seq_len(20),
  KEEP.OUT.ATTRS = FALSE
)
bootstrap_data$value <- seq_len(nrow(bootstrap_data))
set.seed(702)
boot_ids <- fast_group_bootstrap(bootstrap_data, "id", n_boot = 1L)[[1L]]
reference_materialize <- function() {
  ids <- boot_ids
  names(ids)[names(ids) == "original_id"] <- "id"
  key <- as.character(bootstrap_data$id)
  rows <- split(seq_len(nrow(bootstrap_data)), key)
  pieces <- lapply(seq_len(nrow(ids)), function(i) {
    right <- bootstrap_data[rows[[as.character(ids$id[i])]], c("time", "value")]
    cbind(ids[rep(i, nrow(right)), ], right)
  })
  bind_rows_fill(pieces)
}
reference_bootstrap <- timings(reference_materialize)
candidate_bootstrap <- timings(function() {
  materialize_bootstrap_sample(boot_ids, bootstrap_data, "id")
})
add_result(
  "bootstrap materialization",
  candidate_bootstrap,
  reference_bootstrap,
  target = 5
)

brownian <- timings(
  function() {
    sim_trajectories_brownian(
      n_patients = 3000,
      follow_up_time = 60,
      seed = 703
    )
  },
  iterations = 3L
)
add_result("Brownian simulation 3000x60", brownian)

results <- do.call(rbind, results)
utils::write.csv(results, "benchmark-results.csv", row.names = FALSE)
print(results)

failed <- !is.na(results$target) & results$speedup < results$target
if (any(failed)) {
  stop(
    "Serial benchmark gate failed: ",
    paste(results$workflow[failed], collapse = ", ")
  )
}
