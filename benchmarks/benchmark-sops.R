devtools::load_all()
source("tests/testthat/helper-simulation.R")

elapsed <- function(code, iterations = 5L) {
  values <- numeric(iterations)
  for (i in seq_len(iterations)) {
    gc()
    values[i] <- unname(system.time(code())["elapsed"])
  }
  stats::median(values)
}

data <- make_test_data(n_patients = 300, follow_up_time = 15, seed = 701)
baseline <- data[data$time == 1, , drop = FALSE]
time_covariates <- make_time_covariates(
  data,
  "time",
  "time_lin",
  "time_nlin_1"
)
model <- VGAM::vglm(
  ordered(y) ~ (time_lin + time_nlin_1) * tx + yprev,
  family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
  data = data
)

reference <- elapsed(function() {
  soprob_markov_reference(
    model,
    baseline,
    times = 1:15,
    y_levels = 1:6,
    absorb = "6",
    time_covariates = time_covariates
  )
})
optimized <- elapsed(function() {
  soprob_markov(
    model,
    baseline,
    times = 1:15,
    y_levels = 1:6,
    absorb = "6",
    time_covariates = time_covariates
  )
})
averaged <- elapsed(function() {
  avg_sops(
    model,
    baseline,
    variables = list(tx = 0:1),
    times = 1:15,
    y_levels = 1:6,
    absorb = "6",
    time_covariates = time_covariates
  )
})

print(data.frame(
  workflow = c("reference SOP", "optimized SOP", "optimized avg_sops"),
  seconds = c(reference, optimized, averaged),
  reference_speedup = c(1, reference / optimized, NA_real_)
))
