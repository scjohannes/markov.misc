library(markov.misc)
library(VGAM)
library(dplyr)
library(rms)

# Source new files manually
source("R/vgam_helpers.R")
source("R/sops_fast.R")
source("R/mvn_helpers.R") # Needed for get_coef/get_vcov_robust
# Source old files to ensure we have the baseline loaded
if (file.exists("R/sops.R")) {
  source("R/sops.R")
}
if (file.exists("R/sops_new.R")) {
  source("R/sops_new.R")
}

# Setup data
set.seed(237567)
N_PATIENTS <- 300
FU <- 30
markov_data <- sim_trajectories_brownian(
  n_patients = N_PATIENTS,
  follow_up_time = FU,
  treatment_prob = 0.5,
  absorbing_state = 6,
  seed = 237567
)

data <- prepare_markov_data(markov_data)
time_spl_m <- rcs(data$time, 4)
data$time_lin <- as.vector(time_spl_m[, 1])
data$time_nlin_1 <- as.vector(time_spl_m[, 2])
data$time_nlin_2 <- as.vector(time_spl_m[, 3])

t_covs <- data |>
  select(time_lin, time_nlin_1, time_nlin_2) |>
  distinct() |>
  arrange(time_lin) |>
  data.frame()

# Constraints setup (Partial PO)
M <- 5
cons_flexible <- list(
  "(Intercept)" = diag(M),
  "time_lin" = cbind(time_lin_1 = 1, time_lin_5 = 0:(M - 1)), # Non-proportional
  "yprev" = cbind(PO = 1),
  "tx" = cbind(PO = 1),
  "time_nlin_1" = cbind(PO = 1),
  "time_nlin_2" = cbind(PO = 1),
  "time_lin:tx" = cbind(PO = 1),
  "time_nlin_1:tx" = cbind(PO = 1),
  "time_nlin_2:tx" = cbind(PO = 1),
  "time_lin:yprev" = cbind(PO = 1)
)

# Fit model
cat("Fitting model...\n")
m_flex <- vglm(
  ordered(y) ~ (time_lin + time_nlin_1 + time_nlin_2) *
    tx +
    yprev +
    yprev:time_lin,
  family = cumulative(reverse = TRUE, parallel = FALSE ~ time_lin),
  data = data,
  constraints = cons_flexible
)
# Ensure constraints are attached (vglm usually does this, but good safety)
m_flex@call$constraints <- cons_flexible

m_flex <- robcov_vglm(m_flex)

# 1. Benchmark Slow
cat("\nRunning Slow Implementation (avg_sops)...\n")
start_slow <- Sys.time()
res_slow <- avg_sops(
  m_flex,
  newdata = data,
  variables = list(tx = c(0, 1)),
  times = 1:FU,
  ylevels = 1:6,
  absorb = "6",
  t_covs = t_covs,
  id_var = "id"
)
end_slow <- Sys.time()
time_slow <- end_slow - start_slow
print(time_slow)

# 2. Benchmark Fast
cat("\nRunning Fast Implementation (avg_sops_fast)...\n")
start_fast <- Sys.time()
res_fast <- avg_sops_fast(
  m_flex,
  newdata = data,
  variables = list(tx = c(0, 1)),
  times = 1:FU,
  ylevels = 1:6,
  absorb = "6",
  t_covs = t_covs,
  id_var = "id"
)
end_fast <- Sys.time()
time_fast <- end_fast - start_fast
print(time_fast)

cat("\nSpeedup Factor:", as.numeric(time_slow) / as.numeric(time_fast), "x\n")

# 3. Validation
cat("\nChecking for equality...\n")
# Ensure ordering is identical
res_slow_sorted <- res_slow[order(res_slow$time, res_slow$state, res_slow$tx), ]
res_fast_sorted <- res_fast[order(res_fast$time, res_fast$state, res_fast$tx), ]

diffs <- abs(res_slow_sorted$estimate - res_fast_sorted$estimate)
max_diff <- max(diffs)
cat("Max absolute difference:", max_diff, "\n")

if (max_diff < 1e-6) {
  cat("SUCCESS: Results are identical!\n")
} else {
  cat("FAILURE: Results differ.\n")
}

# 4. Benchmark Inference (The Main Goal)
cat("\n--- Benchmarking Inference ---\n")
N_SIM <- 100

set.seed(1234)
# A. Slow Inference
cat("Running Slow Inference (standard avg_sops + inferences)...\n")
start_inf_slow <- Sys.time()
# Note: We need to use res_slow (standard avg_sops object)
# And we need to ensure inferences() calls the standard simulation
res_inf_slow <- inferences(res_slow, method = "simulation", n_sim = N_SIM) |>
  ungroup() |>
  arrange(time, state, tx)
end_inf_slow <- Sys.time()
time_inf_slow <- end_inf_slow - start_inf_slow
print(time_inf_slow)

set.seed(1234)
# B. Fast Inference (Aggregated)
cat("Running Fast Inference (avg_sops_fast + inferences_fast)...\n")
start_inf_fast <- Sys.time()
res_inf_fast <- inferences_fast(
  res_fast,
  method = "simulation",
  n_sim = N_SIM
) |>
  ungroup() |>
  arrange(time, state, tx)
end_inf_fast <- Sys.time()
time_inf_fast <- end_inf_fast - start_inf_fast
print(time_inf_fast)

cat(
  "\nInference Speedup Factor:",
  as.numeric(time_inf_slow) / as.numeric(time_inf_fast),
  "x\n"
)

res_inf_fast_sorted <- res_slow[
  order(res_inf_fast$time, res_inf_fast$state, res_inf_fast$tx),
]
res_inf_slow_sorted <- res_fast[
  order(res_inf_slow$time, res_inf_slow$state, res_inf_slow$tx),
]

diffs <- abs(res_inf_slow_sorted$estimate - res_inf_fast_sorted$estimate)
max_diff <- max(diffs)
cat("Max absolute difference:", max_diff, "\n")

if (max_diff < 1e-6) {
  cat("SUCCESS: Results are identical!\n")
} else {
  cat("FAILURE: Results differ.\n")
}

# 5. TEST INDIVIDUAL SOPs & Fallback
cat("\n--- Testing sops_fast and Auto-Build ---\n")

# A. sops_fast (Fast Path)
cat("Running sops_fast...\n")
ind_fast <- sops_fast(
  m_flex,
  newdata = data[data$time == 1, ], # Baseline only for speed
  times = 1:30,
  ylevels = 1:6,
  absorb = "6",
  t_covs = t_covs
)

# B. Inference on sops_fast (Should use stored components)
cat("Running inferences_fast on sops_fast output...\n")
start_ind_inf <- Sys.time()
ind_inf <- inferences_fast(ind_fast, method = "simulation", n_sim = 20)
end_ind_inf <- Sys.time()

print(end_ind_inf - start_ind_inf)
print(head(ind_inf))
if ("conf.low" %in% names(ind_inf)) {
  cat("SUCCESS: CIs added to individual SOPs.\n")
}

# C. Inference on SLOW sops (Should Auto-Build components)
cat("Running inferences_fast on standard sops output (Auto-Build)...\n")
# Create a "slow" object by removing components if any, or just running standard sops
ind_slow <- sops(
  m_flex,
  newdata = data[data$time == 1, ],
  times = 1:30,
  ylevels = 1:6,
  absorb = "6",
  t_covs = t_covs
)
# Ensure attributes needed for rebuild are there (standard sops should have them)
start_auto <- Sys.time()
ind_auto_inf <- inferences_fast(ind_slow, method = "simulation", n_sim = 20)
end_auto <- Sys.time()
print(end_auto - start_auto)
print(head(ind_auto_inf))

if ("conf.low" %in% names(ind_auto_inf)) {
  cat("SUCCESS: Auto-Build worked for individual SOPs.\n")
}
