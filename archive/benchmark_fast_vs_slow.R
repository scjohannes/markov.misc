# =============================================================================
# Benchmark: Fast vs Slow Simulation-Based Inference
# =============================================================================
#
# This script demonstrates the performance difference between:
# 1. SLOW PATH: Calling VGAM::predict() for each simulation draw
# 2. FAST PATH: Pre-computing design matrices and using matrix multiplication
#
# The fast path achieves 10-100x speedup depending on data size and model complexity.
# =============================================================================

library(markov.misc)
library(VGAM)
library(dplyr)
library(rms)

# Source the archived slow implementation
source("archive/sops_old_slow.R")

# =============================================================================
# Setup Test Data and Model
# =============================================================================

cat("=== Setting up test data and model ===\n\n")

set.seed(237567)
N_PATIENTS <- 100
FU <- 20

markov_data <- sim_trajectories_brownian(
  n_patients = N_PATIENTS,
  follow_up_time = FU,
  treatment_prob = 0.5,
  absorbing_state = 6,
  seed = 237567
)

data <- prepare_markov_data(markov_data)

# Create spline basis for time
time_spl_m <- rcs(data$time, 4)
data$time_lin <- as.vector(time_spl_m[, 1])
data$time_nlin_1 <- as.vector(time_spl_m[, 2])
data$time_nlin_2 <- as.vector(time_spl_m[, 3])

t_covs <- data |>
  select(time_lin, time_nlin_1, time_nlin_2) |>
  distinct() |>
  arrange(time_lin) |>
  data.frame()

# Constraints setup (Partial PO - time has NPO, everything else PO)
M <- 5
cons_flexible <- list(
  "(Intercept)" = diag(M),
  "time_lin" = cbind(time_lin_1 = 1, time_lin_5 = 0:(M - 1)),
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
cat("Fitting partial proportional odds model...\n")
m_flex <- vglm(
  ordered(y) ~ (time_lin + time_nlin_1 + time_nlin_2) *
    tx +
    yprev +
    yprev:time_lin,
  family = cumulative(reverse = TRUE, parallel = FALSE ~ time_lin),
  data = data,
  constraints = cons_flexible
)
m_flex@call$constraints <- cons_flexible

# Add robust covariance
m_flex <- robcov_vglm(m_flex, cluster = data$id)

cat("Model fitted with", length(coef(m_flex)), "coefficients\n\n")

# =============================================================================
# Benchmark Point Estimates (avg_sops)
# =============================================================================

cat("=== Benchmarking Point Estimates ===\n\n")

# Calculate avg_sops (point estimates only - should be fast for both)
cat("Running avg_sops()...\n")
start_point <- Sys.time()
res_base <- avg_sops(
  m_flex,
  newdata = data,
  variables = list(tx = c(0, 1)),
  times = 1:FU,
  ylevels = 1:6,
  absorb = "6",
  t_covs = t_covs,
  id_var = "id"
)
time_point <- Sys.time() - start_point
cat("Time for point estimates:", format(time_point, digits = 3), "\n\n")

# =============================================================================
# Benchmark Inference: SLOW vs FAST
# =============================================================================

N_SIM <- 1000 # Number of simulation draws (increase for more accurate benchmark)

cat("=== Benchmarking Inference (N_SIM =", N_SIM, ") ===\n\n")

# --- SLOW PATH ---
cat("Running SLOW inference (predict() per draw)...\n")
set.seed(1234)
start_slow <- Sys.time()
res_slow <- inferences_simulation_slow(
  res_base,
  n_sim = N_SIM,
  vcov = NULL,
  conf_level = 0.95,
  conf_type = "perc",
  parallel = FALSE,
  workers = NULL,
  return_draws = FALSE
)
time_slow <- Sys.time() - start_slow
cat("SLOW path time:", format(time_slow, digits = 3), "\n")

# --- FAST PATH ---
cat("Running FAST inference (pre-computed matrices)...\n")
set.seed(1234)
start_fast <- Sys.time()
res_fast <- inferences(
  res_base,
  method = "simulation",
  n_sim = N_SIM,
  conf_level = 0.95,
  conf_type = "perc",
  parallel = FALSE,
  return_draws = FALSE
)
time_fast <- Sys.time() - start_fast
cat("FAST path time:", format(time_fast, digits = 3), "\n\n")

# =============================================================================
# Results
# =============================================================================

speedup <- as.numeric(time_slow, units = "secs") /
  as.numeric(time_fast, units = "secs")

cat("=== RESULTS ===\n\n")
cat("SLOW path:", format(time_slow, digits = 3), "\n")
cat("FAST path:", format(time_fast, digits = 3), "\n")
cat("SPEEDUP:  ", round(speedup, 1), "x faster\n\n")

# Validate results are numerically equivalent
res_slow_sorted <- res_slow[order(res_slow$time, res_slow$state, res_slow$tx), ]
res_fast_sorted <- res_fast[order(res_fast$time, res_fast$state, res_fast$tx), ]

# Point estimates should be identical (same seed)
diff_estimate <- max(abs(res_slow_sorted$estimate - res_fast_sorted$estimate))
cat(
  "Max difference in point estimates:",
  format(diff_estimate, scientific = TRUE),
  "\n"
)

# CIs should be very similar (same seed, same algorithm)
diff_ci_low <- max(
  abs(res_slow_sorted$conf.low - res_fast_sorted$conf.low),
  na.rm = TRUE
)
diff_ci_high <- max(
  abs(res_slow_sorted$conf.high - res_fast_sorted$conf.high),
  na.rm = TRUE
)
cat("Max difference in conf.low:", format(diff_ci_low, scientific = TRUE), "\n")
cat(
  "Max difference in conf.high:",
  format(diff_ci_high, scientific = TRUE),
  "\n\n"
)

if (diff_estimate < 1e-10 && diff_ci_low < 0.01 && diff_ci_high < 0.01) {
  cat("SUCCESS: Results are numerically equivalent!\n")
} else {
  cat("WARNING: Results differ more than expected.\n")
}

# =============================================================================
# Scaling Analysis (Optional)
# =============================================================================

cat("\n=== Scaling Analysis ===\n\n")

scaling_results <- data.frame(
  n_sim = integer(),
  time_slow = numeric(),
  time_fast = numeric(),
  speedup = numeric()
)

for (n in c(10, 25, 50)) {
  cat("N_SIM =", n, "... ")

  set.seed(42)
  t1 <- system.time({
    inferences_simulation_slow(
      res_base,
      n_sim = n,
      vcov = NULL,
      conf_level = 0.95,
      conf_type = "perc",
      parallel = FALSE,
      workers = NULL,
      return_draws = FALSE
    )
  })["elapsed"]

  set.seed(42)
  t2 <- system.time({
    inferences(
      res_base,
      method = "simulation",
      n_sim = n,
      conf_level = 0.95,
      conf_type = "perc",
      parallel = FALSE,
      return_draws = FALSE
    )
  })["elapsed"]

  scaling_results <- rbind(
    scaling_results,
    data.frame(
      n_sim = n,
      time_slow = t1,
      time_fast = t2,
      speedup = t1 / t2
    )
  )

  cat(
    "slow:",
    round(t1, 2),
    "s, fast:",
    round(t2, 2),
    "s, speedup:",
    round(t1 / t2, 1),
    "x\n"
  )
}

cat("\nScaling Summary:\n")
print(scaling_results)
