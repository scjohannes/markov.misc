
library(markov.misc)
library(VGAM)
library(dplyr)
library(rms)

# Source dependencies
source("R/vgam_helpers.R")
source("R/sops_fast.R")
source("R/mvn_helpers.R")
# We need bootstrap helpers
source("R/bootstrap.R")
source("R/bootstrap_helpers.R")
# Source standard files needed by bootstrap
if(file.exists("R/sops.R")) source("R/sops.R")
if(file.exists("R/sops_new.R")) source("R/sops_new.R")

# Setup data
set.seed(237567)
N_PATIENTS <- 50 # Small for speed
FU <- 10
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

# Simple Model
fit <- vglm(
  ordered(y) ~ time_lin + tx + yprev, 
  family = cumulative(reverse = TRUE, parallel = TRUE), 
  data = data
)

cat("Running avg_sops_fast...\n")
res_fast <- avg_sops_fast(
  fit, 
  newdata = data, 
  variables = list(tx = c(0, 1)),
  times = 1:FU, 
  ylevels = 1:6, 
  absorb = "6", 
  t_covs = t_covs,
  id_var = "id"
)

cat("Running inferences_fast(method='bootstrap')...\n")
# This should call the standard bootstrap internally
start <- Sys.time()
res_boot <- inferences_fast(
  res_fast, 
  method = "bootstrap", 
  n_boot = 10, # Very small for test
  parallel = FALSE
)
end <- Sys.time()
print(end - start)

print(head(res_boot))

if ("conf.low" %in% names(res_boot)) {
  cat("SUCCESS: Bootstrap CIs added.\n")
} else {
  cat("FAILURE: No CIs found.\n")
}
