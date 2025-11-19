# markov.misc

An R package for first-order Markov transition modelling, trajectory
simulation, visualization and operating-characteristic (power) calculations for
discrete health states. The package was developed to support clinical trial
simulation workflows and includes bootstrap helpers for inference on
longitudinal ordinal outcomes.

## Design & architecture

- States are represented by ordered integers (e.g., 1 = home, 6 = death). One or
  more states can be designated absorbing.
- Data can be generated using a markov model with user-specified baseline data
  and user specified linear predictor. Defaults are provided which are suitable
  for an in-hospital viral respiratory disease setting.

## Installation

Install from local source for development with `remotes` or `devtools`:

```r
# install.packages("remotes")
remotes::install_local(path = ".")
```

## Quick start (example)

This snippet shows a minimal simulation and summary workflow. It uses the
package helpers; adapt arguments (n, tmax, extra_params) for your use case.

```r
library(markov.misc)

ids <- sample(1:250000, size = 200, replace = FALSE)
# simulate 200 patients for 28 days with the  VIOLET linear predictor and violet baseline data
set.seed(123)
sim <- simulate_trajectories(
  baseline_data = violet_baseline[violet_baseline$ID %in% ids, ],
  follow_up_time = 30,
  lp_function = lp_violet,      # See R/lp_violet.R for expected args
)

# convert to a t-test style endpoint and to DRS (Days Returned to Baseline)
ttest_df <- markov_to_ttest(sim, target_state = 1)
drs_df <- markov_to_drs(sim, target_state = 1, follow_up_time = 28)

# Expand example
```

## Core functions and responsibilities


## Simulation-to-analysis pipeline (recommended)


