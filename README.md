# markov.misc

Compute and plot state occupancy probabilities (SOPs) and contrasts from first- and second-order Markov models of ordinal health-state trajectories.

## Installation

Install the development version from GitHub with `pak`:

```r
install.packages("pak")
pak::pak("scjohannes/markov.misc")
```

## ACTT-2-Style SOP Workflow

This example simulates ACTT-2-style ordinal outcomes, fits a proportional-odds
Markov transition model with `orm_markov()`, estimates marginal treatment-arm
SOPs, adds uncertainty, and plots the result.

```r
library(markov.misc)
library(rms)
library(ggplot2)

set.seed(20260526)

trial <- sim_actt2_brownian(
  n_patients = 300,
  follow_up_time = 28,
  treatment_prob = 0.5,
  mu_treatment_effect = -0.03,
  seed = 20260526
)

markov_data <- prepare_markov_data(trial, absorbing_state = 8)

dd <- datadist(markov_data)
options(datadist = "dd")

fit <- orm_markov(
  y ~ rms::rcs(time, 4) + tx + yprev,
  data = markov_data,
  id_var = "id",
  opt_method = "LM",
  scale = TRUE
)

sop <- avg_sops(
  fit,
  variables = list(tx = c(0, 1)),
  times = 1:28,
  y_levels = fit$yunique,
  absorb = "8"
)

sop_ci <- inferences(
  sop,
  method = "mvn",
  n_draws = 200,
  conf_level = 0.95
)

plot_sops(sop_ci, facet_var = "tx") +
  labs(
    x = "Day",
    y = "State occupancy probability"
  )
```

## Analytical Confidence Intervals

The default inference method remains `method = "mvn"`. For a deterministic
first-order delta-method calculation, reuse the same averaged full
proportional-odds SOP object (`sop`, created above by `avg_sops()`) and state the
averaging target explicitly:

```r
sop_empirical <- inferences(
  sop,
  method = "delta",
  target = "empirical"
)

sop_population <- inferences(
  sop,
  method = "delta",
  target = "population"
)

# Materialize only the analytical rows needed downstream.
J <- get_jacobian(sop_empirical, rows = 1:8)
V <- stats::vcov(sop_population, rows = 1:8)
```

The empirical target conditions on the observed standardization profiles. The
population target treats the same stored fitting cohort as sampled and combines
profile variation with fitted-model score variation in a patient-level stacked
influence function. `orm_markov(id_var = "id")` and
`vglm_markov(id_var = "id")` store the row-aligned patient IDs needed for this
calculation. Patient-cluster robustness protects the variance against arbitrary
within-patient score correlation; it does not correct transition-model bias or
Markov/proportional-odds misspecification.

The empirical and population targets apply only to averaged SOP or comparison
objects. For an individual `sops()` result, omit `target` or use
`target = "fixed"`; no other analytical target is accepted.

## Learn More

After installation, see:

```r
vignette("full-po-sops", package = "markov.misc")
vignette("partial-po-vglm-sops", package = "markov.misc")
vignette("second-order-orm-sops", package = "markov.misc")
vignette("many-levels-previous-state-spline", package = "markov.misc")
vignette("factor-time-orm-sops", package = "markov.misc")
vignette("analytical-confidence-intervals", package = "markov.misc")
vignette("bayesian-rmsb-sops", package = "markov.misc")
vignette("mvn-vs-bootstrap-orm-sops", package = "markov.misc")
```

## Data Provenance

The included `violet_baseline` dataset is derived from `Hmisc::simlongord`,
created by Frank Harrell and based on the VIOLET trial. Because that data
provenance is GPL-compatible, `markov.misc` is licensed as `GPL (>= 2)`.
