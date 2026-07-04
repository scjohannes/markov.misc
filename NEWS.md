# markov.misc 0.1.0

- Reworked the vignette set
- Documented `violet_baseline` provenance from `Hmisc::simlongord` and kept
  the package under `GPL (>= 2)` for GPL-compatible redistribution.
- Removed stale `taooh*` archive exports and fixed several edge cases in
  `soprob_markov()`, `states_to_tte()`, `sim_trajectories_markov()`,
  `sim_trajectories_tte()`, `predict_blrm_response_markov()`,
  `apply_to_bootstrap()`, `vglm_markov()`, and `lp_violet()`.
- Breaking change: the old dotted VGAM wrapper name `vglm.markov()` has been
  removed. Use `vglm_markov()` instead.
- `soprob_markov()`, `sops()`, and `avg_sops()` now support second-order
  Markov recursion via `p2varname`, and `rmsb::blrm()` models use sampled
  posterior draws as the native uncertainty path with optional fitted random
  effects from `cluster()`.
- `sim_actt2_brownian()` now uses calibrated Brownian-gap defaults that reduce
  state `1:2` to state `3:7` rehospitalization-like churn, while
  `sim_trajectories_brownian_gap()` supports sampled patient-specific drift
  starts and scalar or threshold-specific time and treatment effects.
- `sim_trajectories_markov()` now supports partial proportional odds data
  generation when `lp_function` returns one linear predictor per threshold, and
  stops when threshold-specific predictors imply negative state probabilities.
- `blrm` SOP prediction now caches fitted random-effect draws once per call and
  vectorizes posterior recursion within chunks; state-wise median summaries are
  documented as not necessarily summing to one across states, while draw-level
  probabilities and mean summaries preserve total probability.
- `avg_sops()`, `sops()`, and related Markov SOP workflows now reject models fit
  with offsets because offsets are not supported by the package prediction
  paths.
- `avg_sops()`, `sops()`, `soprob_markov()`, and `inferences()` now treat
  `rms::orm()` as a first-class proportional-odds backend, including fast-path
  SOP prediction, MVN inference with full `rms::robcov()` covariance matrices,
  refit bootstrap inference, and score-bootstrap inference via
  `inferences(..., engine = "score_bootstrap", cluster = <id>)`.
- Added `orm_markov()` and `blrm_markov()`, and expanded
  `vglm_markov(id_var = ...)`, so wrapper-fitted models store their original
  longitudinal data and ID variable, frequentist wrappers compute
  cluster-robust covariance automatically, and `sops()`/`avg_sops()` can extract
  one prediction row per patient when `newdata = NULL`.
- `orm_markov()` now stores an update-safe `rms::orm()` call so fractional
  weighted bootstrap refits work for grouped and ungrouped `sops()` objects.
- `sops()` and `avg_sops()` now treat user-supplied `newdata` rows as fixed
  prediction or standardization profiles, regenerate their internal `rowid`
  values, and reserve `id_var` for stored-data extraction, refit bootstrap
  clustering, and `blrm` random-effect prediction.
- `inferences()` now supports fractional weighted bootstrap refits via
  `method = "bootstrap", engine = "fwb"`, using mean-one exponential
  patient-level weights for fitting and weighted SOP marginalization.
- FWB and score-bootstrap inference now use draw-specific patient weights for
  every empirical averaging step in `avg_sops()` and grouped `sops()`, with
  weights normalized within `by` strata. Ungrouped stored-data `sops()` draws
  expose `fwb_weight` or `score_weight` via `get_draws()` for manual summaries;
  these column names are reserved when draw weights are attached.
- `inferences()` now treats user-supplied `avg_sops(newdata = ...)`
  standardization profiles as fixed targets for score bootstrap and FWB
  inference: coefficient/refit uncertainty still comes from the original data,
  but profile averaging is unweighted and a warning reports that
  `baseline_weights` were set to `NULL`.
- `inferences()` now accepts `Matrix` package covariance objects returned by
  `rms::robcov()` for `orm_markov()` fits by coercing them to base matrices
  before coefficient/covariance validation.
- `inferences(method = "bootstrap", engine = "fwb")` now supports individual
  `sops()` objects when full refit data and `id_var` metadata are available.
  Ordinary standard bootstrap remains limited to `avg_sops()` because it can
  drop state support needed by fixed individual prediction rows.
- `interpolate_sops()` maps visit-scale SOP output to real elapsed time with
  optional empirical baseline anchoring, and `time_in_state()` can now compute
  trapezoidal real-time AUC while `soprob_markov()`, `sops()`, and `avg_sops()`
  support factor-valued visit indices.
- `plot_sops()` now defaults to line plots with viridis discrete scales and
  supports model-derived SOP summaries from `avg_sops()` and `inferences()`,
  including confidence ribbons for line plots and low-alpha draw overlays for
  stacked bar plots when draws are stored.
- `plot_sops()` now respects the stored `ylevels` order on model-derived SOP
  objects, so character state labels such as `"10"` no longer sort before
  `"2"` in color and fill scales.
- `prepare_markov_data()`, `soprob_markov()`, `avg_sops()`, and `inferences()`
  now support numeric previous-state effects, including nonlinear terms such as
  `rms::rcs(yprev, 6)`, while preserving factor previous-state behavior by
  default.
- `robcov_vglm()` now stores `bread = vcov(fit)` and unscaled score
  crossproducts as `meat`, aligns pre-NA cluster vectors when possible, rejects
  missing or single-valued clusters, and warns that clustered z-test p-values
  may be anti-conservative with few clusters.
- `robcov_vglm()` and `get_vcov_robust()` now apply the G/(G-1) small-sample
  correction by default for clustered robust covariance estimates; set
  `adjust = FALSE` to recover the previous unadjusted behavior.
- `standardize_sops()` now reuses the same sampled `blrm` posterior draws for
  treatment and control counterfactuals, preserving draw-wise pairing for
  downstream contrasts.
