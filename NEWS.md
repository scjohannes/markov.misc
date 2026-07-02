# markov.misc 0.1.0

- Reworked the vignette set
- Documented `violet_baseline` provenance from `Hmisc::simlongord` and kept
  the package under `GPL (>= 2)` for GPL-compatible redistribution.
- Removed stale `taooh*` archive exports and fixed several edge cases in
  `soprob_markov()`, `states_to_tte()`, `sim_trajectories_markov()`,
  `sim_trajectories_tte()`, `predict_blrm_response_markov()`,
  `apply_to_bootstrap()`, `vglm.markov()`, and `lp_violet()`.
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
- `inferences()` now supports fractional weighted bootstrap refits via
  `method = "bootstrap", engine = "fwb"`, using mean-one exponential
  patient-level weights for fitting and weighted SOP marginalization.
- `interpolate_sops()` maps visit-scale SOP output to real elapsed time with
  optional empirical baseline anchoring, and `time_in_state()` can now compute
  trapezoidal real-time AUC while `soprob_markov()`, `sops()`, and `avg_sops()`
  support factor-valued visit indices.
- `plot_sops()` now defaults to line plots with viridis discrete scales and
  supports model-derived SOP summaries from `avg_sops()` and `inferences()`,
  including confidence ribbons for line plots and low-alpha draw overlays for
  stacked bar plots when draws are stored.
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
