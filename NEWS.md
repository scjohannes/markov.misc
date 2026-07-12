# markov.misc 0.1.0

- Markov/SOP workflows now reject non-logit `orm` and cumulative `vglm` links
  before prediction instead of applying logistic algebra to incompatible fits.
- `inferences()` now reuses serial execution plans, evaluates fixed-order draw
  cells, and generates score-bootstrap perturbations in deterministic bounded
  matrix blocks instead of binding one data frame per draw.
- `inferences()` now compiles second-order ORM/VGLM visit-pair designs once and
  replays coefficient draws through a fused proportional-odds joint-state
  kernel, retaining bounded rolling state without transition tensors.
- `interpolate_sops()` now applies one compiled interpolation plan across
  canonical point and stored-draw grids, with the existing generic fallback for
  irregular or duplicate grids.
- `sim_trajectories_brownian()` and the default `sim_trajectories_markov()` path
  now use fused serial categorical sampling and direct long-output construction.
- `soprob_markov()` now streams visit designs into fused PO or general-logit
  propagation kernels and chunks active second-order state pairs, avoiding dense
  first- and second-order transition tensors.
- `sops()` and `avg_sops()` now fuse BLRM posterior cumulative-logit conversion,
  category differencing, clipping, and normalization in one serial native pass.
- `states_to_tte_v2()` now collapses trajectories with a linear indexed run
  scan, and bootstrap samples are materialized from reusable row-index plans.
- `vglm_markov()` now uses an extensible RMS basis registry and ships first-class
  assignment metadata for both `rcs()` and `lsp()` terms.
- `avg_sops()` now marginalizes counterfactual SOP arrays directly instead of
  materializing individual prediction data frames, and posterior `blrm`
  workflows reduce draws natively before constructing public output.
- `soprob_markov()` and `sops()` now use batched design matrices and compiled
  first-order Markov propagation for `vglm` and `orm`, while `blrm` prediction
  caches design matrices across posterior chunks and uses compiled draw updates.
- Breaking change: legacy `vgam` model support has been removed; use `vglm`
  cumulative ordinal models instead.
- Standardized SOP outputs as base data frames with package-specific S3
  classes, canonical estimate/inference columns, and one internal `draws`
  attribute.
- Homogenized public SOP and plotting argument names, added `plot_operchar()`,
  and removed `standardize_sops()` and `plot_results()` without aliases.
- Simplified `inferences()` to the `mvn`, `score_bootstrap`, `bootstrap`, and
  `fwb` methods. Percentile and point-centered Wald intervals are available for
  every frequentist method; tests are added only for explicit null values.
- Reworked the vignette set
- `avg_comparisons()` now computes average SOP, time-in-state, and ordinal
  time-benefit contrasts between counterfactual levels, with uncertainty added
  through the existing `inferences()` workflow.
- `avg_comparisons(metric = "time_benefit")` inference now reuses the shared
  SOP simulation draw engine, honors `workers` for simulation draws, and applies
  the same longitudinal-data validation used by other refit bootstrap paths.
- `avg_comparisons(metric = "time_benefit")` inference now applies stored
  `time_map`/`origin_time` settings to draw-level SOPs before computing
  real-time AUC intervals.
- `bootstrap_standardized_sops()` and `plot_bootstrap_sops()` have been moved
  to `archive/`; use `avg_sops()` with `inferences(method = "bootstrap")` and
  `plot_sops()` for active bootstrap SOP workflows.
- `inferences()` now validates non-null user-supplied `vcov` objects directly,
  so invalid covariance inputs error instead of silently falling back to the
  model covariance.
- `inferences()` now preserves the original SOP row order when attaching
  simulation or bootstrap interval summaries.
- `plot_comparisons()` now plots `avg_comparisons()` output on the difference
  or ratio scale, with uncertainty intervals when available.
- `plot_correlation()` and `plot_variogram()` now compute model-implied
  ordinal correlations from exact pairwise moments instead of simulated paths,
  including draw-specific moment calculations for `blrm` fits before posterior
  averaging. Second-order model correlations now propagate each start time
  forward once and reuse the resulting cross-moments across all later plot
  times. `plot_transitions()`, `plot_correlation()`, and `plot_variogram()` now
  expose `n_draws` for `blrm` diagnostic summaries.
- `plot_lp_difference()` now plots profile-based treatment differences in the
  ordinal-model linear predictor over time, faceted by previous state,
  including posterior-median linear predictor contrasts for `blrm` models.
- `sops()`, `avg_sops()`, and `avg_comparisons()` now require an explicit
  `times` argument instead of inferring a prediction grid, expose Markov
  structure arguments directly, and reserve `refit_data` for refit-bootstrap
  inference.
- Documented `violet_baseline` provenance from `Hmisc::simlongord` and kept
  the package under `GPL (>= 2)` for GPL-compatible redistribution.
- Removed stale `taooh*` archive exports and fixed several edge cases in
  `soprob_markov()`, `states_to_tte()`, `sim_trajectories_markov()`,
  `sim_trajectories_tte()`, `predict_blrm_response_markov()`,
  `apply_to_bootstrap()`, `vglm_markov()`, and `lp_violet()`.
- `sample_from_arrow()` now validates sampling controls up front and correctly
  supports `replace = TRUE` samples larger than the observed arm-specific ID
  counts by materializing duplicate patients with synthetic IDs.
- SOP prediction now renormalizes clipped category probabilities from
  cumulative-logit conversions so each finite positive patient/time
  distribution sums to one.
- Breaking change: the old dotted VGAM wrapper name `vglm.markov()` has been
  removed. Use `vglm_markov()` instead.
- `soprob_markov()`, `sops()`, and `avg_sops()` now support second-order
  Markov recursion via `p2_var`, and `rmsb::blrm()` models use sampled
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
  `inferences(..., method = "score_bootstrap", cluster = <id>)`.
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
  `method = "fwb"`, using mean-one exponential
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
- `inferences(method = "fwb")` now supports individual
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
- `plot_sops()` now respects the stored `y_levels` order on model-derived SOP
  objects, so character state labels such as `"10"` no longer sort before
  `"2"` in color and fill scales.
- `plot_transitions()` now plots empirical or model-based joint transition
  proportions as heatmaps, including model-based treatment-difference
  heatmaps from deterministic counterfactual transition traces, and orders
  numeric-looking time facets numerically. For `blrm` model plots, posterior
  draw-specific transition traces are summarized directly using the existing
  manual `blrm` prediction backend. Model-based plots evaluate the full
  intermediate visit grid before subsetting to requested plot times.
- `plot_variogram()` now plots empirical or model-based correlations
  against absolute time differences with a 0-to-1 coordinate window, preserving
  underlying correlation values for ggplot statistics.
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
