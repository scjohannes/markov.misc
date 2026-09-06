# Analytical Confidence Intervals Implementation Plan

Status date: 2026-07-19
Working branch: `analytical-cis`

This is the durable implementation checklist for analytical confidence
intervals for `sops()`, `avg_sops()`, and supported `avg_comparisons()` results.
It records the approved design, what is currently present in the branch, and
the post-implementation concerns that will be investigated one at a time.

## Approved Design

- Add deterministic inference through
  `inferences(method = "delta", vcov = ..., conf_type = "auto")` without
  changing the existing `method = "mvn"` default.
- Differentiate the first-order full proportional-odds SOP recursion
  analytically on each backend's complete raw coefficient scale. Central finite
  differences are a test oracle only, not a production fallback.
- Use patient-cluster robust coefficient covariance for fixed-profile and
  empirical-cohort propagation unless the user supplies a valid complete
  coefficient covariance.
- For averaged objects, treat `vcov = "conditional"` as conditional on the
  observed standardization profiles. Treat `vcov = "unconditional"` as a
  fitted-cohort superpopulation target and use a stacked patient influence
  function that retains profile/model-score covariance. Averaged results default
  to unconditional inference. Individual `sops()` defaults to conditional and
  rejects unconditional inference. A coefficient matrix selects conditional
  inference; `NULL` selects the class default. Remove the public `target` argument.
- Store factorized analytical state and materialize selected Jacobian or
  covariance blocks with `get_jacobian()` and `vcov()`.
- Support linear differences for SOP and time-in-state estimands, including
  visit-to-real-time mapping, observed-baseline anchoring, linear interpolation,
  and trapezoidal integration.
- Fail explicitly outside the approved scope rather than silently switching to
  numerical differentiation, row-level clusters, or another inferential target.

## Target and Interval Contract

| Result | Default / allowed target | Analytical interval |
| --- | --- | --- |
| Individual `sops()` | `"conditional"` or coefficient matrix | Componentwise logit-delta under `conf_type = "auto"` |
| `avg_sops()` | `"unconditional"`; conditional or matrix also accepted | Componentwise logit-delta under `conf_type = "auto"` |
| Supported `avg_comparisons()` | `"unconditional"`; conditional or matrix also accepted | Identity-scale Wald under `conf_type = "auto"` |
| Stored-cohort averages/comparisons | `"unconditional"` | Fitted-cohort stacked influence covariance |

Patient is always the independent cluster. An explicit row-aligned `cluster`
vector takes precedence; otherwise stored fitting data and stored `id_var`
metadata are required. Cluster robustness protects variance calculations from
arbitrary within-patient score correlation. It does not correct transition-model
bias, informative missingness, or Markov/proportional-odds misspecification.

## Implementation Checklist

### Public Entry and Dispatch

- [x] Add `"delta"` to the `inferences()` method choices.
- [x] Select analytical variance through `vcov` and `conf_type = "auto"` while retaining
  `method = "mvn"` as the default.
- [x] Route SOP objects and average-comparison objects to separate analytical
  handlers.
- [x] Remove stale draw attributes from delta results and attach method, target,
  interval, and covariance-source metadata.

### Analytic SOP Core

- [x] Map ORM and VGLM native coefficients to effective threshold/design
  coefficients with exact names and order.
- [x] Differentiate reverse cumulative-logit category probabilities
  analytically.
- [x] Propagate SOP probabilities and raw-coefficient Jacobians through the
  first-order Markov product rule, including absorbing states.
- [x] Reject unsupported order/link/partial-PO structures and crossed raw
  ordinal probabilities.
- [x] Guard analytical allocations with
  `markov.misc.delta_max_bytes` and the typed
  `markov_misc_delta_too_large` condition.

### Covariance and Superpopulation Targets

- [x] Validate complete named, finite, symmetric, positive-semidefinite custom
  covariance matrices for coefficient-form targets.
- [x] Resolve or compute patient-cluster robust covariance while preserving
  VGLM finite-cluster metadata and the ORM HC0 convention.
- [x] Aggregate raw transition-row scores by patient and construct inverse
  per-patient sensitivity from inverse total information, not from the robust
  sandwich covariance.
- [x] Align designated starting-profile and score IDs exactly. Include only
  fitted patients with at least one usable likelihood transition, reject
  score-only patients whose starting profile is unavailable, and never insert
  profile-only zero-score patients.
- [x] Form vector stacked influences and use
  `stats::cov(influence) / n` with the finite-sample convention reported in
  metadata.
- [x] Reject custom covariance and user-supplied external profiles for the
  superpopulation target.

### Result Access and Comparisons

- [x] Store coefficient-form analytical state as `J` plus coefficient `V`.
- [x] Store superpopulation-form analytical state as average `J` plus the
  patient influence matrix.
- [x] Add row-selective `get_jacobian()` and S3 `vcov()` methods without storing
  a dense all-cell covariance.
- [x] Implement linear operators for SOP and time-in-state differences.
- [x] Reproduce point estimates before propagating a comparison operator.
- [x] Support factor-time real-time differences through a supplied `time_map`,
  observed-baseline anchor, interpolation grid, and trapezoidal weights.

### Tests and Documentation

- [x] Add focused core tests for ORM/VGLM derivatives, absorbing states, factor
  visits, validation conditions, and allocation guards.
- [x] Add focused superpopulation tests for cluster resolution, covariance
  conventions, exact score/profile alignment, profile-only exclusion, and the
  stacked cross term.
- [x] Add focused public/accessor/comparison tests for dispatch, interval
  defaults, selected covariance blocks, and linear comparison operators.
- [x] Update `ARCHITECTURE.md`, `README.md`, `NEWS.md`, and the factor-time
  vignette; add a dedicated analytical-CI vignette and this plan.
- [x] Regenerate roxygen documentation after the public code stabilizes.
- [x] Confirm the focused namespace tests pass with the regenerated
  documentation.
- [x] Run and confirm the complete `devtools::test()` suite on the integrated
  branch, with only the expected installed-package skip.
- [x] Run and confirm the final package check with the repository's Windows R
  and Rtools environment settings: 0 errors, 0 warnings, and 0 notes.

## Supported Scope at Integration

- First-order Markov recursion.
- Full proportional odds with reverse cumulative logit.
- Frequentist `orm`, `vglm`, and `robcov_vglm` backends.
- Individual SOPs, empirical average SOPs, and fitted-cohort superpopulation
  average SOPs.
- Difference comparisons for `estimand = "sop"` and
  `estimand = "time_in_state"`.
- Numeric or factor visit designs supported by the first-order execution plan.
- Componentwise logit-delta intervals for probabilities and Wald intervals for
  comparisons.

The following are not part of this integration: ratios, `time_benefit`,
`by`-stratum delta inference, partial/nonproportional odds, second-order delta
recursion, BLRM delta inference, and simultaneous confidence bands. Random or
independent external target populations are intentionally unsupported and out
of scope. MVN/bootstrap numerical agreement is not a required validation
criterion for the first-order delta method.

## Planned GitHub Follow-ups

The milestone and five issues below track the approved follow-ups after the
current implementation passes the final integrated package check.

- [x] Create the **Analytical CI follow-ups**
  [milestone](https://github.com/scjohannes/markov.misc/milestone/1).
- [x] Create a [ratio-comparison issue](https://github.com/scjohannes/markov.misc/issues/75)
  with an explicit transformation and null testing contract.
- [x] Create a [`time_benefit` derivative/operator issue](https://github.com/scjohannes/markov.misc/issues/79).
- [x] Create a [grouped `by`-stratum target and cluster-alignment issue](https://github.com/scjohannes/markov.misc/issues/76).
- [x] Create a [partial-proportional-odds derivative issue](https://github.com/scjohannes/markov.misc/issues/78).
- [x] Create a [second-order recursion derivative issue](https://github.com/scjohannes/markov.misc/issues/77).

## Final Integration Gate

- [x] Review the integrated diff for accidental API/default changes, including
  preserving the established positional argument order of `inferences()`.
- [x] Confirm generated `NAMESPACE` and `man/` changes after documentation.
- [x] Confirm complete tests and package check pass.
- [x] Confirm the planned GitHub milestone/issues have been created or
  explicitly moved to a later release plan.

## Post-Implementation Concern Register

This section is the source of truth for analytical-CI follow-up review. Keep the
stable `ACI-*` identifiers when updating code, tests, documentation, issues, or
commits. Work on one concern at a time unless two concerns are inseparable.

Allowed statuses are:

- **Open**: recorded but not currently being investigated.
- **Active**: the single concern currently under investigation.
- **Resolved**: evidence supports a code, documentation, or design resolution.
- **Accepted**: the behavior is retained as an explicitly accepted limitation
  or trade-off.
- **Deferred**: intentionally postponed, with the reason and destination
  recorded.

For every concern:

1. Change its status to **Active** before starting work.
2. Record the derivation, experiment, or external contract used to assess it.
3. Record alternatives considered and why the selected resolution was chosen.
4. Add or update focused regression tests when behavior or code changes.
5. Update user documentation and `ARCHITECTURE.md` when the contract changes.
6. Record validation commands and results in the concern's resolution log.
7. Finish with **Resolved**, **Accepted**, or **Deferred**; never delete the
   concern or its history.

### Priority Summary

| ID | Priority | Status | Concern |
| --- | --- | --- | --- |
| ACI-01 | High | Resolved | Superpopulation finite-sample normalization and fitted-cohort contract |
| ACI-02 | High | Resolved | Superpopulation score orientation and cross-term scaling |
| ACI-03 | High | Resolved | Penalized ORM superpopulation inference |
| ACI-04 | High | Resolved | Zero-score profile sensitivity scaling |
| ACI-05 | High | Resolved | Stored ORM robust covariance identity and correction metadata |
| ACI-06 | Medium | Resolved | Unnamed penalized-ORM covariance relabeling |
| ACI-07 | Medium | Resolved | VGLM raw-to-effective constraint mapping |
| ACI-08 | Medium | Resolved | Structural boundary classification |
| ACI-09 | Medium | Resolved | First-follow-up profiles and real-time integration conventions |
| ACI-10 | Medium | Open | Memory accounting versus actual process memory |
| ACI-11 | Medium | Resolved | Grouped native execution row-layout contract |
| ACI-12 | Medium | Resolved | Superpopulation `get_jacobian()` semantics |
| ACI-13 | Medium | Open | `vcov()` dispatch for non-delta result objects |
| ACI-14 | Medium | Resolved | Dense comparison and covariance materialization |
| ACI-15 | Medium | Open | Numerical validation tolerances and portability |
| ACI-16 | Medium | Open | Native C++ maintenance and semantic parity |
| ACI-17 | Low | Open | Performance benchmark generalizability |
| ACI-18 | High | Open | Independent superpopulation-inference validation oracle |
| ACI-19 | Low | Resolved | Public target terminology and formal defaults |
| ACI-20 | Low | Open | Generated native build artifacts in the worktree |
| ACI-21 | Medium | Resolved | Explicit zero ORM penalty list handling |

### ACI-01: Superpopulation finite-sample normalization and fitted-cohort contract

- **Status:** Resolved on 2026-07-19
- **Priority:** High
- **Current decision:** Superpopulation covariance is
  `stats::cov(influence) / n`. Because `stats::cov()` uses denominator `n - 1`,
  this is the manuscript's patient-level sample-covariance correction,
  `n / (n - 1)`, relative to the centered patient-level HC0 expression. Do not
  call this backend HC1. The unreleased `target = "population"` spelling is
  removed rather than deprecated; the public target is
  `vcov = "unconditional"`.
- **Concern:** The covariance normalization must remain distinct from backend
  row-level HC/cadjust conventions. Its patient cohort must also match the
  fitted cohort exactly without introducing profile-only zero-score patients.
- **Resolution approach:** Retain `stats::cov(phi) / n` and test its exact
  centered-matrix identity. Refactor fitting wrappers to retain one complete
  designated-origin prediction profile before response-driven model-frame row
  omission. Include exactly the patients who both have that profile and
  contribute at least one usable likelihood transition anywhere in the fitted
  data. A missing first transition outcome alone is allowed; a fitted patient
  without the required origin profile is an error; a profile-only patient is
  excluded. Remove zero-score insertion and sensitivity-ratio rescaling.
- **Resolution log:** Implemented the user-approved contract:
  every fitted patient must have a complete first-follow-up profile containing
  the prediction variables and previous state, but not the first transition
  response. Numeric `first_followup_time` defaults to 1 and factor/character
  time requires an explicit value. Later likelihood rows do not substitute for
  the designated first-follow-up row. Automatic inference uses model-stored profile
  data; fixed user `newdata` remains empirical-only; `refit_data` remains solely
  refitting infrastructure. The implementation removes profile-only zero scores
  and sensitivity-ratio scaling, rejects penalized or weighted ORM
  superpopulation inference before Jacobian allocation, exposes the covariance
  convention and ignored backend HC settings in metadata, and retains
  `stats::cov(phi) / n`. Work began after checkpoint commit `7b04ff6`.
  `air format .`, focused analytical/profile tests, the complete test suite, and
  `devtools::document()` passed. The analytical-CI, factor-time ORM, and
  many-level previous-state spline vignettes rendered successfully. The final
  `devtools::check()` completed with 0 errors and 0 warnings; its sole note was
  the isolated environment's inability to verify the current network time.

### ACI-02: Superpopulation score orientation and cross-term scaling

- **Status:** Resolved on 2026-07-22
- **Priority:** High
- **Current decision:** Use
  `r_i - mu + G %*% Ainv %*% s_i`, with backend score columns aligned to the raw
  coefficient vector.
- **Concern:** The positive sign and scaling follow the approved derivation, and
  score/bread combinations reproduce coefficient covariance, but the profile
  cross-term has not been validated with an independent estimating-equation
  implementation. Existing hand calculations intentionally encode the same
  formula and therefore are not an independent oracle.
- **Resolution approach:** Numerically perturb patient estimating-equation
  weights and solve or one-step-update coefficients to verify the direction and
  scale of each patient's coefficient influence. Check ORM and VGLM separately
  without requiring agreement with MVN draws or bootstrap percentiles.
- **Resolution log:** Independent backend-specific weighted-refit tests now
  validate the current formula for unpenalized ORM and full-PO VGLM fits. For
  three high-score patients per backend, every likelihood row belonging to the
  patient and the corresponding starting-profile averaging weight are perturbed
  symmetrically. The tests compare the numerical coefficient derivative with
  `Ainv %*% s_i` and the numerical derivative of the complete weighted SOP
  functional with the production stacked influence. They also verify that the
  coefficient-mediated SOP change is materially nonzero relative to the
  profile-only derivative, so the cross-term is exercised. ORM agreed using
  ordinary weighted likelihood refits. VGLM required a `1e-2` symmetric weight
  perturbation and tighter optimizer convergence (`epsilon = 1e-10`) to keep
  refit noise below the validation tolerance. Both backends agree at relative
  tolerance `2e-3`, supporting the positive score sign, `n * bread` conversion
  from inverse total information to inverse per-patient sensitivity, and
  patient-level aggregation of repeated transition rows. The oracle uses full
  backend refits and weighted profile averaging rather than the production
  one-step score/Jacobian formula. No production formula change was required.
  `air format .`, the focused superpopulation tests, and the complete test suite
  passed. The factor-time vignette's comparison object was renamed consistently
  to `cmp_days` after the first check exposed its stale `cmp_days_point` name;
  the vignette then rendered successfully. A final `devtools::check()` including
  two complete vignette builds finished with 0 errors, 0 warnings, and the sole
  environment note that the current time could not be verified.

### ACI-03: Penalized ORM superpopulation inference

- **Status:** Resolved on 2026-07-19
- **Priority:** High
- **Current decision:** Empirical inference supports penalized ORM fits.
  Superpopulation inference deliberately rejects penalized ORM fits before
  allocating the analytical Jacobian. Supporting this combination is not part
  of the current analytical-CI contract.
- **Concern:** The reconstructed patient contributions are unpenalized
  likelihood scores, whereas the inverse sensitivity from a penalized fit may
  include curvature from the penalty. Combining those quantities does not
  establish a coherent stacked patient influence function. It would also omit
  uncertainty introduced by selecting or tuning the penalty.
- **Resolution approach:** Retain an explicit, informative early error for
  `vcov = "unconditional"` with a penalized ORM fit. Preserve empirical
  inference because it propagates the fitted coefficient covariance through
  the SOP Jacobian and does not require reconstructing the superpopulation
  score/sensitivity system.
- **Resolution log:** The rejection is an intentional scope decision, not a
  temporary fallback. Regression tests verify that empirical delta inference
  remains available for penalized ORM fits and that superpopulation inference
  fails before the potentially large Jacobian allocation. The
  `many-levels-previous-state-spline` vignette exercises the supported empirical
  path.

### ACI-04: Zero-score profile sensitivity scaling

- **Status:** Resolved on 2026-07-19 as an inseparable part of ACI-01
- **Priority:** High
- **Current decision:** A patient without any usable likelihood transition is
  excluded from the fitted-cohort estimand. No zero-score row or sensitivity
  rescaling is constructed.
- **Concern:** The previous implementation could include profile-only patients
  with synthetic zero scores and rescale sensitivity by
  `n_profiles / n_score_patients`, changing the fitted estimating-equation
  cohort.
- **Resolution approach:** Define the cohort by the fitted patient IDs, retain
  starting profiles before response-driven row omission, and require exact
  profile/score patient matching. A missing first outcome is allowed when a
  later usable transition exists.
- **Resolution log:** Resolved by the ACI-01 contract. Regression tests cover a
  missing first outcome, later likelihood contribution, profile-only exclusion,
  missing fitted-patient profiles, and exact score/profile matching.

### ACI-05: Stored ORM robust covariance identity and corrections

- **Status:** Resolved on 2026-07-22
- **Priority:** High
- **Current decision:** `orm_markov()` uses the package-internal `robcov_orm()`
  implementation, with the same independent HC1 `(n - 1) / (n - p)` and
  finite-cluster `G / (G - 1)` controls as VGLM. ORM row scores are analytic,
  case-weighted, and aggregated by patient; correction counts exclude zero-
  weight rows and clusters represented only by them. Penalized fits requesting
  sandwich penalty variance use the retained `var.from.info.matrix` inverse
  sensitivity as bread. Stored covariance is reused only
  when package provenance, covariance identity, and any requested cluster
  identity agree. Model-based SOP and diagnostic workflows accept only models
  created by `orm_markov()`, `vglm_markov()`, or `blrm_markov()`.
- **Concern:** Package ownership removes dependence on undocumented
  `rms::robcov()` object semantics, but correctness still depends on exact raw
  coefficient alignment, model-based bread recovery, case-weight handling,
  cluster alignment after omission, and preservation of wrapper provenance
  through internal refits.
- **Resolution approach:** Validate the analytic ORM scores and a hand-assembled
  clustered HC0 sandwich independently; verify exact HC1/cadjust ratios,
  weighting, cluster edge cases, covariance mutation detection, and explicit-
  cluster recomputation. Exercise the provenance gate across SOP and diagnostic
  entry points and verify bootstrap/FWB refits inherit provenance without
  unnecessarily computing robust covariance.
- **Resolution log:** Implementation began on 2026-07-22. Production ORM
  covariance now uses package-owned score, bread, sandwich-assembly, correction,
  and integrity metadata paths; VGLM shares the scalar correction helpers.
  Covariance consumers recompute ORM sandwiches without calling
  `rms::robcov()`, and internal refits propagate wrapper provenance. User-facing
  and architecture documentation has been updated. Independent review caught
  and resolved three material edge cases: the double-sandwich risk for penalized
  `var.penalty = "sandwich"` fits, the semantics of explicit HC/cadjust
  overrides when a stored covariance exists, and exclusion of zero-weight rows
  and zero-weight-only clusters from correction counts. Validation included
  targeted `air` formatting, `devtools::document()`, and focused ORM covariance,
  wrapper-provenance, bootstrap, and spline tests. The complete
  `devtools::test()` suite passed with one expected installed-package `callr`
  skip. The analytical-confidence-intervals, factor-time-orm-sops,
  full-po-sops, and many-levels-previous-state-spline vignettes rendered
  successfully. `devtools::check()` completed with 0 errors, 0 warnings, and 1
  environment-only note that the current time could not be verified.

### ACI-06: Unnamed penalized-ORM covariance relabeling

- **Status:** Resolved (2026-09-05)
- **Priority:** Medium
- **Current decision:** Internally produced ORM covariance matrices with missing
  dimnames are labeled in raw coefficient order. Known `Design$mmcolnames`
  aliases are mapped to displayed coefficient names. Explicit user covariance
  remains strictly named.
- **Concern:** This assumes `rms` preserves raw coefficient order when dimnames
  disappear. The assumption is tested for inline restricted cubic splines but
  not for every interaction, transformation, penalization structure, or future
  `rms` version.
- **Resolution approach:** Build a matrix of ORM fixtures covering supported
  transformations, interactions, and penalties. Compare backend model-matrix
  columns, coefficient order, covariance order, and finite-difference
  Jacobians. Replace inferred relabeling with stronger backend metadata if
  available.
- **Resolution log:** Retained production code after source-contract review of
  installed rms 8.1.1 and independent backend checks. `orm()` selects columns
  by `Design$mmcolnames` and relabels them with `Design$colnames`; `orm.fit()`
  uses `c(iname, xname)` for both coefficient names and information blocks.
  `infoMxop()` can lose labels when undoing predictor scaling, without changing
  that order. The metadata corroborates the current order; no additional
  relabeling heuristic or covariance reconstruction is needed.
  One compact regression in `tests/testthat/test-sops-delta-splines.R` reuses
  the existing fixture for spline-by-treatment interactions and log-time plus
  factor designs, each with `var.penalty = "simple"` and `"sandwich"`.
  It compares the actual normalized inverse sensitivity against a labeled,
  unscaled `rms::orm.fit()` information calculation, independently of the
  normalization helper. Coefficients agree at `1e-12`, confirming evaluation
  at the same fitted point, and covariance agrees at `1e-8`. Existing tests
  retain explicit missing-name/alias cases, strict user-matrix rejection, and
  finite-difference spline SOP Jacobians. The new test does not require rms to
  keep dropping labels in future releases. `air format` and
  `devtools::test(filter = "sops-delta-splines")` passed: 38 assertions,
  no failures, warnings, or skips. Commands used the repository's `MAKEFLAGS`
  and `LC_ALL=C` settings. Updated `ARCHITECTURE.md`; no public behavior or
  production implementation changed. A separately discovered zero-penalty
  input error is recorded under ACI-21.

### ACI-07: VGLM raw-to-effective constraint mapping

- **Status:** Resolved (2026-09-05)
- **Priority:** Medium
- **Current decision:** Construct the raw-to-effective map by concatenating
  VGAM constraint-matrix columns in term order and verify that the map
  reconstructs fitted effective coefficients using VGAM's native
  `coef(model, matrix = TRUE)` expansion.
- **Concern:** The contract is correct for tested full-PO cumulative fits but may
  not cover unusual zero constraints, special terms, offsets, or future VGAM
  coefficient-order changes. Reproduction at the fitted coefficient vector is
  strong but not a formal proof of the full derivative map.
- **Resolution approach:** Verify the map on basis vectors or multiple perturbed
  coefficient vectors, add constrained-formula fixtures, and document the exact
  VGAM constraint contract relied upon. Continue to fail explicitly outside
  verified structures.
- **Resolution log:** Installed VGAM 1.1.14 `coefvlm()` confirms raw coefficients
  follow cumulative constraint-column counts in term order. Kept the map
  algorithm and replaced the fitted reconstruction oracle: the old
  `get_effective_coefs()` repeated our own mapping algorithm, while the new
  check calls VGAM's native expansion. A compact basis-vector regression in
  `tests/testthat/test-sops-delta-vglm-map.R` checks every raw coefficient for
  spline-factor interactions, polynomial/reordered terms, equidistant
  thresholds, and a negative non-unit common-slope constraint basis. It also
  checks rejection of zero-column constraints. Existing offset, partial-PO,
  and unsupported-scope checks remain unchanged. No dependency or mapping
  abstraction was added. Focused map/core/superpopulation tests and the full
  integrated test suite passed.

### ACI-08: Structural boundary classification

- **Status:** Resolved (2026-09-05)
- **Priority:** Medium
- **Current decision:** Exact zero/one SOP estimates return warned `NA` logit
  limits regardless of their numerical standard error. Wald intervals are
  unchanged. No structural certainty is inferred without explicit provenance.
- **Concern:** The previous SE threshold was a numerical heuristic. Logistic saturation can resemble
  a structural boundary, while accumulated floating-point derivative noise can
  obscure a genuinely structural absorbing-state probability.
- **Resolution approach:** Remove the SE heuristic and test numerical saturation
  separately from absorbing transition structure. Require explicit provenance
  if a deterministic-SOP path is introduced later.
- **Resolution log:** Removed the tiny-SE heuristic after tracing native
  initialization and absorbing transitions. The current engine predicts the
  initial distribution from the model; absorbing transition rows carry mass
  forward, but the probability of reaching an absorbing state is estimated.
  These deterministic transition rows do not establish zero/one occupancy
  probabilities. The user's observation convention drops or omits `y` after
  absorption; no such rows need to be inserted. Automatic starting profiles
  are restricted to patients represented in fitted likelihood rows.
  Explicit structural metadata was not added because this inference path has
  no corresponding deterministic-SOP producer. If such a path is introduced,
  it must supply provenance rather than revive an SE heuristic.
  `tests/testthat/test-sops-delta-boundaries.R` uses valid nonabsorbing starting
  profiles, follows predicted absorption, and tests actual saturated individual
  and averaged SOPs alongside zero/tiny/noisy-SE boundary inputs and unchanged
  Wald calculations. Its 18 assertions and existing inference tests passed.
  Updated the analytical-CI vignette, NEWS, and architecture; the full integrated
  test suite passed.

### ACI-09: First-follow-up profiles and real-time integration conventions

- **Status:** Resolved on 2026-07-19
- **Priority:** Medium
- **Current decision:** Separate the first-follow-up profile-selection time from
  the observed baseline-state time and from the requested integration grid.
  The fitting wrappers need only `first_followup_time`; real-time summaries use
  `baseline_time`, `time_map`, and `target_times`. Linear interpolation from the
  observed baseline state to the first model-based SOP is retained, and
  trapezoidal integration begins at `min(target_times)`, not automatically at
  baseline.
- **Public API contract:** Replace `start_time` with
  `first_followup_time = NULL` in `orm_markov()`, `blrm_markov()`, and
  `vglm_markov()`. For numeric time, `NULL` resolves to 1; time 1 must exist and
  values below 1 are rejected. An explicit numeric value must exist and be the
  earliest scheduled follow-up value. Factor or character time requires an
  explicit matching level/value. These checks run when `id_var` enables
  automatic patient profiles; wrappers without `id_var` retain refit data but
  do not require a time column or store starting-profile metadata. Remove
  `origin_time` and `time_map` from the fitting wrappers. Backward-compatible
  aliases are not required, and legacy names captured by `...` must receive an
  informative error rather than leak to a backend fitter.
- **Starting-profile contract:** `first_followup_time` selects exactly one
  pre-response-omission row per fitted patient. ID, `yprev`, modeled predictors,
  and time metadata must be complete; the transition response may be missing.
  Every included patient must contribute at least one usable likelihood row
  somewhere. Patients with no likelihood row are excluded, fitted patients
  without a complete first-follow-up profile are errors, and explicit user
  `newdata` continues to replace automatic profiles. The resolved scalar is
  retained only for audit metadata and error messages; it does not control SOP
  recursion, interpolation, or integration.
- **Real-time contract:** Replace downstream `origin_time` with
  `baseline_time = 0` and remove the `origin` selector. A finite
  `baseline_time` adds the observed `yprev` distribution as the interpolation
  anchor; `baseline_time = NULL` disables it. The baseline time must precede the
  earliest mapped model-based SOP. `time_map` maps modeled factor visits to the
  continuous scale. `target_times` defines the returned interpolation grid and
  the integration interval, so no `integration_start_time` argument is added.
  For example, `baseline_time = 0`, a first mapped SOP at day 7, and
  `target_times = 1:28` interpolate from the observed day-0 distribution to the
  day-7 SOP but integrate only days 1 through 28.
- **Defaults when `target_times` is absent:** `interpolate_sops()` returns the
  baseline anchor and mapped SOP nodes when baseline anchoring is enabled.
  `time_in_state()` and time-in-state comparisons integrate the mapped
  follow-up nodes only, excluding the baseline interval by default. To begin
  integration before the first mapped follow-up, the user supplies the desired
  real-time grid; its minimum becomes the integration start. Including
  `baseline_time` explicitly includes the full baseline-to-follow-up interval
  in the AUC.
- **Preserved conventions:** Continue to average source contributions at
  duplicated mapped times, use linear interpolation only within the support
  bounded by `baseline_time` and the last mapped SOP, normalize probability
  mass as currently documented, and use trapezoidal integration on the sorted
  unique `target_times`.
- **Draw-specific baseline contract:** Coefficient-only simulation and posterior
  draws condition on the observed baseline profiles and therefore reuse the
  point-estimate anchor. Empirical score-bootstrap, fractional weighted
  bootstrap, and refit-bootstrap draws instead use the same patient weights or
  resampled cohort as each draw. Individual draw rows propagate their
  `score_weight` or `fwb_weight` to the baseline node as metadata, never as a
  trajectory identity key. Averaged and grouped draws retain only compact
  per-draw baseline-state distributions rather than full patient-by-draw weight
  matrices.
- **Implementation tranches:**
  1. Refactor wrapper formals, first-follow-up validation, stored attribute
     names, robust-wrapper metadata copying, and automatic profile resolution in
     `R/markov-model-data.R`, `R/vglm_helpers.R`, `R/robcov_vglm.R`, and
     `R/sops-api.R`.
  2. Rename and simplify baseline-anchor helpers and public arguments in
     `R/sops-interpolate.R` and `R/sops-time-in-state.R`; keep point estimates,
     stored draws, posterior draws, and interval recomputation aligned.
  3. Propagate `baseline_time` and the target-grid contract through
     `R/sops-comparisons.R`, comparison setup/reduction/replay helpers, bootstrap
     inference, and stored comparison metadata.
  4. Update the analytical comparison operator in
     `R/sops-delta-comparisons.R` so its baseline-node augmentation,
     interpolation matrix, trapezoidal weights, and point-estimate replay use
     the same grid. The shared baseline anchor must continue to cancel exactly
     in supported treatment differences.
  5. Rename internal `origin` data/metadata helpers to starting-profile terms so
     model metadata cannot be confused with the downstream baseline-time
     anchor.
  6. Regenerate roxygen documentation and update `NEWS.md`, `README.md`,
     `ARCHITECTURE.md`, the analytical-CI vignette, and the factor-time vignette.
- **Validation plan:** Add wrapper tests for numeric default 1, missing time 1,
  values below 1, explicit alternative numeric starts, mandatory factor starts,
  absent factor levels, duplicate/missing patient profiles, missing first
  outcomes with later likelihood rows, and all three automatic SOP/comparison
  consumers. Add hand-calculated interpolation/AUC tests showing a fixed day-0
  baseline anchor with `target_times = 1:28`, default exclusion of the baseline
  interval, explicit inclusion when requested, `baseline_time = NULL`, invalid
  anchor ordering, duplicate mapped times, and fixed `newdata`. Repeat these for
  point estimates, simulation/bootstrap/posterior draws, empirical delta
  inference, and supported superpopulation differences. Require the renamed API
  to reproduce current valid `origin_time = 0` estimates exactly when the same
  `target_times` are explicit; test the intentional no-`target_times` default
  change separately. Then run `air format .`, focused tests, the full suite,
  affected vignette renders, `devtools::document()`, and `devtools::check()`.
- **Resolution log:** Implemented the six tranches. Wrapper metadata now uses
  `markov_starting_profile_*` names; automatic `sops()`, `avg_sops()`, and
  `avg_comparisons()` share the validated profiles. Downstream point, draw,
  bootstrap, posterior, and delta paths use `baseline_time`; target grids control
  AUC bounds, and unused `time_map` entries do not affect baseline ordering or
  analytical weights. Focused wrapper, factor-time, comparison, delta, and
  inference suites passed. The complete test suite passed with one expected
  installed-package skip. The analytical-CI, factor-time ORM, and many-level
  previous-state spline vignettes rendered successfully. `devtools::check()`
  completed with 0 errors, 0 warnings, and one environment note because current
  network time could not be verified. Final review then corrected draw-level
  baseline handling: individual weighted draws now retain their weights through
  pre-follow-up interpolation, while averaged/grouped score-bootstrap, FWB, and
  refit-bootstrap draws store baseline-state anchors computed from the matching
  draw weights or resampled cohort. Direct `time_in_state()` summaries now
  integrate those stored draws and recompute draw-based interval and standard
  error columns on the requested AUC grid.

### ACI-10: Memory accounting versus process memory

- **Status:** Open
- **Priority:** Medium
- **Current decision:** `markov.misc.delta_max_bytes` guards principal numeric
  outputs, rolling recursion workspace, retained state, and requested
  Jacobian/covariance blocks.
- **Concern:** It is not a hard cap on total R process memory. Existing design
  matrices, R copies, object headers, allocator overhead, and C++ container
  capacity are not all included. The estimate may therefore differ from peak
  resident memory.
- **Resolution approach:** Profile peak memory for representative small, full-PO,
  and many-level workflows. Reconcile measured peaks with each accounted
  component, add a safety factor if useful, and keep documentation explicit
  about what the option does and does not guarantee.
- **Resolution log:** Pending.

### ACI-11: Grouped native execution row-layout contract

- **Status:** Resolved (2026-09-05)
- **Priority:** Medium
- **Current decision:** Counterfactual profiles are averaged in contiguous,
  equal-sized scenario blocks produced by `avg_sops()`. The R entry point now
  validates grid order and identical starting-profile order in every block
  before compiling the native plan.
- **Concern:** The optimized native kernel is coupled to an internal row-ordering
  contract. Point-estimate replay alone can miss patient permutations that
  preserve averages but misalign patient influence contributions.
- **Resolution approach:** Validate, document, and test the existing layout
  contract directly, without adding a redundant native grouping index.
- **Resolution log:** Traced `create_counterfactual_data()` through native
  block averaging and patient-influence reduction. Added one inline check of
  scenario columns against repeated grid values and all other columns against
  the first profile block. Row selection preserves column classes/dimensions,
  including matrix-valued metadata. This rejects scenario interleaving, swapped
  blocks, and differing patient permutations before computation, even where
  averages alone cannot reveal incorrect influence alignment. A shared profile
  permutation remains valid. Kept the existing native interface instead of
  adding redundant group indices or reordering silently.
  `tests/testthat/test-sops-delta-layout.R` verifies grouped-versus-ungrouped
  probabilities/Jacobians, three invalid layouts, and common profile reversal,
  including matrix-valued metadata. All six assertions and the full integrated
  test suite passed.

### ACI-12: Superpopulation `get_jacobian()` semantics

- **Status:** Resolved (2026-09-05)
- **Priority:** Medium
- **Current decision:** `get_jacobian()` is an unexported inspection helper.
  For a superpopulation target, it returns `G`,
  the average derivative with respect to model coefficients. Profile-distribution
  variation remains available only through the stored influence representation
  and `vcov()`.
- **Concern:** The name `get_jacobian()` may suggest that the returned matrix
  fully represents superpopulation uncertainty, which it does not.
- **Resolution approach:** Assess likely downstream usage. Either retain the
  contract with stronger naming/documentation and metadata, add a separate
  influence accessor, or reject `get_jacobian()` for superpopulation targets if
  the partial interpretation is too easy to misuse.
- **Resolution log:** Adopted the user's internal-only contract rather than
  extending the public API. Removed the export, public help topic, README
  example, and public help links. Kept the vignette example as
  `markov.misc:::get_jacobian()` with an explicit internal-inspection note and
  explanation that coefficient derivatives alone do not describe unconditional
  uncertainty. `inferences()` already uses analytical state directly and does
  not call this accessor; public `vcov()` behavior and internal tests remain
  unchanged. Updated architecture documentation and regenerated roxygen files.
  `air format .`, focused delta comparison/inference tests (146 assertions,
  including the namespace-export regression), and rendering the analytical-CI
  vignette passed with the repository's R environment settings and the installed
  Quarto Pandoc path.

### ACI-13: `vcov()` dispatch for non-delta results

- **Status:** Open
- **Priority:** Medium
- **Current decision:** S3 `vcov()` methods are registered for SOP and average
  comparison classes and require an attached analytical state.
- **Concern:** Calling `vcov()` on MVN/bootstrap results of the same classes now
  reaches the new method and errors. A delta-specific accessor would have
  avoided changing generic dispatch for non-delta objects.
- **Resolution approach:** Characterize prior and current dispatch behavior.
  Decide whether `vcov()` should support draw-based results, delegate when no
  analytical state exists, or be replaced/supplemented by a dedicated accessor.
  Add compatibility tests for every inference method.
- **Resolution log:** Pending.

### ACI-14: Dense comparison and covariance materialization

- **Status:** Resolved (2026-09-05)
- **Priority:** Medium
- **Current decision:** Store nonzero source indices and weights per comparison
  and apply them directly to Jacobian rows or patient-influence columns.
  Accumulate real-time interpolation/integration weights with the existing
  interpolation plan and base `rowsum()`. `vcov(result, rows = NULL)` still
  explicitly requests a full dense covariance subject to the memory guard;
  selected `rows` avoid the full output allocation.
- **Concern:** Large state/time/comparison grids may hit the guard or allocate
  avoidably large matrices even though the operators are sparse and structured.
- **Resolution approach:** Benchmark realistic large grids and evaluate sparse
  matrices, direct index/weight propagation, chunked covariance blocks, or
  requiring explicit `rows` above a threshold. Preserve identical output.
- **Research and alternatives (2026-09-05):** Confirmed installed
  `marginaleffects` 1.0.0 and `rms` 8.1.1. The installed marginaleffects
  `NEWS.md` documents O(np) `~meandev`/`~meanotherdev` transformations instead
  of dense n-by-n contrasts. Its internal `hypothesis_formula_pullback()`
  applies group means/sums directly to Jacobians; it is not an exported general
  SOP comparison engine. The online NEWS page was behind the installed release,
  so the installed changelog and namespace were used as the version-specific
  evidence. The [rms release notes](https://hbiostat.org/r/rms/) and installed
  `infoMxop()` confirm sparse information matrices and on-demand inverse
  operations, which address coefficient-level calculations rather than our
  result-by-SOP operator. Neither dependency automatically fixes this package's
  allocation. Adopted the same matrix-free principle using base R and existing
  interpolation infrastructure. A Matrix sparse operator would also remove
  zeros, but is unnecessary for these weighted selections and would promote an
  optional dependency to a runtime requirement. Chunking a requested full
  covariance cannot eliminate its quadratic output size; retained the existing
  row-selection API instead of adding another accessor or changing defaults.
- **Benchmark evidence:** `benchmarks/benchmark-delta-comparisons.R` checks
  dense-product agreement at tolerance `1e-12` and compares 20/100/250 visits,
  four states, three treatment scenarios, SOP and real-time time-in-state
  contrasts, 12 coefficient columns, and 100 patient-influence rows. For the
  250-visit SOP case (2,000 results, 3,000 source cells), operator object size
  fell from 45.78 MiB to 0.916 MiB. Median propagation times over three runs
  were 0.07 versus 0.01 seconds for Jacobians and 0.49 versus 0.01 seconds for
  influences on this Windows machine. These exclude fitting, recursion, and
  operator compilation, and are not peak-process-memory measurements. Small
  timings often fall below the timer's resolution. Comparison compilation still
  scans source cells per result; this change removes quadratic operator storage,
  not that setup-time cost.
- **Resolution log:** Implemented and validated. Regression coverage includes both propagation orientations,
  zero-weight rows, single-column inputs, duplicate mapped visits, baseline
  anchoring, a single integration node, and guarded output allocation. A
  4,400-result / 8,800-source grid now compiles and replays under a 2 MiB
  numeric-allocation limit although its old dense operator requires 295.41 MiB.
  Large full covariance requests still fail the guard while selected blocks
  work for both coefficient and influence representations.
  `air format .`, `devtools::document()`, the focused comparison/inference
  tests, the complete `devtools::test(reporter = "summary")` suite, and the
  benchmark passed. The full suite had only the expected installed-package
  worker-invariance skip. `devtools::check(args = "--no-tests")` (tests already
  run separately) passed with 0 errors, 0 warnings, and one environment note
  because current network time could not be verified. Both vignette builds
  and examples passed. The first build attempt could not locate Pandoc;
  supplying the existing `C:/Program Files/Quarto/bin/tools` through
  `RSTUDIO_PANDOC` resolved that environment issue. R commands used the
  repository's required `MAKEFLAGS` and `LC_ALL=C` settings.

### ACI-15: Numerical validation tolerances and portability

- **Status:** Open
- **Priority:** Medium
- **Current decision:** Use strict absolute point-replay tolerance `1e-12`, raw
  crossed-probability tolerance `-1e-10`, derivative/mass tolerances, and
  scale-dependent covariance eigenvalue tolerances.
- **Concern:** The point-replay tolerance may be unnecessarily strict across
  platforms, BLAS implementations, or unusually scaled models. PSD tolerances
  may also become permissive or restrictive as coefficient dimension grows.
- **Resolution approach:** Run the focused suite on supported operating systems
  and R/backend versions. Test extreme but valid scales, then define tolerances
  using explicit absolute-plus-relative formulas tied to the numerical quantity
  being validated.
- **Resolution log:** Pending.

### ACI-16: Native C++ maintenance and semantic parity

- **Status:** Open
- **Priority:** Medium
- **Current decision:** Native cpp11 recursion is the sole production analytical
  path. The slower R recursion is a test-only oracle and never a fallback.
- **Concern:** The native implementation materially improves speed but adds a
  second-language maintenance surface and the possibility that later R and C++
  execution semantics diverge.
- **Resolution approach:** Keep shared deterministic fixtures and native-versus-R
  oracle tests for every supported feature. Document invariants next to the
  native interface, add sanitizer or compiled-code checks where practical, and
  require oracle updates with any recursion change.
- **Resolution log:** Pending.

### ACI-17: Performance benchmark generalizability

- **Status:** Open
- **Priority:** Low
- **Current decision:** Report measured speedups from the exact full-PO vignette
  workflow with 100 MVN draws on the development machine.
- **Concern:** The measured VGLM and ORM speedups do not imply the same advantage
  for small models, very few draws, different coefficient dimensions, or
  workloads dominated by execution-plan compilation.
- **Resolution approach:** Add a reproducible benchmark grid varying patients,
  states, visits, coefficients, and MVN draws. Report absolute times and
  break-even regions rather than a single general speedup claim.
- **Resolution log:** Pending.

### ACI-18: Independent superpopulation-inference validation oracle

- **Status:** Open
- **Priority:** High
- **Current decision:** Validate recursive Jacobians against central finite
  differences and superpopulation covariance against hand calculations of the
  approved stacked influence formula. Do not require close MVN/bootstrap
  agreement.
- **Concern:** Finite differences validate `G`, but the superpopulation hand tests
  reproduce the same score/sensitivity formula as production code. This leaves
  no fully independent oracle for the complete stacked influence calculation.
- **Resolution approach:** Implement a small, deliberately slow test-only
  estimating-equation perturbation or symbolic example that is independent of
  production helpers. Use bootstrap or simulation only for qualitative
  diagnostics, not as a close-agreement acceptance criterion.
- **Resolution log:** Pending.

### ACI-19: Public target terminology and formal defaults

- **Status:** Resolved (2026-09-05)
- **Decision:** Keep `method = "delta"`; replace the public `target` argument
  with `vcov = "conditional"` or `"unconditional"`. Averaged results default to
  unconditional, individual results to conditional. Custom coefficient matrices
  select conditional inference. Internal target metadata labels remain unchanged.
- **Compatibility:** Remove `target` immediately. Existing positional arguments
  retain their order. Unsupported unconditional inference errors without fallback.

### ACI-20: Generated native build artifacts in the worktree

- **Status:** Open
- **Priority:** Low
- **Current decision:** Native object and DLL files are generated locally during
  package loading, documentation, testing, or checking and are not committed.
- **Concern:** `src/cpp11.o`, `src/sops.o`, and `src/markov.misc.dll` are currently
  untracked. They can obscure a clean status or be staged accidentally.
- **Resolution approach:** Confirm the package's intended cleanup and ignore
  policy, remove only verified generated artifacts, and add a repeatable clean
  validation step before commits without hiding source or meaningful outputs.
- **Resolution log:** Pending; the files are not part of commit `76c18c3`.

### ACI-21: Explicit zero ORM penalty list handling

- **Status:** Resolved (2026-09-06)
- **Priority:** Medium
- **Concern:** With an explicit `penalty = 0`, rms stores penalty settings as
  a list while its penalty matrix contains only zeros. `orm_model_bread()` then
  reaches `is.finite(penalty)` and errors because the input is a list. Omitting
  the penalty works; a positive penalty avoids this branch via its nonzero
  matrix. This is separate from covariance ordering.
- **Evidence:** Reproduced during ACI-06 experiments using a fitted
  `rms::orm(ordered(y) ~ rms::rcs(time, 3) * tx + yprev, penalty = 0, ...)`
  and `orm_model_bread()`; the error is `default method not implemented for
  type 'list'`.
- **Resolution approach:** Check every penalty-detection caller, handle the
  backend's penalty representation consistently, and add a small explicit-zero
  regression. Preserve rejection of genuinely penalized unconditional inference.
- **Resolution log:** Flatten stored penalty settings with base `unlist()` in
  both `orm_model_bread()` and `delta_reject_penalized_orm_superpopulation()`.
  A regression fits actual wrapper models with omitted and explicit zero
  penalties under both penalty-variance modes, then compares model bread,
  estimates, and conditional/unconditional SOP covariance. The superpopulation,
  spline, and ORM robust-covariance test files passed, including existing
  rejection of positive penalties. No dependency or helper was added.

## Resolved Implementation Decisions

These original implementation problems were corrected before the initial
analytical-CI commit. They remain here because related residual concerns are
tracked above.

### RID-01: Full individual Jacobian retention and conservative preflight

- **Status:** Resolved on 2026-07-18.
- **Original choice:** Retain every individual Jacobian and count every origin's
  derivative workspace as simultaneous, causing a 300,960,000-byte estimate in
  the many-level empirical workflow.
- **Resolution:** Average empirical Jacobians inside native recursion; retain
  individual probabilities only for the superpopulation target; count rolling
  workspace rather than mutually exclusive origin workspaces.
- **Evidence:** The full many-level spline vignette renders under the default
  256 MiB limit. Grouped results match full individual recursion and the test-only
  R oracle.
- **Residual links:** ACI-10 and ACI-11.

### RID-02: Strict names for backend-owned ORM covariance

- **Status:** Resolved on 2026-07-18.
- **Original choice:** Require every covariance matrix, including internally
  generated penalized-ORM matrices, to carry complete raw-coefficient dimnames.
- **Resolution:** Normalize only matrices produced internally from the exact ORM
  fit, using preserved order and recognized design aliases. Continue to require
  strict names for explicit user covariance.
- **Evidence:** Penalized restricted-cubic-spline fixtures and the full
  many-level vignette pass.
- **Residual link:** ACI-06.

### RID-03: `target` inserted into the existing positional API

- **Status:** Resolved on 2026-07-18.
- **Original choice:** Insert `target` before existing `null` and draw-related
  arguments, which could reinterpret positional calls.
- **Resolution:** Append `target` after all pre-existing formals and add a test
  locking their established order.
- **Evidence:** Focused tests and the complete package check pass.
- **Residual link:** ACI-19.

### RID-04: R recursion as the production path

- **Status:** Resolved on 2026-07-18.
- **Original choice:** Use the analytical R recursion in production; profiling
  showed less improvement over MVN simulation than expected.
- **Resolution:** Implement the recursive Jacobian in native cpp11 code and keep
  the R version as a test-only oracle, with no silent fallback.
- **Evidence:** Native and R Jacobians agree below `1e-12` in focused fixtures.
  On the exact 100-draw full-PO vignette workflow, measured delta speedups were
  approximately 4.1x for VGLM and 5.8x for ORM on the development machine.
- **Residual links:** ACI-16 and ACI-17.

## Concern Work Log

Add dated entries here whenever a concern becomes active or changes status.
Link the stable concern ID to relevant commits, tests, documents, and GitHub
issues.

### 2026-07-19

- Added ACI-01 through ACI-20 after a post-implementation confidence audit.
- Recorded RID-01 through RID-04 so corrected decisions and their evidence are
  retained rather than lost.
- Checkpointed the concern register in commit `7b04ff6`.
- Marked ACI-01 **Active**. Confirmed that the manuscript estimator remains
  `stats::cov(phi) / n`, renamed the public target to `"superpopulation"`
  without a `"population"` compatibility alias, and adopted the fitted-cohort
  and preserved-origin-profile contract recorded under ACI-01.
- Resolved ACI-01 after separating likelihood, refit, starting-profile, and user
  prediction data; enforcing exact fitted-patient profile/score matching; and
  completing focused tests, the full suite, three affected vignette renders,
  roxygen regeneration, and the package check. ACI-04 was resolved as an
  inseparable consequence because the approved cohort contract removes the
  zero-score-patient pathway entirely.
- Resolved ACI-03 by making rejection of penalized ORM superpopulation inference
  an intentional scope decision. The reconstructed likelihood scores and the
  penalty-aware sensitivity do not currently define a verified stacked patient
  influence function, and penalty-selection uncertainty is not represented.
  Empirical delta inference for penalized ORM fits remains supported.
- Marked ACI-09 **Active** and approved the first-follow-up/baseline-time
  separation. `first_followup_time` will only select and validate automatic
  profiles; `baseline_time`, `time_map`, and `target_times` will define observed
  baseline anchoring, real-time interpolation, and the integration grid without
  a separate integration-start argument.
- Resolved ACI-09 after implementing the wrapper/profile rename, observed
  baseline anchoring, target-grid AUC defaults, comparison/delta propagation,
  documentation, and regression coverage. Focused and complete tests, three
  affected vignette renders, roxygen regeneration, and the package check passed.
- Corrected ACI-09 during final review so empirical resampling draws carry
  draw-specific baseline anchors. Individual score/FWB weights remain metadata
  rather than interpolation keys; averaged and grouped paths retain compact
  baseline-state distributions computed from each draw's weights or resampled
  cohort. Direct time-in-state summaries reduce the corrected draws and report
  the resulting AUC uncertainty.

### 2026-09-05

- Marked ACI-14 **Active** and investigated the installed marginaleffects 1.0.0
  and rms 8.1.1 implementations. Replaced dense comparison/interpolation
  operators with direct index/weight propagation, added regression coverage
  and a reproducible benchmark, and documented the retained explicit full
  covariance limit. Resolved ACI-14 after full tests passed with the expected
  installed-package skip and package checking passed with 0 errors, 0 warnings,
  and the sole network-time verification note, including both vignette builds.

- Marked ACI-12 **Active**, then **Resolved** under the user-approved
  internal-only accessor contract. Removed the `get_jacobian()` export and
  public help, retained explicit `:::` inspection in the vignette, and verified
  the namespace boundary, focused delta tests, and analytical-CI vignette render.
- Committed the ACI-12 changes in `7a8a86d`, then marked ACI-06 **Active**.
  Resolved ACI-06 through installed-backend source review and a compact
  independent covariance-order regression, with no production changes.
  Recorded the separately reproduced explicit-zero-penalty error as ACI-21;
  ACI-07 remains the next original open concern.
- Reviewed ACI-07, ACI-08, and ACI-11 in parallel using subagents at the
  user's request, then resolved them with native VGAM coefficient validation,
  removal of the boundary-SE heuristic, and explicit scenario/profile-order
  validation. Incorporated the observed-data convention that rows after
  absorption are omitted. Full tests passed with the expected installed-only skip.
  Package checking (tests run separately) passed with 0 errors, 0 warnings,
  and the sole network-time verification note, including both vignette builds.

### 2026-09-06

- Committed the ACI-06/07/08/11 changes in `af9ed07`, then resolved ACI-21
  with two penalty-list normalization changes and an end-to-end regression.
  Focused tests passed; no vignettes were built.
