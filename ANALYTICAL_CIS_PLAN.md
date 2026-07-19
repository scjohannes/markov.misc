# Analytical Confidence Intervals Implementation Plan

Status date: 2026-07-19
Working branch: `analytical-cis`

This is the durable implementation checklist for analytical confidence
intervals for `sops()`, `avg_sops()`, and supported `avg_comparisons()` results.
It records the approved design, what is currently present in the branch, and
the post-implementation concerns that will be investigated one at a time.

## Approved Design

- Add deterministic inference through
  `inferences(method = "delta", target = ..., conf_type = "auto")` without
  changing the existing `method = "mvn"` default.
- Differentiate the first-order full proportional-odds SOP recursion
  analytically on each backend's complete raw coefficient scale. Central finite
  differences are a test oracle only, not a production fallback.
- Use patient-cluster robust coefficient covariance for fixed-profile and
  empirical-cohort propagation unless the user supplies a valid complete
  coefficient covariance.
- For averaged objects, treat `target = "empirical"` as conditional on the
  observed standardization profiles. Treat `target = "population"` as a
  same-cohort population target and use a stacked patient influence function
  that retains profile/model-score covariance. Individual `sops()` accepts only
  an omitted target or `target = "fixed"`.
- Store factorized analytical state and materialize selected Jacobian or
  covariance blocks with `get_jacobian()` and `vcov()`.
- Support linear differences for SOP and time-in-state estimands, including
  visit-to-real-time mapping, empirical-baseline origin handling, linear
  interpolation, and trapezoidal integration.
- Fail explicitly outside the approved scope rather than silently switching to
  numerical differentiation, row-level clusters, or another inferential target.

## Target and Interval Contract

| Result | Default / allowed target | Analytical interval |
| --- | --- | --- |
| Individual `sops()` | Omitted or `fixed` only | Componentwise logit-delta under `conf_type = "auto"` |
| `avg_sops()` | `empirical` | Componentwise logit-delta under `conf_type = "auto"` |
| Supported `avg_comparisons()` | `empirical` | Identity-scale Wald under `conf_type = "auto"` |
| Stored-cohort averages/comparisons | explicit `population` | Same-cohort stacked influence covariance |

Patient is always the independent cluster. An explicit row-aligned `cluster`
vector takes precedence; otherwise stored fitting data and stored `id_var`
metadata are required. Cluster robustness protects variance calculations from
arbitrary within-patient score correlation. It does not correct transition-model
bias, informative missingness, or Markov/proportional-odds misspecification.

## Implementation Checklist

### Public Entry and Dispatch

- [x] Add `"delta"` to the `inferences()` method choices.
- [x] Add `target` and `conf_type = "auto"` routing while retaining
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

### Covariance and Population Targets

- [x] Validate complete named, finite, symmetric, positive-semidefinite custom
  covariance matrices for coefficient-form targets.
- [x] Resolve or compute patient-cluster robust covariance while preserving
  VGLM finite-cluster metadata and the ORM HC0 convention.
- [x] Aggregate raw transition-row scores by patient and construct inverse
  per-patient sensitivity from inverse total information, not from the robust
  sandwich covariance.
- [x] Align same-cohort profile and score IDs exactly, retain zero-score
  profiles, reject score-only patients, and scale sensitivity to the full
  profile cohort.
- [x] Form vector stacked influences and use
  `stats::cov(influence) / n` with the finite-sample convention reported in
  metadata.
- [x] Reject custom covariance and user-supplied external profiles for the
  population target.

### Result Access and Comparisons

- [x] Store coefficient-form analytical state as `J` plus coefficient `V`.
- [x] Store population-form analytical state as average `J` plus the patient
  influence matrix.
- [x] Add row-selective `get_jacobian()` and S3 `vcov()` methods without storing
  a dense all-cell covariance.
- [x] Implement linear operators for SOP and time-in-state differences.
- [x] Reproduce point estimates before propagating a comparison operator.
- [x] Support factor-time real-time differences through the stored `time_map`,
  origin, interpolation grid, and trapezoidal weights.

### Tests and Documentation

- [x] Add focused core tests for ORM/VGLM derivatives, absorbing states, factor
  visits, validation conditions, and allocation guards.
- [x] Add focused population tests for cluster resolution, covariance
  conventions, score scaling/alignment, zero-score profiles, and the stacked
  cross term. The population test file passed 51 assertions with no warnings in
  the implementation session.
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
- Individual SOPs, empirical average SOPs, and same-cohort population average
  SOPs.
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
| ACI-01 | High | Open | Population finite-sample normalization |
| ACI-02 | High | Open | Population score orientation and cross-term scaling |
| ACI-03 | High | Open | Penalized ORM population inference |
| ACI-04 | High | Open | Zero-score profile sensitivity scaling |
| ACI-05 | High | Open | Stored ORM robust covariance identity and correction metadata |
| ACI-06 | Medium | Open | Unnamed penalized-ORM covariance relabeling |
| ACI-07 | Medium | Open | VGLM raw-to-effective constraint mapping |
| ACI-08 | Medium | Open | Structural boundary classification |
| ACI-09 | Medium | Open | Factor-time interpolation and integration conventions |
| ACI-10 | Medium | Open | Memory accounting versus actual process memory |
| ACI-11 | Medium | Open | Grouped native execution row-layout contract |
| ACI-12 | Medium | Open | Population `get_jacobian()` semantics |
| ACI-13 | Medium | Open | `vcov()` dispatch for non-delta result objects |
| ACI-14 | Medium | Open | Dense comparison and covariance materialization |
| ACI-15 | Medium | Open | Numerical validation tolerances and portability |
| ACI-16 | Medium | Open | Native C++ maintenance and semantic parity |
| ACI-17 | Low | Open | Performance benchmark generalizability |
| ACI-18 | High | Open | Independent population-inference validation oracle |
| ACI-19 | Low | Open | Public target terminology and formal defaults |
| ACI-20 | Low | Open | Generated native build artifacts in the worktree |

### ACI-01: Population finite-sample normalization

- **Status:** Open
- **Priority:** High
- **Current decision:** Population covariance is
  `stats::cov(influence) / n`. Because `stats::cov()` uses denominator `n - 1`,
  this is an HC1-like finite-sample convention.
- **Concern:** This follows the approved formula but is not the raw empirical
  influence covariance based on denominator `n^2`. It also does not literally
  preserve every empirical backend's HC/cadjust convention for the population
  target.
- **Resolution approach:** Derive the target covariance and finite-sample
  correction from the manuscript's estimating equations. Compare
  `crossprod(phi) / n^2`, `stats::cov(phi) / n`, and backend correction factors
  in hand-calculated clustered examples. Decide whether population correction
  should be fixed, backend-specific, or user-selectable.
- **Resolution log:** Pending.

### ACI-02: Population score orientation and cross-term scaling

- **Status:** Open
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
- **Resolution log:** Pending.

### ACI-03: Penalized ORM population inference

- **Status:** Open
- **Priority:** High
- **Current decision:** Empirical inference supports penalized ORM fits.
  Population inference is not currently rejected solely because an ORM fit is
  penalized.
- **Concern:** The inverse sensitivity may reflect the penalty while reconstructed
  patient scores are likelihood scores. Uncertainty from selecting or tuning
  the penalty is not represented. The many-level spline vignette exercises the
  empirical target, not this population contract.
- **Resolution approach:** Derive the influence function with fixed penalty and
  determine whether the penalty contribution cancels after centering. Separately
  decide whether penalty-selection uncertainty is in scope. Until the derivation
  is confirmed, consider explicitly rejecting penalized ORM population
  inference while retaining empirical support.
- **Resolution log:** Pending.

### ACI-04: Zero-score profile sensitivity scaling

- **Status:** Open
- **Priority:** High
- **Current decision:** Baseline-profile patients without likelihood rows are
  retained with zero scores, and inverse per-profile sensitivity is rescaled by
  `n_profiles / n_score_patients`.
- **Concern:** This is appropriate if those patients are cohort members with a
  zero contribution to the study-cohort estimating equation. It may be
  inappropriate if their absent transition rows reflect selection, informative
  missingness, or a different risk-set definition.
- **Resolution approach:** State the sampling and missingness assumptions
  explicitly, derive the estimating equation when patients have no transition
  rows, and test fixtures with known score/profile cohort relationships. Decide
  whether additional metadata or a stricter error is needed.
- **Resolution log:** Pending.

### ACI-05: Stored ORM robust covariance identity and corrections

- **Status:** Open
- **Priority:** High
- **Current decision:** An ORM object containing both `var` and `orig.var` is
  treated as carrying a stored `rms::robcov()` sandwich. It is reported as HC0
  without cluster adjustment and reused when no explicit cluster is supplied.
- **Concern:** This is reliable for `orm_markov()` fits but does not independently
  prove that an arbitrary external ORM object's stored covariance uses the same
  patient IDs or finite-sample convention resolved for analytical inference.
- **Resolution approach:** Audit `rms::robcov()` source and stored metadata across
  supported `rms` versions. Add provenance and cluster-identity checks where
  possible; otherwise restrict automatic reuse to package-created fits and
  recompute or require explicit input for external fits.
- **Resolution log:** Pending.

### ACI-06: Unnamed penalized-ORM covariance relabeling

- **Status:** Open
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
- **Resolution log:** Pending.

### ACI-07: VGLM raw-to-effective constraint mapping

- **Status:** Open
- **Priority:** Medium
- **Current decision:** Construct the raw-to-effective map by concatenating
  VGAM constraint-matrix columns in term order and verify that the map
  reconstructs fitted effective coefficients.
- **Concern:** The contract is correct for tested full-PO cumulative fits but may
  not cover unusual zero constraints, special terms, offsets, or future VGAM
  coefficient-order changes. Reproduction at the fitted coefficient vector is
  strong but not a formal proof of the full derivative map.
- **Resolution approach:** Verify the map on basis vectors or multiple perturbed
  coefficient vectors, add constrained-formula fixtures, and document the exact
  VGAM constraint contract relied upon. Continue to fail explicitly outside
  verified structures.
- **Resolution log:** Pending.

### ACI-08: Structural boundary classification

- **Status:** Open
- **Priority:** Medium
- **Current decision:** An SOP estimate exactly zero or one is treated as
  structural when its delta standard error is at most
  `sqrt(.Machine$double.eps)`; otherwise logit limits are returned as `NA` with a
  warning.
- **Concern:** This is a numerical heuristic. Logistic saturation can resemble
  a structural boundary, while accumulated floating-point derivative noise can
  obscure a genuinely structural absorbing-state probability.
- **Resolution approach:** Propagate explicit structural-state metadata from the
  SOP recursion and use it instead of an SE threshold. Add tests separating
  absorbing-state structure, zero-probability design structure, and numerical
  saturation.
- **Resolution log:** Pending.

### ACI-09: Factor-time interpolation and integration conventions

- **Status:** Open
- **Priority:** Medium
- **Current decision:** Use linear interpolation, average source contributions
  at duplicated mapped times, anchor empirical-origin contrasts at zero, and
  apply trapezoidal integration on requested real times.
- **Concern:** The implementation reproduces current vignette point estimates,
  but duplicate-time collapse, empirical-origin anchoring, and interpolation
  outside explicit SOP nodes are estimand conventions that should be confirmed
  independently of code replay.
- **Resolution approach:** Write the estimand mathematically for regular and
  irregular factor times. Work through small hand examples with duplicated
  mappings and origins between visits. Confirm or revise the existing point
  estimator and analytical operator together.
- **Resolution log:** Pending.

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

- **Status:** Open
- **Priority:** Medium
- **Current decision:** Counterfactual profiles are averaged in contiguous,
  equal-sized scenario blocks produced by the current `avg_sops()` workflow.
- **Concern:** The optimized native kernel is coupled to an internal row-ordering
  contract. Strict point-estimate replay should detect a changed layout, but the
  contract is not represented by an explicit grouping index.
- **Resolution approach:** Document and test the layout contract directly.
  Consider passing an explicit scenario/group index to native code so correct
  grouping no longer depends on contiguity, then benchmark the cost.
- **Resolution log:** Pending.

### ACI-12: Population `get_jacobian()` semantics

- **Status:** Open
- **Priority:** Medium
- **Current decision:** For a population target, `get_jacobian()` returns `G`,
  the average derivative with respect to model coefficients. Profile-distribution
  variation remains available only through the stored influence representation
  and `vcov()`.
- **Concern:** The name `get_jacobian()` may suggest that the returned matrix
  fully represents population uncertainty, which it does not.
- **Resolution approach:** Assess likely downstream usage. Either retain the
  contract with stronger naming/documentation and metadata, add a separate
  influence accessor, or reject `get_jacobian()` for population targets if the
  partial interpretation is too easy to misuse.
- **Resolution log:** Pending.

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

- **Status:** Open
- **Priority:** Medium
- **Current decision:** Analytical comparisons build a dense
  result-by-SOP-cell operator. `vcov(result, rows = NULL)` materializes the full
  selected-result covariance subject to the memory guard.
- **Concern:** Large state/time/comparison grids may hit the guard or allocate
  avoidably large matrices even though the operators are sparse and structured.
- **Resolution approach:** Benchmark realistic large grids and evaluate sparse
  matrices, direct index/weight propagation, chunked covariance blocks, or
  requiring explicit `rows` above a threshold. Preserve identical output.
- **Resolution log:** Pending.

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

### ACI-18: Independent population-inference validation oracle

- **Status:** Open
- **Priority:** High
- **Current decision:** Validate recursive Jacobians against central finite
  differences and population covariance against hand calculations of the
  approved stacked influence formula. Do not require close MVN/bootstrap
  agreement.
- **Concern:** Finite differences validate `G`, but the population hand tests
  reproduce the same score/sensitivity formula as production code. This leaves
  no fully independent oracle for the complete stacked influence calculation.
- **Resolution approach:** Implement a small, deliberately slow test-only
  estimating-equation perturbation or symbolic example that is independent of
  production helpers. Use bootstrap or simulation only for qualitative
  diagnostics, not as a close-agreement acceptance criterion.
- **Resolution log:** Pending.

### ACI-19: Public target terminology and formal defaults

- **Status:** Open
- **Priority:** Low
- **Current decision:** Use `fixed` for individual `sops()` and
  `empirical`/`population` for averaged results. Use `conf_type = "auto"` as the
  formal default while preserving percentile behavior for draw methods.
- **Concern:** Users may naturally describe fixed observed `newdata` as
  empirical. Code that inspects `formals(inferences)` also observes a changed
  default even though runtime draw behavior is preserved.
- **Resolution approach:** Collect representative call patterns, review naming
  against related marginal-estimand APIs, and decide whether aliases or clearer
  error messages are warranted. Preserve existing positional argument order.
- **Resolution log:** Pending.

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
  individual probabilities only for the population target; count rolling
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
- No concern has been selected as **Active** yet.
