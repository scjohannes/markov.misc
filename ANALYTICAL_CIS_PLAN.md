# Analytical Confidence Intervals Implementation Plan

Status date: 2026-07-18
Working branch: `analytical-cis`

This is the durable implementation checklist for analytical confidence
intervals for `sops()`, `avg_sops()`, and supported `avg_comparisons()` results.
It records the approved design, what is currently present in the branch, and
the work that remains before integration.

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
