# Plan: SOP Comparisons and Average Comparisons

## Summary

Add a post-processing `comparisons()` function for `avg_sops()` output and a
convenience `avg_comparisons()` wrapper for model-based workflows.

The main design goal is to compute contrasts from paired SOP draws without
rerunning simulation or bootstrap uncertainty separately for each comparison
level. This preserves the covariance structure between counterfactual levels and
keeps inference computationally efficient.

This plan keeps `avg_sops(variables = NULL)` unsupported.

## Public API

Add:

```r
comparisons(x, ...)
avg_comparisons(model, ...)
```

`comparisons()`:

- Takes a `markov_avg_sops` object.
- Computes contrasts from point estimates and stored draws.
- Supports metrics that are valid from marginal SOPs.
- Does not refit models or rerun inference.

`avg_comparisons()`:

- Calls `avg_sops()` once.
- Optionally calls `inferences()` once.
- Calls `comparisons()` for marginal-SOP metrics.
- Uses a separate patient/profile-level internal path for `time_benefit`.

## Contrast Interface

Follow the `marginaleffects` style where possible:

```r
comparisons(
  x,
  variables = list(tx = c(0, 1)),
  comparison = "difference",
  metric = "time_in_state",
  states = "1"
)
```

Use:

- `variables` to identify counterfactual levels to compare.
- `comparison` to identify the transformation:
  - `"difference"`: comparison level minus reference level
  - `"ratio"`: comparison level divided by reference level
  - later: `"lnratio"` or user-supplied functions
- `metric` to identify the estimand.
- `states` to identify state-specific or lumped-state summaries.

Avoid a separate ad hoc argument such as:

```r
contrast = list("1 - 0" = c("1", "0"))
```

because that is not how `marginaleffects` structures the interface.

## State Selection

Support both single states and lumped state sets:

```r
states = "1"
states = c("1", "2", "3")
states = list(
  state_1 = "1",
  good = c("1", "2", "3"),
  poor = c("6", "7", "8")
)
```

Rules:

- `states = NULL`: return one result per state.
- unnamed character vector: treat as one lumped state set.
- named list: return one result per named state set.
- lumped states are summed before contrast calculation.

## Supported Metrics From `comparisons(markov_avg_sops)`

### `metric = "sop"`

Time-specific SOP contrasts.

Output has one row per:

```r
time, state_set, variable, comparison, by
```

This answers: difference or ratio in occupancy probability at each time point.

### `metric = "time_in_state"`

Expected time in selected state(s).

Visit-scale version:

```r
sum_t P(Y_t in states)
```

Real-time version should reuse existing interpolation or `time_in_state()`
machinery instead of duplicating AUC logic.

This metric is linear in SOPs, so it is valid to compute from already
marginalized `avg_sops()` output:

```r
mean_i(sum_t p_i(t)) == sum_t mean_i(p_i(t))
```

Supports:

- differences
- ratios
- lumped state sets
- `by` strata

## `time_benefit`

Use the generic name `time_benefit`, not `days_benefit`, because the time grid
may represent days, weeks, months, visits, or another unit.

Allow a label-only argument such as:

```r
time_unit = "days"
```

when needed for printing or plotting.

### Estimand

For lower states better:

```r
sum_t [P(Y_tx better than Y_ref) - P(Y_ref better than Y_tx)]
```

Equivalently, at the patient/profile level:

```r
sum_t sum_k0 sum_k1 score(k0, k1) *
  P(Y_ref = k0 | profile_i, t) *
  P(Y_tx  = k1 | profile_i, t)
```

with default scoring:

```r
score(k0, k1) =  1  if k0 > k1
                 0  if k0 == k1
                -1  if k0 < k1
```

when lower states are better. Reverse this logic when higher states are better.

### Why `comparisons(markov_avg_sops)` Cannot Support `time_benefit`

`time_benefit` is nonlinear in the two counterfactual state distributions.

The correct quantity requires products within the same patient/profile:

```r
mean_i[p_i_ref(k0, t) * p_i_tx(k1, t)]
```

A plain `markov_avg_sops` object only contains:

```r
mean_i[p_i_ref(k0, t)]
mean_i[p_i_tx(k1, t)]
```

Multiplying marginal averages is generally incorrect:

```r
mean_i[p_i_ref * p_i_tx] != mean_i[p_i_ref] * mean_i[p_i_tx]
```

Therefore:

- `comparisons(markov_avg_sops, metric = "time_benefit")` should error unless
  the object explicitly carries patient/profile-level counterfactual
  probabilities.
- `avg_comparisons(metric = "time_benefit")` should compute patient/profile-level
  counterfactual probabilities internally, reduce to `time_benefit`, and only
  then average over profiles.

## Draw Handling and Inference

For every supported metric:

- Compute the point estimate from point-estimate data.
- If draws are available, compute the complete estimand inside each `draw_id`.
- Summarize draw-level estimands with:
  - percentile confidence intervals
  - `std.error = sd(draw_estimate)`
- Preserve pairing across counterfactual levels by using the same draw for all
  levels in a contrast.
- Never run uncertainty separately for each contrast level and subtract later.

If `comparisons()` is called on an object without stored draws:

- Return point estimates.
- Either omit uncertainty columns or return `NA` intervals.
- Emit a clear warning that uncertainty requires stored draws.

## `avg_comparisons()` Workflow

For `metric = "sop"` and `metric = "time_in_state"`:

1. Call `avg_sops()` once.
2. If requested, call `inferences()` once.
3. Call `comparisons()` on the resulting object.

For `metric = "time_benefit"`:

1. Build the same counterfactual prediction grid as `avg_sops()`.
2. Predict patient/profile-level SOP arrays for all requested counterfactual
   levels.
3. For each profile, time, draw, and contrast, compute the pairwise
   `time_benefit` score from the two counterfactual state distributions.
4. Sum over time.
5. Average over profiles, using bootstrap weights when applicable.
6. Summarize across paired draws.

## Suggested Arguments

Potential interface:

```r
comparisons(
  x,
  variables,
  metric = c("sop", "time_in_state"),
  states = NULL,
  comparison = c("difference", "ratio"),
  conf_level = 0.95,
  time_map = NULL,
  time_unit = NULL,
  ...
)
```

```r
avg_comparisons(
  model,
  variables,
  metric = c("sop", "time_in_state", "time_benefit"),
  states = NULL,
  comparison = c("difference", "ratio"),
  ordinal_direction = c("lower_better", "higher_better"),
  time_unit = NULL,
  inference = NULL,
  return_draws = FALSE,
  ...
)
```

`inference` should be a list of arguments passed to `inferences()`, for example:

```r
inference = list(method = "bootstrap", engine = "fwb", n_sim = 500)
```

This avoids inventing another inference API.

## Output

Return a data frame with class such as `markov_sops_comparisons`.

Core columns:

```r
metric
state_set
variable
contrast
comparison
estimate
conf.low
conf.high
std.error
```

Additional columns:

- `time` for `metric = "sop"`
- `time_unit` when supplied
- `by` columns from the original `avg_sops()` object
- optional draw storage when `return_draws = TRUE`

## Internal Helpers

Add internal helpers such as:

```r
comparison_grid()
normalize_state_sets()
reduce_sop_metric()
reduce_time_in_state_metric()
reduce_time_benefit_metric()
summarize_comparison_draws()
```

Keep metric reduction separate from model prediction and inference dispatch.

## Tests

Add tests for:

- SOP differences at each time point.
- Time-in-state difference for a single state.
- Time-in-state difference for lumped states.
- Ratio of time in state.
- Multiple state sets in one call.
- Multiple contrasts in one call.
- `by` strata are preserved.
- Draw-level contrasts are paired by `draw_id`.
- `comparisons()` warns or omits uncertainty when no draws are stored.
- `comparisons(metric = "time_benefit")` errors on a plain `markov_avg_sops`
  object without patient/profile-level counterfactual probabilities.
- `avg_comparisons(metric = "time_benefit")` computes from patient/profile-level
  counterfactual probabilities.
- `avg_comparisons()` calls simulation/bootstrap inference only once.

## Documentation

Document:

- `comparisons()` as post-processing for `avg_sops()` outputs.
- `avg_comparisons()` as a convenience wrapper.
- The distinction between linear metrics available from marginal SOPs and
  nonlinear metrics requiring patient/profile-level counterfactual probabilities.
- The generic `time_benefit` name and `time_unit` labeling.
- The rule that uncertainty is computed by applying the estimand inside each
  paired draw.

Update `ARCHITECTURE.md` if implemented, because this adds a new endpoint-summary
and comparison layer on top of SOP inference.
