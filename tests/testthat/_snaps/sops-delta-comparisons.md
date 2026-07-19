# superpopulation comparisons retain transformed stacked influence

    Code
      inferences(point, method = "delta", target = "superpopulation", vcov = case$
        covariance)
    Condition
      Error:
      ! `vcov` cannot be supplied with `target = "superpopulation"`; superpopulation inference uses fitted-model score components and the stacked patient influence function.

# analytical comparison scope rejects deferred estimands

    Code
      inferences(ratio, method = "delta")
    Condition
      Error:
      ! Analytical average comparisons currently support differences only; ratios are not yet supported.

---

    Code
      inferences(time_benefit, method = "delta")
    Condition
      Error:
      ! Analytical average comparisons currently support `estimand = "sop"` or `"time_in_state"`; time benefit is not yet supported.

---

    Code
      inferences(stratified, method = "delta")
    Condition
      Error:
      ! Analytical delta inference with `by` is not yet supported.

---

    Code
      inferences(second_order, method = "delta", vcov = case$covariance)
    Condition
      Error:
      ! Analytical delta inference for second-order Markov models is not yet supported.

# analytical comparisons reject partial proportional odds

    Code
      inferences(partial, method = "delta", vcov = stats::vcov(partial_model))
    Condition
      Error:
      ! Analytical delta inference currently requires a full proportional-odds model; partial proportional odds are not yet supported.
