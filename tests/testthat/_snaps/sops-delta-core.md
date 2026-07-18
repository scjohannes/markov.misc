# analytical SOP validation rejects unsupported model structures

    Code
      run_sop_delta_plan(second_order, case$model)
    Condition
      Error:
      ! Analytical SOP inference currently supports first-order Markov plans only.

---

    Code
      run_sop_delta_plan(unsupported_link, case$model)
    Condition
      Error:
      ! Analytical SOP inference currently supports the reverse cumulative logit link only.

---

    Code
      run_sop_delta_plan(case$plan, unsupported_model)
    Condition
      Error:
      ! Analytical SOP inference supports only `vglm`, `robcov_vglm`, and `orm` models.

---

    Code
      run_sop_delta_plan(partial_plan, partial_model)
    Condition
      Error:
      ! Analytical SOP inference requires full proportional odds; design column `time` has threshold-specific slopes.

# crossed raw ordinal probabilities are rejected without clipping

    Code
      run_sop_delta_plan(case$plan, crossed_model)
    Condition
      Error:
      ! Raw ordinal category probabilities are crossed at the first SOP time point (minimum -0.46151).

# analytical allocation guard is typed and configurable

    Code
      delta_assert_bytes(16, "Test allocation")
    Condition
      Error:
      ! Test allocation requires 16 bytes, above the analytical delta-method limit of 1 bytes. Increase option `markov.misc.delta_max_bytes` only if the required allocation is acceptable.
