# public delta scope enforces fixed targets and patient clustering

    Code
      inferences(fixed, method = "delta", target = "empirical", vcov = case$
      covariance)
    Condition
      Error:
      ! `sops()` delta inference supports only `target = "fixed"`. Empirical and superpopulation targets apply to averaged SOP objects.

---

    Code
      inferences(fixed, method = "delta", target = "fixed")
    Condition
      Error:
      ! Analytical superpopulation inference requires patient clustering. Supply `cluster`, or fit with `orm_markov(..., id_var = ...)` or `vglm_markov(..., id_var = ...)` so row-aligned fitting data and patient-ID metadata are stored. Observation rows are not used as implicit clusters.

---

    Code
      inferences(avg, method = "delta", target = "superpopulation", vcov = superpopulation_case$
        model$var)
    Condition
      Error:
      ! `vcov` cannot be supplied with `target = "superpopulation"`; superpopulation inference uses fitted-model score components and the stacked patient influence function.

---

    Code
      inferences(avg, method = "delta", target = "population")
    Condition
      Error:
      ! Averaged delta inference supports `target = "empirical"` or `target = "superpopulation"`.

---

    Code
      inferences(supplied, method = "delta", target = "superpopulation")
    Condition
      Error:
      ! `target = "superpopulation"` requires the stored fitted-patient cohort. User-supplied `newdata` is treated as a fixed cohort.

# logit delta intervals distinguish structural boundaries

    Code
      nonstructural <- delta_interval_bounds(estimate = c(0, 1), standard_error = c(
        0.1, 0.1), conf_level = 0.95, conf_type = "logit")
    Condition
      Warning:
      Logit-delta limits are undefined for nonstructural boundary SOP estimates; their confidence limits were set to NA.
