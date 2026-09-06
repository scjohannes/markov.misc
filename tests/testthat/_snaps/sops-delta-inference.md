# public delta scope enforces fixed targets and patient clustering

    Code
      inferences(fixed, method = "delta", vcov = "unconditional")
    Condition
      Error:
      ! `sops()` delta inference supports only `vcov = "conditional"` or a coefficient covariance matrix.

---

    Code
      inferences(fixed, method = "delta", vcov = "conditional")
    Condition
      Error:
      ! Analytical unconditional inference requires patient clustering. Supply `cluster`, or fit with `orm_markov(..., id_var = ...)` or `vglm_markov(..., id_var = ...)` so row-aligned fitting data and patient-ID metadata are stored. Observation rows are not used as implicit clusters.

---

    Code
      inferences(avg, method = "delta", vcov = c("conditional", "unconditional"))
    Condition
      Error:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(avg, method = "delta", vcov = "population")
    Condition
      Error:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(supplied, method = "delta", vcov = "unconditional")
    Condition
      Error:
      ! `vcov = "unconditional"` requires the stored fitted-patient cohort. For user-supplied `newdata`, use `vcov = "conditional"`.

# logit delta intervals distinguish structural boundaries

    Code
      nonstructural <- delta_interval_bounds(estimate = c(0, 1), standard_error = c(
        0.1, 0.1), conf_level = 0.95, conf_type = "logit")
    Condition
      Warning:
      Logit-delta limits are undefined for nonstructural boundary SOP estimates; their confidence limits were set to NA.

# delta vcov validates strings and does not silently change targets

    Code
      inferences(object, method = "delta")
    Condition
      Error:
      ! `vcov = "unconditional"` requires the stored fitted-patient cohort. For user-supplied `newdata`, use `vcov = "conditional"`.

---

    Code
      inferences(object, method = "delta", vcov = NA_character_)
    Condition
      Error:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(object, method = "delta", vcov = character())
    Condition
      Error:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(object, method = "delta", vcov = "cond")
    Condition
      Error:
      ! `vcov` must be "conditional", "unconditional", NULL, or a coefficient covariance matrix.

---

    Code
      inferences(object, method = "delta", target = "empirical")
    Condition
      Error in `inferences()`:
      ! unused argument (target = "empirical")

---

    Code
      inferences(object, method = method, vcov = "conditional")
    Condition
      Error in `inferences_impl()`:
      ! Character `vcov` choices are only available with `method = "delta"`.

---

    Code
      inferences(object, method = method, vcov = "conditional")
    Condition
      Error in `inferences_impl()`:
      ! Character `vcov` choices are only available with `method = "delta"`.

---

    Code
      inferences(object, method = method, vcov = "conditional")
    Condition
      Error in `inferences_impl()`:
      ! Character `vcov` choices are only available with `method = "delta"`.

---

    Code
      inferences(object, method = method, vcov = "conditional")
    Condition
      Error in `inferences_impl()`:
      ! Character `vcov` choices are only available with `method = "delta"`.

