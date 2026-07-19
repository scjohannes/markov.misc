# automatic SOP prediction never substitutes a later profile row

    Code
      sops(fit, times = 1:3, y_levels = 1:6, absorb = 6)
    Condition
      Error:
      ! Fitted patient(s) lack a row at the designated starting visit `1`: 1. A later likelihood row is not substituted for the starting profile.

# incomplete and duplicated starting profiles fail clearly

    Code
      avg_sops(incomplete_fit, variables = list(tx = c(0, 1)), times = 1:2, y_levels = 1:
        6, absorb = 6)
    Condition
      Error:
      ! Fitted patient(s) have incomplete starting-profile predictors: 1. The transition response is not required, but ID, modeled predictors, time/origin metadata, and the previous state must be complete.

---

    Code
      avg_sops(duplicated_fit, variables = list(tx = c(0, 1)), times = 1:2, y_levels = 1:
        6, absorb = 6)
    Condition
      Error:
      ! Fitted patient(s) have multiple rows at the designated starting visit `1`: 1. Exactly one starting profile is required per fitted patient.
