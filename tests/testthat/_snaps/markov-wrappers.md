# Markov wrappers validate numeric first follow-up schedules

    Code
      vglm_markov(ordered(y) ~ time_lin + tx + yprev, family = VGAM::cumulative(
        reverse = TRUE, parallel = TRUE), data = no_one, id_var = "id")
    Condition
      Error:
      ! The default `first_followup_time = 1` is not observed in the numeric time variable.

---

    Code
      vglm_markov(ordered(y) ~ time_lin + tx + yprev, family = VGAM::cumulative(
        reverse = TRUE, parallel = TRUE), data = below_one, id_var = "id")
    Condition
      Error:
      ! The numeric time variable contains a value below 1; follow-up time must begin at 1 or later.

---

    Code
      vglm_markov(ordered(y) ~ time_lin + tx + yprev, family = VGAM::cumulative(
        reverse = TRUE, parallel = TRUE), data = data, id_var = "id",
      first_followup_time = 2)
    Condition
      Error:
      ! For numeric time, `first_followup_time` must equal the earliest observed scheduled time (`1`).

# Markov wrappers reject legacy starting-profile arguments

    Code
      vglm_markov(ordered(y) ~ time_lin + tx + yprev, family = VGAM::cumulative(
        reverse = TRUE, parallel = TRUE), data = data, id_var = "id", start_time = 1,
      origin_time = 0, time_map = 1:6)
    Condition
      Error:
      ! vglm_markov() no longer accepts legacy wrapper argument(s): `start_time`, `origin_time`, `time_map`. Use `first_followup_time` to select the stored starting profiles; supply baseline and time-mapping arguments to downstream SOP summary functions.

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
      ! Fitted patient(s) have incomplete starting-profile predictors: 1. The transition response is not required, but ID, modeled predictors, time metadata, and the previous state must be complete.

---

    Code
      avg_sops(duplicated_fit, variables = list(tx = c(0, 1)), times = 1:2, y_levels = 1:
        6, absorb = 6)
    Condition
      Error:
      ! Fitted patient(s) have multiple rows at the designated starting visit `1`: 1. Exactly one starting profile is required per fitted patient.

# factor and character time require explicit first follow-up values

    Code
      vglm_markov(ordered(y) ~ time + tx + yprev, family = VGAM::cumulative(reverse = TRUE,
        parallel = TRUE), data = data, id_var = "id")
    Condition
      Error:
      ! Factor and character time require an explicit `first_followup_time`.

---

    Code
      vglm_markov(ordered(y) ~ time_lin + tx + yprev, family = VGAM::cumulative(
        reverse = TRUE, parallel = TRUE), data = data, id_var = "id", time_var = "visit")
    Condition
      Error:
      ! Factor and character time require an explicit `first_followup_time`.

---

    Code
      vglm_markov(ordered(y) ~ time + tx + yprev, family = VGAM::cumulative(reverse = TRUE,
        parallel = TRUE), data = data, id_var = "id", first_followup_time = "v7")
    Condition
      Error:
      ! `first_followup_time = v7` is not observed in the factor or character time variable.
