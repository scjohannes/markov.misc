# sim_actt1_markov validates wrapper inputs

    Code
      sim_actt1_markov(n_patients = 0)
    Condition
      Error in `sim_actt1_markov()`:
      ! n_patients must be a positive integer

---

    Code
      sim_actt1_markov(treatment_prob = 2)
    Condition
      Error in `sim_actt1_markov()`:
      ! treatment_prob must be a numeric scalar between 0 and 1

---

    Code
      sim_actt1_markov(treatment_effect = Inf)
    Condition
      Error in `sim_actt1_markov()`:
      ! treatment_effect must be a finite numeric scalar

---

    Code
      sim_actt1_markov(treatment_effect_decay = -1)
    Condition
      Error in `sim_actt1_markov()`:
      ! treatment_effect_decay must be a nonnegative finite numeric scalar

# sim_actt2_markov validates wrapper inputs

    Code
      sim_actt2_markov(n_patients = 0)
    Condition
      Error in `sim_actt2_markov()`:
      ! n_patients must be a positive integer

---

    Code
      sim_actt2_markov(treatment_prob = 2)
    Condition
      Error in `sim_actt2_markov()`:
      ! treatment_prob must be a numeric scalar between 0 and 1

---

    Code
      sim_actt2_markov(treatment_effect = Inf)
    Condition
      Error in `sim_actt2_markov()`:
      ! treatment_effect must be a finite numeric scalar

---

    Code
      sim_actt2_markov(treatment_effect_decay = -1)
    Condition
      Error in `sim_actt2_markov()`:
      ! treatment_effect_decay must be a nonnegative finite numeric scalar

# sim_actt2_markov_60day validates wrapper inputs

    Code
      sim_actt2_markov_60day(n_patients = 0)
    Condition
      Error in `sim_actt2_markov_60day()`:
      ! n_patients must be a positive integer

---

    Code
      sim_actt2_markov_60day(treatment_prob = 2)
    Condition
      Error in `sim_actt2_markov_60day()`:
      ! treatment_prob must be a numeric scalar between 0 and 1

---

    Code
      sim_actt2_markov_60day(treatment_effect = Inf)
    Condition
      Error in `sim_actt2_markov_60day()`:
      ! treatment_effect must be a finite numeric scalar

