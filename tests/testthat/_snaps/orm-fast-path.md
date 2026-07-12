# orm fast path validates state support and initial probabilities

    Code
      markov.misc::soprob_markov(fit, baseline, times = 1, y_levels = fit$yunique[
        -length(fit$yunique)])
    Condition
      Error:
      ! `y_levels` defines 5 states, but the fitted model defines 6 states through 5 threshold coefficients.

---

    Code
      markov.misc::soprob_markov(fit, baseline[1L, , drop = FALSE], times = 1,
      y_levels = fit$yunique)
    Condition
      Error:
      ! Model prediction returned missing transition probabilities during the first SOP time point. This often happens when the SOP recursion asks for transitions from a state level that is not represented in the fitted transition model. If that level is absorbing, pass it via `absorb`; otherwise check that `data`, `y_levels`, and fitted factor levels agree.

