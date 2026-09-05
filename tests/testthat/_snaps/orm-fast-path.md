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
      ! Model prediction requires transitions from state level(s) not represented in the fitted transition data: 6. If a level is absorbing, pass it via `absorb`.

