# grouped blrm summaries reject entirely incomplete grouping rows

    Code
      sops(model, newdata = newdata, times = 1:2, by = "grp", n_draws = 2)
    Condition
      Error:
      ! No complete prediction rows remain after omitting missing grouping values.

---

    Code
      avg_sops(model, newdata = newdata, variables = list(tx = c(0, 1)), times = 1:2,
      by = "grp", n_draws = 2)
    Condition
      Error:
      ! No complete prediction rows remain after omitting missing grouping values.

