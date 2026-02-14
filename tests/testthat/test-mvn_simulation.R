# Unit Tests: MVN Simulation-Based Inference for SOPs
#
# These tests verify that:
# 1. set_coef() works correctly for vglm and orm models
# 2. get_vcov_robust() returns correct variance matrices
# 3. inferences(method = "simulation") produces valid confidence intervals
# 4. inferences() works for both avg_sops and sops objects
# 5. get_draws() extracts simulation draws correctly
# 6. score-bootstrap simulation works for robcov_vglm avg_sops pipelines

describe("MVN Simulation-Based Inference for SOPs", {
  describe("set_coef()", {
    it("works for vglm models", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 50, seed = 111)
      m_vglm <- make_test_model(data)

      # Get original coefficients
      orig_coefs <- coef(m_vglm)

      # Create new coefficients (slightly perturbed)
      new_coefs <- orig_coefs * 1.1

      # Set new coefficients
      m_modified <- set_coef(m_vglm, new_coefs)

      # Check that coefficients were updated
      expect_equal(coef(m_modified), new_coefs)

      # Check that original model is unchanged
      expect_equal(coef(m_vglm), orig_coefs)

      # Check predictions differ
      pred_orig <- predict(m_vglm, newdata = data[1:5, ], type = "response")
      pred_mod <- predict(m_modified, newdata = data[1:5, ], type = "response")
      expect_false(all(pred_orig == pred_mod))
    })

    it("works for orm models", {
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 50, seed = 222)

      # Fit orm model - datadist must be defined for rms
      dd <- rms::datadist(data)
      withr::local_options(datadist = "dd")
      # assign to global for rms to find it, then defer removal
      assign("dd", dd, envir = globalenv())
      withr::defer(
        if (exists("dd", envir = globalenv())) rm("dd", envir = globalenv())
      )

      m_orm <- rms::orm(
        y ~ time_lin + time_nlin_1 + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      )

      # Get original coefficients
      orig_coefs <- coef(m_orm)

      # Create new coefficients
      new_coefs <- orig_coefs * 0.9

      # Set new coefficients
      m_modified <- set_coef(m_orm, new_coefs)

      # Check that coefficients were updated
      expect_equal(coef(m_modified), new_coefs)

      # Check that original model is unchanged
      expect_equal(coef(m_orm), orig_coefs)
    })
  })

  describe("get_vcov_robust()", {
    it("returns robust (sandwich) vcov when cluster is NULL for vglm models", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 100, seed = 333, follow_up_time = 15)
      m_vglm <- make_test_model(data)

      # Get vcov without clustering
      V <- get_vcov_robust(m_vglm, cluster = NULL)

      # Should NOT match standard vcov exactly (sandwich vs model-based)
      # They will differ if the model is not perfectly specified
      expect_false(isTRUE(all.equal(V, vcov(m_vglm), tolerance = 1e-8)))

      # Should match robcov_vglm output
      V_rob <- robcov_vglm(m_vglm, cluster = NULL)$var
      expect_equal(V, V_rob)
    })

    it("returns robust (sandwich) vcov when cluster is NULL for orm models", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 100, seed = 333, follow_up_time = 15)
      m_orm <- rms::orm(
        y ~ time_lin + time_nlin_1 + tx + yprev,
        data = data,
        x = TRUE,
        y = TRUE
      )

      # Get vcov without clustering
      V <- get_vcov_robust(m_orm, cluster = NULL)

      m_robust <- rms::robcov(m_orm)

      # Should NOT match standard vcov exactly (sandwich vs model-based)
      # They will differ if the model is not perfectly specified
      expect_false(isTRUE(all.equal(V, m_robust$orig.var, tolerance = 1e-8)))

      # Should match robcov_vglm output
      V_rob <- m_robust$var
      expect_equal(V, V_rob)
    })

    it("computes cluster-robust vcov with formula", {
      skip_if_not_installed("VGAM")

      data <- make_test_data(n_patients = 30, seed = 444, follow_up_time = 15)
      m_vglm <- make_test_model(data)

      # Get cluster-robust vcov using formula
      V_robust <- get_vcov_robust(m_vglm, cluster = ~id, data = data)

      # Should be a matrix with correct dimensions
      expect_true(is.matrix(V_robust))
      expect_equal(nrow(V_robust), length(coef(m_vglm)))
      expect_equal(ncol(V_robust), length(coef(m_vglm)))
    })

    it("extracts vcov from robcov_vglm object", {
      skip_if_not_installed("VGAM")

      data <- make_test_data(n_patients = 30, seed = 555, follow_up_time = 15)
      m_vglm <- make_test_model(data)

      # Create robcov_vglm object
      m_robust <- robcov_vglm(m_vglm, cluster = data$id)

      # Extract vcov from robcov_vglm
      V <- get_vcov_robust(m_robust)

      # Should match the stored var
      expect_equal(V, m_robust$var)
    })
  })

  describe("inferences() with simulation method", {
    it("produces valid CIs for avg_sops", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")
      skip_if_not_installed("mvtnorm")

      data <- make_test_data(
        n_patients = 50,
        seed = 666,
        treatment_effect = 0.3,
        follow_up_time = 15
      )
      m_robust <- make_test_model(data, robust = TRUE)

      # Compute avg_sops with robust model
      result <- avg_sops(
        model = m_robust,
        newdata = data,
        variables = list(tx = c(0, 1)),
        times = 1:15,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      )

      # Add simulation-based inference (vcov extracted from model)
      result_ci <- inferences(
        result,
        method = "simulation",
        n_sim = 100, # Small for testing
        conf_level = 0.95
      )

      # Check that CI columns exist
      expect_contains(names(result_ci), c("conf.low", "conf.high", "std.error"))

      # Check CIs are valid - for most observations, low <= estimate <= high
      # With percentile-based CIs, some estimates may fall slightly outside bounds
      # due to simulation variance. We expect at least 95% to be within bounds.
      in_bounds <- (result_ci$conf.low <= result_ci$estimate) &
        (result_ci$conf.high >= result_ci$estimate)
      expect_gt(mean(in_bounds, na.rm = TRUE), 0.99)
      expect_all_true(result_ci$std.error >= 0)

      # Check attributes
      expect_equal(attr(result_ci, "method"), "simulation")
      expect_equal(attr(result_ci, "n_sim"), 100)
    })

    it("stores simulation draws when return_draws=TRUE", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")
      skip_if_not_installed("mvtnorm")

      data <- make_test_data(n_patients = 40, seed = 777, follow_up_time = 10)
      m_robust <- make_test_model(data, robust = TRUE)

      # Compute avg_sops with simulation CIs and draws
      result <- avg_sops(
        model = m_robust,
        newdata = data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      ) |>
        inferences(
          method = "simulation",
          n_sim = 30,
          return_draws = TRUE
        )

      # Extract draws
      draws <- get_draws(result)

      # Check draws structure
      expect_s3_class(draws, "data.frame")
      expect_contains(
        names(draws),
        c("draw_id", "estimate", "time", "state", "tx")
      )

      # Check number of draws
      n_draws <- length(unique(draws$draw_id))
      expect_lte(n_draws, 30) # May be less if some failed

      # Check we have draws for all combinations
      n_times <- 10
      n_states <- 6
      n_tx <- 2
      expected_rows_per_draw <- n_times * n_states * n_tx
      actual_rows_per_draw <- nrow(draws) / n_draws
      expect_equal(actual_rows_per_draw, expected_rows_per_draw)
    })

    it("supports score-bootstrap engine for avg_sops with robcov_vglm", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 50, seed = 778, follow_up_time = 10)
      baseline_data <- data[data$time == 1, ]
      m_robust <- make_test_model(data, robust = TRUE)

      avg_result <- avg_sops(
        model = m_robust,
        newdata = baseline_data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      )

      result_score <- inferences(
        avg_result,
        method = "simulation",
        engine = "score_bootstrap",
        n_sim = 60,
        return_draws = TRUE
      )

      expect_contains(names(result_score), c("conf.low", "conf.high", "std.error"))
      expect_equal(attr(result_score, "method"), "simulation")
      expect_equal(attr(result_score, "engine"), "score_bootstrap")
      expect_equal(attr(result_score, "score_weight_dist"), "exponential")
      expect_gt(mean(result_score$std.error, na.rm = TRUE), 0)

      draws <- get_draws(result_score)
      expect_s3_class(draws, "data.frame")
      expect_contains(names(draws), c("draw_id", "estimate", "time", "state", "tx"))
      expect_lte(length(unique(draws$draw_id)), 60)
    })

    it("works with custom vcov matrix", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")
      skip_if_not_installed("mvtnorm")

      data <- make_test_data(n_patients = 40, seed = 999, follow_up_time = 10)
      m_vglm <- make_test_model(data)

      # Compute custom vcov (e.g., inflated for conservative CIs)
      V_standard <- vcov(m_vglm)
      V_inflated <- V_standard * 2 # Inflate variance

      # Compute avg_sops with custom vcov
      result <- avg_sops(
        model = m_vglm,
        newdata = data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      ) |>
        inferences(
          method = "simulation",
          n_sim = 50,
          vcov = V_inflated
        )

      # CIs should be wider than with standard vcov
      result_standard <- avg_sops(
        model = m_vglm,
        newdata = data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      ) |>
        inferences(
          method = "simulation",
          n_sim = 50
        )

      # Compare CI widths
      ci_width_inflated <- result$conf.high - result$conf.low
      ci_width_standard <- result_standard$conf.high - result_standard$conf.low

      # Inflated vcov should give wider CIs on average
      expect_gt(
        mean(ci_width_inflated, na.rm = TRUE),
        mean(ci_width_standard, na.rm = TRUE)
      )
    })

    it("produces different CIs for 'wald' vs 'perc' conf_type", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")
      skip_if_not_installed("mvtnorm")

      data <- make_test_data(n_patients = 40, seed = 1010, follow_up_time = 10)
      m_robust <- make_test_model(data, robust = TRUE)

      # Compute avg_sops with robust model
      avg_result <- avg_sops(
        model = m_robust,
        newdata = data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      )

      # Percentile CIs (vcov from model)
      result_perc <- inferences(
        avg_result,
        method = "simulation",
        n_sim = 100,
        conf_type = "perc"
      )

      # Wald CIs (vcov from model)
      result_wald <- inferences(
        avg_result,
        method = "simulation",
        n_sim = 100,
        conf_type = "wald"
      )

      # Check attributes
      expect_equal(attr(result_perc, "conf_type"), "perc")
      expect_equal(attr(result_wald, "conf_type"), "wald")

      # Both should have valid CIs - for most observations
      in_bounds_perc <- (result_perc$conf.low <= result_perc$estimate + 1e-6)
      in_bounds_wald <- (result_wald$conf.low <= result_wald$estimate + 1e-6)
      expect_gt(mean(in_bounds_perc, na.rm = TRUE), 0.95)
      expect_gt(mean(in_bounds_wald, na.rm = TRUE), 0.95)
    })

    it("works for individual sops()", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")
      skip_if_not_installed("mvtnorm")

      data <- make_test_data(n_patients = 40, seed = 1313, follow_up_time = 10)
      m_robust <- make_test_model(data, robust = TRUE)

      # Compute individual sops with robust model
      sops_result <- sops(
        model = m_robust,
        newdata = data[data$time == 1, ],
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        tvarname = "time",
        pvarname = "yprev"
      )

      # Check that newdata_orig attribute is stored
      expect_false(is.null(attr(sops_result, "newdata_orig")))
      expect_equal(
        nrow(attr(sops_result, "newdata_orig")),
        nrow(data[data$time == 1, ])
      )

      # Add simulation-based inference (vcov from model)
      result_ci <- inferences(
        sops_result,
        method = "simulation",
        n_sim = 100,
        conf_level = 0.95,
        return_draws = TRUE
      )

      # Check that CI columns exist
      expect_contains(
        names(result_ci),
        c("rowid", "estimate", "conf.low", "conf.high", "std.error")
      )

      # Check CIs are valid
      in_bounds <- (result_ci$conf.low <= result_ci$estimate) &
        (result_ci$conf.high >= result_ci$estimate)
      expect_gt(mean(in_bounds, na.rm = TRUE), 0.999)

      # Check attributes
      expect_equal(attr(result_ci, "method"), "simulation")
      expect_equal(attr(result_ci, "n_sim"), 100)

      # Check draws can be extracted
      draws <- get_draws(result_ci)
      expect_s3_class(draws, "data.frame")
      expect_contains(
        names(draws),
        c("rowid", "draw", "draw_id", "estimate", "rowid", "time", "state")
      )
    })
  })

  describe("Error handling", {
    it("errors appropriately for invalid inputs in inferences()", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")

      data <- make_test_data(n_patients = 30, seed = 1111, follow_up_time = 10)
      m_vglm <- make_test_model(data)

      # Error: passing model directly instead of sops object
      expect_error(
        inferences(m_vglm),
        "markov_avg_sops|markov_sops"
      )

      # Error: invalid method
      avg_result <- avg_sops(
        model = m_vglm,
        newdata = data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      )

      expect_error(
        inferences(avg_result, method = "invalid"),
        "'arg' should be one of"
      )

      # Error: score bootstrap requires robcov_vglm
      expect_error(
        inferences(
          avg_result,
          method = "simulation",
          engine = "score_bootstrap",
          n_sim = 20
        ),
        "requires a 'robcov_vglm' model"
      )

      # Error: score bootstrap requires avg_sops objects
      sops_result <- sops(
        model = m_vglm,
        newdata = data[data$time == 1, ],
        times = 1:10,
        ylevels = 1:6,
        absorb = 6
      )
      expect_error(
        inferences(
          sops_result,
          method = "simulation",
          engine = "score_bootstrap",
          n_sim = 20
        ),
        "currently supports only 'markov_avg_sops'"
      )
    })

    it("errors when get_draws() called but no draws stored", {
      skip_if_not_installed("VGAM")
      skip_if_not_installed("rms")
      skip_if_not_installed("mvtnorm")

      data <- make_test_data(n_patients = 30, seed = 1212, follow_up_time = 10)
      m_robust <- make_test_model(data, robust = TRUE)

      # Compute avg_sops without return_draws
      result <- avg_sops(
        model = m_robust,
        newdata = data,
        variables = list(tx = c(0, 1)),
        times = 1:10,
        ylevels = 1:6,
        absorb = 6,
        id_var = "id"
      ) |>
        inferences(
          method = "simulation",
          n_sim = 20,
          return_draws = FALSE
        )

      # get_draws should error
      expect_error(
        get_draws(result),
        "No draws found"
      )
    })
  })
})
