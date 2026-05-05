describe("avg_sops() and inferences() pipeline", {
  add_patient_age <- function(data) {
    ids <- sort(unique(data$id))
    age_lookup <- data.frame(
      id = ids,
      age = 45 + ((ids * 11) %% 31)
    )

    data$age <- NULL
    merge(data, age_lookup, by = "id", all.x = TRUE, sort = FALSE)
  }

  add_time_basis <- function(data) {
    spl <- rms::rcs(data$time, 4)
    data$time_lin <- as.vector(spl[, 1])
    data$time_nlin_1 <- as.vector(spl[, 2])
    data$time_nlin_2 <- as.vector(spl[, 3])
    data
  }

  get_time_covariates <- function(data) {
    cols <- c("time", "time_lin", "time_nlin_1", "time_nlin_2")
    out <- unique(data[, cols, drop = FALSE])
    out <- out[order(out$time), , drop = FALSE]
    rownames(out) <- NULL
    out[, c("time_lin", "time_nlin_1", "time_nlin_2"), drop = FALSE]
  }

  build_brownian_pipeline_case <- function() {
    raw_data <- markov.misc::sim_trajectories_brownian(
      n_patients = 70,
      follow_up_time = 8,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 4041,
      mu_treatment_effect = 0.18
    )
    data <- markov.misc::prepare_markov_data(raw_data)
    data <- add_patient_age(data)
    data <- add_time_basis(data)

    list(
      data = data,
      baseline = data[data$time == 1, , drop = FALSE],
      t_covs = get_time_covariates(data),
      ylevels = 1:6,
      absorb = "6"
    )
  }

  fit_pipeline_model <- function(data) {
    suppressWarnings(
      VGAM::vglm(
        ordered(y) ~ (time_lin + time_nlin_1 + time_nlin_2) * tx + yprev + age,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data
      )
    )
  }

  fit_inline_pipeline_model <- function(data) {
    suppressWarnings(
      markov.misc::vglm.markov(
        ordered(y) ~ rms::rcs(time, 4) * tx + yprev + age,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data
      )
    )
  }

  pipeline_signature <- function(result) {
    keep <- result$time %in% c(1, 3, 5, 8) &
      result$state %in% c(1, 3, 6)
    out <- result[keep, c("tx", "time", "state", "estimate", "conf.low", "conf.high", "std.error")]
    out <- out[order(out$tx, out$time, out$state), , drop = FALSE]
    out$state <- as.integer(as.character(out$state))
    numeric_cols <- c("estimate", "conf.low", "conf.high", "std.error")
    out[numeric_cols] <- lapply(out[numeric_cols], round, digits = 8)
    rownames(out) <- NULL
    lapply(out, unname)
  }

  draw_signature <- function(draws) {
    keep <- draws$draw_id %in% c(1L, 2L) &
      draws$time %in% c(1, 4, 8) &
      draws$state %in% c(1, 6)
    out <- draws[keep, c("draw_id", "tx", "time", "state", "draw")]
    out <- out[order(out$draw_id, out$tx, out$time, out$state), , drop = FALSE]
    out$state <- as.integer(as.character(out$state))
    out$draw <- round(out$draw, digits = 8)
    rownames(out) <- NULL
    lapply(out, unname)
  }

  test_that("create_counterfactual_data() stacks one baseline copy per scenario", {
    baseline <- data.frame(
      id = 1:2,
      tx = c(0, 1),
      sex = c("F", "M"),
      marker = c(10, 20)
    )
    # 2x2 grid of counterfactual scenarios
    variables <- list(tx = c(0, 1), sex = c("F", "M"))
    grid <- expand.grid(variables, KEEP.OUT.ATTRS = FALSE)

    out <- markov.misc:::create_counterfactual_data(baseline, grid, variables)

    expect_equal(nrow(out), nrow(baseline) * nrow(grid))
    expect_equal(out$id, rep(baseline$id, times = nrow(grid)))
    expect_equal(out$marker, rep(baseline$marker, times = nrow(grid)))

    # Check that the counterfactual variables are correctly assigned in blocks
    for (i in seq_len(nrow(grid))) {
      rows <- ((i - 1) * nrow(baseline) + 1):(i * nrow(baseline))
      expect_equal(out$tx[rows], rep(grid$tx[i], nrow(baseline)))
      expect_equal(as.character(out$sex[rows]), rep(as.character(grid$sex[i]), nrow(baseline)))
    }
  })

  test_that("marginalize_sops_array() averages counterfactual blocks", {
    sops_array <- array(
      seq_len(4 * 2 * 3) / 100,
      dim = c(4, 2, 3),
      dimnames = list(NULL, NULL, as.character(1:3))
    )
    grid <- data.frame(tx = c(0, 1))

    out <- markov.misc:::marginalize_sops_array(
      sops_array = sops_array,
      grid = grid,
      times = c(1, 2),
      ylevels = 1:3,
      variables = list(tx = c(0, 1)),
      n_cf = 2,
      n_each = 2
    )

    expected_tx0 <- apply(sops_array[1:2, , , drop = FALSE], c(2, 3), mean)
    expected_tx1 <- apply(sops_array[3:4, , , drop = FALSE], c(2, 3), mean)
    expected <- c(as.vector(expected_tx0), as.vector(expected_tx1))

    expect_equal(out$estimate, expected)
    expect_equal(out$tx, rep(c(0, 1), each = 6))
  })

  test_that("marginalize_sops_array() supports patient weights", {
    sops_array <- array(
      c(
        0.1, 0.7, 0.2,
        0.2, 0.6, 0.2,
        0.3, 0.4, 0.3,
        0.4, 0.3, 0.3
      ),
      dim = c(2, 2, 3)
    )
    weights <- c(0.25, 0.75)

    out <- markov.misc:::marginalize_sops_array(
      sops_array = sops_array,
      grid = data.frame(tx = 0),
      times = c(1, 2),
      ylevels = 1:3,
      variables = list(tx = 0),
      n_cf = 1,
      n_each = 2,
      weights = weights
    )

    sops_matrix <- matrix(aperm(sops_array, c(2, 3, 1)), nrow = 6, ncol = 2)
    expected <- as.vector(sops_matrix %*% weights)

    expect_equal(out$estimate, expected)
    expect_error(
      markov.misc:::marginalize_sops_array(
        sops_array = sops_array,
        grid = data.frame(tx = 0),
        times = c(1, 2),
        ylevels = 1:3,
        variables = list(tx = 0),
        n_cf = 1,
        n_each = 2,
        weights = c(1, -1)
      ),
      "non-negative"
    )
  })

  test_that("array_to_df_individual() preserves row identity and optional strata", {
    sops_array <- array(seq_len(3 * 2 * 2) / 20, dim = c(3, 2, 2))
    newdata <- data.frame(
      id = c(101, 102, 103),
      rowid = c(11, 12, 13),
      tx = c(0, 0, 1)
    )

    out <- markov.misc:::array_to_df_individual(
      sops_array = sops_array,
      times = c(1, 2),
      ylevels = 1:2,
      newdata = newdata
    )

    expect_equal(nrow(out), length(sops_array))
    expect_equal(out$rowid[1:3], newdata$rowid)
    expect_equal(out$estimate, as.vector(sops_array))

    stratified <- markov.misc:::array_to_df_individual(
      sops_array = sops_array,
      times = c(1, 2),
      ylevels = 1:2,
      newdata = newdata,
      by = "tx"
    )

    manual <- aggregate(
      estimate ~ time + state + tx,
      data = out,
      FUN = mean
    )
    stratified <- stratified[order(stratified$time, stratified$state, stratified$tx), ]
    manual <- manual[order(manual$time, manual$state, manual$tx), ]
    rownames(stratified) <- NULL
    rownames(manual) <- NULL

    expect_equal(stratified, manual)
  })

  test_that("compute_ci_from_draws() computes percentile and wald intervals", {
    draws <- data.frame(
      draw_id = rep(1:4, times = 2),
      time = rep(1, 8),
      state = rep(c(1, 2), each = 4),
      estimate = c(0.1, 0.2, 0.4, 0.5, 0.5, 0.6, 0.8, 0.9)
    )

    perc <- markov.misc:::compute_ci_from_draws(
      draws_df = draws,
      group_cols = c("time", "state"),
      conf_level = 0.5,
      conf_type = "perc"
    )
    wald <- markov.misc:::compute_ci_from_draws(
      draws_df = draws,
      group_cols = c("time", "state"),
      conf_level = 0.5,
      conf_type = "wald"
    )

    expect_equal(perc$conf.low, c(0.175, 0.575))
    expect_equal(perc$conf.high, c(0.425, 0.825))
    expect_equal(perc$std.error, c(stats::sd(draws$estimate[1:4]), stats::sd(draws$estimate[5:8])))

    critical <- abs(stats::qnorm(0.25))
    expected_wald_low <- c(
      mean(draws$estimate[1:4]) - critical * stats::sd(draws$estimate[1:4]),
      mean(draws$estimate[5:8]) - critical * stats::sd(draws$estimate[5:8])
    )
    expected_wald_high <- c(
      mean(draws$estimate[1:4]) + critical * stats::sd(draws$estimate[1:4]),
      mean(draws$estimate[5:8]) + critical * stats::sd(draws$estimate[5:8])
    )

    expect_equal(wald$conf.low, expected_wald_low)
    expect_equal(wald$conf.high, expected_wald_high)
    expect_error(
      markov.misc:::compute_ci_from_draws(draws, c("time", "state"), conf_type = "other"),
      "conf_type"
    )
  })

  test_that("lp_to_probs() converts cumulative logits into valid category probabilities", {
    eta <- matrix(
      c(
        -1, 0, 1,
        0.5, -0.5, -1
      ),
      nrow = 2,
      byrow = TRUE
    )

    probs <- markov.misc:::lp_to_probs(eta, M = 3)

    expected <- cbind(
      1 - stats::plogis(eta[, 1]),
      stats::plogis(eta[, 1]) - stats::plogis(eta[, 2]),
      stats::plogis(eta[, 2]) - stats::plogis(eta[, 3]),
      stats::plogis(eta[, 3])
    )
    expected[expected < 0] <- 0

    expect_equal(probs, expected)
    expect_true(all(probs >= 0))
  })

  test_that("fast Markov components match soprob_markov() on a Brownian model", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_brownian_pipeline_case()
    model <- fit_pipeline_model(case$data)
    baseline <- case$baseline[1:12, , drop = FALSE]

    components <- markov.misc:::markov_msm_build(
      model = model,
      data = baseline,
      t_covs = case$t_covs,
      times = 1:8,
      ylevels = case$ylevels,
      tvarname = "time_lin",
      pvarname = "yprev"
    )
    gamma <- markov.misc:::compute_Gamma(
      stats::coef(model),
      VGAM::constraints(model)
    )
    fast <- markov.misc:::markov_msm_run(
      components = components,
      Gamma = gamma,
      times = 1:8,
      absorb = case$absorb
    )
    slow <- markov.misc::soprob_markov(
      object = model,
      data = baseline,
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    expect_equal(dim(fast), dim(slow))
    expect_equal(unname(fast), unname(slow), tolerance = 1e-10)
    expect_equal(rowSums(fast[, 8, , drop = FALSE][, 1, ]), rep(1, nrow(baseline)), tolerance = 1e-10)
  })

  test_that("fast Markov components support inline rcs() time terms", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_brownian_pipeline_case()
    model <- fit_inline_pipeline_model(case$data)
    baseline <- case$baseline[1:12, , drop = FALSE]

    components <- markov.misc:::markov_msm_build(
      model = model,
      data = baseline,
      times = 1:8,
      ylevels = case$ylevels,
      pvarname = "yprev"
    )
    gamma <- markov.misc:::compute_Gamma(
      stats::coef(model),
      VGAM::constraints(model)
    )
    fast <- markov.misc:::markov_msm_run(
      components = components,
      Gamma = gamma,
      times = 1:8,
      absorb = case$absorb
    )
    slow <- markov.misc::soprob_markov(
      object = model,
      data = baseline,
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      pvarname = "yprev"
    )

    expect_equal(dim(fast), dim(slow))
    expect_equal(unname(fast), unname(slow), tolerance = 1e-10)
  })

  test_that("seeded Brownian avg_sops() and fast simulation inference remain stable", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    case <- build_brownian_pipeline_case()
    model <- fit_pipeline_model(case$data)

    avg <- markov.misc::avg_sops(
      model = model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )
    expect_equal(nrow(attr(avg, "newdata_orig")), nrow(case$baseline))

    withr::local_seed(4401)
    inferred <- markov.misc::inferences(
      avg,
      method = "simulation",
      engine = "mvn",
      n_sim = 3,
      return_draws = TRUE
    )

    draws <- markov.misc::get_draws(inferred)

    expect_equal(attr(inferred, "method"), "simulation")
    expect_equal(attr(inferred, "engine"), "mvn")
    expect_equal(attr(inferred, "n_successful"), 3L)
    expect_false(anyNA(inferred$conf.low))
    expect_false(anyNA(inferred$conf.high))
    expect_true(all(inferred$conf.low <= inferred$conf.high))
    expect_equal(sort(unique(draws$draw_id)), 1:3)

    expect_snapshot_value(
      list(
        result = pipeline_signature(inferred),
        draws = draw_signature(draws)
      ),
      style = "json2",
      cran = TRUE
    )
  })

  test_that("Brownian avg_sops() supports bootstrap inference", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_brownian_pipeline_case()
    model <- fit_pipeline_model(case$data)

    avg <- markov.misc::avg_sops(
      model = model,
      newdata = case$data,
      variables = "tx",
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs,
      id_var = "id"
    )
    expect_equal(nrow(attr(avg, "newdata_orig")), nrow(case$data))

    withr::local_seed(4402)
    inferred <- markov.misc::inferences(
      avg,
      method = "bootstrap",
      n_sim = 2,
      return_draws = TRUE
    )

    draws <- markov.misc::get_draws(inferred)

    expect_equal(attr(inferred, "method"), "bootstrap")
    expect_equal(attr(inferred, "n_boot"), 2)
    expect_equal(attr(inferred, "n_successful"), 2L)
    expect_false(anyNA(inferred$conf.low))
    expect_false(anyNA(inferred$conf.high))
    expect_true(all(inferred$conf.low <= inferred$conf.high))
    expect_equal(sort(unique(draws$draw_id)), 1:2)

    expect_snapshot_value(
      list(
        result = pipeline_signature(inferred),
        draws = draw_signature(draws)
      ),
      style = "json2",
      cran = TRUE
    )
  })

  test_that("Brownian avg_sops() supports score-bootstrap simulation on the fast path", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_brownian_pipeline_case()
    model <- fit_pipeline_model(case$data)
    robust_model <- markov.misc::robcov_vglm(model, cluster = case$data$id)

    avg <- markov.misc::avg_sops(
      model = robust_model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs,
      id_var = "id"
    )
    expect_equal(nrow(attr(avg, "newdata_orig")), nrow(case$baseline))

    withr::local_seed(4403)
    inferred <- markov.misc::inferences(
      avg,
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 3,
      return_draws = TRUE
    )

    draws <- markov.misc::get_draws(inferred)

    expect_equal(attr(inferred, "method"), "simulation")
    expect_equal(attr(inferred, "engine"), "score_bootstrap")
    expect_equal(attr(inferred, "score_weight_dist"), "exponential")
    expect_equal(attr(inferred, "n_successful"), 3L)
    expect_false(anyNA(inferred$conf.low))
    expect_false(anyNA(inferred$conf.high))
    expect_true(all(inferred$std.error >= 0))
    expect_true(any(inferred$std.error > 0))
    expect_equal(sort(unique(draws$draw_id)), 1:3)

    expect_snapshot_value(
      list(
        result = pipeline_signature(inferred),
        draws = draw_signature(draws)
      ),
      style = "json2",
      cran = TRUE
    )
  })

  test_that("inline rcs() and explicit basis agree for MVN simulation inference", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    case <- build_brownian_pipeline_case()
    explicit_model <- fit_pipeline_model(case$data)
    inline_model <- fit_inline_pipeline_model(case$data)

    explicit_avg <- markov.misc::avg_sops(
      model = explicit_model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )
    inline_avg <- markov.misc::avg_sops(
      model = inline_model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      pvarname = "yprev"
    )

    expect_equal(inline_avg$estimate, explicit_avg$estimate, tolerance = 1e-10)

    withr::local_seed(4411)
    explicit_inferred <- markov.misc::inferences(
      explicit_avg,
      method = "simulation",
      engine = "mvn",
      n_sim = 3,
      return_draws = TRUE
    )
    withr::local_seed(4411)
    inline_inferred <- markov.misc::inferences(
      inline_avg,
      method = "simulation",
      engine = "mvn",
      n_sim = 3,
      return_draws = TRUE
    )

    explicit_draws <- markov.misc::get_draws(explicit_inferred)
    inline_draws <- markov.misc::get_draws(inline_inferred)

    expect_equal(inline_inferred$estimate, explicit_inferred$estimate, tolerance = 1e-10)
    expect_equal(inline_inferred$conf.low, explicit_inferred$conf.low, tolerance = 1e-10)
    expect_equal(inline_inferred$conf.high, explicit_inferred$conf.high, tolerance = 1e-10)
    expect_equal(inline_draws$draw, explicit_draws$draw, tolerance = 1e-10)
  })

  test_that("inline rcs() and explicit basis agree for score-bootstrap simulation inference", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_brownian_pipeline_case()
    explicit_model <- fit_pipeline_model(case$data)
    inline_model <- fit_inline_pipeline_model(case$data)
    explicit_robust <- markov.misc::robcov_vglm(explicit_model, cluster = case$data$id)
    inline_robust <- markov.misc::robcov_vglm(inline_model, cluster = case$data$id)

    explicit_avg <- markov.misc::avg_sops(
      model = explicit_robust,
      newdata = case$baseline,
      variables = "tx",
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs,
      id_var = "id"
    )
    inline_avg <- markov.misc::avg_sops(
      model = inline_robust,
      newdata = case$baseline,
      variables = "tx",
      times = 1:8,
      ylevels = case$ylevels,
      absorb = case$absorb,
      pvarname = "yprev",
      id_var = "id"
    )

    expect_equal(inline_avg$estimate, explicit_avg$estimate, tolerance = 1e-10)

    withr::local_seed(4412)
    explicit_inferred <- markov.misc::inferences(
      explicit_avg,
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 3,
      return_draws = TRUE
    )
    withr::local_seed(4412)
    inline_inferred <- markov.misc::inferences(
      inline_avg,
      method = "simulation",
      engine = "score_bootstrap",
      n_sim = 3,
      return_draws = TRUE
    )

    explicit_draws <- markov.misc::get_draws(explicit_inferred)
    inline_draws <- markov.misc::get_draws(inline_inferred)

    expect_equal(inline_inferred$estimate, explicit_inferred$estimate, tolerance = 1e-10)
    expect_equal(inline_inferred$conf.low, explicit_inferred$conf.low, tolerance = 1e-10)
    expect_equal(inline_inferred$conf.high, explicit_inferred$conf.high, tolerance = 1e-10)
    expect_equal(inline_draws$draw, explicit_draws$draw, tolerance = 1e-10)
  })

  test_that("full-PO vglm.markov() with inline rcs matches VGAM::vglm()", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    raw_data <- markov.misc::sim_trajectories_brownian_gap(
      n_patients = 200,
      follow_up_time = 12,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 2375679,
      mu_treatment_effect = 0
    )
    data <- markov.misc::prepare_markov_data(raw_data, absorbing_state = 6)
    data <- add_patient_age(data)
    form <- ordered(y) ~ rms::rcs(time, 4) + tx + age + yprev

    ordinary_fit <- suppressWarnings(
      VGAM::vglm(
        form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data
      )
    )
    markov_fit <- suppressWarnings(
      markov.misc::vglm.markov(
        form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data
      )
    )

    ordinary_fitted <- stats::predict(ordinary_fit, type = "response")
    markov_fitted <- markov.misc:::predict_vglm_response_markov(markov_fit, data)

    expect_s4_class(markov_fit, "vglm")
    expect_true(isTRUE(attr(markov_fit, "markov_vglm")))
    expect_equal(unname(markov_fitted), unname(ordinary_fitted), tolerance = 1e-10)
  })

  test_that("inline rcs vglm.markov() model can use column-level PPO constraints", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    raw_data <- markov.misc::sim_trajectories_brownian_gap(
      n_patients = 200,
      follow_up_time = 12,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 2375679,
      mu_treatment_effect = 0
    )
    data <- markov.misc::prepare_markov_data(raw_data, absorbing_state = 6)
    data <- add_patient_age(data)
    baseline <- data[data$time == 1, , drop = FALSE]
    form <- ordered(y) ~ rms::rcs(time, 4) + tx + age + yprev

    fit_po <- suppressWarnings(
      markov.misc::vglm.markov(
        form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data
      )
    )

    cons <- fit_po@constraints
    linear_time_term <- grep("^rms::rcs\\(time, 4\\).*time$", names(cons), value = TRUE)[[1]]
    n_thresholds <- nrow(cons[[linear_time_term]])
    yprev_terms <- grep("^yprev", names(cons), value = TRUE)
    cons[yprev_terms] <- NULL
    cons[["yprev"]] <- cbind(PO_effect = rep(1, n_thresholds))
    cons[[linear_time_term]] <- cbind(
      PO_effect = rep(1, n_thresholds),
      linear_deviation = seq_len(n_thresholds) - 1
    )

    fit_ppo <- suppressWarnings(
      markov.misc::vglm.markov(
        form,
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
        data = data,
        constraints = cons
      )
    )
    fit_ppo@call$constraints <- cons

    withr::local_seed(42)
    out <- markov.misc::avg_sops(
      markov.misc::robcov_vglm(fit_ppo, cluster = data$id),
      newdata = baseline,
      variables = "tx",
      times = 1:4,
      ylevels = 1:6,
      absorb = "6",
      id_var = "id"
    ) |>
      markov.misc::inferences(
        method = "simulation",
        engine = "mvn",
        n_sim = 2,
        workers = 1
      )

    expect_s3_class(out, "markov_avg_sops")
    expect_true(all(is.finite(out$estimate)))
  })

  test_that("vglm.markov() constrained PPO matches explicit spline basis", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    raw_data <- markov.misc::sim_trajectories_brownian_gap(
      n_patients = 200,
      follow_up_time = 12,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 2375679,
      mu_treatment_effect = 0
    )
    data <- markov.misc::prepare_markov_data(raw_data, absorbing_state = 6)
    data <- add_patient_age(data)
    data <- add_time_basis(data)

    explicit_basis <- unname(as.matrix(data[, c("time_lin", "time_nlin_1", "time_nlin_2")]))
    inline_basis_raw <- rms::rcs(data$time, 4)
    inline_basis <- matrix(
      as.numeric(inline_basis_raw),
      nrow = nrow(inline_basis_raw),
      ncol = ncol(inline_basis_raw)
    )
    expect_equal(inline_basis, explicit_basis, tolerance = 1e-14)

    baseline <- data[data$time == 1, , drop = FALSE]
    t_covs <- get_time_covariates(data)
    ylevels <- 1:6
    n_thresholds <- length(ylevels) - 1

    explicit_constraints <- list(
      "(Intercept)" = diag(n_thresholds),
      time_lin = cbind(
        PO_effect = rep(1, n_thresholds),
        linear_deviation = seq_len(n_thresholds) - 1
      ),
      time_nlin_1 = cbind(PO_effect = rep(1, n_thresholds)),
      time_nlin_2 = cbind(PO_effect = rep(1, n_thresholds)),
      tx = cbind(PO_effect = rep(1, n_thresholds)),
      age = cbind(PO_effect = rep(1, n_thresholds)),
      yprev = cbind(PO_effect = rep(1, n_thresholds))
    )

    explicit_fit <- suppressWarnings(
      VGAM::vglm(
        ordered(y) ~ time_lin + time_nlin_1 + time_nlin_2 + tx + age + yprev,
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
        data = data,
        constraints = explicit_constraints
      )
    )
    explicit_fit@call$constraints <- explicit_constraints

    inline_form <- ordered(y) ~ rms::rcs(time, 4) + tx + age + yprev
    inline_po <- suppressWarnings(
      markov.misc::vglm.markov(
        inline_form,
        family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
        data = data
      )
    )

    inline_constraints <- inline_po@constraints
    inline_time_terms <- grep("^rms::rcs\\(time, 4\\)", names(inline_constraints), value = TRUE)
    expect_length(inline_time_terms, 3)
    inline_linear_time_term <- inline_time_terms[[1]]
    expect_match(inline_linear_time_term, "time$")
    yprev_terms <- grep("^yprev", names(inline_constraints), value = TRUE)
    inline_constraints[yprev_terms] <- NULL
    inline_constraints[["yprev"]] <- cbind(PO_effect = rep(1, n_thresholds))
    inline_constraints[[inline_linear_time_term]] <- cbind(
      PO_effect = rep(1, n_thresholds),
      linear_deviation = seq_len(n_thresholds) - 1
    )

    inline_fit <- suppressWarnings(
      markov.misc::vglm.markov(
        inline_form,
        family = VGAM::cumulative(reverse = TRUE, parallel = FALSE),
        data = data,
        constraints = inline_constraints
      )
    )
    inline_fit@call$constraints <- inline_constraints

    explicit_fitted <- stats::predict(explicit_fit, type = "response")
    inline_fitted <- markov.misc:::predict_vglm_response_markov(inline_fit, data)
    expect_equal(dim(inline_fitted), dim(explicit_fitted))
    expect_equal(unname(inline_fitted), unname(explicit_fitted), tolerance = 1e-10)

    explicit_sops <- markov.misc::avg_sops(
      explicit_fit,
      newdata = baseline,
      variables = "tx",
      times = 1:12,
      ylevels = ylevels,
      absorb = "6",
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = t_covs,
      id_var = "id"
    )
    inline_sops <- markov.misc::avg_sops(
      inline_fit,
      newdata = baseline,
      variables = "tx",
      times = 1:12,
      ylevels = ylevels,
      absorb = "6",
      pvarname = "yprev",
      id_var = "id"
    )

    order_sops <- function(x) {
      x[order(x$tx, x$time, as.integer(as.character(x$state))), , drop = FALSE]
    }
    explicit_sops <- order_sops(explicit_sops)
    inline_sops <- order_sops(inline_sops)

    expect_equal(inline_sops$tx, explicit_sops$tx)
    expect_equal(inline_sops$time, explicit_sops$time)
    expect_equal(inline_sops$state, explicit_sops$state)
    expect_equal(inline_sops$estimate, explicit_sops$estimate, tolerance = 1e-10)
  })
})
