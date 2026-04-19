describe("SOP regression baselines", {
  add_patient_age <- function(data) {
    data$age <- NULL

    age_lookup <- data.frame(
      id = sort(unique(data$id)),
      age = 50 + ((sort(unique(data$id)) * 7) %% 23)
    )

    merge(data, age_lookup, by = "id", all.x = TRUE, sort = FALSE)
  }

  add_explicit_time_basis <- function(data, nk = 4) {
    time_spl <- rms::rcs(data$time, nk)
    data$time_lin <- as.vector(time_spl[, 1])
    data$time_nlin_1 <- as.vector(time_spl[, 2])
    data$time_nlin_2 <- as.vector(time_spl[, 3])
    data
  }

  get_time_covariates <- function(data) {
    t_covs <- unique(data[, c("time", "time_lin", "time_nlin_1", "time_nlin_2")])
    t_covs <- t_covs[order(t_covs$time), , drop = FALSE]
    as.data.frame(t_covs[, c("time_lin", "time_nlin_1", "time_nlin_2"), drop = FALSE])
  }

  fit_explicit_spline_model <- function(data) {
    VGAM::vglm(
      ordered(y) ~ (time_lin + time_nlin_1 + time_nlin_2) * tx + yprev + age,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  }

  fit_inline_spline_model <- function(data) {
    VGAM::vglm(
      ordered(y) ~ rms::rcs(time, 4) * tx + yprev + age,
      family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
      data = data
    )
  }

  sops_signature <- function(
    sops_array,
    patient_idx = 1:3,
    time_idx = c(1, 2, 4, 6, 10),
    state_idx = c(1, 3, dim(sops_array)[3])
  ) {
    patient_idx <- patient_idx[patient_idx <= dim(sops_array)[1]]
    time_idx <- time_idx[time_idx <= dim(sops_array)[2]]
    state_idx <- state_idx[state_idx <= dim(sops_array)[3]]

    out <- expand.grid(
      patient = patient_idx,
      time = time_idx,
      state = state_idx,
      KEEP.OUT.ATTRS = FALSE
    )
    out$prob <- round(
      mapply(
        function(i, t, s) sops_array[i, t, s],
        out$patient,
        out$time,
        out$state
      ),
      8
    )

    out[order(out$patient, out$time, out$state), , drop = FALSE]
  }

  avg_signature <- function(
    result,
    time_idx = c(1, 2, 4, 6, 10),
    state_idx = c(1, 3, max(as.integer(as.character(result$state))))
  ) {
    out <- subset(result, time %in% time_idx & state %in% state_idx)
    out <- out[order(out$tx, out$time, out$state), c("tx", "time", "state", "estimate")]
    out$state <- as.integer(as.character(out$state))
    out$estimate <- round(out$estimate, 8)
    rownames(out) <- NULL
    out
  }

  draw_signature <- function(
    draws,
    draw_id = 1L,
    time_idx = c(1, 2, 4, 6, 8),
    state_idx = c(1, 3, max(as.integer(as.character(draws$state))))
  ) {
    out <- subset(
      draws,
      draw_id == !!draw_id & time %in% time_idx & state %in% state_idx
    )
    out <- out[order(out$tx, out$time, out$state), c("draw_id", "tx", "time", "state", "draw")]
    out$state <- as.integer(as.character(out$state))
    out$draw <- round(out$draw, 8)
    rownames(out) <- NULL
    out
  }

  build_brownian_case <- function() {
    raw_data <- markov.misc::sim_trajectories_brownian(
      n_patients = 80,
      follow_up_time = 10,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 1001,
      mu_treatment_effect = 0.2
    )

    data <- markov.misc::prepare_markov_data(raw_data)
    data <- add_patient_age(data)
    data <- add_explicit_time_basis(data, nk = 4)

    list(
      data = data,
      baseline = subset(data, time == 1),
      t_covs = get_time_covariates(data),
      ylevels = 1:6,
      absorb = "6"
    )
  }

  build_brownian_gap_case <- function() {
    raw_data <- markov.misc::sim_trajectories_brownian_gap(
      n_patients = 80,
      follow_up_time = 10,
      treatment_prob = 0.5,
      absorbing_state = 6,
      seed = 1002,
      mu_treatment_effect = 0.15
    )

    data <- markov.misc::prepare_markov_data(raw_data)
    data <- add_patient_age(data)
    data <- add_explicit_time_basis(data, nk = 4)

    list(
      data = data,
      baseline = subset(data, time == 1),
      t_covs = get_time_covariates(data),
      ylevels = 1:6,
      absorb = "6"
    )
  }

  build_markov_case <- function() {
    baseline <- markov.misc::violet_baseline[1:80, c("id", "yprev", "tx", "age", "sofa")]

    raw_data <- markov.misc::sim_trajectories_markov(
      baseline_data = baseline,
      follow_up_time = 10,
      parameter = log(0.85),
      seed = 1004
    )

    data <- markov.misc::prepare_markov_data(raw_data)
    data <- add_patient_age(data)
    data <- add_explicit_time_basis(data, nk = 4)

    list(
      data = data,
      baseline = subset(data, time == 1),
      t_covs = get_time_covariates(data),
      ylevels = 1:6,
      absorb = "6"
    )
  }

  test_that("brownian SOP signatures remain stable", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")
    skip_if_not_installed("mvtnorm")

    case <- build_brownian_case()
    model <- suppressWarnings(fit_explicit_spline_model(case$data))

    sops_result <- markov.misc::soprob_markov(
      object = model,
      data = case$baseline,
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    avg_result <- markov.misc::avg_sops(
      model = model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    withr::local_seed(5001)
    sim_result <- markov.misc::inferences(
      avg_result,
      method = "simulation",
      n_sim = 2,
      return_draws = TRUE
    )

    expect_snapshot_value(
      list(
        sops = sops_signature(sops_result),
        avg = avg_signature(avg_result),
        sim_draw = draw_signature(markov.misc::get_draws(sim_result))
      ),
      style = "json2",
      cran = TRUE
    )
  })

  test_that("brownian gap SOP signatures remain stable", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_brownian_gap_case()
    model <- suppressWarnings(fit_explicit_spline_model(case$data))

    sops_result <- markov.misc::soprob_markov(
      object = model,
      data = case$baseline,
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    avg_result <- markov.misc::avg_sops(
      model = model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    expect_snapshot_value(
      list(
        sops = sops_signature(sops_result),
        avg = avg_signature(avg_result)
      ),
      style = "json2",
      cran = TRUE
    )
  })

  test_that("markov SOP signatures remain stable", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_markov_case()
    model <- suppressWarnings(fit_explicit_spline_model(case$data))

    sops_result <- markov.misc::soprob_markov(
      object = model,
      data = case$baseline,
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    avg_result <- markov.misc::avg_sops(
      model = model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    expect_snapshot_value(
      list(
        sops = sops_signature(sops_result),
        avg = avg_signature(avg_result)
      ),
      style = "json2",
      cran = TRUE
    )
  })

  test_that("inline and explicit spline slow-path SOPs agree for interaction models", {
    skip_if_not_installed("VGAM")
    skip_if_not_installed("rms")

    case <- build_brownian_case()
    explicit_model <- fit_explicit_spline_model(case$data)
    inline_model <- fit_inline_spline_model(case$data)

    explicit_sops <- markov.misc::soprob_markov(
      object = explicit_model,
      data = case$baseline[1:10, , drop = FALSE],
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb,
      tvarname = "time_lin",
      pvarname = "yprev",
      t_covs = case$t_covs
    )

    inline_sops <- markov.misc::soprob_markov(
      object = inline_model,
      data = case$baseline[1:10, , drop = FALSE],
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb
    )

    explicit_avg <- markov.misc::avg_sops(
      model = explicit_model,
      newdata = case$baseline,
      variables = "tx",
      times = 1:10,
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
      times = 1:10,
      ylevels = case$ylevels,
      absorb = case$absorb
    )

    expect_equal(inline_sops, explicit_sops, tolerance = 1e-10)
    expect_equal(inline_avg$estimate, explicit_avg$estimate, tolerance = 1e-10)
  })
})
