test_that("analytical averaging checks scenario and starting-profile order", {
  skip_if_not_installed("VGAM")
  skip_if_not_installed("rms")
  data <- make_test_data(n_patients = 48, follow_up_time = 6, seed = 2711)
  model <- make_test_model(data)
  baseline <- data[!duplicated(data$id), , drop = FALSE][1:8, , drop = FALSE]
  baseline$metadata <- I(cbind(first = 1:8, second = 8:1))
  variables <- list(tx = c(0, 1))
  average <- avg_sops(
    model,
    newdata = baseline,
    variables = variables,
    times = 1:3,
    y_levels = 1:6,
    absorb = 6
  )
  inferred <- inferences(average, method = "delta", vcov = stats::vcov(model))
  newdata <- attr(average, "newdata_pred")
  plan <- delta_compile_plan(
    average,
    model,
    newdata,
    times = 1:3,
    output = "average"
  )
  full <- run_sop_delta_plan(plan, model)
  reference <- delta_average_arrays(
    full$probabilities,
    full$jacobian,
    grid = do.call(expand.grid, variables),
    times = 1:3,
    y_levels = 1:6,
    variables = variables,
    n_each = nrow(baseline)
  )
  index <- delta_match_cells(average, reference$cells, c("time", "state", "tx"))
  expect_equal(
    average$estimate,
    reference$cells$estimate[index],
    tolerance = 1e-12
  )
  expect_equal(
    unname(get_jacobian(inferred)),
    unname(reference$jacobian[index, , drop = FALSE]),
    tolerance = 1e-12
  )

  n <- nrow(baseline)
  reordered <- average
  attr(reordered, "newdata_pred") <- newdata[
    c(rev(seq_len(n)), n + rev(seq_len(n))),
    ,
    drop = FALSE
  ]
  expect_equal(
    inferences(
      reordered,
      method = "delta",
      vcov = stats::vcov(model)
    )$std.error,
    inferred$std.error,
    tolerance = 1e-12
  )
  layouts <- list(
    interleaved = as.vector(rbind(seq_len(n), n + seq_len(n))),
    reversed_scenarios = c(n + seq_len(n), seq_len(n)),
    reversed_second_profile_order = c(seq_len(n), n + rev(seq_len(n)))
  )
  for (rows in layouts) {
    malformed <- average
    attr(malformed, "newdata_pred") <- newdata[rows, , drop = FALSE]
    expect_error(
      inferences(malformed, method = "delta", vcov = stats::vcov(model)),
      "same starting-profile order"
    )
  }
})
