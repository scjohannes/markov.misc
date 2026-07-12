test_that("native first-order propagation matches the reference update", {
  initial <- matrix(c(0.6, 0.4, 0, 0.2, 0.7, 0.1), nrow = 2, byrow = TRUE)
  transition <- rbind(
    matrix(c(0.7, 0.2, 0.1, 0.6, 0.3, 0.1), nrow = 2, byrow = TRUE),
    matrix(c(0.2, 0.6, 0.2, 0.1, 0.7, 0.2), nrow = 2, byrow = TRUE)
  )

  actual <- markov.misc:::markov_native_run(
    initial,
    list(transition),
    non_absorb = 1:2,
    absorb = 3L
  )
  expected <- matrix(0, nrow = 2, ncol = 3)
  expected <- expected + transition[1:2, ] * initial[, 1]
  expected <- expected + transition[3:4, ] * initial[, 2]
  expected[, 3] <- expected[, 3] + initial[, 3]

  expect_equal(actual[, 1, ], initial)
  expect_equal(actual[, 2, ], expected, tolerance = 1e-15)
})

test_that("native posterior reduction averages within groups", {
  values <- array(seq_len(2 * 4 * 2 * 2), dim = c(2, 4, 2, 2))
  actual <- markov.misc:::reduce_sops_draw_array(values, c(1, 1, 2, 2), 2)

  expect_equal(actual[, 1, , ], apply(values[, 1:2, , ], c(1, 3, 4), mean))
  expect_equal(actual[, 2, , ], apply(values[, 3:4, , ], c(1, 3, 4), mean))
})

test_that("native logit updates preserve origin-wise normalization", {
  previous <- matrix(c(0.7, 0.3, 0, 0.2, 0.7, 0.1), nrow = 2, byrow = TRUE)
  logits <- rbind(
    c(1.5, -0.5),
    c(0.8, -0.3),
    c(1.0, -1.0),
    c(0.4, -0.8)
  )
  transitions <- markov.misc:::lp_to_probs(logits, 2L)
  expected <- matrix(0, nrow = 2, ncol = 3)
  expected <- expected + transitions[1:2, ] * previous[, 1L]
  expected <- expected + transitions[3:4, ] * previous[, 2L]
  expected[, 3L] <- expected[, 3L] + previous[, 3L]

  actual <- markov.misc:::markov_update_logits_native(
    previous,
    logits,
    non_absorb = 1:2,
    absorb = 3L
  )
  expect_equal(actual, expected, tolerance = 1e-15)
})

test_that("native PO updates match the general logit kernel", {
  previous <- matrix(c(0.7, 0.3, 0, 0.2, 0.7, 0.1), nrow = 2, byrow = TRUE)
  scalar <- c(-0.4, 0.2, 0.8, -0.1)
  cutpoints <- c(1.2, -0.6)
  logits <- outer(scalar, cutpoints, "+")

  expected <- markov.misc:::markov_update_logits_native(
    previous,
    logits,
    non_absorb = 1:2,
    absorb = 3L
  )
  actual <- markov.misc:::markov_update_po_native(
    previous,
    scalar,
    cutpoints,
    non_absorb = 1:2,
    absorb = 3L
  )
  expect_equal(actual, expected, tolerance = 1e-15)
})

test_that("native second-order updates match direct pair accumulation", {
  previous <- array(seq_len(2 * 3 * 3), dim = c(2, 3, 3)) / 100
  older <- c(1L, 2L, 1L)
  current <- c(1L, 1L, 2L)
  transition <- matrix(
    seq_len(2 * length(older) * 3) / 50,
    nrow = 2 * length(older)
  )
  expected <- array(0, dim = c(2, 3, 3))
  for (pair in seq_along(older)) {
    rows <- (pair - 1L) * 2L + seq_len(2L)
    for (target in seq_len(3L)) {
      expected[, current[pair], target] <-
        expected[, current[pair], target] +
        previous[, older[pair], current[pair]] * transition[rows, target]
    }
  }
  expected[, 3L, 3L] <- expected[, 3L, 3L] +
    apply(previous[,, 3L, drop = FALSE], 1L, sum)

  actual <- markov.misc:::markov_update_second_order_native(
    previous,
    transition,
    older,
    current,
    absorb = 3L
  )
  expect_equal(actual, expected, tolerance = 1e-15)
})

test_that("native BLRM probability conversion preserves R semantics", {
  base_eta <- matrix(c(-2, 1, 0.5, -0.2), nrow = 2)
  intercepts <- matrix(c(2, 1, 0, -1), nrow = 2)
  threshold_eta <- matrix(c(0.2, -0.1, 0.4, 0.3), nrow = 2)
  threshold_scale <- c(0.5, -0.25)
  cumulative <- array(NA_real_, dim = c(2, 2, 2))
  for (threshold in seq_len(2L)) {
    cumulative[,, threshold] <- stats::plogis(
      sweep(base_eta, 1L, intercepts[, threshold], "+") +
        threshold_scale[threshold] * threshold_eta
    )
  }
  expected <- array(NA_real_, dim = c(2, 2, 3))
  expected[,, 1L] <- 1 - cumulative[,, 1L]
  expected[,, 2L] <- cumulative[,, 1L] - cumulative[,, 2L]
  expected[,, 3L] <- cumulative[,, 2L]
  expected[expected < 0] <- 0
  expected <- markov.misc:::normalize_probability_array(expected)

  actual <- markov.misc:::blrm_probabilities_native(
    base_eta,
    intercepts,
    threshold_eta,
    threshold_scale
  )
  expect_equal(actual, expected, tolerance = 1e-15)
})

test_that("PO structure detection requires threshold-invariant slopes", {
  X <- cbind("(Intercept)" = 1, x = c(-1, 0, 1))
  po <- markov.misc:::markov_po_structure(
    rbind(c(1, 0.5), c(-1, 0.5)),
    colnames(X),
    X
  )
  non_po <- markov.misc:::markov_po_structure(
    rbind(c(1, 0.5), c(-1, 0.7)),
    colnames(X),
    X
  )

  expect_equal(po$cutpoints, c(1, -1))
  expect_equal(unname(po$beta), c(0, 0.5))
  expect_null(non_po)
})

test_that("vgam fits are no longer supported", {
  model <- structure(list(), class = "vgam")
  expect_snapshot(
    error = TRUE,
    soprob_markov(
      model,
      newdata = data.frame(time = 1, yprev = 1),
      times = 1,
      y_levels = 1:2
    )
  )
})
