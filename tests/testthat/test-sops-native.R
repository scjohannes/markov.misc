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
