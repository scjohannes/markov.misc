make_fake_blrm <- function(draws = NULL, gamma_draws = NULL, pppo = 0L) {
  if (is.null(draws)) {
    draws <- rbind(
      c(0.50, -0.50, 0.25),
      c(1.00, -1.00, -0.25),
      c(0.25, -0.75, 0.10)
    )
    colnames(draws) <- c("y>=2", "y>=3", "tx")
  }

  tau_info <- if (pppo > 0) {
    data.frame(name = "tau_tx")
  } else {
    data.frame(name = character())
  }

  structure(
    list(
      draws = draws,
      non.slopes = 2L,
      pppo = pppo,
      tauInfo = tau_info,
      ylevels = 1:3,
      yname = "y",
      clusterInfo = list(name = "id", cluster = c("a", "b")),
      gamma_draws = gamma_draws
    ),
    class = c("blrm", "orm")
  )
}

fake_blrm_design <- function(model, newdata, second = FALSE) {
  if (second) {
    out <- matrix(newdata$tx, ncol = 1)
    colnames(out) <- "tau_tx"
    return(out)
  }

  out <- matrix(newdata$tx, ncol = 1)
  colnames(out) <- "tx"
  out
}

test_that("soprob_markov handles second-order recursion and absorbing states", {
  model <- structure(list(), class = "vglm")

  transition <- function(h, j) {
    if (j == 3) return(c(0, 0, 1))
    if (h == 1 && j == 1) return(c(0.2, 0.8, 0.0))
    if (h == 1 && j == 2) return(c(0.1, 0.6, 0.3))
    if (h == 2 && j == 1) return(c(0.5, 0.5, 0.0))
    if (h == 2 && j == 2) return(c(0.0, 0.7, 0.3))
    c(0, 0, 1)
  }

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      out <- matrix(NA_real_, nrow = nrow(newdata), ncol = 3)
      for (i in seq_len(nrow(newdata))) {
        h <- as.integer(as.character(newdata$ypprev[i]))
        j <- as.integer(as.character(newdata$yprev[i]))
        out[i, ] <- transition(h, j)
      }
      colnames(out) <- as.character(1:3)
      out
    },
    {
      data <- data.frame(
        id = 1,
        time = 1,
        ypprev = factor(1, levels = 1:3),
        yprev = factor(1, levels = 1:3)
      )

      out <- soprob_markov(
        object = model,
        data = data,
        times = 1:2,
        ylevels = 1:3,
        absorb = 3,
        p2varname = "ypprev"
      )

      absorb_data <- data
      absorb_data$ypprev <- factor(3, levels = 1:3)
      absorb_data$yprev <- factor(3, levels = 1:3)
      absorb_out <- soprob_markov(
        object = model,
        data = absorb_data,
        times = 1:2,
        ylevels = 1:3,
        absorb = 3,
        p2varname = "ypprev"
      )
    }
  )

  expect_equal(out[1, 1, ], c("1" = 0.2, "2" = 0.8, "3" = 0), tolerance = 1e-12)
  expect_equal(out[1, 2, ], c("1" = 0.12, "2" = 0.64, "3" = 0.24), tolerance = 1e-12)
  expect_equal(absorb_out[1, 1, ], c("1" = 0, "2" = 0, "3" = 1), tolerance = 1e-12)
  expect_equal(absorb_out[1, 2, ], c("1" = 0, "2" = 0, "3" = 1), tolerance = 1e-12)
})

test_that("soprob_markov handles single first-order time points", {
  model <- structure(list(), class = "vglm")

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      out <- matrix(c(0.2, 0.3, 0.5), nrow = nrow(newdata), ncol = 3, byrow = TRUE)
      colnames(out) <- as.character(1:3)
      out
    },
    {
      out <- soprob_markov(
        object = model,
        data = data.frame(
          id = 1,
          time = 1,
          yprev = factor(1, levels = 1:3)
        ),
        times = 1,
        ylevels = 1:3,
        absorb = 3
      )
    }
  )

  expect_equal(dim(out), c(1L, 1L, 3L))
  expect_equal(out[1, 1, ], c("1" = 0.2, "2" = 0.3, "3" = 0.5), tolerance = 1e-12)
})

test_that("soprob_markov carries absorbing labels that are not column positions", {
  model <- structure(list(), class = "vglm")

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      out <- matrix(0, nrow = nrow(newdata), ncol = 3)
      colnames(out) <- c("0", "1", "2")

      if (nrow(newdata) == 1L) {
        out[1, ] <- c(0.1, 0.2, 0.7)
        return(out)
      }

      for (i in seq_len(nrow(newdata))) {
        out[i, ] <- switch(
          as.character(newdata$yprev[i]),
          "0" = c(1, 0, 0),
          "1" = c(0, 1, 0),
          "2" = c(0, 0, 1)
        )
      }
      out
    },
    {
      out <- soprob_markov(
        object = model,
        data = data.frame(
          id = 1,
          time = 1,
          yprev = factor(0, levels = 0:2)
        ),
        times = 1:2,
        ylevels = 0:2,
        absorb = 2
      )
    }
  )

  expect_equal(out[1, 2, ], c("0" = 0.1, "1" = 0.2, "2" = 0.7), tolerance = 1e-12)
})

test_that("blrm posterior draw sampling is random, capped, and reproducible", {
  model <- make_fake_blrm(draws = matrix(0, nrow = 150, ncol = 3))

  draw_ids_a <- markov.misc:::select_posterior_draws(model, n_draws = 100L, seed = 11)
  draw_ids_b <- markov.misc:::select_posterior_draws(model, n_draws = 100L, seed = 11)

  expect_length(draw_ids_a, 100)
  expect_identical(draw_ids_a, draw_ids_b)
  expect_false(identical(draw_ids_a, seq_len(100)))
  expect_identical(
    markov.misc:::select_posterior_draws(model, n_draws = NULL),
    seq_len(150)
  )
})

test_that("standardize_sops() reuses blrm posterior draws across arms", {
  model <- make_fake_blrm()
  data <- data.frame(
    id = c("a", "b"),
    time = c(1, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    tx = c(0, 1)
  )
  captured_draws <- list()
  select_calls <- 0L

  with_mocked_bindings(
    select_posterior_draws = function(model, n_draws = 100L, seed = NULL) {
      select_calls <<- select_calls + 1L
      c(3L, 1L)
    },
    soprob_markov = function(...) {
      args <- list(...)
      captured_draws[[length(captured_draws) + 1L]] <<- args$.draw_indices
      n_draws <- length(args$.draw_indices)
      n_pat <- nrow(args$data)
      n_times <- length(args$times)
      n_states <- length(args$ylevels)
      array(
        1 / n_states,
        dim = c(n_draws, n_pat, n_times, n_states)
      )
    },
    {
      out <- standardize_sops(
        model,
        data = data,
        times = 1:2,
        ylevels = 1:3,
        absorb = 3,
        n_draws = 2
      )
    }
  )

  expect_equal(select_calls, 1L)
  expect_length(captured_draws, 2)
  expect_identical(captured_draws[[1]], c(3L, 1L))
  expect_identical(captured_draws[[2]], c(3L, 1L))
  expect_equal(dim(out$sop_tx), c(2L, 2L, 3L))
  expect_equal(dim(out$sop_ctrl), c(2L, 2L, 3L))
})

test_that("manual blrm prediction supports PO, constrained PPO, and random effects", {
  gamma <- rbind(
    c(0.50, -0.50),
    c(1.00, -1.00),
    c(0.25, -0.25)
  )
  colnames(gamma) <- c("a", "b")
  model <- make_fake_blrm(gamma_draws = gamma)
  newdata <- data.frame(
    id = c("a", "b"),
    tx = c(1, 0),
    yprev = factor(c(1, 1), levels = 1:3)
  )

  with_mocked_bindings(
    blrm_design_matrix = fake_blrm_design,
    {
      no_re <- markov.misc:::predict_blrm_response_markov(
        model,
        newdata,
        draw_indices = 1:2
      )
      with_re <- markov.misc:::predict_blrm_response_markov(
        model,
        newdata,
        include_re = TRUE,
        draw_indices = 1:2
      )

      ppo_draws <- cbind(model$draws, tau_tx = c(0.10, -0.20, 0.15))
      ppo_model <- make_fake_blrm(draws = ppo_draws, pppo = 1L)
      ppo_model$cppo <- function(y) as.numeric(y) - 1
      ppo <- markov.misc:::predict_blrm_response_markov(
        ppo_model,
        newdata[1, , drop = FALSE],
        draw_indices = 1
      )

      cppo_string_test <- function(y) as.numeric(y) - 1
      ppo_model$cppo <- "cppo_string_test"
      ppo_string <- markov.misc:::predict_blrm_response_markov(
        ppo_model,
        newdata[1, , drop = FALSE],
        draw_indices = 1
      )

      ppo_model$cppo <- "function(y) as.numeric(y) - 1"
      cppo_error <- try(
        markov.misc:::predict_blrm_response_markov(
          ppo_model,
          newdata[1, , drop = FALSE],
          draw_indices = 1
        ),
        silent = TRUE
      )

      unknown <- newdata[1, , drop = FALSE]
      unknown$id <- "z"
      re_error <- try(
        markov.misc:::predict_blrm_response_markov(
          model,
          unknown,
          include_re = TRUE,
          draw_indices = 1
        ),
        silent = TRUE
      )
    }
  )

  expected_p1 <- 1 - stats::plogis(0.50 + 0.25)
  expected_p2 <- stats::plogis(0.50 + 0.25) - stats::plogis(-0.50 + 0.25)
  expected_p3 <- stats::plogis(-0.50 + 0.25)

  expect_equal(dim(no_re), c(2L, 2L, 3L))
  expect_equal(unname(no_re[1, 1, ]), c(expected_p1, expected_p2, expected_p3), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(no_re, with_re)))
  expect_equal(sum(ppo[1, 1, ]), 1, tolerance = 1e-12)
  expect_equal(ppo_string, ppo, tolerance = 1e-12)
  expect_s3_class(cppo_error, "try-error")
  expect_match(as.character(cppo_error), "Inline expressions are not supported")
  expect_s3_class(re_error, "try-error")
  expect_match(as.character(re_error), "not present")
})

test_that("sops() summarizes blrm posterior SOP draws and stores optional draws", {
  model <- make_fake_blrm()
  newdata <- data.frame(
    id = "a",
    tx = 1,
    yprev = factor(1, levels = 1:3),
    time = 1
  )

  with_mocked_bindings(
    blrm_design_matrix = fake_blrm_design,
    {
      result <- sops(
        model,
        newdata = newdata,
        times = 1:2,
        n_draws = 2,
        seed = 4,
        posterior_summary = "mean",
        return_draws = TRUE
      )
      inferred <- inferences(result, method = "bootstrap", n_sim = 2)
    }
  )

  expect_s3_class(result, "markov_sops")
  expect_contains(names(result), c("estimate", "conf.low", "conf.high", "std.error"))
  expect_equal(attr(result, "n_draws"), 2L)
  expect_length(attr(result, "draw_ids"), 2)
  expect_s3_class(attr(result, "draws"), "data.frame")
  expect_identical(inferred, result)

  draws <- get_draws(result)
  expect_contains(names(draws), c("draw_id", "draw"))

  first_cell <- attr(result, "draws")
  first_cell <- first_cell[first_cell$time == 1 & first_cell$state == 1, ]
  expect_equal(result$estimate[result$time == 1 & result$state == 1], mean(first_cell$estimate))
  expect_s3_class(plot_sops(result, geom = "line"), "ggplot")
})

test_that("avg_sops() streams and summarizes blrm posterior draws", {
  gamma <- rbind(
    c(0.50, -0.50),
    c(1.00, -1.00),
    c(0.25, -0.25)
  )
  colnames(gamma) <- c("a", "b")
  model <- make_fake_blrm(gamma_draws = gamma)
  model$clusterInfo <- list(name = "subject", cluster = c("a", "b"))
  newdata <- data.frame(
    subject = c("a", "b"),
    tx = c(0, 1),
    yprev = factor(c(1, 2), levels = 1:3),
    time = 1
  )

  gamma_calls <- 0L
  with_mocked_bindings(
    blrm_design_matrix = fake_blrm_design,
    get_blrm_gamma_draws = function(model, draw_indices) {
      gamma_calls <<- gamma_calls + 1L
      model$gamma_draws[draw_indices, , drop = FALSE]
    },
    {
      withr::local_options(list(markov.misc.blrm_avg_chunk_size = 1L))
      result <- avg_sops(
        model,
        newdata = newdata,
        variables = list(tx = c(0, 1)),
        times = 1:2,
        include_re = TRUE,
        n_draws = 3,
        seed = 12,
        posterior_summary = "median",
        return_draws = TRUE
      )
    }
  )

  expect_s3_class(result, "markov_avg_sops")
  expect_contains(names(result), c("estimate", "conf.low", "conf.high", "std.error"))
  expect_equal(attr(result, "n_draws"), 3L)
  expect_equal(gamma_calls, 1L)
  expect_s3_class(attr(result, "draws"), "data.frame")
  expect_contains(names(get_draws(result)), c("draw_id", "draw"))
  draw_sums <- stats::aggregate(estimate ~ draw_id + tx + time, attr(result, "draws"), sum)
  expect_equal(draw_sums$estimate, rep(1, nrow(draw_sums)), tolerance = 1e-12)
  expect_s3_class(plot_sops(result, geom = "line", facet_var = "tx"), "ggplot")
})
