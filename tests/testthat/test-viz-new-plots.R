test_that("plot_transitions creates empirical joint-proportion heatmaps", {
  data <- data.frame(
    id = rep(1:4, each = 2),
    week = rep(1:2, times = 4),
    yprev = c(1, 1, 1, 2, 2, 2, 2, 3),
    y = c(1, 2, 2, 2, 2, 3, 3, 3),
    tx = rep(c(0, 1), each = 4)
  )

  plot <- plot_transitions(
    data,
    time_var = "week",
    y_var = "y",
    pvarname = "yprev"
  )

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$layers[[1]]$geom, "GeomTile")
  expect_equal(plot$labels$fill, "Transition proportion")
  expect_equal(
    as.numeric(tapply(plot$data$estimate, plot$data$.time_key, sum)),
    c(1, 1)
  )
})

test_that("plot_transitions orders numeric-looking time facets numerically", {
  data <- data.frame(
    id = rep(1:2, each = 12),
    week = rep(1:12, times = 2),
    yprev = rep(1, 24),
    y = rep(1, 24)
  )

  data$week <- as.character(data$week)
  character_plot <- plot_transitions(
    data,
    time_var = "week",
    y_var = "y",
    pvarname = "yprev"
  )

  data$week <- factor(data$week)
  factor_plot <- plot_transitions(
    data,
    time_var = "week",
    y_var = "y",
    pvarname = "yprev"
  )

  expected <- paste0("Time ", 1:12)
  expect_equal(levels(character_plot$data$.panel), expected)
  expect_equal(levels(factor_plot$data$.panel), expected)
})

test_that("plot_transitions orders sparse time facets chronologically", {
  data <- data.frame(
    id = rep(1:2, each = 3),
    week = rep(c(7, 14, 1), times = 2),
    yprev = rep(1, 6),
    y = rep(1, 6)
  )

  plot <- plot_transitions(
    data,
    times = c(1, 7, 14),
    time_var = "week",
    y_var = "y",
    pvarname = "yprev"
  )

  expect_equal(levels(plot$data$.panel), paste0("Time ", c(1, 7, 14)))
})

test_that("model transition diagnostics simulate full sparse numeric grids", {
  model <- structure(list(), class = "vglm")
  newdata <- data.frame(
    id = 1,
    time = 0,
    yprev = factor(1, levels = 1:2)
  )
  predicted_times <- numeric()

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      predicted_times <<- c(predicted_times, unique(newdata$time))
      out <- matrix(c(1, 0), nrow = nrow(newdata), ncol = 2, byrow = TRUE)
      colnames(out) <- c("1", "2")
      out
    },
    {
      plot <- plot_transitions(
        model,
        newdata = newdata,
        times = c(1, 7, 14),
        ylevels = 1:2,
        n_rep = 1,
        seed = 1,
        show_values = FALSE
      )
    }
  )

  expect_equal(predicted_times, 1:14)
  expect_equal(levels(plot$data$.panel), paste0("Time ", c(1, 7, 14)))
  expect_setequal(as.character(plot$data$.time_key), c("1", "7", "14"))
})

test_that("second-order model diagnostics simulate full sparse numeric grids", {
  model <- structure(list(), class = "vglm")
  newdata <- data.frame(
    id = 1,
    time = 0,
    ypprev = factor(1, levels = 1:2),
    yprev = factor(1, levels = 1:2)
  )
  predicted_times <- numeric()

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      predicted_times <<- c(predicted_times, unique(newdata$time))
      out <- matrix(c(1, 0), nrow = nrow(newdata), ncol = 2, byrow = TRUE)
      colnames(out) <- c("1", "2")
      out
    },
    {
      plot <- plot_transitions(
        model,
        newdata = newdata,
        times = c(1, 7, 14),
        ylevels = 1:2,
        p2varname = "ypprev",
        n_rep = 1,
        seed = 1,
        show_values = FALSE
      )
    }
  )

  expect_equal(predicted_times, 1:14)
  expect_equal(levels(plot$data$.panel), paste0("Time ", c(1, 7, 14)))
  expect_setequal(as.character(plot$data$.time_key), c("1", "7", "14"))
})

test_that("model diagnostics expand sparse factor visit grids", {
  time_info <- list(
    is_factor = TRUE,
    levels = paste0("visit", 1:4),
    ordered = TRUE
  )
  requested <- factor(
    c("visit2", "visit4"),
    levels = time_info$levels,
    ordered = TRUE
  )

  out <- complete_plot_simulation_times(requested, time_info)

  expect_equal(as.character(out), paste0("visit", 1:4))
  expect_equal(levels(out), time_info$levels)
  expect_equal(is.ordered(out), TRUE)
})

test_that("plot_transitions simulates second-order treatment differences", {
  model <- structure(list(), class = "vglm")
  newdata <- data.frame(
    id = 1:2,
    time = 0,
    ypprev = factor(1, levels = 1:2),
    yprev = factor(1, levels = 1:2),
    tx = 0
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      out <- cbind(
        "1" = ifelse(newdata$tx == 0, 1, 0),
        "2" = ifelse(newdata$tx == 0, 0, 1)
      )
      out
    },
    {
      plot <- plot_transitions(
        model,
        newdata = newdata,
        times = 1,
        variables = list(tx = c(0, 1)),
        comparison = "difference",
        ylevels = 1:2,
        p2varname = "ypprev",
        n_rep = 1,
        seed = 1
      )
    }
  )

  diff_12 <- plot$data$estimate[
    plot$data$previous_state == "1" & plot$data$state == "2"
  ]
  diff_11 <- plot$data$estimate[
    plot$data$previous_state == "1" & plot$data$state == "1"
  ]

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$fill, "Difference in transition proportion")
  expect_equal(diff_12, 1)
  expect_equal(diff_11, -1)
})

test_that("plot_transitions carries absorbing states absent from yprev levels", {
  model <- structure(list(), class = "vglm")
  newdata <- data.frame(
    id = 1,
    time = 0,
    yprev = factor(1, levels = 1)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      if (anyNA(newdata$yprev)) {
        stop("missing previous state")
      }
      out <- matrix(c(0, 1), nrow = nrow(newdata), ncol = 2, byrow = TRUE)
      colnames(out) <- c("1", "2")
      out
    },
    {
      plot <- plot_transitions(
        model,
        newdata = newdata,
        times = 1:2,
        ylevels = 1:2,
        absorb = 2,
        n_rep = 1,
        seed = 1
      )
    }
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(sum(plot$data$n), 2)
})

test_that("plot_transitions uses manual blrm predictions averaged over draws", {
  model <- structure(
    list(
      ylevels = 1:2,
      draws = matrix(0, nrow = 3, ncol = 1)
    ),
    class = "blrm"
  )
  newdata <- data.frame(
    time = 0,
    yprev = factor(1, levels = 1:2)
  )
  calls <- list()
  sampled_probabilities <- NULL

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_blrm_response_markov = function(
      object,
      newdata,
      include_re,
      id_var,
      draw_indices,
      gamma_draws
    ) {
      calls[[length(calls) + 1L]] <<- list(
        include_re = include_re,
        id_var = id_var,
        draw_indices = draw_indices,
        gamma_draws = gamma_draws
      )
      out <- array(
        NA_real_,
        dim = c(length(draw_indices), nrow(newdata), 2L),
        dimnames = list(draw_indices, rownames(newdata), c("1", "2"))
      )
      out[,, 1] <- matrix(
        rep(c(1, 0, 0), nrow(newdata)),
        nrow = length(draw_indices)
      )
      out[,, 2] <- matrix(
        rep(c(0, 1, 1), nrow(newdata)),
        nrow = length(draw_indices)
      )
      out
    },
    markov_sample_states = function(probabilities, ylevel_names) {
      sampled_probabilities <<- probabilities
      rep("2", nrow(probabilities))
    },
    {
      plot <- plot_transitions(
        model,
        newdata = newdata,
        times = 1,
        ylevels = 1:2,
        n_rep = 1,
        seed = 1
      )
    }
  )

  estimate <- plot$data$estimate[
    plot$data$previous_state == "1" & plot$data$state == "2"
  ]

  expect_equal(estimate, 1)
  expect_length(calls, 1)
  expect_equal(sampled_probabilities[, "2"], 2 / 3)
  expect_false(calls[[1]]$include_re)
  expect_equal(calls[[1]]$id_var, "id")
  expect_equal(calls[[1]]$draw_indices, 1:3)
  expect_null(calls[[1]]$gamma_draws)
})

test_that("plot_correlation creates ordinal correlation heatmaps", {
  data <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(1:3, times = 4),
    y = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 2, 3, 3)
  )

  plot <- plot_correlation(data, facet_var = NULL)

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$layers[[1]]$geom, "GeomTile")
  expect_equal(nrow(plot$data), 3L)
  expect_equal(levels(plot$data$time_1), c("1", "2", "3"))
})

test_that("plot_variogram plots correlations by time difference", {
  data <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(1:3, times = 4),
    y = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 2, 3, 3)
  )

  plot <- plot_variogram(data, smooth = FALSE)

  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot$layers[[1]]$geom, "GeomPoint")
  expect_equal(plot$labels$x, "Absolute Time Difference")
  expect_equal(plot$scales$get_scales("y")$limits, c(0, 1))
})

test_that("plot_lp_difference plots profile-based contrasts by previous state", {
  model <- structure(list(), class = "vglm")
  profile <- data.frame(
    tx = 0,
    age = 59,
    yprev = factor(1, levels = 1:2)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    markov_linear_predictor_matrix = function(model, newdata) {
      cbind(
        eta1 = as.numeric(newdata$tx) + as.numeric(as.character(newdata$yprev)),
        eta2 = 2 * as.numeric(newdata$tx)
      )
    },
    {
      plot <- plot_lp_difference(
        model,
        profile = profile,
        variables = list(tx = c(0, 1)),
        times = 1:2,
        previous_states = 1:2,
        ylevels = 1:2
      )
    }
  )

  eta1 <- plot$data$estimate[plot$data$threshold == "eta1"]
  eta2 <- plot$data$estimate[plot$data$threshold == "eta2"]

  expect_s3_class(plot, "ggplot")
  expect_equal(unique(eta1), 1)
  expect_equal(unique(eta2), 2)
})

test_that("plot_lp_difference supports blrm posterior median linear predictors", {
  model <- structure(list(ylevels = 1:3), class = "blrm")
  profile <- data.frame(
    tx = 0,
    yprev = factor(1, levels = 1:3)
  )
  calls <- list()

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    plot_blrm_predict = function(model, newdata, type, ...) {
      args <- list(...)
      calls[[length(calls) + 1L]] <<- c(list(type = type), args)
      10 * args$kint + as.numeric(newdata$tx)
    },
    {
      plot <- plot_lp_difference(
        model,
        profile = profile,
        variables = list(tx = c(0, 1)),
        times = 1:2,
        previous_states = 1:2,
        ylevels = 1:3
      )
    }
  )

  eta1 <- plot$data$estimate[plot$data$threshold == "eta1"]
  eta2 <- plot$data$estimate[plot$data$threshold == "eta2"]

  expect_equal(vapply(calls, `[[`, integer(1), "kint"), 1:2)
  expect_equal(vapply(calls, `[[`, character(1), "type"), c("lp", "lp"))
  expect_equal(
    vapply(calls, `[[`, character(1), "posterior.summary"),
    c("median", "median")
  )
  expect_identical(vapply(calls, `[[`, logical(1), "cint"), c(FALSE, FALSE))
  expect_equal(unique(eta1), 1)
  expect_equal(unique(eta2), 1)
})
