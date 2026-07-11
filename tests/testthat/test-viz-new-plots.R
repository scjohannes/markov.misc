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
    p_var = "yprev"
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
    p_var = "yprev"
  )

  data$week <- factor(data$week)
  factor_plot <- plot_transitions(
    data,
    time_var = "week",
    y_var = "y",
    p_var = "yprev"
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
    p_var = "yprev"
  )

  expect_equal(levels(plot$data$.panel), paste0("Time ", c(1, 7, 14)))
})

test_that("model transition diagnostics evaluate full sparse numeric grids", {
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
        y_levels = 1:2,
        seed = 1,
        show_values = FALSE
      )
    }
  )

  expect_equal(predicted_times, 1:14)
  expect_equal(levels(plot$data$.panel), paste0("Time ", c(1, 7, 14)))
  expect_setequal(as.character(plot$data$.time_key), c("1", "7", "14"))
})

test_that("second-order model diagnostics evaluate full sparse numeric grids", {
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
        y_levels = 1:2,
        p2_var = "ypprev",
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

  out <- complete_plot_recursion_times(requested, time_info)

  expect_equal(as.character(out), paste0("visit", 1:4))
  expect_equal(levels(out), time_info$levels)
  expect_equal(is.ordered(out), TRUE)
})

test_that("transition traces preserve first- and second-order SOP marginals", {
  model <- structure(list(), class = "vglm")
  first_order <- data.frame(
    id = 1:2,
    time = 0,
    yprev = factor(c(1, 2), levels = 1:2)
  )
  second_order <- data.frame(
    id = 1:2,
    time = 0,
    ypprev = factor(c(1, 2), levels = 1:2),
    yprev = factor(c(1, 2), levels = 1:2)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      p2 <- ifelse(as.character(newdata$yprev) == "1", 0.25, 0.75)
      cbind("1" = 1 - p2, "2" = p2)
    },
    {
      trace_first <- markov.misc:::markov_transition_trace(
        model,
        data = first_order,
        times = 1:3,
        y_levels = 1:2
      )
      sops_first <- soprob_markov(
        model,
        newdata = first_order,
        times = 1:3,
        y_levels = 1:2
      )
      trace_second <- markov.misc:::markov_transition_trace(
        model,
        data = second_order,
        times = 1:3,
        y_levels = 1:2,
        p2_var = "ypprev"
      )
      sops_second <- soprob_markov(
        model,
        newdata = second_order,
        times = 1:3,
        y_levels = 1:2,
        p2_var = "ypprev"
      )
    }
  )

  expect_equal(trace_first$sops, sops_first)
  expect_equal(
    apply(trace_first$transitions, c(1, 2, 4), sum),
    trace_first$sops
  )
  expect_equal(trace_second$sops, sops_second)
  expect_equal(
    apply(trace_second$transitions, c(1, 2, 4), sum),
    trace_second$sops
  )
})

test_that("plot_transitions computes second-order treatment differences", {
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
        y_levels = 1:2,
        p2_var = "ypprev",
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
        y_levels = 1:2,
        absorb = 2,
        seed = 1
      )
    }
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(sum(plot$data$n), 2)
})

test_that("plot_transitions summarizes blrm transition draws without sampling", {
  model <- structure(
    list(
      y_levels = 1:2,
      draws = matrix(0, nrow = 3, ncol = 1)
    ),
    class = "blrm"
  )
  newdata <- data.frame(
    time = 0,
    yprev = factor(1, levels = 1:2)
  )
  calls <- list()

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
    {
      plot <- plot_transitions(
        model,
        newdata = newdata,
        times = 1,
        y_levels = 1:2,
        seed = 1,
        n_draws = NULL
      )
    }
  )

  estimate <- plot$data$estimate[
    plot$data$previous_state == "1" & plot$data$state == "2"
  ]

  expect_equal(estimate, 2 / 3)
  expect_length(calls, 1)
  expect_false(calls[[1]]$include_re)
  expect_equal(calls[[1]]$id_var, "id")
  expect_equal(calls[[1]]$draw_indices, 1:3)
  expect_null(calls[[1]]$gamma_draws)
})

test_that("plot_transitions passes n_draws to blrm diagnostics", {
  model <- structure(
    list(
      y_levels = 1:2,
      draws = matrix(0, nrow = 3, ncol = 1)
    ),
    class = "blrm"
  )
  newdata <- data.frame(
    time = 0,
    yprev = factor(1, levels = 1:2)
  )
  selected <- NULL
  calls <- list()

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    select_posterior_draws = function(model, n_draws = 100L, seed = NULL) {
      selected <<- list(n_draws = n_draws, seed = seed)
      seq_len(n_draws)
    },
    predict_blrm_response_markov = function(
      object,
      newdata,
      include_re,
      id_var,
      draw_indices,
      gamma_draws
    ) {
      calls[[length(calls) + 1L]] <<- draw_indices
      out <- array(
        0,
        dim = c(length(draw_indices), nrow(newdata), 2L),
        dimnames = list(draw_indices, rownames(newdata), c("1", "2"))
      )
      out[,, 1] <- 1
      out
    },
    {
      plot <- plot_transitions(
        model,
        newdata = newdata,
        times = 1,
        y_levels = 1:2,
        seed = 123,
        n_draws = 1
      )
    }
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(selected$n_draws, 1)
  expect_equal(selected$seed, 123)
  expect_equal(calls, list(1L))
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
  expect_s3_class(plot$coordinates, "CoordCartesian")
  expect_equal(plot$coordinates$limits$y, c(0, 1))
})

test_that("correlation diagnostics preserve negative correlations", {
  data <- data.frame(
    id = rep(1:4, each = 2),
    time = rep(1:2, times = 4),
    y = c(1, 2, 2, 1, 1, 2, 2, 1)
  )

  correlation <- plot_correlation(data, show_values = FALSE)
  correlation_build <- ggplot2::ggplot_build(correlation)

  expect_equal(as.numeric(correlation$data$correlation), -1)
  expect_equal(correlation$scales$get_scales("fill")$limits, c(-1, 1))
  expect_equal(sum(correlation_build$data[[1]]$fill == "grey90"), 0L)

  variogram <- plot_variogram(data, smooth = FALSE)
  variogram_build <- ggplot2::ggplot_build(variogram)

  expect_equal(as.numeric(variogram$data$correlation), -1)
  expect_equal(variogram$coordinates$limits$y, c(0, 1))
  expect_equal(variogram_build$data[[1]]$y, -1)
})

test_that("plot_correlation uses first-order model-implied moments", {
  model <- structure(list(), class = "vglm")
  newdata <- data.frame(
    time = 0,
    yprev = factor(c(1, 2), levels = 1:2)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      yprev <- as.character(newdata$yprev)
      p2 <- ifelse(
        newdata$time == 1,
        ifelse(yprev == "2", 1, 0),
        ifelse(yprev == "1", 0.5, 1)
      )
      cbind("1" = 1 - p2, "2" = p2)
    },
    {
      plot <- plot_correlation(
        model,
        newdata = newdata,
        times = 1:2,
        y_levels = 1:2,
        show_values = FALSE
      )
      variogram <- plot_variogram(
        model,
        newdata = newdata,
        times = 1:2,
        y_levels = 1:2,
        smooth = FALSE
      )
    }
  )

  expected <- 1 / sqrt(3)
  expect_equal(as.numeric(plot$data$correlation), expected, tolerance = 1e-6)
  expect_equal(
    as.numeric(variogram$data$correlation),
    expected,
    tolerance = 1e-6
  )
})

test_that("plot_correlation uses second-order model-implied moments", {
  model <- structure(list(), class = "vglm")
  newdata <- data.frame(
    time = 0,
    ypprev = factor(1, levels = 1:2),
    yprev = factor(1, levels = 1:2)
  )

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    predict_vglm_response_markov = function(object, newdata) {
      p2 <- ifelse(newdata$time == 1, 0.5, as.character(newdata$yprev) == "2")
      cbind("1" = 1 - p2, "2" = p2)
    },
    {
      plot <- plot_correlation(
        model,
        newdata = newdata,
        times = 1:3,
        y_levels = 1:2,
        p2_var = "ypprev",
        show_values = FALSE
      )
    }
  )

  expect_equal(as.numeric(plot$data$correlation), rep(1, 3), tolerance = 1e-6)
})

test_that("plot_correlation summarizes blrm correlations by draw", {
  model <- structure(
    list(
      y_levels = 1:2,
      draws = matrix(0, nrow = 2, ncol = 1)
    ),
    class = "blrm"
  )
  newdata <- data.frame(
    time = 0,
    yprev = factor(1, levels = 1:2)
  )
  calls <- list()

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
      calls[[length(calls) + 1L]] <<- draw_indices
      out <- array(
        NA_real_,
        dim = c(length(draw_indices), nrow(newdata), 2L),
        dimnames = list(draw_indices, rownames(newdata), c("1", "2"))
      )
      yprev <- as.character(newdata$yprev)
      for (i in seq_along(draw_indices)) {
        draw <- draw_indices[i]
        if (all(newdata$time == 1)) {
          p2 <- rep(0.5, nrow(newdata))
        } else if (draw == 1L) {
          p2 <- yprev == "2"
        } else {
          p2 <- yprev == "1"
        }
        out[i, , 1] <- 1 - p2
        out[i, , 2] <- p2
      }
      out
    },
    {
      plot <- plot_correlation(
        model,
        newdata = newdata,
        times = 1:2,
        y_levels = 1:2,
        seed = 1,
        n_draws = NULL,
        show_values = FALSE
      )
    }
  )

  expect_equal(as.numeric(plot$data$correlation), 0, tolerance = 1e-6)
  expect_equal(calls, list(1:2, 1:2))
})

test_that("plot_correlation and plot_variogram pass n_draws to blrm diagnostics", {
  model <- structure(
    list(
      y_levels = 1:2,
      draws = matrix(0, nrow = 4, ncol = 1)
    ),
    class = "blrm"
  )
  newdata <- data.frame(
    time = 0,
    yprev = factor(1, levels = 1:2)
  )
  selected <- list()
  calls <- list()

  with_mocked_bindings(
    validate_markov_model = function(object) NULL,
    select_posterior_draws = function(model, n_draws = 100L, seed = NULL) {
      selected[[length(selected) + 1L]] <<- list(n_draws = n_draws, seed = seed)
      seq_len(n_draws)
    },
    predict_blrm_response_markov = function(
      object,
      newdata,
      include_re,
      id_var,
      draw_indices,
      gamma_draws
    ) {
      calls[[length(calls) + 1L]] <<- draw_indices
      out <- array(
        NA_real_,
        dim = c(length(draw_indices), nrow(newdata), 2L),
        dimnames = list(draw_indices, rownames(newdata), c("1", "2"))
      )
      p2 <- ifelse(newdata$time == 1, 0.5, as.character(newdata$yprev) == "2")
      for (i in seq_along(draw_indices)) {
        out[i, , 1] <- 1 - p2
        out[i, , 2] <- p2
      }
      out
    },
    {
      correlation <- plot_correlation(
        model,
        newdata = newdata,
        times = 1:2,
        y_levels = 1:2,
        seed = 321,
        n_draws = 1,
        show_values = FALSE
      )
      variogram <- plot_variogram(
        model,
        newdata = newdata,
        times = 1:2,
        y_levels = 1:2,
        seed = 654,
        n_draws = 1,
        smooth = FALSE
      )
    }
  )

  expect_s3_class(correlation, "ggplot")
  expect_s3_class(variogram, "ggplot")
  expect_equal(
    selected,
    list(
      list(n_draws = 1, seed = 321),
      list(n_draws = 1, seed = 654)
    )
  )
  expect_true(all(vapply(calls, identical, logical(1), 1L)))
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
        y_levels = 1:2
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
        y_levels = 1:3
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
