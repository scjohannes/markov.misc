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
