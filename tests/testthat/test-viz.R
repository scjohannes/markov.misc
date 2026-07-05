test_that("plot_sops creates bar plots with optional faceting", {
  data <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(1:3, 4),
    y = c(1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 3, 3),
    tx = rep(c(0, 1), each = 6)
  )

  plot <- plot_sops(data, geom = "bar")
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 1)
  expect_equal(unique(built$data[[1]]$PANEL), factor(c(1, 2)))
  expect_equal(
    plot$labels$title,
    "Empirical State Occupancy Probabilities Over Time"
  )
})

test_that("plot_sops creates line plots with and without linetype groups", {
  data <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(1:3, 4),
    y = c(1, 1, 2, 1, 2, 2, 2, 2, 3, 1, 3, 3),
    tx = rep(c(0, 1), each = 6),
    site = rep(c("a", "b"), times = 6)
  )

  default_plot <- plot_sops(data, facet_var = NULL)
  expect_s3_class(default_plot$layers[[1]]$geom, "GeomLine")
  expect_true(default_plot$scales$has_scale("colour"))
  expect_true(default_plot$scales$has_scale("fill"))

  plot <- plot_sops(
    data,
    facet_var = NULL,
    linetype_var = "tx",
    geom = "line"
  )
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 1)
  expect_equal(plot$labels$linetype, "tx")
  expect_true(all(built$data[[1]]$y >= 0 & built$data[[1]]$y <= 1))

  no_linetype <- plot_sops(data, facet_var = NULL, geom = "line")
  expect_null(no_linetype$labels$linetype)

  override <- suppressMessages(
    default_plot +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::scale_fill_brewer(palette = "Dark2")
  )
  expect_s3_class(override, "ggplot")
  expect_true(override$scales$has_scale("colour"))
  expect_true(override$scales$has_scale("fill"))
})

make_avg_sops_plot_data <- function(with_draws = TRUE) {
  data <- expand.grid(
    time = 1:2,
    state = factor(c("1", "2", "3"), levels = c("1", "2", "3")),
    tx = c(0, 1),
    KEEP.OUT.ATTRS = FALSE
  )
  data$estimate <- c("1" = 0.2, "2" = 0.3, "3" = 0.5)[
    as.character(data$state)
  ]
  data$conf.low <- pmax(0, data$estimate - 0.05)
  data$conf.high <- pmin(1, data$estimate + 0.05)
  class(data) <- c("markov_avg_sops", class(data))

  if (with_draws) {
    draws <- expand.grid(
      draw_id = 1:7,
      time = 1:2,
      state = factor(c("1", "2", "3"), levels = c("1", "2", "3")),
      tx = c(0, 1),
      KEEP.OUT.ATTRS = FALSE
    )
    draws$estimate <- c("1" = 0.18, "2" = 0.31, "3" = 0.51)[
      as.character(draws$state)
    ]
    attr(data, "simulation_draws") <- draws
  }

  data
}

test_that("plot_sops creates model-derived line plots with ribbons", {
  data <- make_avg_sops_plot_data()

  plot <- plot_sops(data, geom = "line", facet_var = "tx")
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 2)
  expect_equal(plot$labels$title, "State Occupancy Probabilities Over Time")
  expect_equal(nrow(built$data[[1]]), nrow(data))
})

test_that("plot_sops respects stored ylevels for model-derived state order", {
  state_levels <- as.character(1:12)
  data <- expand.grid(
    time = 1,
    state = sort(state_levels),
    KEEP.OUT.ATTRS = FALSE
  )
  data$estimate <- 1 / length(state_levels)
  class(data) <- c("markov_avg_sops", class(data))
  attr(data, "ylevels") <- state_levels

  plot <- plot_sops(data, geom = "bar", facet_var = NULL, n_draws = 0)
  built <- ggplot2::ggplot_build(plot)
  fill_scale <- built$plot$scales$get_scales("fill")

  expect_equal(fill_scale$get_limits(), state_levels)
})

test_that("plot_sops warns but plots avg_sops bars without stored draws", {
  data <- make_avg_sops_plot_data(with_draws = FALSE)

  expect_warning(
    plot <- plot_sops(data, geom = "bar", facet_var = "tx"),
    "No stored draws found",
    fixed = TRUE
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 1)
})

test_that("plot_sops overlays a deterministic subset of avg_sops draws", {
  data <- make_avg_sops_plot_data()

  plot <- plot_sops(data, geom = "bar", facet_var = "tx", n_draws = 3)
  built <- ggplot2::ggplot_build(plot)
  draw_layer_ids <- unname(vapply(
    plot$layers[-1],
    function(layer) {
      unique(layer$data$draw_id)
    },
    numeric(1)
  ))
  layer_max_y <- vapply(
    built$data,
    function(layer_data) {
      max(layer_data$y)
    },
    numeric(1)
  )
  draw_layer_nrows <- unname(vapply(
    plot$layers[-1],
    function(layer) nrow(layer$data),
    integer(1)
  ))

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 4)
  expect_equal(draw_layer_ids, c(1, 4, 7))
  expect_equal(draw_layer_nrows, rep(12L, 3))
  expect_equal(layer_max_y, rep(1, 4))
})

test_that("plot_sops validates uncertainty controls", {
  data <- make_avg_sops_plot_data()

  expect_error(
    plot_sops(data, geom = "bar", n_draws = Inf),
    "`n_draws`"
  )
  expect_error(
    plot_sops(data, geom = "bar", draw_alpha = NA),
    "`draw_alpha`"
  )
  expect_error(
    plot_sops(data, geom = "bar", show_uncertainty = NA),
    "`show_uncertainty`"
  )
})

test_that("plot_sops rejects ambiguous summary groups", {
  data <- make_avg_sops_plot_data(with_draws = FALSE)
  class(data) <- "data.frame"
  data <- rbind(
    transform(data, scenario = "a"),
    transform(data, scenario = "b")
  )

  line_msg <- tryCatch(
    plot_sops(data, geom = "line", facet_var = "tx"),
    error = conditionMessage
  )
  bar_msg <- tryCatch(
    plot_sops(data, geom = "bar", facet_var = "tx"),
    error = conditionMessage
  )

  expect_match(line_msg, "`facet_var` or `linetype_var`", fixed = TRUE)
  expect_match(bar_msg, "`facet_var`, aggregate", fixed = TRUE)
  expect_no_match(bar_msg, "`linetype_var`", fixed = TRUE)

  plot <- plot_sops(data, geom = "line", facet_var = c("tx", "scenario"))
  expect_s3_class(plot, "ggplot")
})

test_that("plot_sops supports summary data frames and grid facets", {
  data <- make_avg_sops_plot_data()
  attr(data, "simulation_draws") <- NULL
  class(data) <- "data.frame"
  data$source <- rep(c("a", "b"), length.out = nrow(data))

  plot <- plot_sops(data, geom = "line", facet_var = c("tx", "source"))

  expect_s3_class(plot, "ggplot")
  expect_equal(class(plot$facet)[1], "FacetGrid")
})

make_comparison_plot_data <- function(
  comparison = "difference",
  with_time = TRUE
) {
  if (isTRUE(with_time)) {
    data <- expand.grid(
      time = 1:2,
      state_set = c("1", "2"),
      contrast = c("1 - 0", "2 - 0"),
      KEEP.OUT.ATTRS = FALSE
    )
    data$metric <- "sop"
  } else {
    data <- data.frame(
      state_set = c("1", "2"),
      contrast = c("1 / 0", "1 / 0"),
      metric = "time_in_state"
    )
  }
  data$variable <- "tx"
  data$reference_level <- 0
  data$comparison_level <- ifelse(data$contrast %in% c("2 - 0", "2 / 0"), 2, 1)
  data$comparison <- comparison
  data$estimate <- seq_len(nrow(data)) / 10
  if (identical(comparison, "ratio")) {
    data$estimate <- data$estimate + 0.8
  }
  data$conf.low <- data$estimate - 0.05
  data$conf.high <- data$estimate + 0.05
  class(data) <- c("markov_avg_comparisons", class(data))
  data
}

test_that("plot_comparisons plots time-varying differences on contrast scale", {
  data <- make_comparison_plot_data()

  plot <- plot_comparisons(data)
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 3)
  expect_s3_class(plot$layers[[3]]$geom, "GeomLine")
  expect_equal(unique(built$data[[1]]$yintercept), 0)
  expect_equal(plot$labels$y, "Difference")
  expect_equal(plot$labels$colour, "State set")
  expect_equal(plot$labels$linetype, "Contrast")
})

test_that("plot_comparisons plots ratios with point intervals and ratio reference", {
  data <- make_comparison_plot_data(comparison = "ratio", with_time = FALSE)

  plot <- plot_comparisons(data)
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 3)
  expect_s3_class(plot$layers[[3]]$geom, "GeomPoint")
  expect_equal(unique(built$data[[1]]$yintercept), 1)
  expect_equal(plot$labels$x, "State set")
  expect_equal(plot$labels$y, "Ratio")
})

test_that("plot_comparisons accepts custom grouping labels", {
  data <- make_comparison_plot_data(comparison = "ratio", with_time = FALSE)
  data$initial_y <- c("1", "2")

  plot <- plot_comparisons(data, color_var = "initial_y")

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$colour, "initial_y")
  expect_equal(plot$labels$fill, "initial_y")
})

test_that("plot_comparisons orders numeric state sets naturally", {
  data <- data.frame(
    state_set = c("1", "10", "11", "2", "3"),
    contrast = "1 / 0",
    metric = "time_in_state",
    variable = "tx",
    reference_level = 0,
    comparison_level = 1,
    comparison = "ratio",
    estimate = c(1, 1.1, 1.2, 0.9, 0.8)
  )
  data$conf.low <- data$estimate - 0.05
  data$conf.high <- data$estimate + 0.05
  class(data) <- c("markov_avg_comparisons", class(data))
  attr(data, "ylevels") <- as.character(1:11)

  plot <- plot_comparisons(data)
  built <- ggplot2::ggplot_build(plot)

  expect_equal(
    built$layout$panel_params[[1]]$x$get_labels(),
    c("1", "2", "3", "10", "11")
  )
})

test_that("plot_comparisons validates comparison data", {
  expect_error(
    plot_comparisons(data.frame(estimate = 1)),
    "missing required comparison column",
    fixed = TRUE
  )

  data <- make_comparison_plot_data()
  data$comparison[1] <- "ratio"

  expect_error(
    plot_comparisons(data),
    "requires one comparison scale",
    fixed = TRUE
  )
})

test_that("plot_results validates grouping and x variables", {
  data <- data.frame(
    analysis = c("a", "b"),
    sample_size = c(50, 100),
    power = c(0.7, 0.8)
  )

  expect_error(
    plot_results(data, power, x = missing_x, group = analysis),
    "Variable 'analysis' must be a factor",
    fixed = TRUE
  )

  data$analysis <- factor(data$analysis)
  expect_error(
    plot_results(data, power, x = missing_x, group = analysis),
    "Variable 'missing_x' not found in data",
    fixed = TRUE
  )
})

test_that("plot_results returns individual and combined plots", {
  data <- data.frame(
    analysis = factor(rep(c("a", "b"), each = 2)),
    sample_size = rep(c(50, 100), 2),
    power = c(0.6, 0.7, 0.5, 0.8),
    coverage = c(0.9, 0.91, 0.92, 0.93)
  )

  result <- plot_results(
    data,
    power,
    coverage,
    x = sample_size,
    group = analysis
  )

  expect_named(result, c("single_plots", "grid"))
  expect_length(result$single_plots, 2)
  expect_s3_class(result$single_plots[[1]], "ggplot")
  expect_s3_class(result$grid, "patchwork")

  separate <- plot_results(
    data,
    power,
    x = sample_size,
    group = analysis,
    combine = FALSE
  )
  expect_null(separate$grid)
})

test_that("plot_results uses summary geoms for factor x variables", {
  data <- data.frame(
    analysis = factor(rep(c("a", "b"), each = 2)),
    scenario = factor(rep(c("low", "high"), 2)),
    power = c(0.6, 0.7, 0.5, 0.8)
  )

  result <- plot_results(data, power, x = scenario, group = analysis)

  expect_length(result$single_plots[[1]]$layers, 2)
  expect_equal(result$single_plots[[1]]$labels$y, "power")
})
