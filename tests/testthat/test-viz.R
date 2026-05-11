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
  expect_equal(plot$labels$title, "Empirical State Occupancy Probabilities Over Time")
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
  draw_layer_ids <- unname(vapply(plot$layers[-1], function(layer) {
    unique(layer$data$draw_id)
  }, numeric(1)))
  layer_max_y <- vapply(built$data, function(layer_data) {
    max(layer_data$y)
  }, numeric(1))
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

  result <- plot_results(data, power, coverage, x = sample_size, group = analysis)

  expect_named(result, c("single_plots", "grid"))
  expect_length(result$single_plots, 2)
  expect_s3_class(result$single_plots[[1]], "ggplot")
  expect_s3_class(result$grid, "patchwork")

  separate <- plot_results(data, power, x = sample_size, group = analysis, combine = FALSE)
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

test_that("plot_bootstrap_sops creates confidence bands with grouping and faceting", {
  bootstrap_data <- expand.grid(
    boot_id = 1:4,
    time = 1:3,
    tx = c(0, 1),
    site = c("a", "b")
  )
  bootstrap_data$state_1 <- seq(0.2, 0.8, length.out = nrow(bootstrap_data))
  bootstrap_data$state_2 <- 1 - bootstrap_data$state_1

  plot <- plot_bootstrap_sops(
    bootstrap_data,
    facet_var = "site",
    conf_level = 0.90,
    title = "Bootstrap SOPs"
  )
  built <- ggplot2::ggplot_build(plot)

  expect_s3_class(plot, "ggplot")
  expect_equal(length(plot$layers), 2)
  expect_equal(plot$labels$title, "Bootstrap SOPs")
  expect_equal(plot$labels$linetype, "tx")
  expect_true(all(built$data[[1]]$y >= 0 & built$data[[1]]$y <= 1))
})

test_that("plot_bootstrap_sops supports ungrouped default titles", {
  bootstrap_data <- data.frame(
    boot_id = rep(1:4, each = 2),
    time = rep(1:2, 4),
    state_1 = c(0.2, 0.3, 0.25, 0.35, 0.22, 0.32, 0.28, 0.38),
    state_2 = c(0.8, 0.7, 0.75, 0.65, 0.78, 0.68, 0.72, 0.62)
  )

  plot <- plot_bootstrap_sops(bootstrap_data, group_var = NULL)

  expect_s3_class(plot, "ggplot")
  expect_null(plot$labels$linetype)
  expect_match(plot$labels$title, "95% Confidence Bands", fixed = TRUE)
})
