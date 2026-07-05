# Archived on 2026-07-05.
#
# This legacy plot helper consumed the wide state_* output produced by
# bootstrap_standardized_sops(). The active package uses plot_sops() for tidy SOP
# outputs from sops()/avg_sops()/inferences().

plot_bootstrap_sops <- function(
  bootstrap_data,
  time_var = "time",
  group_var = "tx",
  facet_var = NULL,
  conf_level = 0.95,
  title = NULL
) {
  conf_level <- validate_conf_level(conf_level)
  alpha <- 1 - conf_level
  lower_q <- alpha / 2
  upper_q <- 1 - alpha / 2

  sops_long <- pivot_state_columns_long(bootstrap_data)

  group_vars <- c(time_var, "state")
  if (!is.null(group_var)) {
    group_vars <- c(group_var, group_vars)
  }
  if (!is.null(facet_var)) {
    group_vars <- c(facet_var, group_vars)
  }
  group_vars <- unique(group_vars)

  group_key <- do.call(
    interaction,
    c(sops_long[, group_vars, drop = FALSE], drop = TRUE, sep = "\r")
  )
  ci_bands <- bind_rows_fill(lapply(
    split(seq_len(nrow(sops_long)), group_key),
    function(idx) {
      group_data <- sops_long[idx, , drop = FALSE]
      out <- group_data[1, group_vars, drop = FALSE]
      out$lower <- as.numeric(stats::quantile(group_data$probability, lower_q))
      out$median <- stats::median(group_data$probability)
      out$upper <- as.numeric(stats::quantile(group_data$probability, upper_q))
      out
    }
  ))

  aes_mapping <- ggplot2::aes(
    x = .data[[time_var]],
    y = .data[["median"]],
    color = .data[["state"]]
  )

  if (!is.null(group_var)) {
    aes_mapping$linetype <- ggplot2::aes(
      linetype = factor(.data[[group_var]])
    )$linetype
    aes_mapping$group <- ggplot2::aes(
      group = interaction(.data[["state"]], factor(.data[[group_var]]))
    )$group
  } else {
    aes_mapping$group <- ggplot2::aes(group = .data[["state"]])$group
  }

  p <- ggplot2::ggplot(ci_bands, aes_mapping) +
    ggplot2::geom_line()

  if (!is.null(group_var)) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["lower"]],
          ymax = .data[["upper"]],
          fill = .data[["state"]],
          group = interaction(.data[["state"]], factor(.data[[group_var]]))
        ),
        alpha = 0.2,
        color = NA
      )
  } else {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data[["lower"]],
          ymax = .data[["upper"]],
          fill = .data[["state"]]
        ),
        alpha = 0.2,
        color = NA
      )
  }

  if (is.null(title)) {
    ci_pct <- round(conf_level * 100)
    title <- paste0(
      "Standardized State Occupation Probabilities with Bootstrapped ",
      ci_pct,
      "% Confidence Bands"
    )
  }

  p <- p +
    ggplot2::labs(
      x = "Time",
      y = "State Occupation Probability",
      color = "State",
      fill = "State",
      title = title
    )

  if (!is.null(group_var)) {
    p <- p + ggplot2::labs(linetype = group_var)
  }

  if (!is.null(facet_var)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet_var]]))
  }

  p
}
