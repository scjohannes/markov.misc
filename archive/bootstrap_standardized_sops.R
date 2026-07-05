# Archived on 2026-07-05.
#
# This legacy helper produced wide state_* bootstrap SOP tables for
# plot_bootstrap_sops() and time_in_state() examples. The active package now uses
# avg_sops() |> inferences(method = "bootstrap") for bootstrap uncertainty.

bootstrap_standardized_sops <- function(
  model,
  data,
  n_boot,
  workers = NULL,
  parallel = NULL,
  ylevels = factor(1:6),
  absorb = 6,
  include_coefs = TRUE,
  times = NULL,
  update_datadist = TRUE,
  use_coefstart = FALSE,
  varnames = list(tvarname = "time", pvarname = "yprev", id = "id", tx = "tx"),
  t_covs = NULL
) {
  if (!is.null(parallel)) {
    warning(
      "The 'parallel' argument is deprecated. ",
      "Please use 'workers' instead.\n",
      "  - For sequential processing: workers = NULL or workers = 1\n",
      "  - For parallel processing: workers = N (e.g., workers = 8)"
    )
    if (parallel && is.null(workers)) {
      workers <- parallel::detectCores() - 1
    }
  }

  if (!inherits(model, "orm") && !inherits(model, "vglm")) {
    stop(
      "model must be an orm object (rms package) or vglm object (VGAM package)."
    )
  }

  if (is.null(times)) {
    times <- 1:max(data[[varnames$tvarname]])
  }

  n_times <- length(times)
  n_states <- length(ylevels)

  boot_ids <- fast_group_bootstrap(
    data = data,
    id_var = varnames$id,
    n_boot = n_boot
  )

  analysis_fn <- function(boot_data) {
    boot_result <- bootstrap_analysis_wrapper(
      boot_data = boot_data,
      model = model,
      factor_cols = c("y", varnames$pvarname),
      original_data = data,
      ylevels = ylevels,
      absorb = absorb,
      update_datadist = update_datadist,
      use_coefstart = use_coefstart
    )

    m_boot <- boot_result$model
    boot_data <- boot_result$data
    boot_ylevels <- boot_result$ylevels
    boot_absorb <- boot_result$absorb
    missing_states <- boot_result$missing_states

    if (is.null(m_boot)) {
      return(list(
        sop_result = NULL,
        coefs = NULL
      ))
    }

    ylevel_names <- as_state_labels(ylevels)
    states_present <- ylevel_names[
      !ylevel_names %in% as_state_labels(missing_states)
    ]

    sop_result <- tryCatch(
      standardize_sops(
        model = m_boot,
        data = boot_data,
        times = times,
        ylevels = factor(boot_ylevels),
        absorb = boot_absorb,
        varnames = list(
          tvarname = varnames$tvarname,
          pvarname = varnames$pvarname,
          id = "new_id",
          tx = varnames$tx
        ),
        t_covs = t_covs
      ),
      error = function(e) {
        warning("standardize_sops failed: ", e$message)
        NULL
      }
    )

    if (include_coefs == TRUE) {
      coefs <- as.list(coef(m_boot))
    } else {
      coefs <- NULL
    }

    if (!is.null(sop_result) && length(missing_states) > 0) {
      sop_tx_full <- matrix(0, nrow = n_times, ncol = n_states)
      sop_ctrl_full <- matrix(0, nrow = n_times, ncol = n_states)
      colnames(sop_tx_full) <- as.character(ylevels)
      colnames(sop_ctrl_full) <- as.character(ylevels)

      for (i in seq_along(states_present)) {
        original_state <- states_present[i]
        sop_tx_full[, as.character(original_state)] <- sop_result$sop_tx[, i]
        sop_ctrl_full[, as.character(original_state)] <- sop_result$sop_ctrl[, i]
      }

      sop_result <- list(
        sop_tx = sop_tx_full,
        sop_ctrl = sop_ctrl_full
      )
    }

    list(sop_result = sop_result, coefs = coefs)
  }

  boot_results <- apply_to_bootstrap(
    boot_samples = boot_ids,
    analysis_fn = analysis_fn,
    data = data,
    id_var = varnames$id,
    workers = workers,
    packages = c("rms", "VGAM", "Hmisc", "stats"),
    globals = c(
      "model",
      "times",
      "ylevels",
      "absorb",
      "varnames",
      "n_times",
      "n_states",
      "t_covs"
    )
  )

  sop_result_list <- list()

  for (i in seq_along(boot_results)) {
    sop_i <- boot_results[[i]][["sop_result"]]

    if (!is.null(sop_i)) {
      tx_df <- as.data.frame(sop_i$sop_tx, optional = TRUE)
      colnames(tx_df) <- paste0("state_", ylevels)
      tx_df$time <- times
      tx_df$tx <- 1
      tx_df$boot_id <- i

      ctrl_df <- as.data.frame(sop_i$sop_ctrl, optional = TRUE)
      colnames(ctrl_df) <- paste0("state_", ylevels)
      ctrl_df$time <- times
      ctrl_df$tx <- 0
      ctrl_df$boot_id <- i

      sop_result_list[[i]] <- bind_rows_fill(list(tx_df, ctrl_df))
    }
  }

  if (length(sop_result_list) > 0) {
    sop_results <- bind_rows_fill(sop_result_list)
    sop_results <- reorder_columns(sop_results, c("boot_id", "time", "tx"))
  } else {
    state_cols <- stats::setNames(
      replicate(n_states, numeric(0), simplify = FALSE),
      paste0("state_", ylevels)
    )
    sop_results <- data.frame(
      boot_id = integer(0),
      time = times[0],
      tx = integer(0),
      state_cols,
      check.names = FALSE
    )
  }

  coefs <- lapply(boot_results, function(x) x[["coefs"]])
  coefs_results <- named_list_to_wide(coefs, id = seq_len(n_boot))

  list(sops = sop_results, coefs = coefs_results)
}
