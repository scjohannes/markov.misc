# Helper functions for simulating data and fitting models in tests

#' Generate standard test data using Brownian motion simulation
#'
#' @param n_patients Number of patients
#' @param follow_up_time Follow-up time
#' @param seed Random seed
#' @param treatment_effect Treatment effect (mu_treatment_effect)
#'
#' @return A prepared markov data frame
make_test_data <- function(
  n_patients = 50,
  follow_up_time = 20,
  seed = 123,
  treatment_effect = 0
) {
  raw_data <- sim_trajectories_brownian(
    n_patients = n_patients,
    follow_up_time = follow_up_time,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = seed,
    mu_treatment_effect = treatment_effect
  )

  data <- prepare_markov_data(raw_data)

  # Generate manual spline basis
  # We use rcs() to get the basis, then extract columns
  time_spl_m <- rms::rcs(data$time, 3)
  knots <- attr(time_spl_m, "parms")
  data$time_lin <- as.vector(time_spl_m[, 1])
  data$time_nlin_1 <- as.vector(time_spl_m[, 2])

  return(data)
}

make_time_covariates <- function(data, time_col = "time", ...) {
  cols <- c(time_col, ...)
  out <- unique(data[cols])
  out <- out[order(out[[time_col]]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Fit a standard VGLM model for testing
#'
#' @param data Data frame from make_test_data
#' @param robust Logical, whether to return a robust covariance model
#' @param cluster Cluster variable for robust covariance
#'
#' @return A fitted vglm or robcov_vglm object
make_test_model <- function(data, robust = FALSE, cluster = NULL) {
  # fit model
  m <- VGAM::vglm(
    ordered(y) ~ time_lin + time_nlin_1 + tx + yprev,
    family = VGAM::cumulative(reverse = TRUE, parallel = TRUE),
    data = data
  )

  if (robust) {
    if (is.null(cluster)) {
      cluster <- data$id
    }
    return(robcov_vglm(m, cluster = cluster))
  }

  m
}

make_score_bootstrap_case <- function(
  seed,
  n_patients = 60,
  follow_up_time = 10
) {
  data <- make_test_data(
    n_patients = n_patients,
    seed = seed,
    follow_up_time = follow_up_time
  )

  list(
    data = data,
    baseline = data[data$time == 1, , drop = FALSE],
    model = make_test_model(data, robust = TRUE),
    times = seq_len(follow_up_time),
    ylevels = 1:6,
    absorb = 6
  )
}
