expect_trajectory_contract <- function(
  object,
  expected_cols,
  n_patients = NULL,
  follow_up_time = NULL,
  states = NULL,
  check_lag = TRUE
) {
  testthat::expect_s3_class(object, "data.frame")
  testthat::expect_named(object, expected_cols, ignore.order = FALSE)
  testthat::expect_false(
    anyNA(object$id),
    info = "`id` should not contain missing values"
  )
  testthat::expect_false(
    anyNA(object$time),
    info = "`time` should not contain missing values"
  )

  if (!is.null(n_patients) && !is.null(follow_up_time)) {
    testthat::expect_equal(nrow(object), n_patients * follow_up_time)
    testthat::expect_setequal(unique(object$id), seq_len(n_patients))
    testthat::expect_equal(sort(unique(object$time)), seq_len(follow_up_time))
    testthat::expect_equal(
      as.integer(table(object$id)),
      rep(follow_up_time, n_patients),
      info = "Each patient should contribute one row per follow-up time"
    )
  }

  if (!is.null(states)) {
    testthat::expect_true(
      all(object$y[!is.na(object$y)] %in% states),
      info = "`y` should stay inside the declared state support"
    )
    testthat::expect_true(
      all(object$yprev[!is.na(object$yprev)] %in% states),
      info = "`yprev` should stay inside the declared state support"
    )
  }

  if (check_lag) {
    by_id <- split(
      object[order(object$id, object$time), , drop = FALSE],
      object$id
    )
    lag_ok <- vapply(
      by_id,
      function(patient) {
        if (nrow(patient) <= 1) {
          return(TRUE)
        }
        identical(
          unname(patient$yprev[-1]),
          unname(utils::head(patient$y, -1))
        )
      },
      logical(1)
    )

    testthat::expect_equal(
      names(lag_ok)[!lag_ok],
      character(),
      info = "`yprev` should be the previous observed `y` within each patient"
    )
  }

  invisible(object)
}

expect_absorbing_state_sticky <- function(
  object,
  absorbing_state,
  id_col = "id",
  time_col = "time",
  state_col = "y"
) {
  object <- object[order(object[[id_col]], object[[time_col]]), , drop = FALSE]
  state_values <- as.character(object[[state_col]])
  target <- as.character(absorbing_state)

  testthat::expect_true(
    any(state_values == target, na.rm = TRUE),
    info = "The fixture should exercise the absorbing state"
  )

  by_id <- split(object, object[[id_col]])
  sticky <- vapply(
    by_id,
    function(patient) {
      patient_states <- as.character(patient[[state_col]])
      first_absorbed <- match(target, patient_states)
      if (is.na(first_absorbed)) {
        return(TRUE)
      }
      all(
        patient_states[seq.int(first_absorbed, length(patient_states))] ==
          target
      )
    },
    logical(1)
  )

  testthat::expect_equal(
    names(sticky)[!sticky],
    character(),
    info = "Once entered, the absorbing state should be carried forward"
  )

  invisible(object)
}

expect_probability_array <- function(
  object,
  expected_dim = NULL,
  tolerance = 1e-10,
  state_margin = length(dim(object)),
  check_normalized = TRUE
) {
  testthat::expect_true(
    is.array(object),
    info = "Probability output should be an array"
  )
  testthat::expect_true(
    is.numeric(object),
    info = "Probability output should be numeric"
  )

  if (!is.null(expected_dim)) {
    testthat::expect_equal(dim(object), expected_dim)
  }

  testthat::expect_false(
    anyNA(object),
    info = "Probabilities should not be missing"
  )
  testthat::expect_true(
    all(is.finite(object)),
    info = "Probabilities should be finite"
  )
  testthat::expect_true(
    all(object >= -tolerance & object <= 1 + tolerance),
    info = "Probabilities should stay in [0, 1]"
  )

  if (check_normalized) {
    keep_margins <- setdiff(seq_along(dim(object)), state_margin)
    probability_sums <- apply(object, keep_margins, sum)
    testthat::expect_equal(
      as.vector(probability_sums),
      rep(1, length(probability_sums)),
      tolerance = tolerance,
      info = "Probabilities should sum to one over states"
    )
  }

  invisible(object)
}

expect_absorbing_probability_monotone <- function(
  object,
  absorbing_state,
  tolerance = 1e-10,
  time_margin = 2L,
  state_margin = length(dim(object))
) {
  dims <- dim(object)
  testthat::expect_true(
    length(dims) %in% c(3L, 4L),
    info = "Absorbing probability checks support 3D and 4D SOP arrays"
  )

  state_names <- dimnames(object)[[state_margin]]
  state_index <- if (!is.null(state_names)) {
    match(as.character(absorbing_state), as.character(state_names))
  } else {
    as.integer(absorbing_state)
  }

  testthat::expect_false(
    is.na(state_index),
    info = "Absorbing state should be in the state dimension"
  )

  if (length(dims) == 3L) {
    absorbing_probs <- matrix(
      object[,, state_index],
      nrow = dims[1],
      ncol = dims[time_margin]
    )
    increments <- apply(absorbing_probs, 1, diff)
  } else {
    absorbing_probs <- array(
      object[,,, state_index],
      dim = dims[-state_margin]
    )
    increments <- apply(absorbing_probs, c(1, 2), diff)
  }

  testthat::expect_true(
    all(increments >= -tolerance),
    info = "Absorbing-state probability should not decrease over time"
  )

  invisible(object)
}

expect_inference_intervals <- function(
  object,
  require_std_error = FALSE,
  require_positive_std_error = FALSE
) {
  testthat::expect_contains(names(object), c("conf.low", "conf.high"))
  testthat::expect_false(
    anyNA(object$conf.low),
    info = "`conf.low` should be complete"
  )
  testthat::expect_false(
    anyNA(object$conf.high),
    info = "`conf.high` should be complete"
  )
  testthat::expect_true(
    all(object$conf.low <= object$conf.high),
    info = "Interval lower bounds should not exceed upper bounds"
  )

  if (require_std_error || require_positive_std_error) {
    testthat::expect_contains(names(object), "std.error")
    testthat::expect_false(
      anyNA(object$std.error),
      info = "`std.error` should be complete"
    )
    testthat::expect_true(
      all(object$std.error >= 0),
      info = "`std.error` should be non-negative"
    )
  }

  if (require_positive_std_error) {
    testthat::expect_true(
      any(object$std.error > 0),
      info = "At least one standard error should be positive"
    )
  }

  invisible(object)
}
