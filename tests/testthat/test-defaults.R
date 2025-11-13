# Unit Tests: Default LP Function and Dataset
#
# These tests verify that:
# 1. violet_baseline dataset exists and has correct structure
# 2. lp_violet() function works correctly
# 3. sim_trajectories_markov() works with defaults
# 4. Minimal argument calls work as expected

test_that("violet_baseline dataset exists and has correct structure", {
  # Check dataset exists
  expect_true(exists("violet_baseline", envir = asNamespace("markov.misc")))

  # Load the dataset
  data("violet_baseline", package = "markov.misc")

  # Check dimensions
  expect_equal(nrow(violet_baseline), 250000)
  expect_equal(ncol(violet_baseline), 5)

  # Check columns exist
  expected_cols <- c("id", "yprev", "age", "sofa", "tx")
  expect_true(all(expected_cols %in% names(violet_baseline)))
  expect_setequal(names(violet_baseline), expected_cols)

  # Check data types
  expect_type(violet_baseline$id, "integer")
  expect_type(violet_baseline$yprev, "double")
  expect_type(violet_baseline$age, "integer")
  expect_type(violet_baseline$sofa, "integer")
  expect_type(violet_baseline$tx, "double")

  # Check no missing values in required columns
  expect_false(any(is.na(violet_baseline$id)))
  expect_false(any(is.na(violet_baseline$yprev)))
  expect_false(any(is.na(violet_baseline$age)))
  expect_false(any(is.na(violet_baseline$sofa)))
  expect_false(any(is.na(violet_baseline$tx)))

  # Check yprev values are valid states (1-6)
  expect_true(all(violet_baseline$yprev %in% 1:6))

  # Check tx is binary (0 or 1)
  expect_true(all(violet_baseline$tx %in% c(0, 1)))

  # Check age is reasonable (e.g., 18-120)
  expect_true(all(violet_baseline$age >= 18 & violet_baseline$age <= 120))

  # Check sofa is reasonable (0-24)
  expect_true(all(violet_baseline$sofa >= 0 & violet_baseline$sofa <= 24))
})

test_that("lp_violet function works correctly", {
  # Test with sample data
  test_lp <- lp_violet(
    yprev = c(2, 3, 5),
    t = 10,
    age = c(65, 70, 55),
    sofa = c(5, 8, 6),
    tx = c(0, 1, 1),
    parameter = log(0.8),
    extra_params = c(
      "time" = -0.738194,
      "time'" = 0.7464006,
      "age" = 0.010321,
      "sofa" = 0.046901,
      "yprev=1" = -8.518344,
      "yprev=3" = 0,
      "yprev=4" = 1.315332,
      "yprev=5" = 6.576662,
      "yprev=1 * time" = 0,
      "yprev=3 * time" = 0,
      "yprev=4 * time" = 0,
      "yprev=5 * time" = 0
    )
  )

  # Check output length
  expect_length(test_lp, 3)

  # Check no NA values
  expect_false(any(is.na(test_lp)))

  # Check output is numeric
  expect_type(test_lp, "double")

  # Check values are finite
  expect_true(all(is.finite(test_lp)))
})

test_that("lp_violet handles edge cases", {
  # Test with single patient
  single_lp <- lp_violet(
    yprev = 2,
    t = 1,
    age = 65,
    sofa = 5,
    tx = 0,
    parameter = 0,
    extra_params = c(
      "time" = -0.738194,
      "time'" = 0.7464006,
      "age" = 0.010321,
      "sofa" = 0.046901,
      "yprev=1" = -8.518344,
      "yprev=3" = 0,
      "yprev=4" = 1.315332,
      "yprev=5" = 6.576662,
      "yprev=1 * time" = 0,
      "yprev=3 * time" = 0,
      "yprev=4 * time" = 0,
      "yprev=5 * time" = 0
    )
  )

  expect_length(single_lp, 1)
  expect_false(is.na(single_lp))

  # Test with different previous states
  for (state in c(1, 2, 3, 4, 5)) {
    lp_state <- lp_violet(
      yprev = state,
      t = 5,
      age = 65,
      sofa = 5,
      tx = 0,
      parameter = 0,
      extra_params = c(
        "time" = -0.738194,
        "time'" = 0.7464006,
        "age" = 0.010321,
        "sofa" = 0.046901,
        "yprev=1" = -8.518344,
        "yprev=3" = 0,
        "yprev=4" = 1.315332,
        "yprev=5" = 6.576662,
        "yprev=1 * time" = 0,
        "yprev=3 * time" = 0,
        "yprev=4 * time" = 0,
        "yprev=5 * time" = 0
      )
    )
    expect_false(is.na(lp_state))
  }
})

test_that("sim_trajectories_markov works with all defaults", {
  # Load dataset
  data("violet_baseline", package = "markov.misc")

  # Use a small subset for testing
  test_baseline <- violet_baseline[1:100, ]

  # Run simulation with explicit seed
  set.seed(12345)
  trajectories <- sim_trajectories_markov(
    baseline_data = test_baseline,
    follow_up_time = 30,
    parameter = log(0.8), # OR = 0.8
    seed = 12345
  )

  # Check output structure
  expected_rows <- 100 * 30 # 100 patients × 30 days
  expect_equal(nrow(trajectories), expected_rows)

  # Check required columns exist
  required_cols <- c("id", "time", "y", "yprev", "age", "sofa", "tx")
  expect_true(all(required_cols %in% names(trajectories)))

  # Check data types
  expect_type(trajectories$id, "integer")
  expect_type(trajectories$time, "integer")
  expect_type(trajectories$y, "double")
  expect_type(trajectories$yprev, "double")

  # Check time sequences are correct
  expect_equal(unique(trajectories$time), 1:30)

  # Check each patient has 30 observations
  patient_counts <- table(trajectories$id)
  expect_true(all(patient_counts == 30))

  # Check states are valid (1-6)
  expect_true(all(trajectories$y %in% 1:6))
  expect_true(all(trajectories$yprev %in% 1:6))

  # Check death state (6) is absorbing
  deaths <- trajectories[trajectories$y == 6, ]
  if (nrow(deaths) > 0) {
    # For each death, check all subsequent states are also 6
    for (patient_id in unique(deaths$id)) {
      patient_data <- trajectories[trajectories$id == patient_id, ]
      first_death <- min(patient_data$time[patient_data$y == 6])
      subsequent_states <- patient_data$y[patient_data$time >= first_death]
      expect_true(all(subsequent_states == 6))
    }
  }
})

test_that("sim_trajectories_markov works with minimal arguments", {
  # Load dataset
  data("violet_baseline", package = "markov.misc")

  # Run simulation with minimal arguments (uses default follow_up_time = 60)
  set.seed(99999)
  trajectories_minimal <- sim_trajectories_markov(
    baseline_data = violet_baseline[1:50, ],
    parameter = log(0.75),
    seed = 99999
  )

  # Should use default follow_up_time of 60
  expected_rows <- 50 * 60
  expect_equal(nrow(trajectories_minimal), expected_rows)

  # Check output has proper structure
  expect_s3_class(trajectories_minimal, "tbl_df")

  # Check required columns
  required_cols <- c("id", "time", "y", "yprev", "age", "sofa", "tx")
  expect_true(all(required_cols %in% names(trajectories_minimal)))
})

test_that("sim_trajectories_markov produces reproducible results", {
  # Load dataset
  data("violet_baseline", package = "markov.misc")

  test_baseline <- violet_baseline[1:20, ]

  # Run simulation twice with same seed
  traj1 <- sim_trajectories_markov(
    baseline_data = test_baseline,
    follow_up_time = 10,
    parameter = log(0.8),
    seed = 42
  )

  traj2 <- sim_trajectories_markov(
    baseline_data = test_baseline,
    follow_up_time = 10,
    parameter = log(0.8),
    seed = 42
  )

  # Results should be identical
  expect_identical(traj1, traj2)

  # Run with different seed
  traj3 <- sim_trajectories_markov(
    baseline_data = test_baseline,
    follow_up_time = 10,
    parameter = log(0.8),
    seed = 999
  )

  # Results should be different
  expect_false(identical(traj1$y, traj3$y))
})

test_that("sim_trajectories_markov handles different treatment effects", {
  # Load dataset
  data("violet_baseline", package = "markov.misc")

  test_baseline <- violet_baseline[1:50, ]

  # Test with null effect (parameter = 0)
  traj_null <- sim_trajectories_markov(
    baseline_data = test_baseline,
    follow_up_time = 20,
    parameter = 0,
    seed = 123
  )
  expect_equal(nrow(traj_null), 50 * 20)

  # Test with negative effect (worse outcomes)
  traj_neg <- sim_trajectories_markov(
    baseline_data = test_baseline,
    follow_up_time = 20,
    parameter = log(1.2),
    seed = 123
  )
  expect_equal(nrow(traj_neg), 50 * 20)

  # Test with positive effect (better outcomes)
  traj_pos <- sim_trajectories_markov(
    baseline_data = test_baseline,
    follow_up_time = 20,
    parameter = log(0.7),
    seed = 123
  )
  expect_equal(nrow(traj_pos), 50 * 20)
})
