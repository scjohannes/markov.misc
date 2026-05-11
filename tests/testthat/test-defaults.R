# Unit Tests: Default LP Function and Dataset
#
# These tests verify that:
# 1. violet_baseline dataset exists and has correct structure
# 2. lp_violet() function works correctly
# 3. sim_trajectories_markov() works with defaults
# 4. Minimal argument calls work as expected

describe("violet_baseline", {
  it("exists and has correct structure", {
    # Load the dataset
    data("violet_baseline", package = "markov.misc")
    expect_true(exists("violet_baseline"))

    # Check dimensions
    expect_equal(nrow(violet_baseline), 250000)
    expect_equal(ncol(violet_baseline), 5)

    # Check columns exist
    expected_cols <- c("id", "yprev", "age", "sofa", "tx")
    expect_named(violet_baseline, expected_cols, ignore.order = TRUE)

    # Check data types
    expect_type(violet_baseline$id, "integer")
    expect_type(violet_baseline$yprev, "double")
    expect_type(violet_baseline$age, "integer")
    expect_type(violet_baseline$sofa, "integer")
    expect_type(violet_baseline$tx, "double")

    # Check no missing values in required columns
    expect_all_false(is.na(violet_baseline$id))
    expect_all_false(is.na(violet_baseline$yprev))
    expect_all_false(is.na(violet_baseline$age))
    expect_all_false(is.na(violet_baseline$sofa))
    expect_all_false(is.na(violet_baseline$tx))

    # Check yprev values are valid states (1-6)
    expect_all_true(violet_baseline$yprev %in% 1:6)

    # Check tx is binary (0 or 1)
    expect_all_true(violet_baseline$tx %in% c(0, 1))

    # Check age is reasonable (e.g., 18-120)
    expect_gte(min(violet_baseline$age), 18)
    expect_lte(max(violet_baseline$age), 120)

    # Check sofa is reasonable (0-24)
    expect_gte(min(violet_baseline$sofa), 0)
    expect_lte(max(violet_baseline$sofa), 24)
  })
})

describe("lp_violet()", {
  extra_params <- c(
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

  it("works correctly with sample data", {
    test_lp <- lp_violet(
      yprev = c(2, 3, 5),
      t = 10,
      age = c(65, 70, 55),
      sofa = c(5, 8, 6),
      tx = c(0, 1, 1),
      parameter = log(0.8),
      extra_params = extra_params
    )

    expect_length(test_lp, 3)
    expect_all_false(is.na(test_lp))
    expect_type(test_lp, "double")
    expect_all_true(is.finite(test_lp))

    # Snapshot the values to catch regressions
    expect_snapshot(test_lp)
  })

  it("handles edge cases", {
    # Test with single patient
    single_lp <- lp_violet(
      yprev = 2,
      t = 1,
      age = 65,
      sofa = 5,
      tx = 0,
      parameter = 0,
      extra_params = extra_params
    )

    expect_length(single_lp, 1)
    expect_false(is.na(single_lp))

    # Test with different previous states
    for (state in 1:5) {
      lp_state <- lp_violet(
        yprev = state,
        t = 5,
        age = 65,
        sofa = 5,
        tx = 0,
        parameter = 0,
        extra_params = extra_params
      )
      expect_false(is.na(lp_state))
    }
  })
})

describe("sim_trajectories_markov()", {
  it("works with all defaults", {
    data("violet_baseline", package = "markov.misc")
    test_baseline <- violet_baseline[1:100, ]

    trajectories <- sim_trajectories_markov(
      baseline_data = test_baseline,
      follow_up_time = 30,
      parameter = log(0.8),
      seed = 12345
    )

    expect_equal(nrow(trajectories), 100 * 30)
    expect_contains(
      names(trajectories),
      c("id", "time", "y", "yprev", "age", "sofa", "tx")
    )
    expect_type(trajectories$id, "integer")
    expect_type(trajectories$time, "integer")
    expect_type(trajectories$y, "double")
    expect_type(trajectories$yprev, "double")

    expect_equal(unique(trajectories$time), 1:30)
    expect_true(all(table(trajectories$id) == 30))
    expect_all_true(trajectories$y %in% 1:6)
    expect_all_true(trajectories$yprev %in% 1:6)

    # Check death state (6) is absorbing
    deaths <- trajectories[trajectories$y == 6, ]
    if (nrow(deaths) > 0) {
      for (patient_id in unique(deaths$id)) {
        patient_data <- trajectories[trajectories$id == patient_id, ]
        first_death <- min(patient_data$time[patient_data$y == 6])
        subsequent_states <- patient_data$y[patient_data$time >= first_death]
        expect_all_true(subsequent_states == 6)
      }
    }
  })

  it("works with minimal arguments", {
    data("violet_baseline", package = "markov.misc")
    trajectories_minimal <- sim_trajectories_markov(
      baseline_data = violet_baseline[1:50, ],
      parameter = log(0.75),
      seed = 99999
    )

    expect_equal(nrow(trajectories_minimal), 50 * 60) # Default follow_up_time = 60
    expect_s3_class(trajectories_minimal, "data.frame")
    expect_contains(
      names(trajectories_minimal),
      c("id", "time", "y", "yprev", "age", "sofa", "tx")
    )
  })

  it("produces reproducible results", {
    data("violet_baseline", package = "markov.misc")
    test_baseline <- violet_baseline[1:20, ]

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

    expect_identical(traj1, traj2)

    traj3 <- sim_trajectories_markov(
      baseline_data = test_baseline,
      follow_up_time = 10,
      parameter = log(0.8),
      seed = 999
    )

    expect_false(identical(traj1$y, traj3$y))
  })

  it("handles different treatment effects", {
    data("violet_baseline", package = "markov.misc")
    test_baseline <- violet_baseline[1:50, ]

    # Test with null effect
    traj_null <- sim_trajectories_markov(
      baseline_data = test_baseline,
      follow_up_time = 20,
      parameter = 0,
      seed = 123
    )
    expect_equal(nrow(traj_null), 50 * 20)

    # Test with negative effect
    traj_neg <- sim_trajectories_markov(
      baseline_data = test_baseline,
      follow_up_time = 20,
      parameter = log(1.2),
      seed = 123
    )
    expect_equal(nrow(traj_neg), 50 * 20)

    # Test with positive effect
    traj_pos <- sim_trajectories_markov(
      baseline_data = test_baseline,
      follow_up_time = 20,
      parameter = log(0.7),
      seed = 123
    )
    expect_equal(nrow(traj_pos), 50 * 20)
  })
})
