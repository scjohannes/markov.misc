## Code to prepare `violet_baseline` dataset
##
## This script creates the default baseline dataset for markov.misc package.
## The dataset is based on the VIOLET trial data (Hmisc::simlongord).
##
## Credit: Frank Harrell created the simlongord dataset in the Hmisc package.
## This baseline dataset is derived from that source.

library(Hmisc)
library(dplyr)

# Load VIOLET data
getHdata(simlongord)

# Transform to 4-state model first
simlongord <- simlongord |>
  mutate(
    y = case_when(
      y == "Home" ~ 1,
      y == "In Hospital/Facility" ~ 2,
      y == "Vent/ARDS" ~ 3,
      y == "Dead" ~ 4
    ),
    yprev = case_when(
      yprev == "Home" ~ "1",
      yprev == "In Hospital/Facility" ~ "2",
      yprev == "Vent/ARDS" ~ "3"
    ),
    yprev = factor(yprev, levels = c(1:3))
  )

# Extract baseline data (time == 1) from first 10000 patients
X_base_original <- simlongord |>
  select(id, yprev, age, sofa) |>
  # filter(id <= 10000, !is.na(yprev)) |>
  distinct(id, .keep_all = TRUE) |>
  mutate(
    # Map original yprev to expanded 6-state scale
    # Original yprev=2 (Hospital) → new states 2,3,4 with probabilities
    # Original yprev=3 (Vent) → new state 5
    yprev = case_when(
      yprev == "2" ~ sample(
        c(2, 3, 4),
        n(),
        replace = TRUE,
        prob = c(0.5, 0.3, 0.2)
      ),
      yprev == "3" ~ 5,
      TRUE ~ as.numeric(as.character(yprev))
    )
  )

# Reduce proportion of ventilated patients from ~30% to ~5%
# This makes the baseline more representative of typical respiratory illness of
# interest for PROACT
X_base_vent <- X_base_original |> filter(yprev == 5)
X_base_other <- X_base_original |> filter(yprev != 5)

# Calculate target sample sizes
n_vent_original <- nrow(X_base_vent)
n_vent_reduced <- round(n_vent_original * 0.05) # 5% ventilated
n_other <- 250000 - n_vent_reduced # Fill the rest with non-vent

# Resample with replacement to achieve target distribution
set.seed(12345) # For reproducibility of package data
X_base_resampled <- bind_rows(
  X_base_vent |> slice_sample(n = n_vent_reduced, replace = TRUE),
  X_base_other |> slice_sample(n = n_other, replace = TRUE)
)

# Create final baseline dataset
violet_baseline <- X_base_resampled |>
  mutate(
    id = row_number(),
    tx = sample(c(0, 1), n(), replace = TRUE)
  ) |>
  select(id, yprev, age, sofa, tx)

# Save as package data
usethis::use_data(violet_baseline, overwrite = TRUE)
