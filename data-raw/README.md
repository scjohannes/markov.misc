# data-raw

This directory contains R scripts used to prepare package datasets.

## Files

- `violet_baseline.R`: Creates the `violet_baseline` dataset from Hmisc::simlongord data

## Usage

To regenerate the package data:

```r
# Install required packages
install.packages("Hmisc")
install.packages("usethis")

# Run the data preparation script
source("data-raw/violet_baseline.R")
```

This will create/update the `data/violet_baseline.rda` file that is included in the package.

## Attribution

The `violet_baseline` dataset is derived from the `simlongord` dataset created by Frank Harrell, available in the Hmisc package. The original data is based on the VIOLET trial.
