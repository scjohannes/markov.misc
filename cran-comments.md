## Submission Context

This is a first submission of markov.misc to CRAN.

## Test Environments

- local Windows 11 x64 (build 26200), R 4.5.2
- win-builder, R-devel: pending
- macOS builder or GitHub Actions macOS, R release: pending
- Linux builder or GitHub Actions Linux, R release: pending

## R CMD check Results

0 errors | 0 warnings | 0 notes

Local command:

```r
R CMD check markov.misc_0.1.0.tar.gz --no-manual
```

Manual PDF building was skipped locally because the Windows environment used for
this check is not configured with LaTeX. The source package builds with
vignettes, and package vignettes are rebuilt successfully during check.

## Reverse Dependencies

There are no reverse dependencies because this is a first CRAN submission.

## Additional Comments

The package includes `violet_baseline`, a dataset derived from
`Hmisc::simlongord`. `Hmisc` is licensed under GPL (>= 2), and this package is
therefore licensed under GPL (>= 2). The data preparation script is included in
`data-raw/violet_baseline.R`, and attribution is documented in
`?violet_baseline`.
