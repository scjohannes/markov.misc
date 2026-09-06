# Repository Guidelines

## Project Structure & Module Organization
`markov.misc` is an R package. Core source code lives in `R/` (for example, `R/sops.R`, `R/simulate_trajectories.R`, `R/bootstrap.R`, `R/power.R`). Tests are in `tests/testthat/` with a package entrypoint at `tests/testthat.R`. Generated documentation is in `man/`, long-form analyses are in `vignettes/`, and reusable datasets are in `data/` with creation scripts in `data-raw/`. Keep architecture notes in `ARCHITECTURE.md` aligned with code changes.

## Architecture Documentation
When making code changes, update `ARCHITECTURE.md` in the same work so it stays aligned with the package design. Use the `describe-design` skill for architecture updates, especially when changes affect module boundaries, public APIs, data contracts, model backends, inference workflows, simulation behavior, endpoint summaries, or operating-characteristic workflows.

## Build, Test, and Development Commands
Use standard R package workflows from the repository root:

- `R -q -e "devtools::document()"`: regenerate `NAMESPACE` and `man/*.Rd` from roxygen2 comments.
- `R -q -e "devtools::test()"`: run the full testthat suite.
- `R -q -e "devtools::check(build_vignettes = FALSE)"`: run package checks (tests, examples, metadata) without building vignettes.
- `R -q -e "remotes::install_local('.')"`: install the local package for interactive use.
- `R -q -e "testthat::test_file('tests/testthat/test-sops.R')"`: run a focused test file.

By default, skip vignette building when running package checks because it takes too long. Use `devtools::check(build_vignettes = FALSE)`; build vignettes only when explicitly requested.

## Coding Style & Naming Conventions
Follow tidyverse-oriented R style: 2-space indentation, `|>` pipes, and explicit namespaces for non-base calls (for example `dplyr::mutate`, `stats::predict`). Use `snake_case` for functions and variables. Preserve project naming for state-model fields (`id`, `time`, `y`, `yprev`, `tx`). Exported functions require roxygen2 docs with `@export`; add `@importFrom` tags as needed. If you use NSE column names in dplyr/ggplot2, register them in `R/globals.R`.

## Writing for Package Users
- Write help pages and vignettes for readers who know neither the implementation nor prior discussions. Explain what they can calculate, which arguments to use, what the output means, and any restrictions that affect their choices.
- Make explanations standalone. Do not rely on an inaccessible manuscript, internal design notes, or development history to define a method or justify its behavior.
- Use the current public API's terminology consistently. For analytical confidence intervals, use **conditional variance** and **unconditional variance**, selected through `vcov` in `inferences()`. Conditional variance accounts for coefficient estimation while treating the patients' starting states and covariates used for prediction as given. Unconditional variance also accounts for which patients were sampled and their role in estimating the coefficients.
- Apply conditional/unconditional variance terminology throughout the repository, including internal function and variable names, result metadata, errors, tests, filenames, and developer notes. Reserve "superpopulation" for an actual source population in simulation or sampling descriptions, not as a name for a variance estimate.
- Avoid invented labels and unexplained implementation jargon such as "fixed-profile", "empirical-cohort", "profile-distribution uncertainty", "Delta results", "dense covariance", or "Jacobian inspection". Describe the concrete meaning instead: for example, "a covariance for every pair of estimates" or "how each estimate changes with each model coefficient".
- Introduce necessary statistical terms and formulas with a plain-language explanation, define their symbols, and connect them to the user's calculation. Keep storage formats, internal metadata, and execution details in developer documentation unless they help the user make a practical decision.
- Explain examples and output directly. For example, say that the square roots of a covariance matrix's diagonal entries are standard errors, and that requesting fewer rows reduces memory use. Do not merely rename technical concepts with other jargon.

## Testing Guidelines
Testing uses `testthat` edition 3. Name files `test-<feature>.R` and keep scenarios deterministic (`set.seed()`, fixed fixtures). Add regression tests for bug fixes and edge cases (absorbing states, model-class compatibility, bootstrap paths). Snapshot updates in `tests/testthat/_snaps/` should only be committed when behavior changes intentionally.

## Commit & Pull Request Guidelines
Recent history favors short, imperative commit messages (often prefixed with `//`). Keep commits focused and include regenerated docs when interfaces change. Pull requests should include: a brief problem/solution summary, linked issue (if available), affected files/modules, and evidence of validation (`devtools::test()` and, for larger changes, `devtools::check(build_vignettes = FALSE)`).
