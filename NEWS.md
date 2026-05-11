# markov.misc 0.0.2.0000

- `avg_sops()`, `sops()`, and related Markov SOP workflows now reject models fit with offsets because offsets are not supported by the package prediction paths.
- `avg_sops()`, `sops()`, `soprob_markov()`, and `inferences()` now treat `rms::orm()` as a first-class proportional-odds backend, including fast-path SOP prediction, MVN inference with full `rms::robcov()` covariance matrices, refit bootstrap inference, and score-bootstrap inference via `inferences(..., engine = "score_bootstrap", cluster = <id>)`.
- `prepare_markov_data()`, `soprob_markov()`, `avg_sops()`, and `inferences()` now support numeric previous-state effects, including nonlinear terms such as `rms::rcs(yprev, 6)`, while preserving factor previous-state behavior by default.
- `robcov_vglm()` now stores `bread = vcov(fit)` and unscaled score crossproducts as `meat`, aligns pre-NA cluster vectors when possible, rejects missing or single-valued clusters, and warns that clustered z-test p-values may be anti-conservative with few clusters.
