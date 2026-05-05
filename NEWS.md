# markov.misc 0.0.2.0000

- `prepare_markov_data()`, `soprob_markov()`, `avg_sops()`, and `inferences()` now support numeric previous-state effects, including nonlinear terms such as `rms::rcs(yprev, 6)`, while preserving factor previous-state behavior by default.
