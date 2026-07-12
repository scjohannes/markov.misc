markov_native_run <- function(initial, transitions, non_absorb, absorb) {
  cpp_markov_propagate(initial, transitions, non_absorb, absorb)
}

markov_update_logits_native <- function(
  previous,
  logits,
  non_absorb,
  absorb
) {
  cpp_markov_update_logits(previous, logits, non_absorb, absorb)
}

markov_update_po_native <- function(
  previous,
  scalar_predictor,
  cutpoints,
  non_absorb,
  absorb
) {
  cpp_markov_update_po(
    previous,
    scalar_predictor,
    cutpoints,
    non_absorb,
    absorb
  )
}

normalize_probability_array_native <- function(probs) {
  dims <- dim(probs)
  values <- if (is.double(probs)) probs else as.numeric(probs)
  out <- cpp_normalize_probability_array(
    values,
    dims[1],
    dims[2],
    dims[3]
  )
  dimnames(out) <- dimnames(probs)
  out
}

reduce_sops_draw_array <- function(probs, groups, group_count) {
  dims <- dim(probs)
  if (length(dims) != 4L || length(groups) != dims[2]) {
    stop("Posterior SOP array and grouping index are not aligned.")
  }
  values <- if (is.double(probs)) probs else as.numeric(probs)
  cpp_reduce_sops_draws(
    values,
    dims[1],
    dims[2],
    dims[3],
    dims[4],
    as.integer(groups),
    as.integer(group_count)
  )
}

markov_update_draws_native <- function(
  previous,
  transition,
  non_absorb,
  absorb
) {
  dims <- dim(previous)
  previous_values <- if (is.double(previous)) previous else as.numeric(previous)
  transition_values <- if (is.double(transition)) {
    transition
  } else {
    as.numeric(transition)
  }
  cpp_markov_update_draws(
    previous_values,
    transition_values,
    dims[1],
    dims[2],
    dims[3],
    as.integer(non_absorb),
    as.integer(absorb)
  )
}

markov_update_second_order_native <- function(
  previous,
  transition,
  older,
  current,
  absorb
) {
  dims <- dim(previous)
  cpp_markov_update_second_order(
    previous,
    transition,
    dims[1L],
    dims[2L],
    as.integer(older),
    as.integer(current),
    as.integer(absorb)
  )
}

blrm_probabilities_native <- function(
  base_eta,
  intercepts,
  threshold_eta = numeric(),
  threshold_scale = numeric()
) {
  cpp_blrm_probabilities(
    base_eta,
    intercepts,
    threshold_eta,
    threshold_scale,
    nrow(base_eta),
    ncol(base_eta),
    ncol(intercepts)
  )
}

markov_update_second_order_po_native <- function(
  previous,
  scalar_predictor,
  cutpoints,
  older,
  current,
  absorb
) {
  dims <- dim(previous)
  cpp_markov_update_second_order_po(
    previous,
    scalar_predictor,
    cutpoints,
    dims[1L],
    dims[2L],
    as.integer(older),
    as.integer(current),
    as.integer(absorb)
  )
}
