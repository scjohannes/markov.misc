markov_native_run <- function(initial, transitions, non_absorb, absorb) {
  cpp_markov_propagate(initial, transitions, non_absorb, absorb)
}

normalize_probability_array_native <- function(probs) {
  dims <- dim(probs)
  out <- cpp_normalize_probability_array(
    as.numeric(probs),
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
  cpp_reduce_sops_draws(
    as.numeric(probs),
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
  cpp_markov_update_draws(
    as.numeric(previous),
    as.numeric(transition),
    dims[1],
    dims[2],
    dims[3],
    as.integer(non_absorb),
    as.integer(absorb)
  )
}
