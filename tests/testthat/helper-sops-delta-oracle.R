run_sop_delta_r_oracle <- function(plan, model) {
  coef <- sop_delta_raw_coef(model)
  map <- get_effective_coef_map(model)
  validated <- sop_delta_validate_plan(plan, model, map, coef)
  components <- validated$components
  Gamma <- validated$Gamma
  map <- validated$map

  n <- components$n_pat
  n_times <- length(plan$times)
  n_states <- components$n_states
  n_coef <- length(coef)
  absorb <- which(
    as.character(plan$y_levels) %in% as.character(plan$absorb)
  )
  non_absorb <- setdiff(seq_len(n_states), absorb)

  initial <- sop_delta_raw_probabilities(
    components$X_init,
    Gamma,
    map,
    "the first SOP time point"
  )
  probabilities <- array(0, dim = c(n, n_times, n_states))
  jacobian <- array(0, dim = c(n, n_times, n_states, n_coef))
  probabilities[, 1L, ] <- initial$probabilities
  jacobian[, 1L, , ] <- initial$derivative
  current_probability <- initial$probabilities
  current_jacobian <- initial$derivative

  if (n_times >= 2L) {
    for (visit in 2L:n_times) {
      X_visit <- sop_delta_visit_design(components, visit, non_absorb)
      next_probability <- matrix(0, nrow = n, ncol = n_states)
      next_jacobian <- array(0, dim = c(n, n_states, n_coef))

      for (origin in absorb) {
        next_probability[, origin] <- current_probability[, origin]
        next_jacobian[, origin, ] <- current_jacobian[, origin, ]
      }
      for (origin_position in seq_along(non_absorb)) {
        origin <- non_absorb[origin_position]
        rows <- (origin_position - 1L) * n + seq_len(n)
        transition <- sop_delta_raw_probabilities(
          X_visit[rows, , drop = FALSE],
          Gamma,
          map,
          paste0(
            "SOP time point ",
            visit,
            ", origin state ",
            plan$y_levels[origin]
          )
        )
        origin_probability <- current_probability[, origin]
        for (destination in seq_len(n_states)) {
          destination_probability <- transition$probabilities[, destination]
          next_probability[, destination] <-
            next_probability[, destination] +
            origin_probability * destination_probability
          for (coefficient in seq_len(n_coef)) {
            next_jacobian[, destination, coefficient] <-
              next_jacobian[, destination, coefficient] +
              current_jacobian[, origin, coefficient] *
                destination_probability +
              origin_probability *
                transition$derivative[, destination, coefficient]
          }
        }
      }
      probabilities[, visit, ] <- next_probability
      jacobian[, visit, , ] <- next_jacobian
      current_probability <- next_probability
      current_jacobian <- next_jacobian
    }
  }

  dimnames(probabilities) <- list(NULL, NULL, plan$y_levels)
  dimnames(jacobian) <- list(NULL, NULL, plan$y_levels, names(coef))
  list(probabilities = probabilities, jacobian = jacobian)
}
