# Internal serializable execution plans for frequentist first-order SOPs.

compile_sop_execution_plan <- function(
  model,
  newdata,
  times,
  y_levels,
  absorb = NULL,
  time_var = "time",
  p_var = "yprev",
  gap_var = NULL,
  time_covariates = NULL,
  builder = c("batched", "streamed"),
  output = "array"
) {
  builder <- match.arg(builder)
  validate_markov_model(model)
  resolved_times <- resolve_sop_times(
    model,
    newdata,
    times,
    time_var,
    time_covariates = time_covariates,
    default = "fast"
  )$times
  build <- if (identical(builder, "batched")) {
    markov_msm_build_batched
  } else {
    markov_msm_build
  }
  components <- build(
    model = model,
    newdata = newdata,
    time_covariates = time_covariates,
    times = resolved_times,
    y_levels = y_levels,
    absorb = absorb,
    time_var = time_var,
    p_var = p_var,
    gap_var = gap_var
  )
  structure(
    list(
      version = 1L,
      model_class = class(model)[1L],
      link = "logit",
      recursion_order = 1L,
      times = resolved_times,
      y_levels = as_state_labels(y_levels),
      absorb = as_state_labels(absorb),
      time_var = time_var,
      p_var = p_var,
      gap_var = gap_var,
      output = output,
      workspace_bytes = 256 * 1024^2,
      basis_terms = attr(model, "markov_basis_terms", exact = TRUE),
      components = components
    ),
    class = "markov_sop_exec_plan"
  )
}

run_sop_execution_plan <- function(plan, Gamma) {
  if (!inherits(plan, "markov_sop_exec_plan") || plan$version != 1L) {
    stop("Invalid or unsupported SOP execution plan.")
  }
  markov_msm_run(
    plan$components,
    Gamma,
    times = plan$times,
    absorb = plan$absorb
  )
}
