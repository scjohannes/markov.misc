# Patient-cluster covariance and superpopulation delta-method helpers.

delta_superpopulation_coefficients <- function(model) {
  beta <- get_coef(model)
  if (
    !is.numeric(beta) ||
      length(beta) < 1L ||
      any(!is.finite(beta)) ||
      is.null(names(beta)) ||
      anyNA(names(beta)) ||
      any(names(beta) == "") ||
      anyDuplicated(names(beta))
  ) {
    stop(
      "Analytical inference requires a finite coefficient vector with ",
      "unique, non-missing names.",
      call. = FALSE
    )
  }
  beta
}

delta_superpopulation_model <- function(model) {
  if (inherits(model, "robcov_vglm")) {
    fit <- model$vglm_fit
    if (is.null(fit) || !inherits(fit, "vglm")) {
      stop(
        "The `robcov_vglm` object does not contain its fitted vglm model.",
        call. = FALSE
      )
    }
    return(fit)
  }
  model
}

delta_validate_full_po_model <- function(model) {
  fit <- delta_superpopulation_model(model)

  if (inherits(fit, "blrm")) {
    stop(
      "Analytical superpopulation inference supports only frequentist full ",
      "proportional-odds orm and vglm models; blrm models are not supported.",
      call. = FALSE
    )
  }
  if (!inherits(fit, c("orm", "vglm"))) {
    stop(
      "Analytical superpopulation inference supports only full proportional-odds ",
      "orm and vglm models.",
      call. = FALSE
    )
  }

  validate_markov_model(fit)

  if (inherits(fit, "vglm")) {
    constraints <- VGAM::constraints(fit)
    constraint_names <- names(constraints)
    slope_index <- if (is.null(constraint_names)) {
      seq_along(constraints)[-1L]
    } else {
      which(constraint_names != "(Intercept)")
    }
    is_common_slope <- vapply(
      constraints[slope_index],
      function(constraint) {
        constraint <- as.matrix(constraint)
        ncol(constraint) == 1L &&
          all(is.finite(constraint)) &&
          max(abs(constraint[, 1L] - constraint[1L, 1L])) <= 1e-12
      },
      logical(1)
    )
    if (length(is_common_slope) && !all(is_common_slope)) {
      stop(
        "Analytical inference currently requires a full proportional-odds ",
        "vglm model (`parallel = TRUE`); partial or nonproportional-odds ",
        "constraints are not supported.",
        call. = FALSE
      )
    }
  }

  invisible(fit)
}

delta_validate_named_matrix <- function(
  value,
  names,
  label,
  positive_semidefinite = FALSE,
  positive_definite = FALSE
) {
  if (methods::is(value, "Matrix")) {
    value <- as.matrix(value)
  }
  q <- length(names)
  if (
    !is.matrix(value) || !is.numeric(value) || !identical(dim(value), c(q, q))
  ) {
    stop(
      "`",
      label,
      "` must be a numeric [",
      q,
      " x ",
      q,
      "] matrix with one row and column per raw model coefficient.",
      call. = FALSE
    )
  }
  row_names <- rownames(value)
  column_names <- colnames(value)
  valid_names <- function(x) {
    !is.null(x) &&
      length(x) == q &&
      !anyNA(x) &&
      !any(x == "") &&
      !anyDuplicated(x)
  }
  if (
    !valid_names(row_names) ||
      !valid_names(column_names) ||
      !setequal(row_names, names) ||
      !setequal(column_names, names)
  ) {
    stop(
      "Rows and columns of `",
      label,
      "` must be uniquely named and match the complete raw coefficient vector.",
      call. = FALSE
    )
  }
  value <- value[names, names, drop = FALSE]
  if (any(!is.finite(value))) {
    stop("`", label, "` contains non-finite values.", call. = FALSE)
  }

  scale <- max(1, max(abs(value)))
  symmetry_error <- max(abs(value - t(value))) / scale
  if (!is.finite(symmetry_error) || symmetry_error > 1e-8) {
    stop("`", label, "` is not numerically symmetric.", call. = FALSE)
  }
  value <- (value + t(value)) / 2

  if (positive_semidefinite || positive_definite) {
    eigenvalues <- eigen(value, symmetric = TRUE, only.values = TRUE)$values
    eigen_scale <- max(1, max(abs(eigenvalues))) * q
    semidefinite_tolerance <- sqrt(.Machine$double.eps) * eigen_scale
    definite_tolerance <- 100 * .Machine$double.eps * eigen_scale
    if (
      positive_semidefinite &&
        min(eigenvalues) < -semidefinite_tolerance
    ) {
      stop(
        "`",
        label,
        "` is not numerically positive semidefinite.",
        call. = FALSE
      )
    }
    if (positive_definite && min(eigenvalues) <= definite_tolerance) {
      stop(
        "`",
        label,
        "` is not numerically positive definite.",
        call. = FALSE
      )
    }
  }

  value
}

delta_fitting_row_count <- function(model) {
  if (inherits(model, "robcov_vglm")) {
    if (!is.matrix(model$scores)) {
      stop(
        "The `robcov_vglm` object does not contain row-aligned score ",
        "contributions. Refit it with `robcov_vglm()`.",
        call. = FALSE
      )
    }
    return(nrow(model$scores))
  }

  fit <- delta_superpopulation_model(model)
  if (inherits(fit, "vglm")) {
    return(as.integer(stats::nobs(fit, type = "lm")))
  }
  if (inherits(fit, "orm")) {
    if (!is.null(fit$x)) {
      return(nrow(fit$x))
    }
    if (!is.null(fit$y)) {
      return(length(fit$y))
    }
    n <- tryCatch(stats::nobs(fit), error = function(e) NULL)
    if (!is.null(n) && length(n) == 1L && is.finite(n)) {
      return(as.integer(n))
    }
    stop(
      "Could not determine the fitted orm row count. Refit with ",
      "`x = TRUE, y = TRUE`.",
      call. = FALSE
    )
  }

  stop(
    "Analytical superpopulation inference supports only full proportional-odds ",
    "orm and vglm models.",
    call. = FALSE
  )
}

delta_align_explicit_cluster <- function(model, cluster, n_rows) {
  fit <- delta_superpopulation_model(model)
  if (inherits(fit, "vglm")) {
    return(align_cluster_vglm(fit, cluster, n_rows))
  }
  align_cluster_orm(fit, cluster, n_rows)
}

delta_validate_cluster <- function(cluster, n_rows) {
  if (length(cluster) != n_rows) {
    stop(
      "The patient cluster vector has length ",
      length(cluster),
      ", but the fitted model has ",
      n_rows,
      " score rows.",
      call. = FALSE
    )
  }
  if (is.list(cluster) || is.data.frame(cluster) || is.matrix(cluster)) {
    stop(
      "`cluster` must be a one-dimensional patient ID vector.",
      call. = FALSE
    )
  }
  if (anyNA(cluster)) {
    stop(
      "The patient cluster vector contains missing values after alignment ",
      "with the fitted rows.",
      call. = FALSE
    )
  }
  cluster_character <- as.character(cluster)
  if (anyNA(cluster_character) || any(cluster_character == "")) {
    stop(
      "Patient cluster IDs must be non-missing and non-empty.",
      call. = FALSE
    )
  }
  ids <- unique(cluster_character)
  if (length(ids) < 2L) {
    stop(
      "Patient-cluster inference requires at least two independent patients.",
      call. = FALSE
    )
  }
  list(cluster = cluster, ids = ids)
}

# Resolve patient clusters for analytical superpopulation inference. Explicit
# row-level IDs take precedence; fitting rows are never implicit clusters.
resolve_delta_cluster <- function(model, cluster = NULL) {
  delta_validate_full_po_model(model)
  n_rows <- delta_fitting_row_count(model)
  source <- NULL
  id_var <- NULL

  if (!is.null(cluster)) {
    if (inherits(cluster, "formula")) {
      cluster_vars <- all.vars(cluster)
      if (length(cluster_vars) != 1L) {
        stop(
          "A `cluster` formula must select exactly one patient ID column, ",
          "for example `~ id`.",
          call. = FALSE
        )
      }
      data <- markov_model_data(model)
      id_var <- cluster_vars[[1L]]
      if (is.null(data) || !is.data.frame(data) || !id_var %in% names(data)) {
        stop(
          "Explicit cluster formula `~ ",
          id_var,
          "` requires that column in the model's stored fitting data.",
          call. = FALSE
        )
      }
      fit <- delta_superpopulation_model(model)
      data <- markov_align_model_data(data, fit)
      cluster <- data[[id_var]]
      source <- "explicit_formula"
    } else {
      cluster <- delta_align_explicit_cluster(model, cluster, n_rows)
      source <- "explicit_vector"
    }
  } else {
    id_var <- markov_model_id_var(model)
    data <- markov_model_data(model)
    if (
      is.null(id_var) ||
        !is.character(id_var) ||
        length(id_var) != 1L ||
        is.na(id_var) ||
        is.null(data) ||
        !is.data.frame(data)
    ) {
      stop(
        "Analytical superpopulation inference requires patient clustering. ",
        "Supply ",
        "`cluster`, or fit with `orm_markov(..., id_var = ...)` or ",
        "`vglm_markov(..., id_var = ...)` so row-aligned fitting data and ",
        "patient-ID metadata are stored. Observation rows are not used as ",
        "implicit clusters.",
        call. = FALSE
      )
    }
    if (!id_var %in% names(data)) {
      stop(
        "Stored patient ID column `",
        id_var,
        "` is absent from the stored fitting data.",
        call. = FALSE
      )
    }
    fit <- delta_superpopulation_model(model)
    data <- markov_align_model_data(data, fit)
    cluster <- data[[id_var]]
    source <- "stored_id_var"
  }

  checked <- delta_validate_cluster(cluster, n_rows)
  list(
    cluster = checked$cluster,
    ids = checked$ids,
    metadata = list(
      source = source,
      id_var = id_var,
      cluster_level = "patient",
      n_rows = n_rows,
      n_clusters = length(checked$ids)
    )
  )
}

delta_vglm_covariance_metadata <- function(robust, source) {
  list(
    source = source,
    covariance = "patient_cluster_robust",
    backend = "vglm",
    type = robust$type,
    cadjust = robust$cadjust,
    adjustment_factor = robust$adjustment_factor,
    bread_type = robust$bread_type,
    n_clusters = robust$n_clusters
  )
}

delta_orm_backend_aliases <- function(model, coefficient_names) {
  threshold_count <- as.integer(model$non.slopes %||% 0L)
  slope_count <- length(coefficient_names) - threshold_count
  if (
    threshold_count < 1L ||
      slope_count < 0L ||
      threshold_count + slope_count != length(coefficient_names)
  ) {
    return(list(coefficient_names))
  }

  threshold_names <- coefficient_names[seq_len(threshold_count)]
  design <- model$Design
  slope_aliases <- list(
    design$colnames,
    design$mmcolnames,
    colnames(model$x)
  )
  aliases <- c(
    list(coefficient_names),
    lapply(slope_aliases, function(value) {
      if (is.null(value) || length(value) != slope_count) {
        return(NULL)
      }
      c(threshold_names, as.character(value))
    })
  )
  Filter(
    function(value) {
      length(value) == length(coefficient_names) &&
        !anyNA(value) &&
        !any(value == "") &&
        !anyDuplicated(value)
    },
    aliases
  )
}

# rms::robcov() constructs covariance matrices in coefficient order, but
# penalized fits can lose dimnames and inline rms transforms can use model-matrix
# labels instead of the display labels in coef(). Normalize only matrices that
# were produced internally from this exact orm fit.
delta_normalize_orm_backend_matrix <- function(
  value,
  model,
  coefficient_names
) {
  if (methods::is(value, "Matrix")) {
    value <- as.matrix(value)
  }
  q <- length(coefficient_names)
  if (!is.matrix(value) || !identical(dim(value), c(q, q))) {
    return(value)
  }

  aliases <- delta_orm_backend_aliases(model, coefficient_names)
  normalize_axis <- function(axis_names) {
    if (is.null(axis_names)) {
      return(coefficient_names)
    }
    for (candidate in aliases) {
      index <- match(axis_names, candidate)
      if (!anyNA(index) && !anyDuplicated(index)) {
        return(coefficient_names[index])
      }
    }
    axis_names
  }

  rownames(value) <- normalize_axis(rownames(value))
  colnames(value) <- normalize_axis(colnames(value))
  value
}

delta_orm_model_bread <- function(model, coefficient_names) {
  bread <- model$orig.var
  if (is.null(bread)) {
    if (!requireNamespace("rms", quietly = TRUE)) {
      stop("Package 'rms' is required for orm superpopulation inference.")
    }
    robust <- rms::robcov(model)
    bread <- robust$orig.var
  }
  if (is.null(bread)) {
    stop(
      "Could not obtain the full model-based covariance for orm superpopulation ",
      "inference. Refit with `x = TRUE, y = TRUE`.",
      call. = FALSE
    )
  }

  delta_normalize_orm_backend_matrix(
    bread,
    model,
    coefficient_names
  )
}

delta_reject_weighted_orm_superpopulation <- function(model) {
  weights <- model$weights
  if (
    is.null(weights) ||
      length(weights) == 0L ||
      all(is.finite(weights) & weights == 1)
  ) {
    return(invisible(NULL))
  }

  stop(
    "Superpopulation delta inference does not currently support orm fits with ",
    "non-unit case weights because the required row score contributions ",
    "must incorporate those weights. Use `target = \"empirical\"` or refit ",
    "without case weights.",
    call. = FALSE
  )
}

delta_reject_penalized_orm_superpopulation <- function(model) {
  penalty_matrix <- model$penalty.matrix
  penalty <- model$penalty
  penalized <-
    (!is.null(penalty_matrix) &&
      length(penalty_matrix) > 0L &&
      any(is.finite(penalty_matrix) & penalty_matrix != 0)) ||
    (!is.null(penalty) &&
      length(penalty) > 0L &&
      any(is.finite(penalty) & penalty != 0))
  if (!penalized) {
    return(invisible(NULL))
  }

  stop(
    "Superpopulation delta inference does not currently support penalized ",
    "orm likelihoods because the patient scores and sensitivity have not ",
    "been shown to incorporate the penalty consistently. Ordinary spline ",
    "transformations remain supported with `target = \"empirical\"`.",
    call. = FALSE
  )
}

# Resolve a complete raw-coefficient covariance, preferring an explicit matrix
# and otherwise using the backend's patient-cluster robust convention.
get_delta_cluster_vcov <- function(model, cluster = NULL, vcov = NULL) {
  delta_validate_full_po_model(model)
  beta <- delta_superpopulation_coefficients(model)

  if (!is.null(vcov)) {
    covariance <- delta_validate_named_matrix(
      vcov,
      names(beta),
      "vcov",
      positive_semidefinite = TRUE
    )
    return(list(
      vcov = covariance,
      metadata = list(
        source = "explicit",
        covariance = "user_supplied",
        backend = class(delta_superpopulation_model(model))[[1L]],
        type = NA_character_,
        cadjust = NA,
        adjustment_factor = NA_real_,
        bread_type = NA_character_,
        n_clusters = NA_integer_
      ),
      cluster_info = NULL
    ))
  }

  cluster_info <- resolve_delta_cluster(model, cluster = cluster)
  fit <- delta_superpopulation_model(model)

  if (inherits(fit, "vglm")) {
    stored_cluster_matches <- inherits(model, "robcov_vglm") &&
      !is.null(model$cluster) &&
      identical(
        as.character(model$cluster),
        as.character(cluster_info$cluster)
      )
    use_stored <- inherits(model, "robcov_vglm") &&
      is.null(cluster) &&
      !is.null(model$var) &&
      stored_cluster_matches
    robust <- if (use_stored) {
      model
    } else {
      robcov_vglm(
        fit,
        cluster = cluster_info$cluster,
        bread = if (inherits(model, "robcov_vglm")) {
          model$bread_type %||% "observed"
        } else {
          "observed"
        },
        type = if (inherits(model, "robcov_vglm")) {
          model$type %||% "HC0"
        } else {
          "HC0"
        },
        cadjust = if (inherits(model, "robcov_vglm")) {
          model$cadjust %||% TRUE
        } else {
          TRUE
        }
      )
    }
    covariance <- delta_validate_named_matrix(
      robust$var,
      names(beta),
      "patient-cluster robust covariance",
      positive_semidefinite = TRUE
    )
    metadata <- delta_vglm_covariance_metadata(
      robust,
      if (use_stored) "stored_robcov_vglm" else "computed_robcov_vglm"
    )
  } else {
    use_stored <- is.null(cluster) &&
      !is.null(model$orig.var) &&
      !is.null(model$var)
    robust <- if (use_stored) {
      model
    } else {
      if (!requireNamespace("rms", quietly = TRUE)) {
        stop("Package 'rms' is required for orm robust covariance inference.")
      }
      # rms::robcov() starts from stats::vcov(fit). For an already robust orm
      # fit that accessor returns the sandwich covariance, so restore the full
      # inverse-information covariance before applying an explicit clustering.
      covariance_fit <- fit
      if (!is.null(covariance_fit$orig.var)) {
        covariance_fit$var <- covariance_fit$orig.var
        covariance_fit$orig.var <- NULL
        covariance_fit$clusterInfo <- NULL
      }
      rms::robcov(covariance_fit, cluster = cluster_info$cluster)
    }
    backend_covariance <- delta_normalize_orm_backend_matrix(
      robust$var,
      fit,
      names(beta)
    )
    covariance <- delta_validate_named_matrix(
      backend_covariance,
      names(beta),
      "patient-cluster robust covariance",
      positive_semidefinite = TRUE
    )
    metadata <- list(
      source = if (use_stored) "stored_rms_robcov" else "computed_rms_robcov",
      covariance = "patient_cluster_robust",
      backend = "orm",
      type = "HC0",
      cadjust = FALSE,
      adjustment_factor = 1,
      bread_type = "rms_model",
      n_clusters = length(cluster_info$ids)
    )
  }

  metadata$cluster_source <- cluster_info$metadata$source
  metadata$cluster_level <- "patient"
  list(
    vcov = covariance,
    metadata = metadata,
    cluster_info = cluster_info
  )
}

delta_validate_scores <- function(scores, coefficient_names, n_rows) {
  scores <- as.matrix(scores)
  if (
    !is.numeric(scores) ||
      !identical(dim(scores), c(n_rows, length(coefficient_names)))
  ) {
    stop(
      "Row score contributions must be a numeric [",
      n_rows,
      " x ",
      length(coefficient_names),
      "] matrix.",
      call. = FALSE
    )
  }
  score_names <- colnames(scores)
  if (
    is.null(score_names) ||
      anyNA(score_names) ||
      any(score_names == "") ||
      anyDuplicated(score_names) ||
      !setequal(score_names, coefficient_names)
  ) {
    stop(
      "Score columns must uniquely match the complete raw coefficient vector.",
      call. = FALSE
    )
  }
  scores <- scores[, coefficient_names, drop = FALSE]
  if (any(!is.finite(scores))) {
    stop("Row score contributions contain non-finite values.", call. = FALSE)
  }
  scores
}

# Aggregate raw likelihood scores by patient and construct inverse
# per-participant sensitivity for the stacked influence function.
get_delta_score_components <- function(model, cluster = NULL) {
  delta_validate_full_po_model(model)
  beta <- delta_superpopulation_coefficients(model)
  fit <- delta_superpopulation_model(model)
  if (inherits(fit, "orm")) {
    delta_reject_weighted_orm_superpopulation(fit)
    delta_reject_penalized_orm_superpopulation(fit)
  }
  cluster_info <- resolve_delta_cluster(model, cluster = cluster)
  n_rows <- cluster_info$metadata$n_rows

  if (inherits(model, "robcov_vglm")) {
    row_scores <- model$scores
    bread <- model$bread
    component_metadata <- list(
      bread_source = paste0("robcov_vglm_", model$bread_type %||% "observed"),
      backend_hc_type_ignored = model$type %||% NA_character_,
      backend_cadjust_ignored = model$cadjust %||% NA
    )
  } else if (inherits(fit, "vglm")) {
    row_scores <- compute_scores_vglm(fit)
    bread <- compute_observed_bread_vglm(fit)$bread
    component_metadata <- list(
      bread_source = "vglm_observed_information",
      backend_hc_type_ignored = NA_character_,
      backend_cadjust_ignored = NA
    )
  } else {
    row_scores <- compute_scores_orm(fit)
    bread <- delta_orm_model_bread(fit, names(beta))
    component_metadata <- list(
      bread_source = "rms_inverse_total_information",
      backend_hc_type_ignored = NA_character_,
      backend_cadjust_ignored = NA
    )
  }

  row_scores <- delta_validate_scores(row_scores, names(beta), n_rows)
  bread <- delta_validate_named_matrix(
    bread,
    names(beta),
    "model-based bread",
    positive_definite = TRUE
  )

  cluster_character <- as.character(cluster_info$cluster)
  cluster_factor <- factor(cluster_character, levels = cluster_info$ids)
  scores <- rowsum(row_scores, cluster_factor, reorder = FALSE)
  scores <- as.matrix(scores)
  rownames(scores) <- cluster_info$ids
  colnames(scores) <- names(beta)

  n <- nrow(scores)
  Ainv <- n * bread
  dimnames(Ainv) <- list(names(beta), names(beta))

  component_metadata$cluster_source <- cluster_info$metadata$source
  component_metadata$n_patients <- n
  list(
    ids = cluster_info$ids,
    scores = scores,
    Ainv = Ainv,
    n = n,
    metadata = component_metadata
  )
}

delta_validate_profile_matrix <- function(value, label) {
  value <- as.matrix(value)
  if (!is.numeric(value) || length(value) == 0L || any(!is.finite(value))) {
    stop(
      "`",
      label,
      "` must be a non-empty finite numeric matrix.",
      call. = FALSE
    )
  }
  value
}

# Combine centered profile functionals and fitted-coefficient influence for the
# fitted-cohort superpopulation target.
delta_stacked_influence <- function(
  individual_values,
  average_jacobian,
  profile_ids,
  score_components
) {
  values <- delta_validate_profile_matrix(
    individual_values,
    "individual_values"
  )
  G <- delta_validate_profile_matrix(average_jacobian, "average_jacobian")
  n <- nrow(values)
  m <- ncol(values)

  if (n < 2L) {
    stop(
      "Stacked superpopulation inference requires at least two profiles.",
      call. = FALSE
    )
  }
  if (length(profile_ids) != n || anyNA(profile_ids)) {
    stop(
      "`profile_ids` must contain one non-missing patient ID per profile row.",
      call. = FALSE
    )
  }
  profile_ids <- as.character(profile_ids)
  if (any(profile_ids == "") || anyDuplicated(profile_ids)) {
    stop("`profile_ids` must be non-empty and unique.", call. = FALSE)
  }
  if (!is.list(score_components)) {
    stop("`score_components` must come from `get_delta_score_components()`.")
  }
  required <- c("ids", "scores", "Ainv", "n")
  missing_required <- required[vapply(
    required,
    function(name) is.null(score_components[[name]]),
    logical(1)
  )]
  if (length(missing_required)) {
    stop(
      "`score_components` is missing: ",
      paste(missing_required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  score_ids <- as.character(score_components$ids)
  scores <- as.matrix(score_components$scores)
  if (
    anyNA(score_ids) ||
      any(score_ids == "") ||
      anyDuplicated(score_ids) ||
      nrow(scores) != length(score_ids)
  ) {
    stop("Score-component IDs must be non-missing, unique, and row-aligned.")
  }
  if (
    length(score_components$n) != 1L ||
      !is.numeric(score_components$n) ||
      is.na(score_components$n) ||
      score_components$n != length(score_ids)
  ) {
    stop("`score_components$n` does not match its patient score rows.")
  }

  unmatched_scores <- setdiff(score_ids, profile_ids)
  unmatched_profiles <- setdiff(profile_ids, score_ids)
  if (length(unmatched_scores)) {
    stop(
      "Likelihood-score patients are absent from the target profiles: ",
      paste(utils::head(unmatched_scores, 5L), collapse = ", "),
      if (length(unmatched_scores) > 5L) " ..." else "",
      ". Fitted-cohort superpopulation inference requires every score patient to ",
      "appear exactly once in `profile_ids`.",
      call. = FALSE
    )
  }
  if (length(unmatched_profiles)) {
    stop(
      "Target-profile patients have no usable likelihood transition: ",
      paste(utils::head(unmatched_profiles, 5L), collapse = ", "),
      if (length(unmatched_profiles) > 5L) " ..." else "",
      ". Superpopulation inference includes only fitted patients and does not ",
      "insert zero patient-score rows.",
      call. = FALSE
    )
  }
  if (score_components$n != n) {
    stop(
      "The starting-profile and patient-score cohorts have different sizes.",
      call. = FALSE
    )
  }

  cell_names <- colnames(values)
  if (
    is.null(cell_names) ||
      anyNA(cell_names) ||
      any(cell_names == "") ||
      anyDuplicated(cell_names)
  ) {
    stop(
      "`individual_values` columns must have unique, non-missing cell names."
    )
  }
  jacobian_rows <- rownames(G)
  if (
    is.null(jacobian_rows) ||
      anyNA(jacobian_rows) ||
      any(jacobian_rows == "") ||
      anyDuplicated(jacobian_rows) ||
      !setequal(jacobian_rows, cell_names)
  ) {
    stop(
      "Jacobian row names must uniquely match all `individual_values` cells.",
      call. = FALSE
    )
  }
  if (nrow(G) != m) {
    stop("`average_jacobian` must have one row per estimand cell.")
  }
  G <- G[cell_names, , drop = FALSE]

  coefficient_names <- colnames(scores)
  if (
    is.null(coefficient_names) ||
      anyNA(coefficient_names) ||
      any(coefficient_names == "") ||
      anyDuplicated(coefficient_names)
  ) {
    stop("Score columns must have unique, non-missing coefficient names.")
  }
  jacobian_columns <- colnames(G)
  if (
    is.null(jacobian_columns) ||
      anyNA(jacobian_columns) ||
      any(jacobian_columns == "") ||
      anyDuplicated(jacobian_columns) ||
      !setequal(jacobian_columns, coefficient_names)
  ) {
    stop("Score and Jacobian coefficient names do not match.", call. = FALSE)
  }
  G <- G[, coefficient_names, drop = FALSE]
  score_Ainv <- delta_validate_named_matrix(
    score_components$Ainv,
    coefficient_names,
    "score_components$Ainv",
    positive_definite = TRUE
  )
  if (any(!is.finite(scores))) {
    stop("Patient score contributions contain non-finite values.")
  }

  score_order <- match(profile_ids, score_ids)
  aligned_scores <- scores[
    score_order,
    coefficient_names,
    drop = FALSE
  ]
  rownames(aligned_scores) <- profile_ids
  Ainv <- score_Ainv

  estimate <- colMeans(values)
  profile_component <- sweep(values, 2L, estimate, "-")
  coefficient_influence <- aligned_scores %*% t(Ainv) %*% t(G)
  influence <- profile_component + coefficient_influence
  dimnames(influence) <- list(profile_ids, cell_names)

  covariance <- stats::cov(influence) / n
  covariance <- (covariance + t(covariance)) / 2
  dimnames(covariance) <- list(cell_names, cell_names)
  profile_covariance <- stats::cov(values) / n
  dimnames(profile_covariance) <- list(cell_names, cell_names)

  list(
    influence = influence,
    covariance = covariance,
    std.error = sqrt(pmax(diag(covariance), 0)),
    estimate = estimate,
    profile_component = profile_component,
    profile_covariance = profile_covariance,
    coefficient_influence = coefficient_influence,
    Ainv = Ainv,
    scores = aligned_scores,
    ids = profile_ids,
    n = n,
    metadata = c(
      list(
        source = "stacked_patient_influence",
        covariance_source = "stacked_patient_influence",
        target = "superpopulation",
        cohort = "fitted_patients_with_origin_profiles",
        covariance_convention = "sample_covariance_of_influence_divided_by_n",
        finite_sample_method = "patient_level_sample_covariance",
        finite_sample_factor = n / (n - 1),
        n_patients = n,
        n_score_clusters = score_components$n
      ),
      score_components$metadata %||% list()
    )
  )
}

# Propagate a raw-coefficient covariance through an estimand Jacobian.
delta_empirical_covariance <- function(jacobian, covariance) {
  jacobian <- delta_validate_profile_matrix(jacobian, "jacobian")
  coefficient_names <- colnames(jacobian)
  if (
    is.null(coefficient_names) ||
      anyNA(coefficient_names) ||
      any(coefficient_names == "") ||
      anyDuplicated(coefficient_names)
  ) {
    stop("`jacobian` columns must have unique raw-coefficient names.")
  }
  covariance <- delta_validate_named_matrix(
    covariance,
    coefficient_names,
    "covariance",
    positive_semidefinite = TRUE
  )
  answer <- jacobian %*% covariance %*% t(jacobian)
  answer <- (answer + t(answer)) / 2
  row_names <- rownames(jacobian)
  dimnames(answer) <- list(row_names, row_names)
  answer
}
