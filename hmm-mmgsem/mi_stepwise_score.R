# ============================================================================== 
# Stepwise Score-Based Measurement Invariance for Repeated Cross-Sectional Data
# ==============================================================================
#
# This version implements the classic sequential strategy using lavTestScore:
#   1) Identify equality constraints that are most violated under metric invariance
#   2) Free the strongest one (largest score statistic / smallest p-value)
#   3) Refit the model
#   4) Re-check fit
#   5) Repeat until the metric model is acceptable
#
# This is done separately for each time point (across groups) and each group
# (across time points), and the resulting non-invariances are gathered into the
# final omnibus model across all G x T cells.
#
# NOTE: This implementation follows the general/global strategy, i.e., we use
# group.partial = c("factor =~ item", ...), not cell-specific labels.
# ==============================================================================

compute_mi_long_stepwise_score <- function(data,
                                           model,
                                           group_var,
                                           time_var,
                                           alpha = 0.05,
                                           d_cfi_threshold = 0.010,
                                           max_iter = 20,
                                           fit_indices = c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                                           ...) {

  # Ensure the grouping/time variables are character-valued for stable ordering
  data <- data[!is.na(data[[group_var]]) & !is.na(data[[time_var]]), , drop = FALSE]
  data$group_time <- interaction(data[[group_var]], data[[time_var]], sep = ".Time", drop = TRUE)

  unique_times  <- unique(as.character(data[[time_var]]))
  unique_groups <- unique(as.character(data[[group_var]]))
  all_cells     <- levels(data$group_time)

  if (length(unique_times) == 0 || length(unique_groups) == 0) {
    stop("No valid groups or time points found in the supplied data.")
  }

  # ---------------------------------------------------------------------------
  # Helper functions
  # ---------------------------------------------------------------------------
  build_metric_model <- function(sub_data,
                                 split_var,
                                 partial_constraints = character(0),
                                 ...) {
    metric_syn <- semTools::measEq.syntax(
      configural.model = model,
      data = sub_data,
      group = split_var,
      group.equal = "loadings",
      group.partial = if (length(partial_constraints) > 0) partial_constraints else NULL,
      std.lv = TRUE
    )

    metric_fit <- lavaan::cfa(
      model = as.character(metric_syn),
      data = sub_data,
      group = split_var,
      ...
    )

    return(metric_fit)
  }

  # Extract the most violated invariance constraint based on how often it is
  # implicated across the group-vs-reference comparisons evaluated by lavTestScore.
  candidate_from_score <- function(fit, alpha = 0.05) {
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      return(NULL)
    }

    ts <- tryCatch(lavTestScore(fit), error = function(e) NULL)
    if (is.null(ts) || is.null(ts$uni)) {
      return(NULL)
    }

    uni_tab <- ts$uni                                      # Extract the univariate tests
    uni_tab <- uni_tab[uni_tab$op == "==", , drop = FALSE] # Extract only the equality constraints
    if (nrow(uni_tab) == 0) {
      return(NULL)
    }

    uni_tab <- uni_tab[!is.na(uni_tab$p.value), , drop = FALSE]
    if (nrow(uni_tab) == 0) {
      return(NULL)
    }

    # Bonferroni-adjusted significance threshold for deciding whether a given
    # comparison is a violation.
    alpha_adj <- alpha / nrow(uni_tab)
    uni_tab$constraint_key <- vapply(
      seq_len(nrow(uni_tab)),
      function(i) constraint_key_from_score_row(fit, uni_tab$lhs[i], uni_tab$rhs[i]),
      character(1)
    )
    uni_tab$violated <- as.logical(uni_tab$p.value < alpha_adj)

    if (!any(uni_tab$violated)) {
      return(NULL)
    }

    # Primary ordering criterion: how many times the same loading constraint is
    # flagged across the group comparisons. Tie-break using the largest score
    # statistic, then smallest p-value.
    agg <- aggregate(
      violated ~ constraint_key,
      data = uni_tab,
      FUN = sum,
      na.action = na.omit
    )
    names(agg)[2] <- "n_violations"

    agg <- agg[order(-agg$n_violations), , drop = FALSE]
    chosen_key <- agg$constraint_key[1]

    candidate_rows <- uni_tab[uni_tab$constraint_key == chosen_key, , drop = FALSE]

    return(candidate_rows[1, , drop = FALSE])
  }

  parameter_name_from_constraint <- function(fit, lhs_label, rhs_label) {
    pt <- lavaan::parTable(fit)

    # lavTestScore equality rows are often parameter labels, so we map back to
    # a loading expression such as "eta1 =~ y1" using the partable.
    row1 <- pt[pt$plabel == lhs_label | pt$label == lhs_label, , drop = FALSE]
    row2 <- pt[pt$plabel == rhs_label | pt$label == rhs_label, , drop = FALSE] # In case we add the group-specific version later

    if (nrow(row1) == 0 || nrow(row2) == 0) {
      return(NULL)
    }

    # Keep the row that corresponds to a loading constraint
    idx <- which(row1$op == "=~")
    if (length(idx) == 0) {
      return(NULL)
    }

    r1 <- row1[idx[1], ]
    return(paste(r1$lhs, "=~", r1$rhs))
  }

  constraint_key_from_score_row <- function(fit, lhs_label, rhs_label) {
    # lavTestScore often returns group-specific labels (e.g., a parameter label that
    # differs by group) even when they represent the same underlying loading.
    # We want the stable semantic key: the item loading itself, not the group-specific
    # label used in that contrast.
    pt <- lavaan::parTable(fit)

    row1 <- pt[pt$label == lhs_label | pt$plabel == lhs_label, , drop = FALSE]
    row2 <- pt[pt$label == rhs_label | pt$plabel == rhs_label, , drop = FALSE]

    if (nrow(row1) == 0 || nrow(row2) == 0) {
      return(NULL)
    }

    load_rows <- rbind(row1[row1$op == "=~", , drop = FALSE], row2[row2$op == "=~", , drop = FALSE])
    if (nrow(load_rows) == 0) {
      return(NULL)
    }

    r1 <- load_rows[1, ]
    paste(r1$lhs, "=~", r1$rhs)
  }

  stepwise_slice <- function(sub_data,
                            split_var,
                            label,
                            ...) {
    config_syn <- semTools::measEq.syntax(
      configural.model = model,
      data = sub_data,
      group = split_var,
      std.lv = TRUE
    )
    config_fit <- lavaan::cfa(
      model = as.character(config_syn),
      data = sub_data,
      group = split_var,
      ...
    )

    config_cfi <- lavaan::fitMeasures(config_fit, "cfi")
    if (is.null(config_cfi) || is.na(config_cfi)) {
      config_cfi <- NA_real_
    }

    partial_constraints <- character(0)
    history <- list()
    current_fit <- build_metric_model(sub_data, split_var, partial_constraints, ...)

    iter <- 0
    delta_cfi <- Inf

    while (iter < max_iter && (is.na(delta_cfi) || delta_cfi > d_cfi_threshold)) {
      iter <- iter + 1

      metric_cfi <- lavaan::fitMeasures(current_fit, "cfi")
      if (is.null(metric_cfi) || is.na(metric_cfi)) {
        break
      }

      candidate <- candidate_from_score(current_fit, alpha = alpha)
      if (is.null(candidate)) {
        break
      }

      param_name <- parameter_name_from_constraint(current_fit, candidate$lhs, candidate$rhs)
      if (is.null(param_name) || param_name %in% partial_constraints) {
        break
      }

      delta_cfi <- if (is.na(config_cfi)) NA_real_ else config_cfi - metric_cfi

      partial_constraints <- c(partial_constraints, param_name)
      history[[length(history) + 1]] <- list(
        iter = iter,
        parameter = param_name,
        stat = candidate$stat,
        p_value = candidate$p.value
      )

      current_fit <- build_metric_model(sub_data, split_var, partial_constraints, ...)
      metric_cfi <- lavaan::fitMeasures(current_fit, "cfi")
      delta_cfi <- if (is.na(config_cfi)) NA_real_ else config_cfi - metric_cfi
    }

    final_metric_cfi <- lavaan::fitMeasures(current_fit, "cfi")
    final_delta_cfi <- if (is.na(config_cfi) || is.na(final_metric_cfi)) NA_real_ else config_cfi - final_metric_cfi

    return(list(
      label = label,
      config_fit = config_fit,
      metric_fit = current_fit,
      partial_constraints = partial_constraints,
      history = history,
      config_cfi = config_cfi,
      metric_cfi = final_metric_cfi,
      delta_cfi = final_delta_cfi,
      acceptable = if (is.na(final_delta_cfi)) FALSE else final_delta_cfi <= d_cfi_threshold
    ))
  }

  # ---------------------------------------------------------------------------
  # Step 1: Across groups within each time point
  # ---------------------------------------------------------------------------
  message("Executing Step 1: cross-sectional metric invariance, freeing the strongest violations sequentially...")
  cross_sectional_results <- vector("list", length(unique_times))
  names(cross_sectional_results) <- paste0("Time_", unique_times)

  pb_1 <- txtProgressBar(min = 0, max = length(unique_times), style = 3)
  for (t in seq_along(unique_times)) {
    sub_data <- data[data[[time_var]] == unique_times[t], , drop = FALSE]
    result <- stepwise_slice(sub_data, split_var = group_var, label = unique_times[t], ...)
    cross_sectional_results[[t]] <- result
    setTxtProgressBar(pb_1, t)
  }
  close(pb_1)

  # ---------------------------------------------------------------------------
  # Step 2: Across time points within each group
  # ---------------------------------------------------------------------------
  message("Executing Step 2: longitudinal metric invariance, freeing the strongest violations sequentially...")
  longitudinal_results <- vector("list", length(unique_groups))
  names(longitudinal_results) <- paste0("Group_", unique_groups)

  pb_2 <- txtProgressBar(min = 0, max = length(unique_groups), style = 3)
  for (g in seq_along(unique_groups)) {
    sub_data <- data[data[[group_var]] == unique_groups[g], , drop = FALSE]
    result <- stepwise_slice(sub_data, split_var = time_var, label = unique_groups[g], ...)
    longitudinal_results[[g]] <- result
    setTxtProgressBar(pb_2, g)
  }
  close(pb_2)

  # ---------------------------------------------------------------------------
  # Collect the final partial invariance set (global/general strategy)
  # ---------------------------------------------------------------------------
  final_partials <- character(0)

  for (res in cross_sectional_results) {
    final_partials <- c(final_partials, res$partial_constraints)
  }
  for (res in longitudinal_results) {
    final_partials <- c(final_partials, res$partial_constraints)
  }

  final_partials <- unique(final_partials)

  # ---------------------------------------------------------------------------
  # Final omnibus model across all G x T cells
  # ---------------------------------------------------------------------------
  message("Estimating final omnibus model across all group-by-time cells...")
#   final_syntax <- semTools::measEq.syntax(
#     configural.model = model,
#     data = data,
#     group = "group_time",
#     group.equal = "loadings",
#     group.partial = if (length(final_partials) > 0) final_partials else NULL,
#     std.lv = FALSE
#   )

#   final_fit <- lavaan::cfa(
#     model = as.character(final_syntax),
#     data = data,
#     group = "group_time",
#     ...
#   )

  final_fit <- lavaan::cfa(
    model = model,
    data = data,
    group = "group_time",
    group.equal = "loadings",
    group.partial = if (length(final_partials) > 0) final_partials else NULL,
    std.lv = FALSE
  )

  phi_matrices <- lavaan::lavInspect(final_fit, "cov.lv")

  fit_measures_list <- list(
    step1_across_groups = lapply(cross_sectional_results, function(x) {
      rbind(
        configural = lavaan::fitMeasures(x$config_fit, fit_indices),
        metric = lavaan::fitMeasures(x$metric_fit, fit_indices)
      )
    }),
    step2_across_time = lapply(longitudinal_results, function(x) {
      rbind(
        configural = lavaan::fitMeasures(x$config_fit, fit_indices),
        metric = lavaan::fitMeasures(x$metric_fit, fit_indices)
      )
    }),
    final_model = lavaan::fitMeasures(final_fit, fit_indices)
  )

  return(list(
    step1_across_groups = cross_sectional_results,
    step2_across_time = longitudinal_results,
    final_fit = final_fit,
    phi_matrices = phi_matrices,
    final_partial_constraints = final_partials,
    fit_measures = fit_measures_list #,
    # final_syntax = as.character(final_syntax)
  ))
}
