# ==============================================================================
# Two-Step Measurement Invariance for Repeated Cross-Sectional Designs
# ==============================================================================
#' @title Two-Step Measurement Invariance in Repeated Cross-Sectional Data
#' 
#' @description
#' Evaluates measurement invariance (MI) in repeated cross-sectional SEM designs 
#' through a two-step sequential testing procedure:
#' \enumerate{
#'   \item \strong{Step 1 (Cross-Sectional / Across Groups):} Tests measurement invariance 
#'         across groups separately at each time point \eqn{t \in \{1, \dots, T\}}.
#'   \item \strong{Step 2 (Longitudinal / Across Time):} Tests measurement invariance 
#'         across time points separately for each group \eqn{g \in \{1, \dots, G\}}.
#' }
#'
#' @details
#' In repeated cross-sectional multi-group SEM designs (e.g., European Social Survey), establishing comparability 
#' requires verifying that measurement parameters (e.g., factor loadings)
#' remain invariant across both group and temporal dimensions.
#'
#' Typical levels of measurement invariance evaluated at each step include:
#' \itemize{
#'   \item \strong{Configural Invariance:} Same factor structure / pattern of zero and non-zero loadings.
#'   \item \strong{Metric (Weak) Invariance:} Equal factor loadings (\eqn{\Lambda}).
#'   \item \strong{Scalar (Strong) Invariance:} Equal factor loadings (\eqn{\Lambda}) and item intercepts (\eqn{\nu} / \eqn{\tau}).
#'   \item \strong{Strict (Residual) Invariance:} Equal factor loadings, item intercepts, and residual variances (\eqn{\Theta}).
#' }
#'
#' WE WILL FOCUS ONLY ON CONFIGURAL AND METRIC INVARIANCE (FOR NOW).
#'
#' @param data A `data.frame` containing the observed data, a group identifier, and a time identifier.
#' @param model A character string specifying the measurement model syntax (i.e., `lavaan` model syntax).
#' @param group Character. Name of the column denoting group membership (e.g., country, cluster).
#' @param time Character. Name of the column denoting time point / wave.
#' @param invariance_levels Character vector of MI levels to test sequentially. Default is `c("configural", "metric")`.
#' @param fit_measures Character vector of fit measures to extract for model comparison (e.g., `c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr", "aic", "bic")`).
#' @param d_cfi_threshold Numeric. Cutoff value for \eqn{\Delta \text{CFI}} indicating violation of invariance (e.g., `0.010` following Cheung & Rensvold, 2002).
#' @param d_rmsea_threshold Numeric. Cutoff value for \eqn{\Delta \text{RMSEA}} (e.g., `0.015` following Chen, 2007).
#' @param ... Additional arguments passed to the underlying SEM estimation engine (e.g., `lavaan::cfa` or `lavaan::sem`).
#'
#' @return A list of class `"mi_rcs"` containing:
#' \item{step1_across_groups}{A list of invariance testing results (fit tables, model comparisons, parameter estimates) across groups for each time point \eqn{t}.}
#' \item{step2_across_time}{A list of invariance testing results across time points for each group \eqn{g}.}
#' \item{summary_table}{A consolidated summary table detailing fit indices, \eqn{\Delta} fit metrics, and invariance conclusions across both steps.}
#' \item{model_syntax}{The specified measurement model syntax.}
#' \item{call}{The matched function call.}
#'
#' @references
#' \itemize{
#'   \item Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes for testing measurement invariance. \emph{Structural Equation Modeling}, 9(2), 233-255.
#'   \item Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of measurement invariance. \emph{Structural Equation Modeling}, 14(3), 464-504.
#'   \item Davidov, E., Meuleman, B., Cieciuch, J., Schmidt, P., & Billiet, J. (2014). Measurement equivalence in cross-national research. \emph{Annual Review of Sociology}, 40, 55-75.
#' }
#'
#' @examples
#' \dontrun{
#' # Example measurement syntax
#' cfa_model <- '
#'   Eta1 =~ y1 + y2 + y3
#'   Eta2 =~ y4 + y5 + y6
#' '
#'
#' # Run two-step measurement invariance
#' mi_results <- compute_mi_long(
#'   data        = ess_data,
#'   model       = cfa_model,
#'   group       = "country",
#'   time        = "wave",
#'   estimator   = "MLR"
#' )
#'
#' # Print summary
#' summary(mi_results)
#' }
#'
#' @export
compute_mi_long <- function(data,
                            model,
                            group_var,
                            time_var,
                            alpha = 0.05,
                            fit_indices = c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                            d_cfi_threshold = 0.010,
                            method = c("score", "alignment"),
                            specific_noninv = FALSE,
                            ...) {
  
  method <- match.arg(method)
  
  # Setup interaction group_time for the final model (G x T cells)
  data$group_time <- interaction(data[[group_var]], data[[time_var]], sep = ".Time", drop = TRUE)
  
  unique_times  <- unique(as.character(data[[time_var]]))
  unique_groups <- unique(as.character(data[[group_var]]))
  all_cells     <- levels(data$group_time)
  
  n_times  <- length(unique_times)
  n_groups <- length(unique_groups)
  n_cells  <- length(all_cells)
  
  flagged_global <- character(0)
  
  # Track pairwise/cell-specific violations: lhs, rhs, cell1, cell2
  flagged_pairs <- data.frame(
    lhs   = character(0),
    rhs   = character(0),
    cell1 = character(0),
    cell2 = character(0),
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------------------------
  # Helper: Extract Exact Pairwise Non-Invariances from lavTestScore
  # ----------------------------------------------------------------------------
  extract_score_pairs <- function(fit, fixed_val, mode = c("wave", "group"), alpha = 0.05) {
    if (!lavInspect(fit, "converged")) return(list(global = character(0), pairs = NULL))
    
    ts <- tryCatch(lavTestScore(fit), error = function(e) NULL)
    if (is.null(ts) || is.null(ts$uni)) return(list(global = character(0), pairs = NULL))
    
    uni_tab <- ts$uni
    # Look only at equality constraints (op == "==")
    uni_tab <- uni_tab[uni_tab$op == "==", ]
    if (nrow(uni_tab) == 0) return(list(global = character(0), pairs = NULL))
    
    # Bonferroni correction
    corrected_alpha <- alpha / nrow(uni_tab)
    sig_tests <- uni_tab[!is.na(uni_tab$p.value) & uni_tab$p.value < corrected_alpha, ]
    if (nrow(sig_tests) == 0) return(list(global = character(0), pairs = NULL))
    # browser()
    pt <- parTable(fit)
    grp_labels <- lavInspect(fit, "group.label")
    
    glob_out <- character(0)
    pair_out <- list()
    
    for (i in seq_len(nrow(sig_tests))) {
      p1 <- sig_tests$lhs[i]
      p2 <- sig_tests$rhs[i]
      
      row1 <- pt[pt$plabel == p1 | pt$label == p1, ]
      row2 <- pt[pt$plabel == p2 | pt$label == p2, ]
      
      # We only care about factor loading constraints (=~)
      if (nrow(row1) > 0 && row1$op[1] == "=~") {
        l_fac  <- row1$lhs[1]
        l_item <- row1$rhs[1]
        glob_out <- c(glob_out, paste(l_fac, "=~", l_item))
        
        g1_name <- grp_labels[row1$group[1]]
        g2_name <- grp_labels[row2$group[1]]
        
        # Build group_time cell identifier (Group.TimeX)
        if (mode == "wave") { # It will results in something like group1.Time1 vs group2.Time1 (always same time)
          c1 <- paste0(g1_name, ".Time", fixed_val)
          c2 <- paste0(g2_name, ".Time", fixed_val)
        } else { # Result: group1.Time1 vs group1.Time2 (always same group; g_name becomes the time point in lavaan partable)
          c1 <- paste0(fixed_val, ".Time", g1_name)
          c2 <- paste0(fixed_val, ".Time", g2_name)
        }
        
        pair_out[[length(pair_out) + 1]] <- data.frame(
          lhs = l_fac, rhs = l_item, cell1 = c1, cell2 = c2, stringsAsFactors = FALSE
        )
      }
    }
    
    p_df <- if (length(pair_out) > 0) do.call(rbind, pair_out) else NULL
    return(list(global = unique(glob_out), pairs = p_df))
  }
  
  # ============================================================================
  # METHOD: SCORE TEST (Two-Step Across Waves & Groups)
  # ============================================================================
  if (method == "score") {
    
    # --------------------------------------------------------------------------
    # Step 1: Across Groups at Each Time Point
    # --------------------------------------------------------------------------
    message("Executing Step 1: Cross-Sectional Invariance per Time Point...")
    pb_1 <- txtProgressBar(min = 0, max = n_times, style = 3)
    cross_sectional_results <- vector(length = n_times, mode = "list")
    names(cross_sectional_results) <- paste0("Time_", unique_times)
    
    for (t in seq_len(n_times)) {
      sub_data <- data[data[[time_var]] == unique_times[t], ]
      
      config_syn <- semTools::measEq.syntax(configural.model = model, data = sub_data, group = group_var, std.lv = TRUE)
      metric_syn <- semTools::measEq.syntax(configural.model = model, data = sub_data, group = group_var, group.equal = "loadings", std.lv = TRUE)
      
      config_fit <- lavaan::cfa(model = as.character(config_syn), data = sub_data, group = group_var, ...)
      metric_fit <- lavaan::cfa(model = as.character(metric_syn), data = sub_data, group = group_var, ...)
      
      config_cfi <- lavaan::fitmeasures(config_fit, "cfi")
      metric_cfi <- lavaan::fitmeasures(metric_fit, "cfi")
      
      if (!is.na(metric_cfi) && !is.na(config_cfi) && (metric_cfi < config_cfi - d_cfi_threshold)) {
        res <- extract_score_pairs(metric_fit, fixed_val = unique_times[t], mode = "wave")
        flagged_global <- c(flagged_global, res$global)
        if (!is.null(res$pairs)) flagged_pairs <- rbind(flagged_pairs, res$pairs)
      }
      
      cross_sectional_results[[t]] <- list(configural_fit = config_fit, metric_fit = metric_fit)
      setTxtProgressBar(pb_1, t)
    }
    close(pb_1)
    
    # --------------------------------------------------------------------------
    # Step 2: Across Time Points for Each Group
    # --------------------------------------------------------------------------
    message("Executing Step 2: Longitudinal Invariance per Group...")
    pb_2 <- txtProgressBar(min = 0, max = n_groups, style = 3)
    longitudinal_results <- vector(length = n_groups, mode = "list")
    names(longitudinal_results) <- paste0("Group_", unique_groups)
    
    for (g in seq_len(n_groups)) {
      sub_data <- data[data[[group_var]] == unique_groups[g], ]
      
      config_syn <- semTools::measEq.syntax(configural.model = model, data = sub_data, group = time_var, std.lv = TRUE)
      metric_syn <- semTools::measEq.syntax(configural.model = model, data = sub_data, group = time_var, group.equal = "loadings", std.lv = TRUE)
      
      config_fit <- lavaan::cfa(model = as.character(config_syn), data = sub_data, group = time_var, ...)
      metric_fit <- lavaan::cfa(model = as.character(metric_syn), data = sub_data, group = time_var, ...)
      
      config_cfi <- lavaan::fitmeasures(config_fit, "cfi")
      metric_cfi <- lavaan::fitmeasures(metric_fit, "cfi")
      
      if (!is.na(metric_cfi) && !is.na(config_cfi) && (metric_cfi < config_cfi - d_cfi_threshold)) {
        res <- extract_score_pairs(metric_fit, fixed_val = unique_groups[g], mode = "group")
        flagged_global <- c(flagged_global, res$global)
        if (!is.null(res$pairs)) flagged_pairs <- rbind(flagged_pairs, res$pairs)
      }
      
      longitudinal_results[[g]] <- list(configural_fit = config_fit, metric_fit = metric_fit)
      setTxtProgressBar(pb_2, g)
    }
    close(pb_2)
  }
  
  # ============================================================================
  # METHOD: ALIGNMENT OPTIMIZATION (Alternative)
  # ============================================================================
  if (method == "alignment") {
    browser()
    if (!requireNamespace("sirt", quietly = TRUE)) {
      stop("Package 'sirt' must be installed to use method = 'alignment'.")
    }
    message("Running Alignment Optimization across all G x T cells...")
    
    fit_config_all <- lavaan::cfa(model, data = data, group = "group_time", std.lv = TRUE, ...)
    pars <- lavInspect(fit_config_all, "est")
    parsed_m <- lavaanify(model)
    load_pairs <- parsed_m[parsed_m$op == "=~", c("lhs", "rhs")]
    items <- unique(load_pairs$rhs)
    
    lambda_mat <- matrix(NA, nrow = n_cells, ncol = length(items), dimnames = list(all_cells, items))
    nu_mat     <- matrix(0,  nrow = n_cells, ncol = length(items), dimnames = list(all_cells, items))
    
    for (cl in all_cells) { # Go over every group-time combination
      cl_lambda <- pars[[cl]]$lambda
      for (itm in items) {
        fac <- load_pairs$lhs[load_pairs$rhs == itm]
        lambda_mat[cl, itm] <- cl_lambda[itm, fac]
      }
    }
    
    align_res <- sirt::invariance.alignment(lambda = lambda_mat, nu = nu_mat)
    noninv_loadings <- align_res$par_invariance$loadings
    
    for (itm in items) {
      fac <- load_pairs$lhs[load_pairs$rhs == itm]
      bad_cells <- rownames(noninv_loadings)[noninv_loadings[, itm] == 1]
      if (length(bad_cells) > 0) {
        flagged_global <- c(flagged_global, paste(fac, "=~", itm))
        for (bc in bad_cells) {
          flagged_pairs <- rbind(flagged_pairs, data.frame(
            lhs = fac, rhs = itm, cell1 = bc, cell2 = bc, stringsAsFactors = FALSE
          ))
        }
      }
    }
    cross_sectional_results <- list()
    longitudinal_results    <- list()
  }
  
  flagged_global <- unique(flagged_global)
  if (nrow(flagged_pairs) > 0) flagged_pairs <- unique(flagged_pairs)
  
  # ============================================================================
  # Step 3: Build Final Partial Invariance Syntax
  # ============================================================================
  message("\nEstimating Final Omnibus Model across G x T cells...")
  
  if (!specific_noninv) {
    # Strategy 1: Global Freeing Across All Cells
    final_syntax <- semTools::measEq.syntax(
      configural.model = model,
      data             = data,
      group            = "group_time",
      group.equal      = "loadings",
      std.lv           = TRUE,
      group.partial    = if (length(flagged_global) > 0) flagged_global else NULL
    )
    final_fit <- lavaan::cfa(model = as.character(final_syntax), data = data, group = "group_time", ...)
    
  } else {
    # Strategy 2: Specific Non-Invariance Syntax with c() Labels
    parsed_model <- lavaanify(model)
    load_rows    <- parsed_model[parsed_model$op == "=~", ]
    syntax_lines <- character(0)
    
    for (fac in unique(load_rows$lhs)) {
      fac_items <- load_rows$rhs[load_rows$lhs == fac] # Select the items of the corresponding factor
      fac_elements <- character(0)
      
      for (i in seq_along(fac_items)) {
        itm <- fac_items[i]
        
        # First item fixed to 1 for scale identification if not std.lv
        if (i == 1) {
          fac_elements <- c(fac_elements, itm)
          next
        }
        
        # Check if item has flagged pairwise violations
        item_flags <- flagged_pairs[flagged_pairs$lhs == fac & flagged_pairs$rhs == itm, ]
        
        lbl_vec <- character(n_cells)
        for (c_idx in seq_along(all_cells)) {
          c_name <- all_cells[c_idx]
          
          # Is this cell implicated in any failed pair comparison?
          in_violation <- any(item_flags$cell1 == c_name | item_flags$cell2 == c_name)
          
          if (in_violation) {
            # Give this specific cell a unique parameter label
            lbl_vec[c_idx] <- paste0("l_", fac, "_", itm, "_", gsub("[^[:alnum:]]", "_", c_name))
          } else {
            # Invariant cells retain the shared common label
            lbl_vec[c_idx] <- paste0("l_", fac, "_", itm, "_common")
          }
        }
        
        lbl_string <- paste0("c(", paste(lbl_vec, collapse = ", "), ")")
        fac_elements <- c(fac_elements, paste0(lbl_string, "*", itm))
      }
      syntax_lines <- c(syntax_lines, paste0(fac, " =~ ", paste(fac_elements, collapse = " + ")))
    }
    
    final_syntax <- paste(syntax_lines, collapse = "\n")
    final_fit <- lavaan::cfa(model = final_syntax, data = data, group = "group_time", ...)
  }
  
  phi_matrices <- lavaan::lavInspect(final_fit, "cov.lv")
  
  # ----------------------------------------------------------------------------
  # Fit Measures Summary
  # ----------------------------------------------------------------------------
  fit_measures_list <- list(
    step1_across_groups = lapply(cross_sectional_results, function(x) {
      rbind(configural = lavaan::fitmeasures(x$configural_fit)[fit_indices],
            metric     = lavaan::fitmeasures(x$metric_fit)[fit_indices])
    }),
    step2_across_time = lapply(longitudinal_results, function(x) {
      rbind(configural = lavaan::fitmeasures(x$configural_fit)[fit_indices],
            metric     = lavaan::fitmeasures(x$metric_fit)[fit_indices])
    }),
    final_model = lavaan::fitmeasures(final_fit)[fit_indices]
  )
  
  return(list(
    step1_across_groups = cross_sectional_results,
    step2_across_time   = longitudinal_results,
    final_fit           = final_fit,
    phi_matrices        = phi_matrices,
    flagged_global      = flagged_global,
    flagged_pairs       = flagged_pairs,
    fit_measures        = fit_measures_list,
    final_syntax        = if (specific_noninv) final_syntax else as.character(final_syntax)
  ))
}