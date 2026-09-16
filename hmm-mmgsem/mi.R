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
                            invariance_levels = c("configural", "metric"),
                            fit_indices = c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
                            d_cfi_threshold = 0.020,
                            d_rmsea_threshold = 0.015,
                            ...) {
  
  # We will use the 'lavaan' and 'semTools' packages to automatize the MI testing process.
  # Do some variable preparations
  unique_times  <- unique(data[time_var])
  unique_groups <- unique(data[group_var])

  n_times  <- length(unique_times) 
  n_groups <- length(unique_groups) 

  flagged_params <- character(0)  # Initialize a vector to store flagged parameters for invariance violations

  # Helper to safely extract non-invariant loading strings out of lavaan fit objects
  extract_flagged <- function(fit) {
    if (!lavInspect(fit, "converged")) return(character(0))
    ts <- tryCatch(lavTestScore(fit), error = function(e) NULL)
    if (is.null(ts) || is.null(ts$uni)) return(character(0))
    
    uni_tests <- ts$uni
    corrected_alpha <- alpha / nrow(uni_tests)  # Bonferroni correction for multiple comparisons
    sig_tests <- uni_tests[uni_tests$p.value < corrected_alpha & !is.na(uni_tests$p.value), ]
    if (nrow(sig_tests) == 0) return(character(0))
    
    # Extract lhs, op, rhs to construct 'Factor =~ Item' syntax
    par_table <- parTable(fit)
    flagged <- character(0)
    for (plabel in sig_tests$lhs) {
      row_match <- parTable[parTable$plabel == plabel | parTable$label == plabel, ]
      if (nrow(row_match) > 0 && any(row_match$op == "=~")) {
        load_row <- row_match[row_match$op == "=~", ][1, ]
        flagged <- c(flagged, paste(load_row$lhs, load_row$op, load_row$rhs))
      }
    }
    return(unique(flagged))
  }

  # ----------------------------------------------------------------------------
  # Step 1: Measurement Invariance Across Groups at Each Time Point
  # ----------------------------------------------------------------------------
  # Set up a progress bar for Step 1
  message("Executing Step 1: Cross-Sectional Invariance per Time Point...")
  pb_1 <- progress_bar$new(
    format = "  Step 1 [:bar] :percent in :elapsed | Wave :current/:total",
    total = n_times, clear = FALSE, width = 70
  )

  cross_sectional_results <- vector(length = n_times, mode = "list")
  names(cross_sectional_results) <- paste0(Time:, unique_times)

  for(t in 1:n_times){
    # For progress bar
    pb_1$tick()

    # Subset the data with all groups and a single time point
    sub_data <- data[data[time_var] == t, ]

    # Generate model syntax using semTools
    config_syntax <- semTools::measEq.syntax(configural.model = model, 
                                             data             = sub_data, 
                                             group            = group_var)

    metric_syntax <- semTools::measEq.syntax(configural.model = model, 
                                             data             = sub_data, 
                                             group            = group_var,
                                             group.equal      = "loadings")

    # Fit the models
    config_fit <- lavaan::cfa(model = config_syntax)
    metric_fit <- lavaan::cfa(model = metric_syntax)

    config_fit_idxs <- lavaan::fitmeasures(config_fit)[fit_indices]
    metric_fit_idxs <- lavaan::fitmeasures(metric_fit)[fit_indices] 

    if(metric_fit_idxs["cfi"] < config_fit_idxs["cfi"] - d_cfi_threshold){
       flagged_params <- c(flagged_params, extract_flagged(metric_fit))
    }

    longitudinal_results[[g]] <- list(
        configural_fit <- config_fit,
        metric_fit     <- metric_fit
    )
  }

  # ----------------------------------------------------------------------------
  # Step 2: Measurement Invariance Across Time Points for Each Group
  # ----------------------------------------------------------------------------
  # Set up a progress bar for Step 2
  message("Executing Step 2: Longitudinal Invariance per Group...")
  pb_2 <- progress_bar$new(
    format = "  Step 2 [:bar] :percent in :elapsed | Wave :current/:total",
    total = n_groups, clear = FALSE, width = 70
  )

  longitudinal_results <- vector(length = n_groups, mode = "list")
  names(longitudinal_results) <- paste0(Group:, unique_group)

  for(g in 1:n_groups){
    # For progress bar
    pb_2$tick()

    # Subset the data with all groups and a single time point
    sub_data <- data[data[group_var] == g, ]

    # Generate model syntax using semTools
    config_syntax <- semTools::measEq.syntax(configural.model = model, 
                                             data             = sub_data, 
                                             group            = group_var)

    metric_syntax <- semTools::measEq.syntax(configural.model = model, 
                                             data             = sub_data, 
                                             group            = group_var,
                                             group.equal      = "loadings")

    # Fit the models
    config_fit <- lavaan::cfa(model = config_syntax)
    metric_fit <- lavaan::cfa(model = metric_syntax)

    config_fit_idxs <- lavaan::fitmeasures(config_fit)[fit_indices]
    metric_fit_idxs <- lavaan::fitmeasures(metric_fit)[fit_indices] 

    if(metric_fit_idxs["cfi"] < config_fit_idxs["cfi"] - d_cfi_threshold){
       flagged_params <- c(flagged_params, extract_flagged(metric_fit))
    }

    longitudinal_results[[g]] <- list(
        configural_fit <- config_fit,
        metric_fit     <- metric_fit
    )
  }

  # -------------------------------------------------------------
  # Step 3: Build final syntax
  # -------------------------------------------------------------
  # message("\nBuilding Final Model Syntax with Flagged Parameters for Partial Invariance...")
  
  # pb_3 <- progress_bar$new(
  #   format = "  Building Model Syntax [:bar] :percent | Elapsed: :elapsed",
  #   total = 1, clear = FALSE, width = 70
  # )

  # -------------------------------------------------------------
  # Final Step: Simultaneous Omnibus MG-CFA Across G x T Cells
  # -------------------------------------------------------------
  message("\nEstimating Final (Partial) Invariance Model across G x T cells...")
  data$group_time <- interaction(data[[group_var]], data[[time_var]], sep = ".Time", drop = TRUE)
  
  pb_4 <- progress_bar$new(
    format = "  Final Estimation [:bar] :percent | Elapsed: :elapsed",
    total = 1, clear = FALSE, width = 70
  )

  final_syntax <- semTools::measEq.syntax(configural.model = model, 
                                          data             = data, 
                                          group            = "group_time",
                                          group.equal      = c("loadings"),
                                          group.partial    = if (length(flagged_params) > 0) flagged_params else NULL) # Add flagged parameters (if there are any)

  pb_4$tick()

  final_fit <- lavaan::cfa(model = final_syntax, data = data, group = "group_time", ...)

  pb_4$tick()


  # Extract (in advance) the factor covariance matrices for Step 2 of LMMG-SEM
  phi_matrices <- lavaan::lavInspect(object = final_fit, what = "cov.lv")

  # ----------------------------------------------------------------------------
  # Model Comparisons & Invariance Decision Summary
  # ----------------------------------------------------------------------------
  # Extract relevant fit indices for all models
  # Empty matrix to store the fit indices
  step1_fits <- as.data.frame(matrix(data = NA, nrow = 2, ncol = length(fit_indices), dimnames = list(c("configural", "metric"), fit_indices)))
  step2_fits <- as.data.frame(matrix(data = NA, nrow = 2, ncol = length(fit_indices), dimnames = list(c("configural", "metric"), fit_indices)))

  fit_measures_list <- list(
    step1_across_groups = lapply(cross_sectional_results, function(x) {
      fit_table <- rbind(lavaan::fitmeasures(x$configural_fit)[fit_indices], lavaan::fitmeasures(x$metric_fit)[fit_indices])
    }),
    step2_across_time = lapply(longitudinal_results, function(x) {
      fit_table <- rbind(lavaan::fitmeasures(x$configural_fit)[fit_indices], lavaan::fitmeasures(x$metric_fit)[fit_indices])
    }),
    final_model = lavaan::fitmeasures(final_fit)[fit_indices]
  )

  return(
    list(
      step1_across_groups = cross_sectional_results,
      step2_across_time   = longitudinal_results,
      final_fit           = final_fit,
      phi_matrices        = phi_matrices,
      flagged_params      = unique(flagged_params),
      fit_measures        = fit_measures_list
    )
  )
}
