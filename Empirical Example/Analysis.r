# 2026-10-05
# Proper ESS analysis reproducing the Pokropek paper
library("dplyr")
library("lavaan")

# Source necessary functions
invisible(lapply(list.files("hmm-mmgsem", pattern = "\\.R$", full.names = TRUE), source))

# # Y - Opposition to migration (reverse coded)
# imsmetn
# imdfetn
# impcntr
# 
# # X - Education years
# eduyrs
# 
# # M - Feelings of immigration threat (reverse coded)
# imbgeco
# imueclt
# imwbcnt

################################################################################
# Data -------------------------------------------------------------------------
################################################################################

# Load individual data ---------------------------------------------------------
data_ESS5 <- read.csv(file = "Data/Data_ESS5.csv") |> 
  filter(cntry %in% c("BE", "CH", "CZ", "DE",
                      "EE", "ES", "FI", "FR",
                      "GB", "HU", "IE", "IL",
                      "LT", "NL", "NO", "PL",
                      "PT", "SI")) |> 
  select(essround, idno, agea, cntry, eduyrs, brncntr,
         imsmetn, imdfetn, impcntr,
         imbgeco, imueclt, imwbcnt)

data_ESS6 <- read.csv(file = "Data/Data_ESS6.csv") |> 
  filter(cntry %in% c("BE", "CH", "CZ", "DE",
                      "EE", "ES", "FI", "FR",
                      "GB", "HU", "IE", "IL",
                      "LT", "NL", "NO", "PL",
                      "PT", "SI")) |> 
  select(essround, idno, agea, cntry, eduyrs, brncntr, 
         imsmetn, imdfetn, impcntr,
         imbgeco, imueclt, imwbcnt)

data_ESS7 <- read.csv(file = "Data/Data_ESS7.csv") |> 
  filter(cntry %in% c("BE", "CH", "CZ", "DE",
                      "EE", "ES", "FI", "FR",
                      "GB", "HU", "IE", "IL",
                      "LT", "NL", "NO", "PL",
                      "PT", "SI")) |> 
  select(essround, idno, agea, cntry, eduyrs, brncntr, 
         imsmetn, imdfetn, impcntr,
         imbgeco, imueclt, imwbcnt)

data_ESS8 <- read.csv(file = "Data/Data_ESS8.csv") |> 
  filter(cntry %in% c("BE", "CH", "CZ", "DE",
                      "EE", "ES", "FI", "FR",
                      "GB", "HU", "IE", "IL",
                      "LT", "NL", "NO", "PL",
                      "PT", "SI")) |> 
  select(essround, idno, agea, cntry, eduyrs, brncntr,
         imsmetn, imdfetn, impcntr,
         imbgeco, imueclt, imwbcnt)


# Bind into a single data set
data <- rbind(data_ESS5, data_ESS6, data_ESS7, data_ESS8)

# Remove NAs depending on the variable type
# Threat + education years
threat_cols <- c("imbgeco", "imueclt", "imwbcnt", "eduyrs")
data[threat_cols] <- lapply(data[threat_cols], function(x) {
  x[x %in% c(77, 88, 99)] <- NA
  return(x)
})
rm(threat_cols)

# Opposition to migration has different NA codes
# + Born in country
opposition_cols <- c("imsmetn", "imdfetn", "impcntr", "brncntr")
data[opposition_cols] <- lapply(data[opposition_cols], function(x) {
  x[x %in% c(7, 8, 9)] <- NA
  return(x)
})
rm(opposition_cols)

# Age has different NA codes
age <- c(as.matrix(data["agea"])) # Removed odd R class
age[age == 999] <- NA
data$agea <- age
rm(age)

# Follow the same exclusion criteria of the original paper
data <- data[data$agea >= 25, ]      # Remove people below 25 yo (as in the original study)
data <- data[data$brncntr == 1, ]    # Remove people not born in the country (as in the original study)
data <- na.omit(data)                # Remove data with missing values (as in the original study)

# Final sample of 98167 individuals. Slightly different from the original paper (101351). No apparent reason for the difference.

# Reverse code the items that were initially coded in the opposite direction
data <- data %>% mutate(
  imbgeco = 11 - imbgeco,
  imueclt = 11 - imueclt,
  imwbcnt = 11 - imwbcnt
)

# data$cntry_time <- paste0(data$cntry, "_", data$essround)
# data$eduyrs_z   <- scale(data$eduyrs)

################################################################################
# Estimate measurement models --------------------------------------------------
################################################################################

# Define the model
S1_model <- '
  # Opposition to migration
  mig_opp =~ imsmetn + imdfetn + impcntr
  
  # Feelings of threat
  threath =~ imbgeco + imueclt + imwbcnt
  
  # Education years
  edu_yrs =~ eduyrs
'

# Measurement invariace testing across country-time combinatios using custom function
mm_results <- compute_mi_long_stepwise_score(model           = S1_model, 
                                             data            = data, 
                                             group_var       = "cntry", 
                                             time_var        = "essround",
                                             method          = "test_statistic", 
                                             alpha           = 0.05, 
                                             d_cfi_threshold = 0.02)

# Check results
mm_results$final_syntax
mm_results$final_partial_constraints # Fully invariant, no additional non-invariances
mm_results$fit_measures
metric_fit <- mm_results$final_fit


################################################################################
# Estimate the LMMG-SEM --------------------------------------------------------
################################################################################
# Source all files in the hmm-mmgsem folder
# setwd("C:\\Users\\User\\OneDrive - KU Leuven\\0. Postdoc Leuven\\1. Papers\\Paper 2\\hmm-mmgsem")
# files.sources = list.files()
# sapply(files.sources, source)

# Run the LMMG-SEM
N_gs    <- lavaan::lavInspect(metric_fit, "nobs")
ngroups <- length(unique(data[data$essround == 5, ]$cntry))
ntimes  <- length(unique(data$essround))
group_names <- unique(data$cntry)

# Define the model
S2_model <- '
  mig_opp ~ threath + edu_yrs
  threath ~ edu_yrs
'

# # For debugging
# debug <- lmmgsem(nstates  = 1,
#                  cov_eta  = cov_eta,
#                  S2       = S2_model,
#                  ngroups  = ngroups,
#                  g_names  = group_names,
#                  ntimes   = ntimes,
#                  N_gs     = N_gs,
#                  nstarts  = 20,
#                  max_it   = 10000,
#                  s2_fast  = T,
#                  printing = F)

# Estimate the model
fits <- model_selection(nclus    = 1:8,
                        cov_eta  = mm_results$phi_matrices,
                        S2       = S2_model,
                        ngroups  = ngroups,
                        g_names  = group_names,
                        ntimes   = ntimes,
                        N_gs     = N_gs,
                        nstarts  = 100,
                        max_it   = 10000,
                        s2_fast  = T,
                        printing = F)

# Check scree plots
plot(fits$overview$K, fits$overview$loglik); lines(fits$overview$K, fits$overview$loglik)
plot(fits$overview$K, fits$overview$BIC_N);  lines(fits$overview$K, fits$overview$BIC_N)
plot(fits$overview$K, fits$overview$BIC_G);  lines(fits$overview$K, fits$overview$BIC_G)
plot(fits$overview$K, fits$overview$AIC_3);  lines(fits$overview$K, fits$overview$AIC_3)

# Check, by raw numbers, the selected model per measure
sapply(c("BIC_N", "BIC_G", "AIC_3"), FUN = \(x) fits$overview$K[which.min(fits$overview[[x]])])
fits$overview$K[which.max(fits$overview$CHull)]

# BICs  = 4
# CHull = 4
# AIC_3 = 7

# Extract the "best" model
best_fit_2 <- fits$model_fit[[2]]
best_fit_2$param$beta_ks
best_fit_2$A

best_fit_3 <- fits$model_fit[[3]]
best_fit_3$param$beta_ks

best_fit_4 <- fits$model_fit[[4]]
best_fit_4$param$beta_ks
round(best_fit_4$A, 3)

#############################################################################
# Plots --------------------------------------------------------------------
#############################################################################
# Source the function to generate the alluvial plot and the Europe map transitions
source("R/plotting_function.R")

# Extract country codes
country_codes <- unique(data$cntry)

# Generate plots
plots <- plot_lmm_results(lmm_fit          = best_fit_4, 
                          country_codes    = country_codes, 
                          custom_palette   = c("#2B5C8F", "#D95F02", "#1B9E77", "#7570B3"),
                          sem_model_syntax = S2_model)


plots$alluvial_plot
plots$map_plot
plots$path_plots
plots$path_plots_html

# Save plots
# --- Plot 1: Map (Portrait Full-Width) ---
ggsave(
  filename = "Figures/geographic_clusters.png",
  plot = plots$map_plot,
  dpi = 600,
  width = 165,
  height = 100,
  units = "mm"
)

# --- Plot 2: Alluvial Diagram (Portrait Full-Width) ---
ggsave(
  filename = "Figures/cluster_transitions.png",
  plot = plots$alluvial_plot,
  dpi = 600,
  width = 240,
  height = 150,
  units = "mm"
)


# lavaan::fitmeasures(metric_fit_marker)[c("cfi", "tli", "rmsea", "srmr")]
# cov_eta <- lapply(lavaan::lavInspect(metric_fit_marker, what = "est"), "[[", "psi")
