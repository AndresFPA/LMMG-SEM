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





# # Install required packages
# install.packages("ggalluvial")
# install.packages("tidyverse")
# 
# library(ggalluvial)
# library(tidyverse)
# 
# 
# # Alluvial plot ---------------------------------------------------------------
# # Make an alluvial plot to showcase the cluster changing across time points
# # Create a clustering wide data set that works for the alluvial plotting
# posteriors <- round(best_fit_4$post, 1)
# 
# # List per time point with the cluster to which each country belongs to
# cluster_assignments <- lapply(1:ntimes, function(t) {
#   apply(posteriors[, , t], 1, which.max)
# })
# names(cluster_assignments) <- paste0("Time_", 1:ntimes)
# cluster_assignments
# 
# alluvial_data           <- matrix(NA, nrow = ngroups, ncol = ntimes + 1) %>% as.data.frame()
# colnames(alluvial_data) <- c("Country", paste0("Time_", 1:ntimes))
# alluvial_data$Country   <- unique(data$cntry)
# 
# for(col in 2:(ntimes + 1)){
#   alluvial_data[, col] <- as.factor(cluster_assignments[[col - 1]])
# }
# 
# # Create alluvial plot
# # First, reshape to lodes (long) format explicitly to silence the warning
# data_plot <- to_lodes_form(
#   data.frame(alluvial_data),
#   axes = c("Time_1", "Time_2", "Time_3", "Time_4"),
#   id = "Country_id",
#   key = "Wave",
#   value = "Cluster",
#   discern = FALSE
# )
# 
# # Add baseline state (Time_1) for each country to color the streams
# data_plot <- data_plot %>%
#   group_by(Country) %>%
#   mutate(Baseline_Cluster = Cluster[Wave == "Time_1"]) %>%
#   ungroup()
# 
# # Extract country label positions at Time 1
# country_labels_t1 <- data_plot %>%
#   filter(Wave == "Time_1") %>%
#   arrange(Cluster, Country) %>%
#   mutate(y_pos = ngroups - (row_number() - 1) - 0.5)
# 
# # Some specific colors for each cluster
# cluster_colors <- c("#2B5C8F", "#D95F02", "#1B9E77", "#7570B3")
# 
# p_alluvial <- ggplot(data_plot, aes(x = Wave, stratum = Cluster, alluvium = Country, y = 1)) +
#   geom_alluvium(aes(fill = Baseline_Cluster), width = 0.22, alpha = 0.72, knot.pos = 0.35, color = "white", linewidth = 0.35) +
#   geom_text(stat = "alluvium", aes(label = Country), size = 3, fontface = "plain", color = "white", check_overlap = FALSE, nudge_x = 0.16) +
#   # geom_text(data = ~ filter(.x, Wave %in% c("Time_2", "Time_3")), stat = "flow", aes(label = Country), size = 3, fontface = "plain", color = "white", check_overlap = FALSE, nudge_x = 0.16) +
#   # geom_text(data = ~ filter(.x, Wave %in% c("Time_3", "Time_4")), stat = "flow", aes(label = Country), size = 3, fontface = "plain", color = "white", check_overlap = FALSE, nudge_x = 0.16) +
#   geom_stratum(width = 0.22, fill = "white", color = "grey35", linewidth = 0.6) +
#   geom_text(stat = "stratum", aes(label = paste0("Cluster ", after_stat(stratum))), size = 3, fontface = "bold", color = "grey20") +
#   # geom_text(data = country_labels_t1, aes(x = 1 + (0.22 / 2) + 0.05, y = y_pos, label = Country), size = 2.8, fontface = "bold", inherit.aes = FALSE, color = "white") +
#   scale_x_discrete(
#     limits = c("Time_1", "Time_2", "Time_3", "Time_4"),
#     labels = c("Time 1", "Time 2", "Time 3", "Time 4"),
#     expand = c(0.12, 0.08)
#   ) +
#   scale_fill_manual(
#     values = cluster_colors,
#     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
#     name = "Baseline Cluster (T1)"
#   ) +
#   labs(
#     title = "Cluster Transitions (G = 18 Countries)",
#     # subtitle = "Alluvial flows colored by initial baseline state at Time 1",
#     x = "Wave",
#     y = "Number of Countries"
#   ) +
#   theme_minimal(base_size = 11) +
#   theme(
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom",
#     plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 13)),
#     axis.text.y = element_text(color = "grey30")
#   )
# 
# # Open a fixed-size graphics window before printing
# dev.new(width = 8, height = 5, noRStudioGD = TRUE)
# print(p_alluvial)
# 
# # Save the plot
# ggsave(
#   filename = "Figures/cluster_transitions.png",
#   plot = p_alluvial,
#   width = 3200,
#   height = 2300,
#   units = "px",
#   dpi = 300
# )
# 
# # Geographical plot ---------------------------------------------------------------
# # Required packages
# install.packages("sf")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# 
# library("sf")
# library("rnaturalearth")
# library("rnaturalearthdata")
# 
# # 1. Load world map data as sf object
# world <- ne_countries(scale = "medium", returnclass = "sf")
# 
# # 2. Prepare long-format cluster data (using your exact df)
# alluvial_data
# 
# df_map_long <- alluvial_data %>%
#   pivot_longer(
#     cols = starts_with("Time_"),
#     names_to = "Wave",
#     values_to = "Cluster"
#   ) %>%
#   mutate(
#     Wave = factor(
#       Wave,
#       levels = c("Time_1", "Time_2", "Time_3", "Time_4"),
#       labels = c("Time 1", "Time 2", "Time 3", "Time 4")
#     ),
#     Cluster = factor(Cluster, levels = c("1", "2", "3", "4"))
#   )
# 
# # 2. Get World Spatial Data & Standardize ESS ISO Codes
# # Using `iso_a2_eh` or `gu_a3` resolves Natural Earth's missing ISO bugs (PT, FR, NO)
# world <- ne_countries(scale = "medium", returnclass = "sf") %>%
#   # Fallback logic to guarantee clean 2-letter codes for all European countries
#   mutate(
#     iso2 = case_when(
#       admin == "Portugal" ~ "PT",
#       admin == "France" ~ "FR",
#       admin == "Norway" ~ "NO",
#       admin == "United Kingdom" ~ "GB",
#       admin == "Israel" ~ "IL",
#       TRUE ~ coalesce(iso_a2_eh, iso_a2)
#     )
#   ) %>%
#   select(iso2, name = admin, geometry)
# 
# # 3. Crop map to broader Europe / Mediterranean frame
# europe_base <- world %>%
#   st_crop(xmin = -11, xmax = 36, ymin = 31, ymax = 71)
# 
# # 4. Expand grid across 4 waves so unmeasured countries appear gray in all panels
# map_all_waves <- crossing(
#   Wave = factor(c("Time 1", "Time 2", "Time 3", "Time 4")),
#   iso2 = europe_base$iso2
# ) %>%
#   left_join(europe_base, by = "iso2") %>%
#   left_join(df_map_long, by = c("iso2" = "Country", "Wave" = "Wave")) %>%
#   st_as_sf()
# 
# # 5. Color Palette (consistent across all your figures)
# palette_4 <- c(
#   "1" = "#2B5C8F",
#   "2" = "#D95F02",
#   "3" = "#1B9E77",
#   "4" = "#7570B3"
# )
# 
# # 6. Faceted Map Plot
# p_map <- ggplot(map_all_waves) +
#   geom_sf(aes(fill = Cluster), color = "white", linewidth = 0.25) +
#   scale_fill_manual(
#     values = palette_4,
#     labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
#     na.value = "#ECECEC", # Light neutral gray for non-sample countries
#     name = "Latent Cluster",
#     drop = FALSE
#   ) +
#   # Map coordinates framed for all 18 countries (PT to IL, NO/FI to IL)
#   coord_sf(xlim = c(-10.5, 35.5), ylim = c(31.5, 68), expand = FALSE) +
#   facet_wrap(~Wave, nrow = 1) +
#   labs(
#     title = "Geographic Distribution of Structural Clusters Across Waves",
#     caption = "Gray indicates countries not included in the sample."
#   ) +
#   theme_void(base_size = 11) +
#   theme(
#     strip.text = element_text(face = "bold", size = 11, margin = margin(b = 6)),
#     plot.title = element_text(face = "bold", hjust = 0.5, size = 13, margin = margin(b = 10)),
#     plot.caption = element_text(size = 8, color = "grey40", hjust = 0.5, margin = margin(t = 8)),
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold", size = 10),
#     panel.spacing = unit(0.5, "lines")
#   )
# 
# p_map
# 
# # 6. Save with publication dimensions (Full A4 width, landscape strip)
# ggsave("cluster_spatial_distribution_A4.png", plot = p_map, width = 7.2, height = 3.6, dpi = 300)





# # Full SEM model
# model <- '
#   # Opposition to migration
#   mig_opp =~ NA*imsmetn + 1*imdfetn + impcntr
  
#   # Feelings of threat
#   threath =~ NA*imbgeco + 1*imueclt + imwbcnt
  
#   # Education years
#   edu_yrs =~ eduyrs

#   mig_opp ~ threath + edu_yrs
#   threath ~ edu_yrs
# '

# sem_fit <- lavaan::sam(model = S1_model, 
#                       data  = data[data$essround == "5", ], 
#                       group = "cntry", 
#                       group.equal = c("loadings"),
#                       std.lv = TRUE
#                       )

# summary(sem_fit)



# Configural model -------------------------------------------------------------
# config_fit <- lavaan::cfa(model = S1_model, 
#                           data  = data, 
#                           group = "cntry_time",
#                           missing = "fiml",
#                           std.lv = TRUE
#                           )

# lavaan::fitmeasures(config_fit)[c("cfi", "tli", "rmsea", "srmr")]
# lapply(lavaan::lavInspect(config_fit, what = "est"), "[[", "psi")

# # Metric model -----------------------------------------------------------------
# metric_fit <- lavaan::cfa(model = S1_model, 
#                           data  = data, 
#                           group = "cntry_time", 
#                           group.equal = c("loadings"),
#                           missing = "fiml",
#                           std.lv = TRUE
#                           )

# lavaan::fitmeasures(metric_fit)[c("cfi", "tli", "rmsea", "srmr")]
# lapply(lavaan::lavInspect(metric_fit, what = "est"), "[[", "lambda")

# Get the largest loading per group-time combination
# sapply(lapply(lapply(lavaan::lavInspect(metric_fit, what = "est"), "[[", "lambda"), FUN = \(x) apply(x, 2, which.max)), "[[", "mig_opp")
# sapply(lapply(lapply(lavaan::lavInspect(metric_fit, what = "est"), "[[", "lambda"), FUN = \(x) apply(x, 2, which.max)), "[[", "threath")

# Fit it with the marker variable approach
# Change marker variable to the largest loading
# S1_model <- '
#   # Opposition to migration
#   mig_opp =~ NA*imsmetn + 1*imdfetn + impcntr
  
#   # Feelings of threat
#   threath =~ NA*imbgeco + 1*imueclt + imwbcnt
  
#   # Education years
#   edu_yrs =~ eduyrs
# '

# metric_fit_marker <- lavaan::cfa(model = S1_model, 
#                                  data  = data, # [data$essround == "5", ], 
#                                  group = "cntry_time", 
#                                  group.equal = c("loadings"),
#                                  missing = "fiml"
#                                  )

# lavaan::fitmeasures(metric_fit_marker)[c("cfi", "tli", "rmsea", "srmr")]
# cov_eta <- lapply(lavaan::lavInspect(metric_fit_marker, what = "est"), "[[", "psi")
