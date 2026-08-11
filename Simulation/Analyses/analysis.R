# Analyze the measurement model estimates
# Set the working directory
setwd("C:/Users/User/OneDrive - KU Leuven/0. Postdoc Leuven/1. Papers/Paper 1/R/Simulation")
setwd("C:/Users/u0159267/OneDrive - KU Leuven/0. Postdoc Leuven/1. Papers/Paper 1/R/Simulation")

# Load libraries
library(dplyr)
library(clue)
library(ggplot2)
library(tidyr)

# Create necessary objects 
cond <- 432
K    <- 50
data <- numeric(cond*K)

# Load the design matrix
load("C:/Users/User/OneDrive - KU Leuven/0. Postdoc Leuven/1. Papers/Paper 1/R/Simulation/Fit/design.Rdata")
load("C:/Users/u0159267/OneDrive - KU Leuven/0. Postdoc Leuven/1. Papers/Paper 1/R/Simulation/Fit/design.Rdata")

exp_design           <- design[rep(seq_len(nrow(design)), each = K), ] # Expanded design matrix
exp_design$Condition <- rep(1:cond, each = K)
exp_design$Repe      <- rep(1:K, times = cond)
exp_design$model     <- NULL # remove the model column
last_idx             <- ncol(exp_design)
exp_design           <- exp_design[, c((last_idx - 1), last_idx, (1:(last_idx - 2)))] # Reorder columns

# Add the columns to store the final results
exp_design$ARI        <- NA
exp_design$RMSE       <- NA
rownames(exp_design) <- NULL

# Create functions necessary for the evaluation
# computation of adjusted rand index
adjrandindex <- function(part1, part2){
  part1 <- as.numeric(part1)
  part2 <- as.numeric(part2)
  
  IM1 <- diag(max(part1))
  IM2 <- diag(max(part2))
  A <- IM1[part1,]
  B <- IM2[part2,]
  
  T <- t(A)%*%B
  N <- sum(T)
  Tc <- apply(T,2,sum)
  Tr <- apply(T,1,sum)
  a <- (sum(T^2) - N)/2
  b <- (sum(Tr^2) - sum(T^2))/2
  c <- (sum(Tc^2) - sum(T^2))/2
  d <- (sum(T^2) + N^2 - sum(Tr^2) - sum(Tc^2))/2
  ARI <- (choose(N,2)*(a + d) - ((a+b)*(a+c)+(c+d)*(b+d)))/(choose(N,2)^2 - ((a+b)*(a+c)+(c+d)*(b+d)))
  
  return(ARI)
}

eval_model <- function(beta_list, z_gks, original_clus, nclus, coeff, psi_gks, tol = 1e-8) {
  # 1. Hard Clustering & True Labels
  clusVec_est  <- c(apply(X = z_gks, MARGIN = 3, FUN = max.col))
  clusVec_true <- unlist(original_clus)
  
  if(any(is.na(clusVec_est))){
    return(list(
      ARI         = NA,
      #CorrectClus = correct_clus,
      RMSE        = NA
      #exo_mean    = mean(sapply(psi_gks[, 1], '[[', 1)),
      #cov_mean    = mean(sapply(psi_gks[, 1], '[[', 2))
    ))
  }
  
  # 2. Find the Optimal Mapping
  # Create a confusion matrix (rows = estimated, cols = true)
  con_mat <- table(clusVec_est, clusVec_true)
  
  # solve_LSAP finds the permutation that maximizes the sum of the diagonal (good classification)
  # In other words, we are minimizing the missclassification error rate
  # We use 'max(con_mat) - con_mat' because LSAP minimizes cost
  mapping <- as.vector(solve_LSAP(con_mat, maximum = TRUE))
  
  # 'mapping' tells us: Estimated Cluster [i] corresponds to True Cluster [mapping[i]]
  # We need to reorder our estimated objects to match the "True" order (1, 2, 3...)
  reorder_idx <- order(mapping) 
  
  # 3. Parameter Recovery (Now correctly mapped)
  # Reorder beta_list so that the first element matches the first true cluster
  beta_list_mapped <- beta_list[reorder_idx]
  
  # Extract estimated non-zero betas into a matrix
  # Assumes 5 regression paths per cluster
  est_betas <- t(sapply(beta_list_mapped, function(b) {
    v <- as.vector(b)
    v[v != 0] 
  }))
  # browser()
  # 4. Define True Betas (Mapped to the 1, 2, 3... order)
  true_betas <- matrix(coeff, nrow = nclus, ncol = 5)
  # Next, we replace one regression per cluster to 0
  # When removing the diagonal, we are saying: in cluster 1, beta1 = 0; in cluster 2, beta2 = 0, etc. 
  diag(true_betas) <- 0 
  true_betas <- true_betas[1:nclus,] 
  
  # 5. RMSE and Bias (Floating-point safe)
  sq_err <- (true_betas - est_betas)^2
  rmse_betas <- sqrt(mean(sq_err))
  
  # 6. Cluster Recovery (ARI is permutation invariant, so this is safe)
  # clusVec_est <- max.col(z_gks[reorder_idx]) # Re-order the matrix to get the correct values
  ARI_res <- adjrandindex(clusVec_est, clusVec_true)
  # correct_clus <- all(clusVec_est == clusVec_true) # Simplified check
  
  return(list(
    ARI         = ARI_res,
    #CorrectClus = correct_clus,
    RMSE        = rmse_betas
    #exo_mean    = mean(sapply(psi_gks[, 1], '[[', 1)),
    #cov_mean    = mean(sapply(psi_gks[, 1], '[[', 2))
  ))
}



# Create new loading function
new_load <- function(file) {
  # Create a temporary environment to load the data into
  temp_env <- new.env()
  
  # Load the file into that specific environment
  loaded_names <- load(file, envir = temp_env)
  
  # Check if the file contained more than one object
  if (length(loaded_names) > 1) {
    warning("The file contains multiple objects. Returning a list of all objects.")
    return(as.list(temp_env))
  }
  
  # Return the single object
  return(temp_env[[loaded_names]])
}

# Make it safe to not stop when finding an error
load_safe <- purrr::safely(new_load)

# Load and analyze the data
ik <- 0

# Pre-define the metrics to collect
metrics <- c("ARI", "RMSE") 
# models  <- c("ML", "DWLS", "ML.ign", "DWLS.ign")

for(i in 1:cond) {
  # Pre-calculate condition-specific data if possible
  for(k in 1:K) {
    ik <- ik + 1
    res_path  <- paste0("Fit/ResultsRow", i, "Rep", k, ".Rdata")
    data_path <- paste0("Data/DataRow", i, "Rep", k, ".Rdata")
    
    # Load data once per iteration
    fit  <- load_safe(res_path) 
    data <- load_safe(data_path)
    
    # Evaluate results
    if(!is.null(fit$result)) { # Check that there are actual results in the model (avoid non-converged data)
      fit   <- fit$result
      
      eval <- eval_model(beta_list     = fit$param$beta_ks, 
                         z_gks         = fit$post,  
                         original_clus = data$result$clusters, # or the original cluster mat
                         nclus         = exp_design[ik, "nclus"], 
                         coeff         = exp_design[ik, "coeff"]) #,
                         # psi_gks       = fit$param$psi)
      
      # Assignment using dynamic column names
      exp_design[ik, "ARI"]  <- eval$ARI
      exp_design[ik, "RMSE"] <- eval$RMSE
    } else {
      exp_design[ik, "ARI"]  <- NA
      exp_design[ik, "RMSE"] <- NA
    }
  }
}

################################################################################
# FINAL ANALYSIS
################################################################################
# Preliminary (How many have run until now?)
exp_design %>% group_by(nclus)   %>% summarise(across(ARI:RMSE, \(x) sum(!is.na(x))))
exp_design %>% group_by(coeff)   %>% summarise(across(ARI:RMSE, \(x) sum(!is.na(x))))
exp_design %>% group_by(N_g)     %>% summarise(across(ARI:RMSE, \(x) sum(!is.na(x))))
exp_design %>% group_by(S_size)  %>% summarise(across(ARI:RMSE, \(x) sum(!is.na(x))))
exp_design %>% group_by(A)       %>% summarise(across(ARI:RMSE, \(x) sum(!is.na(x))))
exp_design %>% group_by(ntimes)  %>% summarise(across(ARI:RMSE, \(x) sum(!is.na(x))))
exp_design %>% group_by(ngroups) %>% summarise(across(ARI:RMSE, \(x) sum(!is.na(x))))

# Cluster recovery -----------------------------------------------------------------------------------------------------
# Add first a perfect recovery measure
exp_design$CorrectClus <- exp_design$ARI == 1

# Check the results per simulation factor
exp_design %>% group_by(nclus)   %>% summarise(across(ARI:CorrectClus, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE, show_n = "never")))
exp_design %>% group_by(coeff)   %>% summarise(across(ARI:CorrectClus, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE, show_n = "never"))) 
exp_design %>% group_by(N_g)     %>% summarise(across(ARI:CorrectClus, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE, show_n = "never"))) 
exp_design %>% group_by(S_size)  %>% summarise(across(ARI:CorrectClus, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE, show_n = "never"))) 
exp_design %>% group_by(A)       %>% summarise(across(ARI:CorrectClus, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE, show_n = "never"))) 
exp_design %>% group_by(ntimes)  %>% summarise(across(ARI:CorrectClus, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE, show_n = "never"))) 
exp_design %>% group_by(ngroups) %>% summarise(across(ARI:CorrectClus, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE, show_n = "never"))) 

exp_design %>% select(ARI, RMSE, CorrectClus) %>% apply(., 2, qwraps2::mean_sd, denote_sd = "paren", digits = 3, na_rm = TRUE) 



library(ggplot2)
library(tidyr)
library(dplyr)

# 1. Aggregate the data first to handle NAs and calculate means
plot_summary <- exp_design %>%
  group_by(coeff, nclus, N_g) %>%
  summarise(
    ARI  = mean(ARI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 2. Pivot for ggplot
  pivot_longer(
    cols = c(ARI),
    names_to = "Measure",
    values_to = "Value"
  ) # %>%
  # mutate(Measure = ifelse(Measure == "ARI", "ARI ML", "ARI DWLS"))

# 3. Plot the aggregated data
ggplot(plot_summary, aes(x = N_g, y = Value)) +
  geom_line(linewidth = 1, color = "#FF0000") +
  geom_point(size = 2 , color = "#FF0000") +
  facet_grid(coeff ~ nclus, labeller = labeller(
    coeff = function(x) paste("Beta:", x), 
    nclus = function(x) paste("Number of Clusters:", x)
  )) +
  # scale_color_manual(values = c("ARI" = "#FF0000", "RMSE" = "#008B8B"), 
  #                    labels = c("ARI", "RMSE")) + 
  labs(
    x = "Within-group Sample Size",
    y = "Cluster Recovery"# ,
    # color = "Measure"
  ) +
  scale_y_continuous(limits = c(0, 1.00), breaks = seq(0, 1.00, 0.25)) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  filename = "cluster_recovery_plot.png",
  width    = 16.5,
  height   = 12.5,
  units    = "cm",
  dpi      = 300
)

# Secondary effects plot
# 1. Aggregate the data first to handle NAs and calculate means
plot_summary <- exp_design %>%
  group_by(S_size, A) %>%
  summarise(
    ARI = mean(ARI, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Plot the aggregated data as a heatmap
# 2. Plot the aggregated data as a heatmap
ggplot(plot_summary, aes(
  x = factor(A), 
  y = factor(S_size, levels = c("random", "bal"), labels = c("Unbalanced", "Balanced")), 
  fill = ARI
)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", ARI)), color = "black", size = 5.5) +
  scale_fill_gradient(
    low = "lightgreen", 
    high = "#006837", 
    # limits = c(0.7, 1.0),
    name = "ARI"
  ) +
  labs(
    x = "Transition Probability",
    y = "Initial Cluster Size"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  filename = "cluster_recovery_heatmap.png",
  width    = 16.5,
  height   = 12.5,
  units    = "cm",
  dpi      = 300
)

# Regression recovery --------------------------------------------------------------------------------------------------
exp_design %>% group_by(nclus)   %>% summarise(across(ARI:RMSE, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) %>% select(RMSE)#, RMSE_ML, RMSE_DWLS)
exp_design %>% group_by(coeff)   %>% summarise(across(ARI:RMSE, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) %>% select(RMSE)#, RMSE_ML, RMSE_DWLS)
exp_design %>% group_by(N_g)     %>% summarise(across(ARI:RMSE, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) %>% select(RMSE)#, RMSE_ML, RMSE_DWLS)
exp_design %>% group_by(S_size)  %>% summarise(across(ARI:RMSE, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) %>% select(RMSE)#, RMSE_ML, RMSE_DWLS)
exp_design %>% group_by(A)       %>% summarise(across(ARI:RMSE, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) %>% select(RMSE)#, RMSE_ML, RMSE_DWLS)
exp_design %>% group_by(ntimes)  %>% summarise(across(ARI:RMSE, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) %>% select(RMSE)#, RMSE_ML, RMSE_DWLS)
exp_design %>% group_by(ngroups) %>% summarise(across(ARI:RMSE, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) %>% select(RMSE)#, RMSE_ML, RMSE_DWLS)
exp_design %>% select(RMSE) %>% apply(., 2, qwraps2::mean_sd, denote_sd = "paren", digits = 3, na_rm = TRUE) 

# 1. Aggregate the data first to handle NAs and calculate means
plot_summary <- exp_design %>%
  group_by(S_size, nclus, N_g) %>%
  summarise(
    RMSE  = mean(RMSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 2. Pivot for ggplot
  pivot_longer(
    cols = c(RMSE),
    names_to = "Measure",
    values_to = "Value"
  ) # %>%
# mutate(Measure = ifelse(Measure == "ARI", "ARI ML", "ARI DWLS"))

# 3. Plot the aggregated data
ggplot(plot_summary, aes(x = N_g, y = Value)) +
  geom_line(linewidth = 1, color = "#008B8B") +
  geom_point(size = 2, color = "#008B8B") +
  facet_grid(S_size ~ nclus, labeller = labeller(
    nclus = function(x) paste("Number of Clusters:", x), 
    S_size = as_labeller(c(
     "bal"    = "Initial Cluster Size: Balanced",
     "random" = "Initial Cluster Size: Unbalanced"
    ))
  )) +
  # scale_color_manual(values = c("ARI" = "#FF0000", "RMSE" = "#008B8B"), 
  #                    labels = c("ARI", "RMSE")) + 
  labs(
    x = "Within-group Sample Size",
    y = "Parameter Recovery"# ,
    # color = "Measure"
  ) +
  # scale_y_continuous(limits = c(0, 1.00), breaks = seq(0, 1.00, 0.25)) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  filename = "regression_recovery_plot.png",
  width    = 16.5,
  height   = 12.5,
  units    = "cm",
  dpi      = 300
)


























# JOINT PLOT -------------------------------------------------------------------
# 1. Aggregate the data first to handle NAs and calculate means
plot_summary <- exp_design %>%
  group_by(coeff, A, N_g) %>%
  summarise(
    ARI  = mean(ARI, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 2. Pivot for ggplot
  pivot_longer(
    cols = c(ARI, RMSE),
    names_to = "Measure",
    values_to = "Value"
  ) # %>%
# mutate(Measure = ifelse(Measure == "ARI", "ARI ML", "ARI DWLS"))

# 3. Plot the aggregated data
ggplot(plot_summary, aes(x = N_g, y = Value, color = Measure)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(A ~ coeff, labeller = labeller(
    coeff = function(x) paste("Beta:", x), 
    A = function(x) paste("Transitions (A):", x)
  )) +
  scale_color_manual(values = c("ARI" = "#FF0000", "RMSE" = "#008B8B"), 
                     labels = c("ARI", "RMSE")) + 
  labs(
    x = "Within-group Sample Size",
    y = "Cluster Recovery",
    color = "Measure"
  ) +
  scale_y_continuous(limits = c(0, 1.00), breaks = seq(0, 1.00, 0.25)) +
  theme_bw() +
  theme(legend.position = "bottom")