# Analysis script
# Placeholder for analysis code.
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
exp_design$RMSE_Lambda <- NA
exp_design$RMSE_Theta  <- NA
rownames(exp_design) <- NULL

# Create functions necessary for the evaluation
eval_model <- function(lambda_ori, theta_ori, lambda_est, theta_est) {
  # Extract estimated non-zero lambdas and turn them into a vector
  lambda_ori_vec  <- unlist(lambda_ori)[unlist(lambda_ori) != 0]
  lambda_est_vec  <- unlist(lambda_est)[unlist(lambda_est) != 0]
  
  theta_ori_vec   <- unlist(theta_ori)[unlist(theta_ori) != 0]
  theta_est_vec   <- unlist(theta_est)[unlist(theta_est) != 0]

  # Compute RMSE for Lambda and Theta
  RMSE_Lambda <- sqrt(mean((lambda_ori_vec - lambda_est_vec)^2))
  RMSE_Theta <- sqrt(mean((theta_ori_vec - theta_est_vec)^2))
  
  return(list(
    RMSE_Lambda = RMSE_Lambda,
    RMSE_Theta  = RMSE_Theta
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
metrics <- c("RMSE_Lambda", "RMSE_Theta") 

for (i in 1:cond) {
  # Pre-calculate condition-specific data if possible
  for (k in 1:K) {
    ik <- ik + 1
    lambdas_path  <- paste0("Lambdas/lambdasRow", i, "Rep", k, ".Rdata")
    thetas_path   <- paste0("Thetas/thetasRow", i, "Rep", k, ".Rdata")
    data_path     <- paste0("Data/DataRow", i, "Rep", k, ".Rdata")
    
    # Load data once per iteration
    lambdas_est <- load_safe(lambdas_path) 
    theta_est   <- load_safe(thetas_path) 
    data        <- load_safe(data_path)
    
    # Evaluate results
    if(!is.null(theta_est$result)) { # Check that there are actual results in the model (avoid non-converged data)
      
      eval <- eval_model(lambda_ori = data$result$lambdas, 
                         theta_ori  = data$result$thetas,
                         lambda_est = lambdas_est,
                         theta_est  = theta_est) 

      # Assignment using dynamic column names
      exp_design[ik, "RMSE_Lambda"] <- eval$RMSE_Lambda
      exp_design[ik, "RMSE_Theta"]  <- eval$RMSE_Theta
    } else {
      exp_design[ik, "RMSE_Lambda"]  <- NA
      exp_design[ik, "RMSE_Theta"] <- NA
    }
  }
}

################################################################################
# FINAL ANALYSIS
################################################################################
# Preliminary (How many have run until now?)
exp_design %>% group_by(nclus)   %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) sum(!is.na(x))))
exp_design %>% group_by(coeff)   %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) sum(!is.na(x))))
exp_design %>% group_by(N_g)     %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) sum(!is.na(x))))
exp_design %>% group_by(S_size)  %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) sum(!is.na(x))))
exp_design %>% group_by(A)       %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) sum(!is.na(x))))
exp_design %>% group_by(ntimes)  %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) sum(!is.na(x))))
exp_design %>% group_by(ngroups) %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) sum(!is.na(x))))

# Performance per simulation factor -----------------------------------------------------------------------------------------------------
exp_design %>% group_by(nclus)   %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) 
exp_design %>% group_by(coeff)   %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) 
exp_design %>% group_by(N_g)     %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) 
exp_design %>% group_by(S_size)  %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) 
exp_design %>% group_by(A)       %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) 
exp_design %>% group_by(ntimes)  %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) 
exp_design %>% group_by(ngroups) %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE)))

# General performance
exp_design %>% summarise(across(RMSE_Lambda:RMSE_Theta, \(x) qwraps2::mean_sd(x, denote_sd = "paren", digits = 3, na_rm = TRUE))) 
