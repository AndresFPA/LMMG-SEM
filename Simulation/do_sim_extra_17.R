# Source necessary functions
invisible(lapply(list.files("hmm-mmgsem", pattern = "\\.R$", full.names = TRUE), source))

# Load final analysis evironment
load("Simulation/Analyses/Analysis_env.Rdata")

# Identify 17 non-converged data sets
non_converged <- row.names(exp_design[is.na(exp_design$ARI), ])

# Prepare model syntax
# Define the models
S1 <- '
  F1 =~ x1 + x2 + x3 + x4 + x5
  F2 =~ x6 + x7 + x8 + x9 + x10
  F3 =~ x11 + x12 + x13 + x14 + x15
  F4 =~ x16 + x17 + x18 + x19 + x20
'

S2 <- '
  F2 ~ F1
  F3 ~ F1
  F4 ~ F1 + F2 + F3
'

NonInv <- c("F1 =~ x2", "F1 =~ x3",
            "F2 =~ x7", "F2 =~ x8",
            "F3 =~ x12", "F3 =~ x13",
            "F4 =~ x17", "F4 =~ x18")

# Prepare simulation functions
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

# Main simulation function
do_sim <- function(Condition){
  # Ensure it is numeric
  Condition <- as.numeric(Condition)
  
  # Get the correct RowDesign and replication number
  RowDesign <- ceiling(Condition/K) # Get the row
  k         <- Condition %% K       # Get correct replication using the remainder function
  if(k == 0){k <- K}                # If the remainde is 0, we know that k = 50
  cat("\n", "Condition", RowDesign, "out of", nrow(design), "\n")
  
  # Load the simulated data
  SimData <- new_load(paste0("C:/Users/User/OneDrive - KU Leuven/0. Postdoc Leuven/1. Papers/Paper 1/R/Simulation/DataRow", RowDesign, "Rep", k, ".Rdata"))

  # Check that there are no empty categories
  data <- as.data.frame(SimData$data)
  
  # First, run the measurement model -------------------------------------------------------
  cfa.fit <- lavaan::cfa(model         = S1, 
                         data          = data, 
                         group         = "group_time", 
                         group.equal   = "loadings", 
                         group.partial = NonInv, 
                         h1            = FALSE,
                         baseline      = FALSE, 
                         se            = "none",
                         test          = "none",
                         check.gradient = FALSE,
                         implied       = FALSE)
  
  # Extract the covariance matrix 
  cov_eta <- lapply(lavaan::lavInspect(cfa.fit, what = "est"), "[[", "psi")
  N_gs    <- lavaan::lavInspect(cfa.fit, "nobs")
  
  # Fit model, including non-invariances
  ctime <- system.time(
    fit <- hmm_mmgsem(cov_eta = cov_eta, 
                      S2      = S2, 
                      ngroups = design[RowDesign, "ngroups"], 
                      nstates = design[RowDesign, "nclus"], 
                      ntimes  = design[RowDesign, "ntimes"], 
                      N_gs    = N_gs, 
                      max_it  = 1000, 
                      s2_fast = T, 
                      nstarts = 20)
    )

  # Save computation times
  save(ctime, file = paste("Times/Time", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  # Save results
  save(fit, file = paste("Fit/Results", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  results <- results$result
}