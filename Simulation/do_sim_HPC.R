# 2026-03-31
# First draft: Simulation script

# Load libraries
library(lavaan)
library(purrr)

# Source necessary functions
source("hmm-mmgsem/hmm_mmgsem.R")
source("hmm-mmgsem/E-Step.R")
source("hmm-mmgsem/M-Step.R")
source("hmm-mmgsem/sem_est_fast.R")
source("hmm-mmgsem/Helper.R")

# Set second working directory
source("DatGen.R")

# Simulation Design
# Which factors are going to be tested? For now:
nclus   <- c(2, 4, 5)      # Number of clusters
ntimes  <- c(2, 3, 4)      # Time points
ngroups <- c(20, 40)       # Number of groups
coeff   <- c(0.3, 0.4)     # Initial regression parameters
N_g     <- c(50, 150, 500) # Sample size per groups
A       <- c(0.7, 0.9)     # Transition probabilities
S_size  <- c("bal", "unb") # Cluster size

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

# Get design matrix
design <- expand.grid(nclus, ntimes, ngroups, coeff, N_g, A, S_size) 
colnames(design) <- c("nclus", "ntimes", "ngroups", "coeff", "N_g", "A", "S_size")

rownames(design) <- NULL
rm(nclus, ntimes, ngroups, coeff, N_g, A, S_size) 

# Functions for the simulation
# First, to avoid stopping due to errors, create a function with data generation and MMGSEM
genDat_analysis <- function(seed, RowDesign, k, NonInv){
  # browser()
  # Set seed per design condition (row) and replication (K)
  set.seed(seed)
  # Generate data
  #SimData <- do.call(what = DataGeneration, args = design[RowDesign, ])$SimData
  SimData <- dat_gen(S      = design[RowDesign, "nclus"], 
                     Time   = design[RowDesign, "ntimes"], 
                     G      = design[RowDesign, "ngroups"], 
                     Beta   = design[RowDesign, "coeff"], 
                     N_g    = design[RowDesign, "N_g"], 
                     a      = design[RowDesign, "A"],
                     S_size = design[RowDesign, "S_size"])
  
  # Check that there are no empty categories
  data <- as.data.frame(SimData$data)
  
  # Save data
  save(SimData, file = paste("Data/Data", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
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
  
  # Extract the measurement model parameters 
  lambdas <- lapply(lavaan::lavInspect(cfa.fit, what = "est"), "[[", "lambda")
  save(lambdas, file = paste("Lambdas/lambdas", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
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
  
  # If everything goes right, return results
  return(fit)
}

# Create SAFE function (in case of errors)
genDat_analysis <- purrr::safely(.f = genDat_analysis, otherwise = NULL)

# Main simulation function
do_sim <- function(Condition){
  # Ensure it is numeric
  Condition <- as.numeric(Condition)
  
  # Get the correct RowDesign and replication number
  RowDesign <- ceiling(Condition/K) # Get the row
  k         <- Condition %% K       # Get correct replication using the remainder function
  if(k == 0){k <- K}                # If the remainde is 0, we know that k = 50
  cat("\n", "Condition", RowDesign, "out of", nrow(design), "\n")
  
  # Code to re-sample in case there is some error in the process (e.g., invalid covariance matrices)
  attempts <- 1
  for(j in 1:attempts){
    print(paste("Attempt #", j))
    # Seed will change if there is an error
    results <- genDat_analysis(seed = (RowDesign * k * j), RowDesign = RowDesign, k = k, NonInv = NonInv)
    test    <- results$result 
    if(!is.null(test)){
      # If there was no error, break the loop and continue
      print("all good")
      break
    }
  }
  
  results <- results$result
}

# ###################################################################### #
# ######################## START PARALLELIZATION ####################### #
# ###################################################################### #
# Define number of replications
K <- 50 # Number of replications per condition

# Parallel execution over all possible RowDesign*K conditions
Cond <- nrow(design) * K

Args <- commandArgs(TRUE)
do_sim(Args[1])
