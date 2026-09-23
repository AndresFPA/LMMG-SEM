# ==============================================================================
# Example: Data Simulation and Measurement Invariance Analysis (LMMG-SEM)
# Using repository functions from Simulation/ and hmm-mmgsem/
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries and Source Repository Functions
# ------------------------------------------------------------------------------
library(lavaan)
library(MASS)
library(semTools)
library(sirt)

# Source data generation function
source("Simulation/DatGen.R")

# Source HMM-MMGSEM helper and estimation functions
source("hmm-mmgsem/Helper.R")
source("hmm-mmgsem/E-Step.R")
source("hmm-mmgsem/M-Step.R")
source("hmm-mmgsem/sem_est_fast.R")
source("hmm-mmgsem/hmm_mmgsem.R")

# Source two-step longitudinal measurement invariance testing function
source("hmm-mmgsem/mi.R")
source("hmm-mmgsem/mi_stepwise_score.R")

# Set seed for reproducibility
set.seed(1)

# ------------------------------------------------------------------------------
# 2. Simulate Data Using dat_gen()
# ------------------------------------------------------------------------------
# Simulation conditions
n_groups   <- 12     # Number of groups (G)
n_times    <- 2      # Number of time points (T)
n_states   <- 2      # Number of latent states/clusters (S)
n_per_grp  <- 150    # Sample size per group (N_g)
trans_prob <- 0.8    # Probability of staying in the same state (a)
reg_coeff  <- 0.4    # Structural regression baseline coefficient (Beta)
s_size     <- "bal"  # Balanced state size distribution ("bal" or "unb")

sim_res <- dat_gen(
  N_g    = n_per_grp,
  G      = n_groups,
  S      = n_states,
  S_size = s_size,
  a      = trans_prob,
  Time   = n_times,
  Beta   = reg_coeff
)

# Extract simulated dataset and components
sim_data     <- sim_res$data
true_clus    <- sim_res$clusters
true_lambdas <- sim_res$lambdas

# ------------------------------------------------------------------------------
# 3. Define Measurement and Structural Model Syntax
# ------------------------------------------------------------------------------
# S1: Measurement model (4 latent factors, 20 indicator items)
S1 <- '
  F1 =~ x1 + x2 + x3 + x4 + x5
  F2 =~ x6 + x7 + x8 + x9 + x10
  F3 =~ x11 + x12 + x13 + x14 + x15
  F4 =~ x16 + x17 + x18 + x19 + x20
'

# S2: Structural model (regressions between latent factors across states)
S2 <- '
  F2 ~ F1
  F3 ~ F1
  F4 ~ F1 + F2 + F3
'

# ------------------------------------------------------------------------------
# 4. Step 1: Sequential Measurement Invariance Analysis (compute_mi_long)
# ------------------------------------------------------------------------------
mi_results <- compute_mi_long_stepwise_score(
  data              = sim_data,
  model             = S1,
  group_var         = "Group",
  time_var          = "Timepoint",
  alpha             = 0.01,
  fit_indices       = c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
  d_cfi_threshold   = 0.020, 
  max_iter          = 1000, 
  method            = "test_statistic"
)

# mi_results <- compute_mi_long(
#   data              = sim_data,
#   model             = S1,
#   group_var         = "Group",
#   time_var          = "Timepoint",
#   alpha             = 0.01,
#   fit_indices       = c("chisq", "df", "pvalue", "cfi", "rmsea", "srmr"),
#   d_cfi_threshold   = 0.010,
#   method            = "alignment",
#   specific_noninv   = TRUE
# )

# ------------------------------------------------------------------------------
# 5. Review Measurement Invariance Results
# ------------------------------------------------------------------------------
print(mi_results$fit_measures)
print(sort(mi_results$final_partial_constraints))
print(round(mi_results$fit_measures$final_model, 4))
summary(mi_results$final_fit)

# ------------------------------------------------------------------------------
# 6. Step 2: Longitudinal Mixture Multi-Group SEM (hmm_mmgsem)
# ------------------------------------------------------------------------------
cat("\nRunning Step 2: HMM-MMGSEM on estimated latent factor covariance matrices...\n")

# Extract factor covariance matrices (phi) across G x T cells from the MI step
cov_eta <- mi_results$phi_matrices
n_obs   <- rep(n_per_grp, n_groups * n_times)

# Fit HMM-MMGSEM model
fit_hmm <- hmm_mmgsem(
  cov_eta  = cov_eta,
  S2       = S2,
  ngroups  = n_groups,
  nstates  = n_states,
  ntimes   = n_times,
  N_gs     = n_obs,
  max_it   = 100,
  s2_fast  = TRUE,
  nstarts  = 10,
  printing = FALSE
)

# ------------------------------------------------------------------------------
# 7. Model Results Summary
# ------------------------------------------------------------------------------
cat("\n================ HMM-MMGSEM Results ================\n")
cat("Estimated Transition Matrix (A):\n")
print(round(fit_hmm$A, 4))

cat("\nEstimated Posterior State Probabilities (dim: Groups x States x Time):\n")
print(round(fit_hmm$post[, , 1], 4)) # Time 1 posteriors for each group

cat("\nScript completed successfully using repository simulation and analysis functions.\n")