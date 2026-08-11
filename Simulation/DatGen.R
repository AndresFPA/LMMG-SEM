# 2026-03-25
# Generating repeated cross-sectional data with 'states' depending on regression parameters

dat_gen <- function(N_g,    # Sample size per group
                    G,      # Number of groups
                    S,      # Number of states
                    S_size, # Cluster size
                    a,      # Transition probability
                    Time,   # Number of time points
                    Beta,   # Regression difference
                    Lambda, # Factor loadings
                    Theta   # Residuals 
){

  ######################################################################################################################
  # Define transition matrix
  ######################################################################################################################
  
  # Transition probabilities
  # a = probability of staying in the same cluster
  others_a <- 1 - a
  others_a <- others_a/(S-1)
  
  # Create empty A matrix
  A <- matrix(data = 0, 
              nrow = S, 
              ncol = S)
  
  # Diagonals are the probablity of staying in the same cluster
  diag(A) <- a
  
  # Off-diagonal is the probability of moving to another cluster
  A[col(A) != row(A)] <- others_a
  
  ######################################################################################################################
  # Define cluster/state memberships
  ######################################################################################################################
  
  # Define who is in which cluster at each time point
  s_memberships_t <- vector(mode = "list", length = Time)
  
  # Initial state - Timepoint 1
  if(S_size == "bal"){
    s_memberships_t[[1]] <- rep(x = 1:S, each = G/S) # Balanced clustering at time point 1
  } else if (S_size == "unb"){
    # First cluster contains 75% of the groups
    K1_idx <- (G*0.75)
    s_memberships_t[[1]][1:K1_idx] <- 1
    
    # Remaining clusters contain 25% of the groups. They are assigned randomly
    # Temporal object to ensure all clusters contain at least one group
    tmp_test <- 0
    while(tmp_test < (S-1)){ # While loop to ensure that cluster contains at least one group
      if(S == 2){
        tmp_memb <- rep(x = 2, times = (G*0.25))
      } else {
        tmp_memb <- sort(sample(x = c(2:S), size = (G*0.25), replace = T)) # Random clustering at time point 1
      }
      tmp_test <- length(unique(tmp_memb))
      s_memberships_t[[1]][(K1_idx+1):G] <- tmp_memb
    }
  }
  
  
  # Define cluster membership in posterior timepoints based on the transition matrix
  for (t in 2:Time) {
    # We sample the next state based on the probabilities of the current state
    s_memberships_t[[t]] <- sapply(X = 1:G,                                                   # We find cluster membership for each group individually
                                   FUN = \(x) sample(x = 1:S,                                 # Sample from all possible states
                                                     size = 1,                                # Only one cluster membership (per group)
                                                     prob = A[s_memberships_t[[t - 1]][x], ]  # Probabilities are given based on the membership of the same group at the previous time point
                                                     ) # End sample
                                   ) # End sapply
  }
  
  ######################################################################################################################
  # Define structural model (which defines the clusters/states)
  ######################################################################################################################

  # Betas ----------------------------------------------------------------------
  # Important:
  # In the beta vectors, the positions of each refer to the following regressions
  # beta[1]: F1 -> F2
  # beta[2]: F1 -> F3
  # beta[3]: F1 -> F4
  # beta[4]: F2 -> F4
  # beta[5]: F3 -> F4
  
  # Initialize empty betas matrix and vector for each cluster
  betas_mat_s <- vector(mode = "list", length = S)
  betas_vec_s <- vector(mode = "list", length = S)
  
  # Create empty beta matrix
  beta_mat <- matrix(data = 0, nrow = 4, ncol = 4)
  colnames(beta_mat) <- rownames(beta_mat) <- c("F1", "F2", "F3", "F4")
  
  # Define (general) betas vector
  betas <- rep(Beta, 5)
  
  # Fill in the different beta matrices depending on the cluster membership
  for(s in 1:S){
    # Which beta is 0 depends on the cluster
    betas_vec_s[[s]]    <- betas # General betas input
    betas_vec_s[[s]][s] <- 0     # Replace one beta with 0 depending on the cluster

    # Input the regressions in their correct positions in the matrix
    betas_mat_s[[s]] <- beta_mat
    betas_mat_s[[s]][c(2, 3, 4, 8, 12)] <- betas_vec_s[[s]]
  }
  
  # Psi ------------------------------------------------------------------------
  # Define the psi matrix PER TIME POINT!
  # Initialize necessary objects
  psi_gs_t          <- vector(mode = "list", length = Time) # For the psi matrices
  psi_mat           <- matrix(data = 0, nrow = 4, ncol = 4)
  colnames(psi_mat) <- rownames(psi_mat) <- c("F1", "F2", "F3", "F4")
  
  # Define the cov_eta (phi) per time point!
  cov_eta_t    <- vector(mode = "list", length = Time)
  
  # Identity matrix
  I <- diag(4)
  
  # Empty objects for the individual variances
  F1_var_t     <- vector(mode = "list", length = Time) # F1 variances 
  F2_var_t     <- vector(mode = "list", length = Time) # F2 (total) variances 
  F3_var_t     <- vector(mode = "list", length = Time) # F3 (total) variances
  F4_var_t     <- vector(mode = "list", length = Time) # F4 (total) variances
  
  F2_res_var_t <- vector(mode = "list", length = Time) # F2 residual variances
  F3_res_var_t <- vector(mode = "list", length = Time) # F3 residual variances
  F4_res_var_t <- vector(mode = "list", length = Time) # F4 residual variances
  
  # Input necessary objects within each time point
  # Compute residual variances
  for(t in 1:Time){
    psi_gs_t[[t]]  <- vector(mode = "list", length = G) # Within each time point, enough lists for all groups
    cov_eta_t[[t]] <- vector(mode = "list", length = G) # Within each time point, enough lists for all groups
    
    F1_var_t[[t]] <- runif(n = G, min = 0.8, max = 1.2) # Total variances around 1
    F2_var_t[[t]] <- runif(n = G, min = 0.8, max = 1.2) # Total variances around 1
    F3_var_t[[t]] <- runif(n = G, min = 0.8, max = 1.2) # Total variances around 1
    F4_var_t[[t]] <- runif(n = G, min = 0.8, max = 1.2) # Total variances around 1
    
    # Prepare residual variance empty objects
    F2_res_var_t[[t]] <- numeric(G)
    F3_res_var_t[[t]] <- numeric(G)
    F4_res_var_t[[t]] <- numeric(G)
    
    # Extract cluster memberships current time point
    current_membs <- s_memberships_t[[t]]
    
    # compute residual variances
    for(g in 1:G){
      memb_this_g  <- current_membs[[g]]
      betas_this_g <- betas_vec_s[[memb_this_g]]
      F2_res_var_t[[t]][g] <- F2_var_t[[t]][g] - (betas_this_g[1]^2 * F1_var_t[[t]][g])
      F3_res_var_t[[t]][g] <- F3_var_t[[t]][g] - (betas_this_g[2]^2 * F1_var_t[[t]][g])
      F4_res_var_t[[t]][g] <- F4_var_t[[t]][g] - (
        (betas_this_g[3]^2 * F1_var_t[[t]][g]) +     # Direct effect from F1 to F4
        (betas_this_g[4]^2 * F2_res_var_t[[t]][g]) + # Direct effect from F2 to F4
        (betas_this_g[5]^2 * F3_res_var_t[[t]][g]) + # Direct effect from F3 to F4
        (2 * betas_this_g[3] * betas_this_g[4] * (betas_this_g[1]^2 * F1_var_t[[t]][g])) +                   # Adjustment: cov between F1 and F2
        (2 * betas_this_g[3] * betas_this_g[5] * (betas_this_g[2]^2 * F1_var_t[[t]][g])) +                   # Adjustment: cov between F1 and F3
        (2 * betas_this_g[4] * betas_this_g[5] * (betas_this_g[1]^2 * betas_this_g[2]^2 * F1_var_t[[t]][g])) # Adjustment: cov between F2 and F3
      )
      
      # Input the final residual covariance matrix into psi
      psi_gs_t[[t]][[g]] <- psi_mat
      psi_gs_t[[t]][[g]]["F1", "F1"] <- F1_var_t[[t]][g]
      psi_gs_t[[t]][[g]]["F2", "F2"] <- F2_res_var_t[[t]][g]
      psi_gs_t[[t]][[g]]["F3", "F3"] <- F3_res_var_t[[t]][g]
      psi_gs_t[[t]][[g]]["F4", "F4"] <- F4_res_var_t[[t]][g]
      
      # Compute the complete covariance matrix (cov_eta)
      betas_mat_this_g <- betas_mat_s[[memb_this_g]]
      cov_eta_t[[t]][[g]] <- solve(I - betas_mat_this_g) %*% psi_gs_t[[t]][[g]] %*% solve(t(I - betas_mat_this_g))
    }
  }
  
  ######################################################################################################################
  # Define measurement model
  ######################################################################################################################
  # Non-invariant Lambda
  lambda_vec <- rep(x = c(1, rep(sqrt(0.6), 4), rep(0, 20)), times = 4)[1:(20*4)]
  lambda     <- matrix(data = lambda_vec,
                       nrow = 20, # Items
                       ncol = 4)  # Factors
  
  # Initialize group-specific lambda
  lambda_gs_t <- vector(mode = "list", length = Time)
  
  for(t in 1:Time){
    lambda_gs_t[[t]] <- vector(mode = "list", length = G)
    for(g in 1:G){
      # Non-invariance - How many groups?
      # Sample random non-invariances
      NonInvariantLoadings <- sample(x = c(runif(100, min = (sqrt(0.6) - 0.4) - .1, max = (sqrt(0.6) - 0.4) + .1), 
                                           runif(100, min = (sqrt(0.6) + 0.4) - .1, max = (sqrt(0.6) + 0.4) + .1)),
                                     size = 2*4)
      
      # Create a non-invariant lambda matrix
      lambda_non_inv <- lambda
      lambda_non_inv[2:3,1] <- NonInvariantLoadings[1:2]
      lambda_non_inv[7:8,2] <- NonInvariantLoadings[3:4]
      lambda_non_inv[12:13,3] <- NonInvariantLoadings[5:6]
      lambda_non_inv[17:18,4] <- NonInvariantLoadings[7:8]
      
      # Save lambda for future evaluation
      lambda_gs_t[[t]][[g]] <- lambda_non_inv 
    }
  }
  
  
  # Theta
  theta_gs_t <- vector(mode = "list", length = Time)
  for(t in 1:Time){
    theta_gs_t[[t]] <- lapply(1:G, \(x) diag(runif(n = 20, min = 0.3, max = 0.5)))
  }
  
  ######################################################################################################################
  # Generate observed covariance matrix
  ######################################################################################################################
  # Generate the observed covariance matrix
  Sigma_t <- vector(mode = "list", length = Time)
  for(t in 1:Time){
    Sigma_t[[t]] <- lapply(1:G, \(x) lambda_gs_t[[t]][[x]] %*% cov_eta_t[[t]][[x]] %*% t(lambda_gs_t[[t]][[x]]) + theta_gs_t[[t]][[x]])
  }
  
  # Create empty data matrix
  data <- as.data.frame(matrix(data = NA, 
                               nrow = N_g*G*Time,
                               ncol = 22))
  colnames(data) <- c("Group", "Timepoint", paste0("x", 1:20))
  
  # Initialize counter
  gt <- 0

  # Fill the matrix
  for(g in 1:G){
    for(t in 1:Time){
      gt <- gt + 1
      start_idx <- ((N_g*(gt - 1)) + 1)
      end_idx   <- (N_g*gt)
      data[start_idx:end_idx, 3:22] <- MASS::mvrnorm(n = N_g, 
                                                     mu = rep(0, 20),
                                                     Sigma = Sigma_t[[t]][[g]],
                                                     empirical = F)
      data$Timepoint[start_idx:end_idx] <- t
      data$Group[start_idx:end_idx]     <- g
    }
  }
  
  # Reorganize data (must show all groups per timepoint)
  idx  <- sort(data$Timepoint, index.return = T)$ix
  data <- data[idx, ]
  
  # Add extra column for all group-time combinations
  data$group_time <- rep(1:(g*t), each = N_g)
  
  return(list(data     = data, 
              clusters = s_memberships_t, 
              lambdas  = lambda_gs_t,
              thetas   = theta_gs_t)
         )
}
