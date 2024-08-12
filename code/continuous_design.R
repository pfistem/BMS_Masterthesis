# Drop the loser for MAMS | Continuous Endpoint
# Bristol Myers Squibb
# Contact: manuel.pfister@bms.com

# -------------------------------------------------------------------------

library(mvtnorm)

# Function to calculate the integral of the Type I error
integral_typeIerror <- function(c, requiredtypeIerror, mean, var, combinations, K) {
  lower <- rep(0, length(mean))
  lower[length(lower)] <- c
  upper <- rep(Inf, length(mean))
  int <- pmvnorm(lower = lower, upper = upper, mean = mean, sigma = var)
  return(as.double(int) * factorial(K) - requiredtypeIerror)
}

# Function to calculate the integral of the power
integral_power <- function(n, requiredpower, c, A, var, combinations, J, K, delta1, delta0, cumgroupsizes) {
  # Get mu under LFC (assuming standardized data)
  mu <- rep(0, length(A[1,]))
  for (i in 1:(K - 1)) {
    mu[((i - 1) * J + 1):(i * J)] <- sqrt(n * cumgroupsizes / 2) * delta0
  }
  mu[((K - 1) * J + 1):(K * J)] <- sqrt(n * cumgroupsizes / 2) * delta1
  mean <- as.double(A %*% mu)
  lower <- rep(0, length(mean))
  lower[length(lower)] <- c
  upper <- rep(Inf, length(mean))
  int <- pmvnorm(lower = lower, upper = upper, mean = mean, sigma = var)
  return(as.double(int) * factorial(K - 1) - requiredpower)
}

# Function to find the design for continuous endpoints with specified family-wise error rate and power
finddesign_continuous <- function(J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas) {
  if (J < 2 || J > 10) stop("J must be an integer between 2 and 10")
  if (K < 2 || K > 100) stop("K must be an integer between 2 and 100")
  if (length(ns) != J) stop("ns must be a vector of length J")
  if (ns[1] != K || ns[length(ns)] != 1 || min(ns[-length(ns)] - ns[-1]) < 1)
    stop("ns must have first entry K, last entry 1 and each entry be strictly less than the previous one")
  if (requiredfwer < 1e-8 || requiredfwer > 1 - 1e-8) stop("requiredfwer must be strictly between 0 and 1")
  if (requiredpower < 1e-8 || requiredpower > 1 - 1e-8) stop("requiredpower must be strictly between 0 and 1")
  
  # Cumulative group sizes
  cumgroupsizes <- cumsum(groupsizes)
  cumgroupsizes <- cumgroupsizes / cumgroupsizes[1]
  
  # Initialize the covariance matrix
  sigma <- matrix(0, J * K, J * K)
  
  # Fill in the sigma blocks for continuous outcomes
  for (i in 1:J) {
    for (j in 1:J) {
      # Continuous outcome-specific covariance calculation
      sigma[i, j] <- sqrt(min(cumgroupsizes[i], cumgroupsizes[j]) / max(cumgroupsizes[i], cumgroupsizes[j]))
    }
  }
  
  # Populate the full covariance matrix
  for (i in 1:K) {
    # Diagonal blocks (covariance within the same treatment across stages)
    sigma[((i - 1) * J + 1):(i * J), ((i - 1) * J + 1):(i * J)] <- sigma[1:J, 1:J]
    
    # Off-diagonal blocks (covariance between different treatments)
    for (j in (1:K)[-i]) {
      sigma[((i - 1) * J + 1):(i * J), ((j - 1) * J + 1):(j * J)] <- sigma[1:J, 1:J] * 
        sqrt((treatmentsigmas[1]^4) / ((treatmentsigmas[i + 1]^2 + treatmentsigmas[1]^2) * (treatmentsigmas[j + 1]^2 + treatmentsigmas[1]^2)))
    }
  }
  
  # Construct the design matrix A
  rowsofA <- ns - 1
  rowsofA[length(rowsofA)] <- 1
  A <- matrix(0, sum(rowsofA), length(sigma[, 1]))
  
  treatments <- 1:K
  whichstagedropped <- rep(1:J, times = c(ns[1:(J - 1)] - ns[-1], 1))
  
  tempint <- 0
  for (i in 1:(J - 1)) {
    treatments_thisstage <- treatments[which(whichstagedropped == i)]
    treatments_futurestage <- treatments[which(whichstagedropped > i)]
    
    # Each treatment in future stages must beat the last entry in current stage treatments
    for (k1 in treatments_futurestage) {
      tempint <- tempint + 1
      A[tempint, (treatments_thisstage[length(treatments_thisstage)] - 1) * J + i] <- -1
      A[tempint, (k1 - 1) * J + i] <- 1
    }
    
    # Each treatment in the current stage must beat the one below it
    if (length(treatments_thisstage) > 1) {
      for (k2 in 1:(length(treatments_thisstage) - 1)) {
        tempint <- tempint + 1
        A[tempint, (treatments_thisstage[k2] - 1) * J + i] <- -1
        A[tempint, (treatments_thisstage[k2 + 1] - 1) * J + i] <- 1
      }
    }
  }
  
  # The last entry in A corresponds to the final treatment surviving all stages
  A[length(A[, 1]), length(A[1, ])] <- 1
  
  # Get the mean vector under the null hypothesis
  mu <- rep(0, length(sigma[, 1]))
  mean <- as.double(A %*% mu)
  var <- A %*% sigma %*% t(A)
  
  # Find the critical value c for the required FWER
  c <- uniroot(integral_typeIerror, lower = -2, upper = 5, requiredtypeIerror = requiredfwer, mean = mean, var = var, combinations = combinations, K = K)$root
  
  # Find the sample size for the given power
  n <- uniroot(integral_power, lower = 0, upper = 2000, requiredpower = requiredpower, c = c, A = A, var = var, combinations = combinations, J = J, K = K, delta1 = delta1, delta0 = delta0, cumgroupsizes = cumgroupsizes)$root
  
  # Print the covariance matrix for verification
  print(sigma)
  
  return(list(n = n, c = c, totalSS = n * sum((ns + 1) * groupsizes / groupsizes[1])))
}

# -------------------------------------------------------------------------
