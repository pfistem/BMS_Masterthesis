# Drop the loser for MAMS (Binary Outcome)
# Bristol Myers Squibb
# Contact: manuel.pfister@bms.com
# Date: 2024-07-17

# -------------------------------------------------------------------------
library(mvtnorm)

# integral_typeIerror is a function called by finddesign that finds the difference 
# between the FWER of a given design and the required FWER for binary outcomes.
integral_typeIerror = function(c, requiredtypeIerror, mean, var, combinations, K) {
  lower = rep(0, length(mean))
  lower[length(lower)] = c
  
  upper = rep(Inf, length(mean))
  
  int = pmvnorm(lower=lower, upper=upper, mean=mean, sigma=var)
  
  return(as.double(int) * factorial(K) - requiredtypeIerror)
}

# -------------------------------------------------------------------------
# integral_power is a function called by finddesign that finds the difference between 
# the power under the LFC of a given design and the required power for binary outcomes.
integral_power = function(n, requiredpower, c, A, var, combinations, J, K, p_control, p_treatment) {
  # get mu under LFC
  mu = rep(0, length(A[1, ]))
  for (i in 1:(K-1)) {
    mu[((i-1)*J+1):(i*J)] = sqrt(n * p_control * (1 - p_control)) * log(p_treatment[i] / p_control)
  }
  
  mu[((K-1)*J+1):(K*J)] = sqrt(n * p_control * (1 - p_control)) * log(p_treatment[K] / p_control)
  
  mean = as.double(A %*% mu)
  
  lower = rep(0, length(mean))
  lower[length(lower)] = c
  
  upper = rep(Inf, length(mean))
  
  int = pmvnorm(lower=lower, upper=upper, mean=mean, sigma=var)
  
  return(as.double(int) * factorial(K-1) - requiredpower)
}

# finddesign returns a design with specified family-wise error rate and power. 
# Arguments:
# J - number of stages
# K - number of experimental arms
# ns - vector of length J stating the number of treatments at each stage - 
#      first entry must be K and last entry must be 1, and entries must be strictly decreasing
# requiredfwer - required family-wise error rate at the global null
# requiredpower - required power at the least favorable configuration
# groupsizes - relative proportion of patients recruited per arm at each stage 
# p_control - probability of success in the control group
# p_treatment - vector of probabilities of success for each treatment arm

# -------------------------------------------------------------------------
# output is a list with entries n, c and totalSS. 
# n is the required group-size, c is the required critical value and totalSS is the total sample size used by the design.
finddesign = function(J, K, ns, requiredfwer, requiredpower, groupsizes, p_control, p_treatment) {
  
  if (J < 2 || J > 10) {
    stop("J must be an integer between 2 and 10")
  }
  
  if (K < 2 || K > 100) {
    stop("K must be an integer between 2 and 100")
  }
  
  # test ns:
  if (length(ns) != J) {
    stop("ns must be a vector of length J")
  }
  
  if (ns[1] != K || ns[length(ns)] != 1 || min(ns[-length(ns)] - ns[-1]) < 1) {
    stop("ns must have first entry K, last entry 1 and each entry be strictly less than the previous one")
  }
  if (requiredfwer < 1e-8 || requiredfwer > 1-1e-8) {
    stop("requiredfwer must be strictly between 0 and 1")
  }
  
  if (requiredpower < 1e-8 || requiredpower > 1-1e-8) {
    stop("requiredpower must be strictly between 0 and 1")
  }
  
  cumgroupsizes = cumsum(groupsizes)
  
  # normalize so that first stage is 1:
  cumgroupsizes = cumgroupsizes / cumgroupsizes[1]
  
  sigma = matrix(0, J * K, J * K)
  
  # fill in sigma blocks for binary outcomes
  for (i in 1:K) {
    for (j in 1:K) {
      if (i == j) {
        for (m in 1:J) {
          for (n in 1:J) {
            sigma[(i-1)*J + m, (j-1)*J + n] = min(cumgroupsizes[m], cumgroupsizes[n]) / sqrt(cumgroupsizes[m] * cumgroupsizes[n])
          }
        }
      } else {
        p1 = p_treatment[i]
        p2 = p_treatment[j]
        cov_factor = sqrt((p_control * (1 - p_control)) / ((p1 * (1 - p1)) * (p2 * (1 - p2))))
        for (m in 1:J) {
          for (n in 1:J) {
            sigma[(i-1)*J + m, (j-1)*J + n] = min(cumgroupsizes[m], cumgroupsizes[n]) / sqrt(cumgroupsizes[m] * cumgroupsizes[n]) * cov_factor
          }
        }
      }
    }
  }
  
  print(sigma)
  
  rowsofA = ns - 1
  rowsofA[length(rowsofA)] = 1
  
  A = matrix(0, sum(rowsofA), length(sigma[, 1]))
  
  treatments = 1:K
  whichstagedropped = rep(1:J, times=c(ns[1:(J-1)] - ns[-1], 1))
  
  # go through each stage; each treatment that drops out should have a row in A together with each treatment dropping out in a subsequent stage
  tempint = 0
  
  for (i in 1:(J-1)) {
    treatments_thisstage = treatments[which(whichstagedropped == i)]
    treatments_futurestage = treatments[which(whichstagedropped > i)]
    
    # each treatment in treatments_futurestage must beat last entry in treatments_thistage
    for (k1 in treatments_futurestage) {
      tempint = tempint + 1
      A[tempint, (treatments_thisstage[length(treatments_thisstage)]-1)*J+i] = -1
      A[tempint, (k1-1)*J+i] = 1
    }
    
    # each treatment in treatments_thisstage must beat the one below it:
    # check if more than one treatment is dropped
    if (length(treatments_thisstage) > 1) {
      for (k2 in 1:(length(treatments_thisstage)-1)) {
        tempint = tempint + 1
        A[tempint, (treatments_thisstage[k2]-1)*J+i] = -1
        A[tempint, (treatments_thisstage[k2+1]-1)*J+i] = 1
      }
    }
  }
  
  A[length(A[, 1]), length(A[1, ])] = 1
  
  # get mean vector under HG:
  mu = rep(0, length(sigma[, 1]))
  mean = as.double(A %*% mu)
  
  var = A %*% sigma %*% t(A)
  
  c = uniroot(integral_typeIerror, lower=-2, upper=5, requiredtypeIerror=requiredfwer, mean=mean, var=var, combinations=combinations, K=K)$root
  
  # find sample size for given power power
  n = uniroot(integral_power, lower=0, upper=2000, requiredpower=requiredpower, c=c, A=A, var=var, combinations=combinations, J=J, K=K, p_control=p_control, p_treatment=p_treatment)$root
  
  return(list(n=n, c=c, totalSS=n * sum((ns + 1) * groupsizes / groupsizes[1])))
}

# -------------------------------------------------------------------------

# Example usage:
finddesign(J=3, K=4, ns=c(4, 2, 1), requiredfwer=0.05, requiredpower=0.9, 
           groupsizes=c(1, 1, 1), p_control=0.5, p_treatment=c(0.5, 0.6, 0.7, 0.8))









