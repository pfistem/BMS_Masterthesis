# Drop the loser for MAMS | Survival Endpoint
# Bristol Myers Squibb
# Contact: manuel.pfister@bms.com
# Date: 2024-12-09

#' Design Calculation for Multi-Arm Multi-Stage (MAMS) Survival Trials
#'
#' This function calculates the required design parameters (sample size, critical value)
#' for multi-arm multi-stage (MAMS) trials with survival endpoints. The design
#' controls the family-wise error rate (FWER) and achieves the desired statistical power.
#'
#' @param J Integer. Number of stages in the trial (must be between 2 and 10).
#' @param K Integer. Number of treatment arms (must be between 2 and 100).
#' @param ns Integer vector. Number of treatments remaining at each stage. The first entry
#'   must equal `K`, the last entry must be `1`, and entries must decrease strictly.
#' @param delta1 Numeric. Expected effect size (log hazard ratio) for the effective treatment
#'   under the alternative hypothesis.
#' @param delta0 Numeric. Expected effect size (log hazard ratio) for ineffective treatments
#'   under the null hypothesis.
#' @param requiredfwer Numeric. Desired family-wise error rate (between 0 and 1).
#' @param requiredpower Numeric. Desired statistical power (between 0 and 1).
#' @param groupsizes Numeric vector. Relative sizes of treatment groups at each stage.
#' @param events Matrix. A matrix of cumulative event counts for each arm at each stage.
#'   The matrix should have `K` rows (one per treatment arm) and `J` columns (one per stage).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{n}{Numeric. Required number of events for the study.}
#'   \item{c}{Numeric. Critical value for the test statistics.}
#'   \item{totalSS}{Numeric. Total sample size across all stages and arms.}
#' }
#'
#' @details
#' The function first validates the input parameters, then constructs a covariance matrix
#' to model the correlations between event counts across stages and treatment arms.
#' It uses a design matrix to enforce constraints at each stage, ensuring proper
#' elimination of inferior treatments while retaining the best-performing option.
#'
#' The critical value (`c`) for the FWER is calculated using numerical integration
#' of the multivariate normal distribution. Similarly, the required number of events
#' for achieving the desired power is determined through root-finding techniques.
#'
#' The function is specific to survival endpoints, incorporating survival-specific
#' covariance calculations based on cumulative event counts.
#'
#' @examples
#' events <- matrix(c(
#'   20, 40, 60,  # Events for treatment arm 1 at stages 1, 2, and 3
#'   30, 50, 100, # Events for treatment arm 2
#'   30, 50, 100, # Events for treatment arm 3
#'   30, 50, 100  # Events for treatment arm 4
#' ), ncol = 3, byrow = TRUE)
#'
#' finddesign_survival(
#'   J = 3,
#'   K = 4,
#'   ns = c(4, 2, 1),
#'   delta1 = -log(0.7),
#'   delta0 = 0,
#'   requiredfwer = 0.05,
#'   requiredpower = 0.8,
#'   groupsizes = c(1, 2/3, 5/3),
#'   events = events
#' )
#'
#' @seealso \code{\link[mvtnorm]{pmvnorm}}, \code{\link[stats]{uniroot}}
#'
#' @export

library(mvtnorm)

# Function to calculate the integral of the Type I error
# Arguments:
# - c: critical value to be determined
# - requiredtypeIerror: desired family-wise error rate
# - mean: mean vector under the null hypothesis
# - var: covariance matrix
# - K: number of treatment arms
integral_typeIerror <- function(c, requiredtypeIerror, mean, var, K, ...) {
  lower <- rep(0, length(mean))
  lower[length(lower)] <- c
  upper <- rep(Inf, length(mean))
  int <- pmvnorm(lower = lower, upper = upper, mean = mean, sigma = var)
  return(as.double(int) * factorial(K) - requiredtypeIerror)
}

integral_power <- function(max_events, events, requiredpower, c, A, var, J, K, delta1, delta0, cumgroupsizes, ...) {
  # max_events is the number of events planned at FA for Control vs 1 experimental arm

  # Get mu under the Least Favorable Configuration (LFC)
  mu <- rep(0, length(A[1, ]))

  # Modify mu to depend on the number of events instead of sample size  ### Update subscripts
  for (i in 2:K) {
    mu[((i - 1) * J + 1):(i * J)] <- sqrt(max_events * events[1,]/max(events[1,]) / 4) * delta0  ### Update 2 to 4
  }
  mu[1:J] <- sqrt(max_events * events[1,]/max(events[1,]) / 4) * delta1  ### Update events to events[1,]

  # Adjust mean based on design matrix A
  mean <- as.double(A %*% mu)

  # Define integration bounds
  lower <- rep(0, length(mean))
  lower[length(lower)] <- c
  upper <- rep(Inf, length(mean))

  # Compute the integral using pmvnorm
  int <- pmvnorm(lower = lower, upper = upper, mean = mean, sigma = var)

  # Return the result adjusted by the number of treatment arms
  return(as.double(int) * factorial(K - 1) - requiredpower)
}

finddesign_survival <- function(J, K, ns, delta1, delta0, requiredfwer,
                                requiredpower, groupsizes, events, ...) {
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

  # Fill in the sigma blocks using the survival covariance formula
  for (k in 1:K){
    for (i in 1:J) {
      for (j in 1:J) {
        # Survival-specific covariance calculation
        sigma[i+J*(k-1), j+J*(k-1)] <- sqrt(ifelse(i < j, events[1, i] / events[1, j], events[1, j] / events[1, i]))
      }
    }
  }

  # Correlation between different arms
  for (k1 in 1:K){
    for(k2 in 1:K){
      if(k1 != k2){
        for (i in 1:J) {
          for (j in 1:J) {
            sigma[i+J*(k1-1), j+J*(k2-1)] <- 0.5 * sqrt(ifelse(i < j, events[1, i] / events[1, j], events[1, j] / events[1, i]))
          }
        }
      }
    }
  }

  # Print the covariance matrix for verification
  print("Covariance Matrix (sigma):")
  print(sigma)

  # Construct the design matrix A
  rowsofA <- ns - 1
  rowsofA[length(rowsofA)] <- 1
  A <- matrix(0, sum(rowsofA), length(sigma[, 1]))

  treatments <- 1:K
  whichstagedropped <- c(3, 2, 1, 1)  # Reverted order as suggested for clarity

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
  A[length(A[, 1]), J] <- 1  # Adjusted as per the comment



  # Get the mean vector under the null hypothesis
  mu <- rep(0, length(sigma[, 1]))
  mean <- as.double(A %*% mu)
  var <- A %*% sigma %*% t(A)

  print(var)
  print(A)
  print(J)
  print(K)
  print(cumgroupsizes)
  print(delta1)
  print(delta0)
  print(requiredpower)

  # Find the critical value c for the required FWER
  c <- uniroot(integral_typeIerror, lower = -2, upper = 5, requiredtypeIerror = requiredfwer, mean = mean, var = var, K = K, tol = 1e-6)$root
  print(c)

  # Find the number of events for the given power
  n <- uniroot(integral_power, lower = 0, upper = 3000, events = events, requiredpower = requiredpower, c = c, A = A, var = var, J = J, K = K, delta1 = delta1, delta0 = delta0, cumgroupsizes = cumgroupsizes, tol = 1e-6)$root

  # Calculate the total number of events using the control randomization weight instead of hardcoded 1
  total_sample_size <- n * sum(((ns + sqrt(K)) * groupsizes) / groupsizes[1])

  # Final stage
  total_sample_size <- n
  # Interim stages
  for(i in (length(ns)-1):1){

    total_sample_size <- total_sample_size + (ns[i]-ns[i+1])*(n * events[1,i]/max(events[1,]))/2

  }

  return(list(n = n, c = c, totalSS = total_sample_size))
}

# -------------------------------------------------------------------------


events <- matrix(c(
  20, 40, 60,   # Events for treatment arm 1 at stages 1, 2, and 3
  30, 50, 100,   # Events for treatment arm 2
  30, 50, 100,   # Events for treatment arm 3
  30, 50, 100    # Events for treatment arm 4
), ncol = 3, byrow = TRUE)

finddesign_survival(J = 3, K = 4, ns = c(4, 2, 1), delta1 = -log(0.7), delta0 = 0,
                    requiredfwer = 0.05, requiredpower = 0.8,
                    groupsizes = c(1, 2/3, 5/3), events = events)






