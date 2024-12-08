# Drop the loser for MAMS | Continuous Endpoint
# Bristol Myers Squibb
# Contact: manuel.pfister@bms.com
# Date: 2024-12-09

# -------------------------------------------------------------------------

#' Find Design for MAMS with Continuous Endpoints
#'
#' This function calculates the required sample size and critical thresholds for
#' a multi-arm, multi-stage (MAMS) clinical trial with continuous endpoints. The
#' design ensures control of the family-wise error rate (FWER) while achieving
#' the desired statistical power.
#'
#' @param J Integer. The number of stages in the trial (must be between 2 and 10).
#' @param K Integer. The number of treatment arms (must be between 2 and 100).
#' @param ns Integer vector. A vector of length `J` specifying the number of
#'   arms remaining at each stage. The first entry must be `K`, the last entry must
#'   be `1`, and each entry must be strictly less than the previous one.
#' @param delta1 Numeric. The treatment effect for the best treatment.
#' @param delta0 Numeric. The treatment effect under the null hypothesis.
#' @param requiredfwer Numeric. The family-wise error rate to control (between 0 and 1).
#' @param requiredpower Numeric. The desired statistical power (between 0 and 1).
#' @param groupsizes Numeric vector. The sizes of the groups in each stage.
#' @param treatmentsigmas Numeric vector. The standard deviations of the treatment effects
#'   for each arm.
#'
#' @return A list containing:
#'   \item{n}{The required sample size per group.}
#'   \item{c}{The critical value for the test statistic.}
#'   \item{totalSS}{The total sample size across all groups and stages.}
#'
#' @details
#' This function constructs a covariance matrix for the continuous outcomes, generates
#' a design matrix (`A`) to enforce constraints at each stage, and calculates the
#' critical value to control the FWER and the sample size needed to achieve the desired power.
#'
#' The function utilizes numerical integration via the \code{\link[mvtnorm]{pmvnorm}} function
#' to compute probabilities under the multivariate normal distribution and root-finding
#' methods (\code{\link[stats]{uniroot}}) to solve for the critical value and sample size.
#'
#' @examples
#' finddesign_continuous(
#'   J = 3,
#'   K = 4,
#'   ns = c(4, 2, 1),
#'   delta1 = 0.545,
#'   delta0 = 0.178,
#'   requiredfwer = 0.05,
#'   requiredpower = 0.9,
#'   treatmentsigmas = c(1, 1, 1, 100, 1),
#'   groupsizes = c(20, 40, 60)
#' )
#'
#' @seealso \code{\link[mvtnorm]{pmvnorm}}, \code{\link[stats]{uniroot}}
#'
#' @export


library(mvtnorm)

# Function to calculate the integral of the Type I error
integral_typeIerror <- function(c, requiredtypeIerror, mean, var, K) {
  lower <- rep(0, length(mean))
  lower[length(lower)] <- c
  upper <- rep(Inf, length(mean))
  int <- pmvnorm(lower = lower, upper = upper, mean = mean, sigma = var)
  return(as.double(int) * factorial(K) - requiredtypeIerror)
}

# Function to calculate the integral of the power
integral_power <- function(n, requiredpower, c, A, var, J, K, delta1, delta0, cumgroupsizes) {
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
  print(A)

  # Find the critical value c for the required FWER
  c <- uniroot(integral_typeIerror, lower = -2, upper = 5, requiredtypeIerror = requiredfwer, mean = mean, var = var, K = K)$root

  # Find the sample size for the given power
  n <- uniroot(integral_power, lower = 0, upper = 2000, requiredpower = requiredpower, c = c, A = A, var = var, J = J, K = K, delta1 = delta1, delta0 = delta0, cumgroupsizes = cumgroupsizes)$root

  # Print the covariance matrix for verification
  print(sigma)
  print(var)
  return(list(n = n, c = c, totalSS = n * sum((ns + 1) * groupsizes / groupsizes[1])))
}

finddesign_continuous(J = 3,
           K = 4,
           ns = c(4, 2, 1),
           delta1 = 0.545,
           delta0 = 0.178,
           requiredfwer = 0.05,
           requiredpower = 0.9,
           treatmentsigmas = c(1, 1, 1, 100, 1),
           groupsizes = c(20, 40, 60))



