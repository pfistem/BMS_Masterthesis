# Drop the loser for MAMS | MAIN FUNCTION
# Bristol Myers Squibb
# Contact: manuel.pfister@bms.com

# -------------------------------------------------------------------------

# finddesign_main.R

finddesign <- function(endpoint_type, J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas, events = NULL) {
  # Load the appropriate design function based on endpoint_type
  if (endpoint_type == "continuous") {
    source("continuous_design.R")
    return(finddesign_continuous(J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas))
  } else if (endpoint_type == "binary") {
    source("binary_design.R")
    return(finddesign_binary(J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas))
  } else if (endpoint_type == "survival") {
    source("survival_design.R")
    return(finddesign_survival(J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas, events))
  } else {
    stop("Invalid endpoint_type. Choose from 'continuous', 'binary', or 'survival'.")
  }
}

# Example usage:
# Continuous:
finddesign("continuous", J = 3, K = 4, ns = c(4, 2, 1), delta1 = 0.545, delta0 = 0.178, requiredfwer = 0.05, requiredpower = 0.9, treatmentsigmas = c(1, 1, 1, 100, 1), groupsizes = c(1, 1, 1))

# Survival:
events <- matrix(c(30, 20, 50, 30, 40, 25), ncol = 2, byrow = TRUE)
finddesign("survival", J = 3, K = 4, ns = c(4, 2, 1), delta1 = 0.545, delta0 = 0.178, requiredfwer = 0.05, requiredpower = 0.9, treatmentsigmas = c(1, 1, 1, 100, 1), groupsizes = c(1, 1, 1), events = events)

# Binary (assuming similar structure to continuous):
# finddesign("binary", J = 3, K = 4, ns = c(4, 2, 1), delta1 = 0.545, delta0 = 0.178, requiredfwer = 0.05, requiredpower = 0.9, treatmentsigmas = c(1, 1, 1, 100, 1), groupsizes = c(1, 1, 1))




