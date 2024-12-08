# Drop the Loser for Multi-Arm Multi-Stage (MAMS) Designs | MAIN FUNCTION
# Bristol Myers Squibb
# Contact: manuel.pfister@bms.com
# Date: 2024-12-09

# This script provides a main function 'finddesign' to compute the design parameters for MAMS trials
# with continuous or survival endpoints.
# It utilizes the 'rpact' package for group sequential design calculations.

# -------------------------------------------------------------------------

if (!requireNamespace("rpact", quietly = TRUE)) {
  install.packages("rpact")
}
library(rpact)

# Main function to find the design based on the type of endpoint
finddesign <- function(endpoint_type, J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas, events = NULL) {

  # Depending on the endpoint type, load the appropriate design function
  if (endpoint_type == "continuous") {
    # Source the function for continuous endpoints
    source("code/continuous_design.R")
    # Call the function to find the design for continuous endpoints
    return(finddesign_continuous(J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, treatmentsigmas))

  } else if (endpoint_type == "survival") {
    # Source the function for survival endpoints
    source("code/survival_design.R")
    # Call the function to find the design for survival endpoints
    return(finddesign_survival(J, K, ns, delta1, delta0, requiredfwer, requiredpower, groupsizes, events))

  } else {
    # If an invalid endpoint type is provided, stop with an error message
    stop("Invalid endpoint_type. Choose from 'continuous' or 'survival'.")
  }
}


# Example usage for Survival endpoint:
events <- matrix(c(
  30, 50, 100,   # Events for treatment arm 1 at stages 1, 2, and 3
  30, 50, 100,   # Events for treatment arm 2
  30, 50, 100,   # Events for treatment arm 3
  30, 50, 100    # Events for treatment arm 4
), ncol = 3, byrow = TRUE)



# Compute design parameters for a survival endpoint
finddesign("survival",
           J = 3,                      # Number of stages
           K = 4,                      # Number of treatment arms
           ns = c(4, 2, 1),            # Sample sizes per arm at each stage
           delta1 = -log(0.7),         # Effect size under the alternative hypothesis
           delta0 = 0,                 # Effect size under the null hypothesis
           requiredfwer = 0.025,       # Alpha spending at each stage
           requiredpower = 0.9,        # Required power
           groupsizes = c(1, 2/3, 5/3),# Group sizes at each stage
           events = events)            # Matrix of events per arm per stage

# Example usage for Continuous endpoint:
finddesign("continuous",
           J = 3,
           K = 4,
           ns = c(4, 2, 1),
           delta1 = 0.545,
           delta0 = 0.178,
           requiredfwer = 0.05,       # Overall family-wise error rate
           requiredpower = 0.9,
           treatmentsigmas = c(1, 1, 1, 100, 1),
           groupsizes = c(1, 1, 1))

# -------------------------------------------------------------------------
# Extension Group Sequential Design

# Define the Pocock spending function
# This function returns the alpha spending at each stage based on the Pocock type of group sequential design
spending_function <- function(J, alpha, type, events) {
  # Get the group sequential design using rpact
  design <- getDesignGroupSequential(typeOfDesign = type, alpha = alpha, kMax = J, informationRates = events[1,] / events[1, 3])
  # Return the alpha spent at each stage
  return(design$alphaSpent)
}

# Obtain the alpha spending function using Pocock design for J=3 stages
alpha_spending <- spending_function(J = 3, alpha = 0.025, type = "asOF", events = events)



