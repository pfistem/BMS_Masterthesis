#' Survival Design Simulation for Drop-the-Losers Trials
#'
#' This script simulates oncology trials using a drop-the-losers design with survival endpoints.
#' It calculates trial outcomes based on various scenarios for progression-free survival (PFS)
#' and tests hypotheses using log-rank tests. The simulations allow for the assessment of design
#' parameters, including recruitment rates, dropout, and effect sizes.
#'
#' @title Survival Design Simulation
#' @subtitle Drop-the-losers Design in Oncology Trials
#' @author Manuel Pfister (\email{manuel.pfister@bms.com})
#' @date Last modified on 12 Nov 2024
#'
#' @section Design Parameters:
#' - **Censoring Settings**: Dropout rate is set to 10% (\code{dropout_PFS = 0.10}).
#' - **Recruitment**: Recruitment rate of 10 patients per month (\code{RecruitmentRate = 10}).
#' - **Block Randomization**: Block size is set to 6 (\code{block_size = 6}).
#' - **Median PFS**: Control arm median PFS is 8 months (\code{medianPFS_control = 8}).
#'
#' @section Simulation Parameters:
#' - **Scenarios**: Six scenarios based on hazard ratio configurations:
#'   - Null Hypothesis: \code{HR = c(1, 1, 1, 1)}.
#'   - Hybrid 1: \code{HR = c(0.67, 1, 1, 1)}.
#'   - Hybrid 2: \code{HR = c(0.67, 0.67, 1, 1)}.
#'   - Hybrid 1 + 1: \code{HR = c(0.67, 0.80, 1, 1)}.
#'   - Harmful Drug: \code{HR = c(1.2, 1, 1, 1)}.
#'   - All Effective: \code{HR = c(0.67, 0.67, 0.67, 0.67)}.
#' - **Events**: Stages defined by cumulative event thresholds:
#'   - Stage 1: 30 events.
#'   - Stage 2: 50 events.
#'   - Stage 3: 100 events.
#'
#' @section Simulation Details:
#' The simulation uses parallel processing via \code{foreach} and \code{doParallel} to execute
#' up to \code{Nsim = 200000} iterations across six predefined scenarios. Each iteration models
#' patient recruitment, treatment randomization, and PFS using survival analysis techniques.
#'
#' @section Outputs:
#' - **Results_sim**: Data frame containing simulation results for each scenario and iteration.
#' - **Saved Output**: Simulation results are saved as an \code{RData} file in the \code{Simulation/} directory
#'   with a timestamped filename.
#'
#' @section Examples:
#' To run the script, source it in R after ensuring required packages are installed:
#' \dontrun{
#' # Run the script
#' source("survival_design_simulation.R")
#'
#' # Load results
#' load("Simulation/Sim_2024-11-12_12-00.RData")
#' }
#'
#' @section Parallel Processing:
#' The script utilizes parallel processing for simulations:
#' - Number of cores: All available cores minus one (\code{detectCores() - 1}).
#' - Cluster setup: \code{makeCluster} and \code{registerDoParallel}.
#' - Simulation execution: \code{foreach} loop with \code{%dopar%}.
#' - Cluster cleanup: \code{stopCluster}.
#'
#' @note Ensure the external function file \code{simulation_function_with_futility_and_log_rank.R}
#' is available in the specified path.
#'
#' @dependencies
#' - \code{rpact}
#' - \code{survival}
#' - \code{mvtnorm}
#' - \code{doParallel}
#' - \code{foreach}

set.seed(7)

### Packages
library(rpact)
library(survival)
library(mvtnorm)
library(doParallel)
library(foreach)


# Source external function files
#source("finddesign_main.R")
source("code/simulation_function_with_futility_and_log_rank.R")


### Design parameters

# Censoring settings
dropout_PFS <- 0.10  # 10% dropout rate

# Recruitment rate and block size for randomization
RecruitmentRate <- 10  # Recruitment rate per patient per month
block_size <- 6  # Block size for treatment randomization


### PFS endpoints

# Clinical assumptions based on protocol scenarios
medianPFS_control <- 8    # Control arm: 8 months

## Not used: see Hazard Ratio Scenarios
# medianPFS_not_meaningful <- 10  # Treatment: 10 months
# medianPFS_good_effect <- 12    # Treatment: 12 months (good effect)


### Design

# Drop-the-losers design with critical threshold for final analysis
events <- matrix(c(30, 50,100,    # Stage I: 50 events, Stage II: 70 events, Stage III: 90 events
                   30, 50,100,    # Treatment 2
                   30, 50,100,    # Treatment 3
                   30, 50,100),   # Treatment 4
                 ncol = 3, byrow = TRUE)

### Simulations

# Number of simulated studies per scenario
Nsim <- 200000
K = 4

# Scenarios: Simplified based on the four cases outlined in the protocol
Scenarios <- matrix(
  c(1, 1, 1, 1,                # Null Hypothesis: No effect
    0.67, 1, 1, 1,             # Hybrid 1: One effective drug
    0.67, 0.67, 1, 1,          # Hybrid 2: Two effective drugs
    0.67, 0.80, 1, 1,          # Hybrid 1 + 1: One effective, one half-effect drug
    1.2, 1, 1, 1,              # Harmful drug
    0.67, 0.67, 0.67, 0.67),   # All treatments are effect
  ncol = 4,
  byrow = TRUE,
  dimnames = list(
    c("Null hypothesis", "Hybrid 1", "Hybrid 2", "Hybrid 1 + 1", "Harmful", "All effectiv"),
    NULL
  )
)

# Results table: Record every simulation
Results <- expand.grid(
  Scenario = 1:nrow(Scenarios),
  Sim = 1,
  stringsAsFactors = FALSE
)

Results$HR_1 <- Scenarios[Results$Scenario, 1]
Results$HR_2 <- Scenarios[Results$Scenario, 2]
Results$HR_3 <- Scenarios[Results$Scenario, 3]
Results$HR_4 <- Scenarios[Results$Scenario, 4]

### Parallel Simulation

# Register parallel backend
numCores <- detectCores() - 1  # Use one less than the available cores
cl <- makeCluster(numCores)
registerDoParallel(cl)

#Parallel simulation loop using foreach
Results_sim <- foreach(i = 1:Nsim, .combine = 'rbind', .packages = c('survival'), .inorder = FALSE) %dopar% {
  Results$Sim <- i
  simulate_trial(i, thresholds = c(30, 50, 100),
                 Results = Results, RecruitmentRate = 15, patients_stage1 = 150, patients_stage2 = 150,
                 block_size = 6, K = 4,
                 medianPFS_control = 8, dropout_PFS = 0.1)
}

# For Log Rank Test
Results_sim <- foreach(i = 1:Nsim, .combine = 'rbind', .packages = c('survival'), .inorder = FALSE) %dopar% {
  Results$Sim <- i
  simulate_trial_log_rank(i, thresholds = c(30, 50, 100),
                 Results = Results, RecruitmentRate = 15, patients_stage1 = 150, patients_stage2 = 150,
                 block_size = 6, K = 4,
                 medianPFS_control = 8, dropout_PFS = 0.1)
}

# Stop the cluster after the simulation
stopCluster(cl)

# Store results back into Results dataframe
# final_results <- merge(Results_sim, Results, by = "Scenario", all.x = TRUE)

# Save the simulation results with date and time up to minutes
save(Results_sim, file = paste0("Simulation/Sim_", format(Sys.time(), "%Y-%m-%d_%H-%M"), ".RData"))

