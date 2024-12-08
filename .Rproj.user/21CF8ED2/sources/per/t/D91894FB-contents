####################################################################
### Title: Survival Design Simulation: Drop-the-losers Design in Oncology Trials
### Created by: Manuel Pfister (manuel.pfister@bms.com)
### Last modified on: 21 Nov 2024
####################################################################


#' Simulate a Clinical Trial with Futility and Drop-the-Loser Analyses
#'
#' This function simulates clinical trial scenarios to evaluate progression-free survival (PFS)
#' outcomes with interim and final analyses. The analyses include futility and "drop-the-loser"
#' strategies to optimize trial design.
#'
#' @param i Simulation iteration number.
#' @param thresholds Vector of event thresholds for interim and final analyses.
#' @param patients_stage1 Number of patients recruited in stage 1.
#' @param patients_stage2 Number of patients recruited in stage 2.
#' @param Results Data frame containing hazard ratios (HR) for each treatment arm in different scenarios.
#' @param RecruitmentRate Recruitment rate for patients in the trial.
#' @param block_size Block size for randomizing patients into treatment arms.
#' @param K Number of treatment arms (including the control arm).
#' @param medianPFS_control Median progression-free survival (PFS) time in the control group.
#' @param dropout_PFS Monthly dropout probability for PFS analysis.
#' @return A data frame containing futility analysis results for all simulated scenarios.
#' @examples
#' simulate_trial(1, c(50, 100), 100, 200, Results, 0.1, 4, 3, 6, 0.05)
#' @export
simulate_trial <- function(i,
                           thresholds,
                           patients_stage1,
                           patients_stage2,
                           Results,
                           RecruitmentRate,
                           block_size,
                           K,
                           medianPFS_control,
                           dropout_PFS) {
  # Check if the input parameters are valid, otherwise throw an error
  if (!is.numeric(i) ||
      i <= 0)
    stop("Simulation iteration 'i' must be a positive integer.")
  if (!is.numeric(thresholds) ||
      any(thresholds <= 0))
    stop("Thresholds must be a numeric vector with positive values.")
  if (!is.numeric(patients_stage1) ||
      patients_stage1 <= 0)
    stop("patients_stage1 must be a positive integer.")
  if (!is.numeric(patients_stage2) ||
      patients_stage2 <= 0)
    stop("patients_stage2 must be a positive integer.")
  if (!is.data.frame(Results))
    stop("Results must be a data frame.")
  if (!is.numeric(RecruitmentRate) ||
      RecruitmentRate <= 0)
    stop("RecruitmentRate must be a positive number.")
  if (!is.numeric(block_size) ||
      block_size <= 0)
    stop("block_size must be a positive integer.")
  if (!is.numeric(K) ||
      K <= 0)
    stop("K must be a positive integer.")
  if (!is.numeric(medianPFS_control) ||
      medianPFS_control <= 0)
    stop("medianPFS_control must be a positive number.")
  if (!is.numeric(dropout_PFS) ||
      dropout_PFS < 0 ||
      dropout_PFS > 1)
    stop("dropout_PFS must be a number between 0 and 1.")

  # Initialize a data frame to store results of the trial
  futility_results <- data.frame(
    Scenario = integer(),
    # Scenario identifier
    TreatmentArm = integer(),
    # ID of the treatment arm being analyzed
    Simulation = integer(),
    # Simulation run number
    AnalysisType = character(),
    # Type of analysis (interim or final)
    Time_Futility = numeric(),
    # Time when futility analysis occurs
    Initial_Status = character(),
    # Status of the arm before analysis (e.g., "Continue")
    Futility_Outcome = character(),
    # Result of the analysis (e.g., "Dropped due to Futility")
    HR = numeric(),
    # Hazard ratio for PFS comparison
    TestStatistic = numeric(),
    # Test statistic from Cox regression
    PValue = numeric(),
    # P-value from the test
    Dropped = logical(),
    # Indicates if this arm was dropped
    stringsAsFactors = FALSE
  )

  # Loop over all scenarios in the Results data frame
  for (z in 1:nrow(Results)) {
    # Reset cumulative data for recruitment, treatment, PFS, and censoring
    cumulative_recruit <- c()
    cumulative_trt <- c()
    cumulative_PFS <- c()
    cumulative_cens_PFS <- c()
    futility_drops <- c()  # List of treatment arms dropped due to futility
    current_threshold_index <- 1  # Start with the first threshold for interim analysis

    # Perform interim analyses at each threshold
    while (current_threshold_index <= length(thresholds)) {
      required_events <- thresholds[current_threshold_index]  # Number of events required at this stage

      # Generate new patient recruitment times
      if (current_threshold_index == 1) {
        # Recruitment for stage 1
        new_recruit <- cumsum(rexp(patients_stage1 * 2, rate = RecruitmentRate))
      } else {
        # Recruitment for stage 2, starting at the time of IA1 plus an optional delay
        delay_months <- 2 # Default 2 months
        time_of_IA1 <- futility_results$Time_Futility[futility_results$AnalysisType == "Interim 1" &
                                                        futility_results$Simulation == i]

        if (length(time_of_IA1) > 0) {
          start_time_stage2 <- max(time_of_IA1) + delay_months
        } else {
          stop("Time of IA1 not found. Ensure IA1 results are correctly calculated.")
        }

        new_recruit <- cumsum(rexp(patients_stage2 * 2, rate = RecruitmentRate)) + start_time_stage2
      }

      # Combine new recruitment data with cumulative recruitment data
      cumulative_recruit <- c(cumulative_recruit, new_recruit)

      # Randomly assign patients to treatment arms, ensuring balance
      new_trt <- sample(rep(0:K, length.out = length(new_recruit)))
      cumulative_trt <- c(cumulative_trt, new_trt)

      # Generate PFS times for control and treatment arms
      new_PFS <- ifelse(
        new_trt == 0,
        rexp(length(new_recruit), log(2) / medianPFS_control),
        rexp(
          length(new_recruit),
          log(2) / medianPFS_control * as.numeric(Results[z, new_trt + 2])
        )
      )
      cumulative_PFS <- c(cumulative_PFS, new_PFS)

      # Generate censoring times based on dropout rate
      new_cens_PFS <- rexp(length(new_recruit), rate = -log(1 - dropout_PFS) / 12)
      cumulative_cens_PFS <- c(cumulative_cens_PFS, new_cens_PFS)

      # Prepare to analyze only the treatment arms that have not been dropped
      remaining_arms <- setdiff(1:K, futility_drops)

      # If all arms have been dropped, record the result and stop further analysis for this scenario
      if (length(remaining_arms) == 0) {
        futility_results <- rbind(
          futility_results,
          data.frame(
            Scenario = z,
            TreatmentArm = NA,
            Simulation = i,
            AnalysisType = ifelse(
              current_threshold_index == length(thresholds),
              "Final Analysis",
              paste("Interim", current_threshold_index)
            ),
            Time_Futility = NA,
            Initial_Status = "NA",
            Futility_Outcome = "No arms available due to Futility",
            HR = NA,
            TestStatistic = NA,
            PValue = NA,
            Dropped = NA
          )
        )
        break
      }

      # Data frame to temporarily store results for this threshold
      temp_results <- data.frame()

      # Loop through each remaining treatment arm
      for (k in remaining_arms) {
        # Filter data for the control and current treatment arm
        ad_recruit <- cumulative_recruit[cumulative_trt == 0 |
                                           cumulative_trt == k]
        ad_trt <- cumulative_trt[cumulative_trt == 0 |
                                   cumulative_trt == k]
        ad_PFS <- cumulative_PFS[cumulative_trt == 0 |
                                   cumulative_trt == k]
        ad_cens_PFS <- cumulative_cens_PFS[cumulative_trt == 0 |
                                             cumulative_trt == k]

        # Sort data by recruitment time plus event time for accurate time-to-event analysis
        ordering <- order(ad_recruit + ad_PFS)
        ad_recruit <- ad_recruit[ordering]
        ad_trt <- ad_trt[ordering]
        ad_PFS <- ad_PFS[ordering]
        ad_cens_PFS <- ad_cens_PFS[ordering]

        # Count cumulative events to determine if the required number of events is reached
        cum_events <- cumsum(ad_PFS <= ad_cens_PFS)

        # Find the time of the last required event
        if (any(cum_events >= required_events)) {
          Last_event <- min(which(cum_events >= required_events))
          Time_Futility <- (ad_recruit + ad_PFS)[Last_event]
        } else {
          Last_event <- length(ad_PFS)
          Time_Futility <- max(ad_recruit + ad_PFS)
        }

        # Determine if patients experienced a PFS event or were censored
        event_PFS <- ifelse(ad_PFS <= ad_cens_PFS &
                              ad_recruit + ad_PFS <= Time_Futility,
                            1,
                            0)

        # Adjust PFS times for censored patients
        ad_PFS[event_PFS == 0] <- ifelse(
          ad_recruit[event_PFS == 0] + ad_cens_PFS[event_PFS == 0] <= Time_Futility,
          ad_cens_PFS[event_PFS == 0],
          Time_Futility - ad_recruit[event_PFS == 0]
        )

        # Perform a Cox regression to calculate HR, test statistic, and p-value
        fit_Cox <- coxph(Surv(ad_PFS, event_PFS) ~ ad_trt)
        HR <- exp(fit_Cox$coefficients)
        TestStatistic <- summary(fit_Cox)$coef[1, "z"]
        PValue <- summary(fit_Cox)$coef[1, "Pr(>|z|)"]

        # Add results for this arm to the temporary data frame
        temp_results <- rbind(
          temp_results,
          data.frame(
            Scenario = z,
            TreatmentArm = k,
            Simulation = i,
            AnalysisType = ifelse(
              current_threshold_index == length(thresholds),
              "Final Analysis",
              paste("Interim", current_threshold_index)
            ),
            Time_Futility = Time_Futility,
            Initial_Status = "Continue",
            Futility_Outcome = "Continue",
            HR = HR,
            TestStatistic = TestStatistic,
            PValue = PValue,
            Dropped = FALSE
          )
        )
      }

      # Identify arms to drop based on futility (HR > 1)
      temp_results$Dropped[temp_results$HR > 1] <- TRUE
      temp_results$Futility_Outcome[temp_results$HR > 1] <- "Dropped due to Futility"

      # Apply "drop-the-loser" rules to drop additional arms if needed
      if (current_threshold_index == 1) {
        needed_drops <- max(2 - sum(temp_results$Dropped), 0)
        if (needed_drops > 0) {
          drop_candidates <- temp_results[!temp_results$Dropped, ]
          drop_candidates <- drop_candidates[order(abs(drop_candidates$TestStatistic), decreasing = FALSE), ]
          to_drop <- drop_candidates$TreatmentArm[1:needed_drops]
          temp_results$Dropped[temp_results$TreatmentArm %in% to_drop] <- TRUE
          temp_results$Futility_Outcome[temp_results$TreatmentArm %in% to_drop] <- "Dropped due to Drop-the-Loser Design"
        }
      } else if (current_threshold_index == 2) {
        if (sum(temp_results$Dropped) == 0 && length(remaining_arms) > 1) {
          drop_candidates <- temp_results[!temp_results$Dropped, ]
          drop_candidates <- drop_candidates[order(abs(drop_candidates$TestStatistic), decreasing = FALSE), ]
          to_drop <- drop_candidates$TreatmentArm[1]
          temp_results$Dropped[temp_results$TreatmentArm == to_drop] <- TRUE
          temp_results$Futility_Outcome[temp_results$TreatmentArm == to_drop] <- "Dropped due to Drop-the-Loser Design"
        }
      }

      # Add temporary results to the main results data frame
      futility_results <- rbind(futility_results, temp_results)

      # Update the list of dropped arms
      futility_drops <- unique(c(futility_drops, temp_results$TreatmentArm[temp_results$Dropped]))

      # Move to the next threshold for analysis
      current_threshold_index <- current_threshold_index + 1
    }
  }

  # Return the complete futility results for all scenarios and analyses
  return(futility_results)
}


# # Scenarios: Simplified based on the four cases outlined in the protocol
# Scenarios <- matrix(
#   c(1, 1, 1, 1,                # Null Hypothesis: No effect
#     0.67, 1, 1, 1,             # Hybrid 1: One effective drug
#     0.67, 0.67, 1, 1,          # Hybrid 2: Two effective drugs
#     0.67, 0.80, 1, 1,          # Hybrid 1 + 1: One effective, one half-effect drug
#     1.2, 1, 1, 1,              # Harmful drug
#     0.67, 0.67, 0.67, 0.67),   # All treatments are effect
#   ncol = 4,
#   byrow = TRUE,
#   dimnames = list(
#     c("Null hypothesis", "Hybrid 1", "Hybrid 2", "Hybrid 1 + 1", "Harmful", "All effectiv"),
#     NULL
#   )
# )
#
# # Results table: Record every simulation
# Results <- expand.grid(
#   Scenario = 1:nrow(Scenarios),
#   Sim = 1,
#   stringsAsFactors = FALSE
# )
#
# Results$HR_1 <- Scenarios[Results$Scenario, 1]
# Results$HR_2 <- Scenarios[Results$Scenario, 2]
# Results$HR_3 <- Scenarios[Results$Scenario, 3]
# Results$HR_4 <- Scenarios[Results$Scenario, 4]
#
#
#
# # Example
# results <- simulate_trial(
#   i = 1,
#   thresholds = c(30, 50, 100),
#   patients_stage1 = 100,
#   patients_stage2 = 50,
#   Results = Results,
#   RecruitmentRate = 15,
#   block_size = 6,
#   K = 4,
#   medianPFS_control = 8,
#   dropout_PFS = 0.1
# )


# -------------------------------------------------------------------------

#' Simulate a Clinical Trial with Futility, Drop-the-Loser, and Log-Rank Test Analyses
#'
#' This function simulates clinical trial scenarios to evaluate progression-free survival (PFS)
#' outcomes with interim and final analyses. The analyses include futility using Cox regression,
#' "drop-the-loser" strategy using log-rank tests, and a final log-rank test for efficacy.
#'
#' @param i Simulation iteration number.
#' @param thresholds Vector of event thresholds for interim and final analyses.
#' @param patients_stage1 Number of patients recruited in stage 1.
#' @param patients_stage2 Number of patients recruited in stage 2.
#' @param Results Data frame containing hazard ratios (HR) for each treatment arm in different scenarios.
#' @param RecruitmentRate Recruitment rate for patients in the trial.
#' @param block_size Block size for randomizing patients into treatment arms.
#' @param K Number of treatment arms (including the control arm).
#' @param medianPFS_control Median progression-free survival (PFS) time in the control group.
#' @param dropout_PFS Monthly dropout probability for PFS analysis.
#' @return A data frame containing futility analysis results for all simulated scenarios.
#' @examples
#' simulate_trial_log_rank(1, c(50, 100), 100, 200, Results, 0.1, 4, 3, 6, 0.05)
#' @export


simulate_trial_log_rank <- function(i,
                                    thresholds,
                                    patients_stage1,
                                    patients_stage2,
                                    Results,
                                    RecruitmentRate,
                                    block_size,
                                    K,
                                    medianPFS_control,
                                    dropout_PFS) {
  # Validate inputs
  if (!is.numeric(i) ||
      i <= 0)
    stop("Simulation iteration 'i' must be a positive integer.")
  if (!is.numeric(thresholds) ||
      any(thresholds <= 0))
    stop("Thresholds must be positive.")
  if (!is.numeric(patients_stage1) ||
      patients_stage1 <= 0)
    stop("patients_stage1 must be positive.")
  if (!is.numeric(patients_stage2) ||
      patients_stage2 <= 0)
    stop("patients_stage2 must be positive.")
  if (!is.data.frame(Results))
    stop("Results must be a data frame.")
  if (!is.numeric(RecruitmentRate) ||
      RecruitmentRate <= 0)
    stop("RecruitmentRate must be positive.")
  if (!is.numeric(block_size) ||
      block_size <= 0)
    stop("block_size must be positive.")
  if (!is.numeric(K) ||
      K <= 0)
    stop("K must be a positive integer.")
  if (!is.numeric(medianPFS_control) ||
      medianPFS_control <= 0)
    stop("medianPFS_control must be positive.")
  if (!is.numeric(dropout_PFS) ||
      dropout_PFS < 0 ||
      dropout_PFS > 1)
    stop("dropout_PFS must be between 0 and 1.")

  # Results storage
  futility_results <- data.frame(
    Scenario = integer(),
    TreatmentArm = integer(),
    Simulation = integer(),
    AnalysisType = character(),
    Time_Analysis = numeric(),
    Outcome = character(),
    HR = numeric(),
    TestStatistic_Cox = numeric(),
    TestStatistic_LogRank = numeric(),
    Dropped = logical(),
    stringsAsFactors = FALSE
  )

  # Loop through each scenario in the Results data frame
  for (z in 1:nrow(Results)) {
    cumulative_recruit <- c()
    cumulative_trt <- c()
    cumulative_PFS <- c()
    cumulative_cens_PFS <- c()
    futility_drops <- c()
    current_threshold_index <- 1

    # Process each threshold (interim and final analyses)
    while (current_threshold_index <= length(thresholds)) {
      required_events <- thresholds[current_threshold_index]

      if (current_threshold_index == 1) {
        # Stage 1 recruitment: exponential distribution for recruitment times
        new_recruit <- cumsum(rexp(patients_stage1 * 2, rate = RecruitmentRate))

        # Stage 2 recruitment: add delay between stages and adjust recruitment times
      } else {
        delay_months <- 2
        previous_analysis_time <- max(futility_results$Time_Analysis[futility_results$AnalysisType == paste("Interim", current_threshold_index - 1)], na.rm = TRUE)
        start_time_stage2 <- previous_analysis_time + delay_months
        new_recruit <- cumsum(rexp(patients_stage2 * 2, rate = RecruitmentRate)) + start_time_stage2
      }

      # Update cumulative recruitment and assign treatments
      cumulative_recruit <- c(cumulative_recruit, new_recruit)
      new_trt <- sample(rep(0:K, length.out = length(new_recruit)))
      cumulative_trt <- c(cumulative_trt, new_trt)

      # Simulate PFS times and censoring times
      new_PFS <- ifelse(
        new_trt == 0,
        rexp(length(new_recruit), log(2) / medianPFS_control),
        rexp(
          length(new_recruit),
          log(2) / medianPFS_control * as.numeric(Results[z, new_trt + 2])
        )
      )
      cumulative_PFS <- c(cumulative_PFS, new_PFS)

      # Censoring times follow exponential distribution
      new_cens_PFS <- rexp(length(new_recruit), rate = -log(1 - dropout_PFS) / 12)
      cumulative_cens_PFS <- c(cumulative_cens_PFS, new_cens_PFS)

      # Identify remaining treatment arms not dropped due to futility
      remaining_arms <- unique(setdiff(1:K, futility_drops))
      if (length(remaining_arms) == 0) {
        # If no arms remain, add a "No arms available" entry and stop analysis for this scenario
        futility_results <- rbind(
          futility_results,
          data.frame(
            Scenario = z,
            TreatmentArm = NA,
            Simulation = i,
            AnalysisType = ifelse(
              current_threshold_index == length(thresholds),
              "Final Analysis",
              paste("Interim", current_threshold_index)
            ),
            Time_Analysis = NA,
            Outcome = "No arms available due to Futility",
            HR = NA,
            TestStatistic_Cox = NA,
            TestStatistic_LogRank = NA,
            Dropped = NA
          )
        )
        break
      }

      # Data frame to temporarily store results for this threshold
      temp_results <- data.frame()

      # Analyze each remaining treatment arm
      for (k in remaining_arms) {
        ad_recruit <- cumulative_recruit[cumulative_trt == 0 |
                                           cumulative_trt == k]
        ad_trt <- cumulative_trt[cumulative_trt == 0 |
                                   cumulative_trt == k]
        ad_PFS <- cumulative_PFS[cumulative_trt == 0 |
                                   cumulative_trt == k]
        ad_cens_PFS <- cumulative_cens_PFS[cumulative_trt == 0 |
                                             cumulative_trt == k]

        ordering <- order(ad_recruit + ad_PFS)
        ad_recruit <- ad_recruit[ordering]
        ad_trt <- ad_trt[ordering]
        ad_PFS <- ad_PFS[ordering]
        ad_cens_PFS <- ad_cens_PFS[ordering]

        cum_events <- cumsum(ad_PFS <= ad_cens_PFS)
        if (any(cum_events >= required_events)) {
          last_event <- min(which(cum_events >= required_events))
          time_analysis <- (ad_recruit + ad_PFS)[last_event]
        } else {
          time_analysis <- max(ad_recruit + ad_PFS)
        }

        event_PFS <- ifelse(ad_PFS <= ad_cens_PFS &
                              ad_recruit + ad_PFS <= time_analysis,
                            1,
                            0)
        ad_PFS[event_PFS == 0] <- ifelse(
          ad_recruit[event_PFS == 0] + ad_cens_PFS[event_PFS == 0] <= time_analysis,
          ad_cens_PFS[event_PFS == 0],
          time_analysis - ad_recruit[event_PFS == 0]
        )

        # Futility Analysis (Cox Regression)
        fit_Cox <- coxph(Surv(ad_PFS, event_PFS) ~ ad_trt)
        HR <- exp(fit_Cox$coefficients)
        TestStatistic_Cox <- summary(fit_Cox)$coef[1, "z"]
        Dropped <- ifelse(HR > 1, TRUE, FALSE)

        # Log-Rank Test for Drop-the-Loser
        fit_LogRank <- survdiff(Surv(ad_PFS, event_PFS) ~ ad_trt)
        TestStatistic_LogRank <- fit_LogRank$chisq

        temp_results <- rbind(
          temp_results,
          data.frame(
            Scenario = z,
            TreatmentArm = k,
            Simulation = i,
            AnalysisType = ifelse(
              current_threshold_index == length(thresholds),
              "Final Analysis",
              paste("Interim", current_threshold_index)
            ),
            Time_Analysis = time_analysis,
            Outcome = ifelse(Dropped, "Dropped due to Futility", "Continue"),
            HR = HR,
            TestStatistic_Cox = TestStatistic_Cox,
            TestStatistic_LogRank = TestStatistic_LogRank,
            Dropped = Dropped
          )
        )
      }

      # Apply Drop-the-Loser Logic for Interim Analyses Only
      if (current_threshold_index < length(thresholds)) {
        num_drops_needed <- ifelse(current_threshold_index == 1, 2, 1)
        futile_arms <- temp_results[temp_results$Dropped, ]
        num_futile <- nrow(futile_arms)

        if (num_futile < num_drops_needed) {
          # Drop additional arms based on Log-Rank Test statistics
          remaining_to_drop <- num_drops_needed - num_futile
          drop_candidates <- temp_results[!temp_results$Dropped, ]
          drop_candidates <- drop_candidates[order(drop_candidates$TestStatistic_LogRank,
                                                   decreasing = FALSE), ]
          to_drop <- head(drop_candidates$TreatmentArm, remaining_to_drop)
          temp_results$Dropped[temp_results$TreatmentArm %in% to_drop] <- TRUE
          temp_results$Outcome[temp_results$TreatmentArm %in% to_drop] <- "Dropped due to Drop-the-Loser Design"
        }
      }

      # Update overall results
      futility_results <- rbind(futility_results, temp_results, row.names = NULL)
      futility_drops <- unique(c(futility_drops, temp_results$TreatmentArm[temp_results$Dropped]))

      # Finalize thresholds
      current_threshold_index <- current_threshold_index + 1
    }
  }

  return(futility_results)
}


  # # Scenarios: Simplified based on the four cases outlined in the protocol
# Scenarios <- matrix(
#   c(1, 1, 1, 1,                # Null Hypothesis: No effect
#     0.67, 1, 1, 1,             # Hybrid 1: One effective drug
#     0.67, 0.67, 1, 1,          # Hybrid 2: Two effective drugs
#     0.67, 0.80, 1, 1,          # Hybrid 1 + 1: One effective, one half-effect drug
#     1.2, 1, 1, 1,              # Harmful drug
#     0.67, 0.67, 0.67, 0.67),   # All treatments are effect
#   ncol = 4,
#   byrow = TRUE,
#   dimnames = list(
#     c("Null hypothesis", "Hybrid 1", "Hybrid 2", "Hybrid 1 + 1", "Harmful", "All effectiv"),
#     NULL
#   )
# )
#
# # Results table: Record every simulation
# Results <- expand.grid(
#   Scenario = 1:nrow(Scenarios),
#   Sim = 1,
#   stringsAsFactors = FALSE
# )
#
# Results$HR_1 <- Scenarios[Results$Scenario, 1]
# Results$HR_2 <- Scenarios[Results$Scenario, 2]
# Results$HR_3 <- Scenarios[Results$Scenario, 3]
# Results$HR_4 <- Scenarios[Results$Scenario, 4]
#
# Example
# results <- simulate_trial_log_rank(
#   i = 1,
#   thresholds = c(30, 50, 100),
#   patients_stage1 = 100,
#   patients_stage2 = 100,
#   Results = Results,
#   RecruitmentRate = 15,
#   block_size = 6,
#   K = 4,
#   medianPFS_control = 8,
#   dropout_PFS = 0.1
# )


# -------------------------------------------------------------------------


