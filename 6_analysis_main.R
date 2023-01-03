# ---------------------------------------------------------------------------- #
# Accuracy of the physicians’ intuitive risk estimation in the diagnostic 
# management of pulmonary embolism
# Script author: Kim Luijken
#
# Run the current script to perform the main analysis: RR, sensitivity, and 
# specificity. Focused on results for one subgroup (female versus male)
# ---------------------------------------------------------------------------- #

# pre-processing: stratify on sex
fem_data <- simulated_data[ simulated_data$sex == 1, ]
male_data <- simulated_data[ simulated_data$sex == 0, ]


# RR ----

# ---------------------------------------------------------------------------- #
# Overall RR
# ---------------------------------------------------------------------------- #

# estimate log RRs and ses over all studies
overall_input_RR <- estimating_log_RR_overall( input_data = simulated_data)

# overall RR from two-stage regression
ma_2stage_overall <- rma( yi = overall_input_RR$pooled_log_RRs_per_study,
                          sei = overall_input_RR$pooled_se_log_RRs_per_study,
                          method = "REML",
                          test = "knha") # For CI estimation, 
                                         # corrected for assumed fixed variance

# store predicted output of mixed model
out_2stage_overall <- predict( ma_2stage_overall)

# create table
overall_table_RR <- gt( data.frame( "RR" = round( out_2stage_overall$pred,
                                              digits = 2),
                                "RR_pi" = paste0( round( out_2stage_overall$pi.lb,
                                                         digits = 2),
                                                  " to ",
                                                  round( out_2stage_overall$pi.ub,
                                                         digits = 2)))) %>%
  cols_label(
    RR = "RR",
    RR_pi = "95% Prediction interval"
  ) %>%
  tab_header( "Overall RR for positive versus negative gestalt on risk of 
   pulmonary embolism diagnosis (in simulated data)") %>%
  opt_align_table_header(align = "left")

# save table in tables/ folder
gtsave( overall_table_RR, 
        "./tables/overall_RR.docx")

# ---------------------------------------------------------------------------- #
# RR in females
# ---------------------------------------------------------------------------- #

# estimate log RRs and ses over all studies
female_input_RR <- estimating_log_RR_overall( input_data = fem_data)

# overall RR from two-stage regression
ma_2stage_fem <- rma( yi = female_input_RR$pooled_log_RRs_per_study,
                          sei = female_input_RR$pooled_se_log_RRs_per_study,
                          method = "REML",
                          test = "knha") # For CI estimation, 
# corrected for assumed fixed variance

# store predicted output of mixed model
out_2stage_fem <- predict( ma_2stage_fem)

# create table
female_table_RR <- gt( data.frame( "RR" = round( out_2stage_fem$pred, digits = 2),
                                 "RR_pi" = paste0( round( out_2stage_fem$pi.lb,
                                                          digits = 2),
                                                   " to ",
                                                   round( out_2stage_fem$pi.ub,
                                                          digits = 2)))) %>%
  cols_label(
    RR = "RR",
    RR_pi = "95% Prediction interval"
  ) %>%
  tab_header( "RR for positive versus negative gestalt on risk of 
   pulmonary embolism diagnosis in females (in simulated data)") %>%
  opt_align_table_header(align = "left")

# save table in tables/ folder
gtsave( female_table_RR, 
        "./tables/female_RR.docx")

# ---------------------------------------------------------------------------- #
# RR in males
# ---------------------------------------------------------------------------- #

# estimate log RRs and ses over all studies
male_input_RR <- estimating_log_RR_overall( input_data = male_data)

# overall RR from two-stage regression
ma_2stage_male <- rma( yi = male_input_RR$pooled_log_RRs_per_study,
                      sei = male_input_RR$pooled_se_log_RRs_per_study,
                      method = "REML",
                      test = "knha") # For CI estimation, 
# corrected for assumed fixed variance

# store predicted output of mixed model
out_2stage_male <- predict( ma_2stage_male)

# create table
male_table_RR <- gt( data.frame( "RR" = round( out_2stage_male$pred, digits = 2),
                               "RR_pi" = paste0( round( out_2stage_male$pi.lb, 
                                                        digits = 2),
                                                 " to ",
                                                 round( out_2stage_male$pi.ub,
                                                        digits = 2)))) %>%
  cols_label(
    RR = "RR",
    RR_pi = "95% Prediction interval"
  ) %>%
  tab_header( "RR for positive versus negative gestalt on risk of 
   pulmonary embolism diagnosis in males (in simulated data)") %>%
  opt_align_table_header(align = "left")

# save table in tables/ folder
gtsave( male_table_RR, 
        "./tables/male_RR.docx")

# sensitivity and specificity ----

# ---------------------------------------------------------------------------- #
# overall sensitivity and specificity
# ---------------------------------------------------------------------------- #

overall_input_ss <- extract_subset_input_sensspec( subset_data_level = simulated_data)

# create table
overall_table_ss <- gt( data.frame( "am" = c( "Sensitivity", "Specificity"), 
                                    "pam_ci" = c( overall_input_ss$sensitivity..95..CI.,
                                                  overall_input_ss$specificity..95..CI.))) %>%
  cols_label(
    am = " ",
    pam_ci = "Pooled accuracy measure (95% CI)"
  ) %>%
  tab_header( "Overall sensitivity and specificity of gestalt for diagnosing 
   pulmonary embolism (in simulated data)") %>%
  opt_align_table_header(align = "left")

# save table in tables/ folder
gtsave( overall_table_ss, 
        "./tables/overall_sensspec.docx")

# ---------------------------------------------------------------------------- #
# Sensitivity and specificity for females
# ---------------------------------------------------------------------------- #

fem_input_ss <- extract_subset_input_sensspec( subset_data_level = fem_data)

# create table
female_table_ss <- gt( data.frame( "am" = c( "Sensitivity", "Specificity"),
                                    "pam_ci" = c( fem_input_ss$sensitivity..95..CI.,
                                                  fem_input_ss$specificity..95..CI.))) %>%
  cols_label(
    am = " ",
    pam_ci = "Pooled accuracy measure (95% CI)"
  ) %>%
  tab_header( "Sensitivity and specificity of gestalt for diagnosing 
   pulmonary embolism in females (in simulated data)") %>%
  opt_align_table_header(align = "left")

# save table in tables/ folder
gtsave( female_table_ss, 
        "./tables/female_sensspec.docx")

# ---------------------------------------------------------------------------- #
# Sensitivity and specificity for males
# ---------------------------------------------------------------------------- #

male_input_ss <- extract_subset_input_sensspec( subset_data_level = male_data)

### singularity warning -- alternative estimation
# there was a singular fit of the bivariate generalized linear mixed-effects 
# model with zero variance of the logit-sensitivity parameter.  
# estimate the model again using a Bayesian estimation approach  
male_input_ssB <- extract_subset_input_sensspec_bayesian( subset_data_level = 
                                                             male_data)

# also estimate the model again using fixed effects only  
male_input_ssF <- extract_subset_input_sensspec_fixed( subset_data_level = 
                                                         male_data)

# create table
male_table_ss <- gt( data.frame( "am" = c( "Sensitivity",
                                           "Specificity",
                                           "Sensitivity (Bayesian estimation)",
                                           "Specificity (Bayesian estimation)",
                                           "Sensitivity (fixed effects)",
                                           "Specificity (fixed effects"),
                                 "pam_ci" = c( male_input_ss$sensitivity..95..CI.,
                                               male_input_ss$specificity..95..CI.,
                                               male_input_ssB$sensitivity..95..CI.,
                                               male_input_ssB$specificity..95..CI.,
                                               male_input_ssF$sensitivity..95..CI.,
                                               male_input_ssF$specificity..95..CI.))) %>%
  cols_label(
    am = " ",
    pam_ci = "Pooled accuracy measure (95% CI)"
  ) %>%
  tab_header( "Sensitivity and specificity of gestalt for diagnosing 
   pulmonary embolism in males (in simulated data)") %>%
  opt_align_table_header(align = "left")

# save table in tables/ folder
gtsave( male_table_ss, 
        "./tables/male_sensspec.docx")