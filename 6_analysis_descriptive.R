# ---------------------------------------------------------------------------- #
# Accuracy of the physicians’ intuitive risk estimation in the diagnostic 
# management of pulmonary embolism
# Script author: Kim Luijken
#
# Run the current script to perform descriptive analyses 
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# missingness ----
# ---------------------------------------------------------------------------- #
# describe missing data for original dataset (before imputation)
missing <- prop_miss_table( dataframe = simulated_data_incomplete)

# save table in tables/ folder
gtsave( gt( missing), 
        "./tables/descriptives_missing.docx")

# ---------------------------------------------------------------------------- #
# sample size ----
# ---------------------------------------------------------------------------- #

# total n
n_tot_temp <- c( by( simulated_data[ simulated_data$.imp == 1,],
                     simulated_data[ simulated_data$.imp == 1, "study"],
                     nrow))
n_tot <- data.frame( "study" = names( n_tot_temp),
                     "n_tot" = n_tot_temp)

# individuals with PE
n_finalVTE_temp <- c( by( simulated_data[ simulated_data$PEdiagn == 1,],
                          simulated_data[ simulated_data$PEdiagn == 1, "study"],
                          function(x) nrow(x) / length( unique( simulated_data$.imp))))
n_finalVTE <- data.frame( "study" = names( n_tot_temp),
                          "n_finalVTE" = paste0( round( n_finalVTE_temp),
                                                 " (",
                                                 round( n_finalVTE_temp / n_tot_temp,
                                                        digits = 2),
                                                 ")"
                                                 ))

# individuals with gestalt
n_gestalt_temp <- c( by( simulated_data[ simulated_data$gestalt == 1,],
                              simulated_data[ simulated_data$gestalt == 1, "study"],
                              function(x) nrow(x) / length( unique( simulated_data$.imp))))
n_gestalt <- data.frame( "study" = names( n_tot_temp),
                              "n_gestalt" = paste0( round( n_gestalt_temp),
                                                    " (", 
                                                    round( n_gestalt_temp / n_tot_temp, 
                                                           digits = 2),
                                                    ")"))

overall_n <- data.frame( "Study" = unique( simulated_data$study),
                         "Overall_sample_size" = n_tot$n_tot,
                         "FinalPE" = n_finalVTE$n_finalVTE,
                         "Gestalt_pos" = n_gestalt$n_gestalt)

table_overall_n <- gt( overall_n) %>%
  cols_label( Study = "Study",
              Overall_sample_size = "Overall sample size",
              FinalPE = "Individuals with final PE (prop)",
              Gestalt_pos = "Gestalt-positive individuals")

# save table in tables/ folder
gtsave( table_overall_n, 
        "./tables/descriptives_n.docx")

# ---------------------------------------------------------------------------- #
# RR per individual study ----
# ---------------------------------------------------------------------------- #

# estimate log RRs and ses over all studies
overall_input_RR <- estimating_log_RR_overall( input_data = simulated_data)

# make forest plot on RR scale
# save figure in figures/ folder
pdf( "./figures/overall_RR_descriptive.pdf")
metafor::forest( x = overall_input_RR$pooled_log_RRs_per_study,
                 sei = overall_input_RR$pooled_se_log_RRs_per_study,
                 atransf = exp,
                 at = log(c(0.6, 1, 2, 4)),
                 refline = 0,
                 xlab = "Relative Risk on VTE for PEgestalte (log scale)",
                 main = "RRs per study (descriptive)",
                 ilab.xpos = -1,
                 header="study",
                 xlim = c(-2, 4),
                 psize = 1)
dev.off()


# ---------------------------------------------------------------------------- #
# sensitivity per individual study ----
# ---------------------------------------------------------------------------- #

# compute sens input per study
raw_per_study_sens <- lapply( unique( simulated_data$study),
                              function(x) estimating_sens_per_study( study_index = x,
                                                                     input_data = simulated_data))

raw_input_per_study_sens <- data.frame( matrix( unlist( do.call( rbind, raw_per_study_sens)), ncol = 4))
colnames( raw_input_per_study_sens) <- unique( names( unlist( raw_per_study_sens)))

# compute confidence intervals
cilb_per_study_sens <- raw_input_per_study_sens$pooled_sens - 1.96 * sqrt( raw_input_per_study_sens$pooled_var_sens)
ciub_per_study_sens <- raw_input_per_study_sens$pooled_sens + 1.96 * sqrt( raw_input_per_study_sens$pooled_var_sens)

# cap off confidence intervals at 0 - 1 bound
cilb_per_study_sens <- ifelse( cilb_per_study_sens < 0, 0, cilb_per_study_sens) 
ciub_per_study_sens <- ifelse( ciub_per_study_sens > 1, 1, ciub_per_study_sens)

# save figure in figures/ folder
pdf( "./figures/overall_sens_descriptive.pdf")
metafor::forest( x = raw_input_per_study_sens$pooled_sens,
                 ci.lb = cilb_per_study_sens,
                 ci.ub = ciub_per_study_sens,
                 at = seq( 0, 1, by = 0.1),
                 refline = NA,
                 xlab = "Sensitivity",
                 ilab = paste0( raw_input_per_study_sens$TP, "/", raw_input_per_study_sens$TP + raw_input_per_study_sens$FN),
                 ilab.xpos = 0,
                 ilab.pos = 2,
                 header= F,
                 xlim = c(-1, 1.5),
                 psize = 1)
text( -1, 9, "study", font = 2, pos = 4)
text( 0, 9, "TP / positives", font = 2, pos = 2, cex = 0.8)
text( 1.5, 9, "Sensitivity [ 95% CI]", font = 2, cex = 0.8, pos = 2)
dev.off()

# ---------------------------------------------------------------------------- #
# specificity per individual study ----
# ---------------------------------------------------------------------------- #

# compute spec input per study
raw_per_study_spec <- lapply( unique( simulated_data$study),
                              function(x) estimating_spec_per_study( study_index = x,
                                                                     input_data = simulated_data))

raw_input_per_study_spec <- data.frame( matrix( unlist( do.call( rbind, raw_per_study_spec)), ncol = 4))
colnames( raw_input_per_study_spec) <- unique( names( unlist( raw_per_study_spec)))

# compute confidence intervals
cilb_per_study_spec <- raw_input_per_study_spec$pooled_spec - 1.96 * sqrt( raw_input_per_study_spec$pooled_var_spec)
ciub_per_study_spec <- raw_input_per_study_spec$pooled_spec + 1.96 * sqrt( raw_input_per_study_spec$pooled_var_spec)

# cap off confidence intervals at 0 - 1 bound
cilb_per_study_spec <- ifelse( cilb_per_study_spec < 0, 0, cilb_per_study_spec) 
ciub_per_study_spec <- ifelse( ciub_per_study_spec > 1, 1, ciub_per_study_spec)

# save figure in figures/ folder
pdf( "./figures/overall_spec_descriptive.pdf")
metafor::forest( x = raw_input_per_study_spec$pooled_spec,
                 ci.lb = cilb_per_study_spec,
                 ci.ub = ciub_per_study_spec,
                 at = seq( 0, 1, by = 0.1),
                 refline = NA,
                 xlab = "Specificity",
                 ilab = paste0( raw_input_per_study_spec$TN, "/", raw_input_per_study_spec$TN + raw_input_per_study_spec$FP),
                 ilab.xpos = 0,
                 ilab.pos = 2,
                 header= F,
                 xlim = c(-1, 1.5),
                 psize = 1)
text( -1, 9, "study", font = 2, pos = 4)
text( 0, 9, "TN / negatives", font = 2, pos = 2, cex = 0.8)
text( 1.5, 9, "Specificity [ 95% CI]", font = 2, cex = 0.8, pos = 2)
dev.off()




