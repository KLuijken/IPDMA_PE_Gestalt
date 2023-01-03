# ---------------------------------------------------------------------------- #
# Accuracy of the physicians’ intuitive risk estimation in the diagnostic 
# management of pulmonary embolism
# Script author: Kim Luijken
#
# The current script contains helper functions to perform the analysis
# ---------------------------------------------------------------------------- #

# This script contains the following functions:
### logit: to perform a logit transformation
### rubins_rules_var: to estimate the variance of a point estimate when pooling
###                   across imputed datasets
### estimating_log_RR_per_study: to estimate the log RR PE diagnosis - gestalt 
###                   in an individual study pooled over imputed datasets
### estimating_log_RR_overall: to estimate the log RR PE diagnosis - gestalt and
###                   ses for all studies
### extract_subset_input: to pool the RRs over studies using a 2-stage meta 
###                   analysis ( random effect for gestalt)
### create_subset_forest: to visualise the RR in levels of a specific subgroup
###                   in a small forest plot ( not presented in main text)
### create_input_forest: to rearrange output of extract_subset_input() such that
###                   it can be easily used in the metafor::forest() function
### estimating_sens_per_study: to estimate the sensitivity of gestalt in an
###                   individual study pooled over imputed datasets
### estimating_spec_per_study: to estimate the specificity of gestalt in an
###                   individual study pooled over imputed datasets
### estimating_sensspec_input_per_study: to provide input for the two-stage  
###                   bivariate model estimating the overall sensitivity and 
###                   specificity of gestalt ( true/false positive and true/
###                   false negative per study pooled over imputed datasets)
### extract_subset_input_sensspec: to estimate the overall sensitivity and 
###                   specificity of gestalt in a subgroup 
### extract_subset_input_sensspec_bayesian: to estimate the overall sensitivity and 
###                   specificity of gestalt in a subgroup using a Bayesian
###                   estimation procedure
### extract_subset_input_sensspec_fixed: to estimate the overall sensitivity and 
###                   specificity of gestalt in a subgroup using fixed effects 
###                   only
### prop_miss_table: to describe missingness, written by Toshihiko Takada


# logit ----
logit <- function(x) {log( x/( 1-x))}

# rubins_rules_var ----
rubins_rules_var <- function( estimates,
                              ses,
                              n_imputed_sets = length( unique( simulated_data$.imp))){
  # Within study variance
  within_var <- mean( ses^2)
  
  # Between study variance
  theta_bar <- mean( estimates)
  
  between_var <- sum( ( estimates - theta_bar)^2)/n_imputed_sets
  
  # total variance
  total_var <- within_var + between_var + between_var/n_imputed_sets
  
  return( total_var)
}


# estimating_log_RR_per_study ----
# subset data per study and compute log RR for all imputed data sets + pool
# using Rubin's rules
estimating_log_RR_per_study <- function( study_index,
                                         input_data){
  # take study subset
  single_study_data <- input_data[input_data$study == study_index, ]
  
  # fit log-binomial model in study subset for each imputed data set
  fit_single_study <-  by( single_study_data,
                           # subset by imputed data set
                           INDICES = as.factor( single_study_data$.imp),
                           # log-binomial model
                           function(y) {glm( PEdiagn ~ gestalt, 
                                             data = y,
                                             family = binomial( link = "log"))
                           })
  
  # pool using Rubin's rules
  pooled_fit_single_study <- pool( fit_single_study)
  
  return( pooled_fit_single_study)
}

# estimating_log_RR_overall ----
# combine output from estimating_log_RR_per_study() on individual studies to 
# arrive at overall log RRs and ses
estimating_log_RR_overall <- function( input_data){
  # estimate log RR per study pooled across imputed data sets
  ### use helper function estimating_log_RR_per_study
  MI_log_RRs_per_study <- lapply( unique( input_data$study),
                                  function(x) estimating_log_RR_per_study( 
                                    study_index = x,
                                    input_data = input_data))
  
  # extract relevant results
  pooled_log_results_per_study <- lapply( MI_log_RRs_per_study, 
                                          function(x) summary( x, conf.int = TRUE))
  
  # store overall log RR of gestalt-PEdiagn per study
  pooled_log_RRs_per_study <- unlist( 
    lapply( 1:length(unique( input_data$study)), 
            function(x) pooled_log_results_per_study[[ x]]$estimate[ 
              pooled_log_results_per_study[[ x]]$term == "gestalt"]))
  
  # store standard error of overall log RR of gestalt-PEdiagn per study
  pooled_se_log_RRs_per_study <- unlist( 
    lapply( 1:length(unique( input_data$study)), 
            function(x) pooled_log_results_per_study[[ x]]$std.error[ 
              pooled_log_results_per_study[[ x]]$term == "gestalt"]))
  
  return( list( pooled_log_RRs_per_study = pooled_log_RRs_per_study,
                pooled_se_log_RRs_per_study = pooled_se_log_RRs_per_study)
          )
}

# extract_subset_input ----
# define function to subset data per study and per subset
# compute log RR for all imputed data sets + pool using Rubin's rules
# perform a 2-stage meta analysis
extract_subset_input <- function( stratify_value,
                                  overall_data,
                                  stratify_var,
                                  n_imp = length( unique( simulated_data$.imp))){
  # subset the data
  complete_data <- overall_data[ complete.cases( overall_data[,  stratify_var]), ]
  subset_data_level <- complete_data[ which( complete_data[,  stratify_var] == stratify_value), ]
  
  # estimate log RR per study pooled across imputed data sets in data subset
  ### use function 1 above
  MI_log_RRs_per_study_per_subset <- lapply( unique( subset_data_level$study),
                                             function(x) estimating_log_RR_per_study( study_index = x,
                                                                                      input_data = subset_data_level))
  # extract relevant results
  pooled_log_results_per_study_per_subset <- lapply( MI_log_RRs_per_study_per_subset, 
                                                     function(x) summary( x, conf.int = TRUE))
  
  # store log RR of gestalt-PEdiagn in subset per study
  pooled_log_RRs_per_study_per_subset <- unlist( lapply( 1:length(unique( subset_data_level$study)), 
                                                         function(x)
                                                           pooled_log_results_per_study_per_subset[[x]]$estimate[pooled_log_results_per_study_per_subset[[x]]$term == "gestalt"]))
  
  # store standard error of log RR of gestalt-PEdiagn in subset per study
  pooled_se_log_RRs_per_study_per_subset <- unlist( lapply( 1:length(unique( subset_data_level$study)), 
                                                            function(x)
                                                              pooled_log_results_per_study_per_subset[[x]]$std.error[pooled_log_results_per_study_per_subset[[x]]$term == "gestalt"]))
  
  
  # meta-analysis to obtain RR across all studies in subset
  ma_2stage_subset <- rma( yi = pooled_log_RRs_per_study_per_subset,
                           sei = pooled_se_log_RRs_per_study_per_subset,
                           method = "REML",
                           test = "knha") # For CI estimation, 
  # corrected for assumed fix variance
  
  out_2stage_subset <- predict( ma_2stage_subset)
  
  return( list( out_2stage_subset = out_2stage_subset, n = (round( nrow( subset_data_level) / n_imp))))
  
}

# create_subset_forest ----
# based on input from function extract_subset_input()
create_subset_forest <- function( overall_data,
                                  stratify_var, # character
                                  levels_stratify_var # character vector
){ 
  stratify_var_values <- na.omit( sort( unique( overall_data[ , stratify_var]),
                                        decreasing = T))
  
  stratified_output <- lapply( stratify_var_values,
                               function( x){ extract_subset_input( stratify_value = x,
                                                                   overall_data = overall_data,
                                                                   stratify_var = stratify_var)})
  stratified_log_RRs <- sapply( 1:length( stratify_var_values),
                                function( x) stratified_output[[x]]$out_2stage_subset$pred)
  
  stratified_lb_pi <- sapply( 1:length( stratify_var_values),
                              function( x) stratified_output[[x]]$out_2stage_subset$pi.lb)
  
  stratified_ub_pi <- sapply( 1:length( stratify_var_values),
                              function( x) stratified_output[[x]]$out_2stage_subset$pi.ub)
  
  stratified_n <- sapply( 1:length( stratify_var_values),
                          function( x) stratified_output[[x]]$n)
  
  metafor::forest( x = stratified_log_RRs,
                   ci.lb = stratified_lb_pi,
                   ci.ub = stratified_ub_pi,
                   atransf = exp,
                   at = log(c(0.6, 1, 2, 4)),
                   refline = 0,
                   xlab = "Relative Risk on PEdiagnosis for gestalt (log scale)",
                   main = paste0( "RRs per ", stratify_var),
                   slab = levels_stratify_var,
                   ilab = stratified_n,
                   ilab.xpos = -2,
                   psize = 1)
}



# create_input_forest ----
# based on input from function extract_subset_input()
create_input_forest <- function( overall_data,
                                 stratify_var, # character
                                 levels_stratify_var # character vector
){ 
  stratify_var_values <- na.omit( sort( unique( overall_data[ , stratify_var]),
                                        decreasing = T))
  
  stratified_output <- lapply( stratify_var_values,
                               function( x){ extract_subset_input( stratify_value = x,
                                                                   overall_data = overall_data,
                                                                   stratify_var = stratify_var)})
  stratified_log_RRs <- sapply( 1:length( stratify_var_values),
                                function( x) stratified_output[[x]]$out_2stage_subset$pred)
  
  stratified_lb_pi <- sapply( 1:length( stratify_var_values),
                              function( x) stratified_output[[x]]$out_2stage_subset$pi.lb)
  
  stratified_ub_pi <- sapply( 1:length( stratify_var_values),
                              function( x) stratified_output[[x]]$out_2stage_subset$pi.ub)
  
  stratified_n <- sapply( 1:length( stratify_var_values),
                          function( x) stratified_output[[x]]$n)
  
  return( list( stratified_log_RRs = stratified_log_RRs,
                stratified_lb_pi = stratified_lb_pi,
                stratified_ub_pi = stratified_ub_pi,
                stratified_n = stratified_n))
}

# estimating_sens_per_study ----
# define function to subset data per study and compute sensitivity for all  
# imputed data sets + pool using Rubin's rules
estimating_sens_per_study <- function( study_index,
                                       input_data){
  
  # take study subset
  single_study_data <- input_data[ which( input_data$study == study_index), ]
  
  # true positives per study across all imputations
  true_positives <-  c( by( single_study_data,
                            # subset by imputed data set
                            INDICES = as.factor( single_study_data$.imp),
                            # log-binomial model
                            function(y) { sum( y[ , "PEdiagn"] == 1 & y[ , "gestalt"] == 1)
                            }))
  
  # false negatives per study across all imputations
  false_negatives <-  c( by( single_study_data,
                             # subset by imputed data set
                             INDICES = as.factor( single_study_data$.imp),
                             # log-binomial model
                             function(y) { sum( y[ , "PEdiagn"] == 1 & y[ , "gestalt"] == 0)
                             }))
  
  sens <- true_positives / ( true_positives + false_negatives)
  
  se_sens <- sqrt( (sens * ( 1- sens)) / ( true_positives + false_negatives))
  
  pooled_sens <- mean( sens)
  
  pooled_var_sens <- rubins_rules_var( estimates = sens,
                                       ses = se_sens)
  
  return( list( pooled_sens = pooled_sens,
                pooled_var_sens = pooled_var_sens,
                TP = round( mean( true_positives)),
                FN = round( mean( false_negatives))))
}

# estimating_spec_per_study ----
# define function to subset data per study and compute specificity for all  
# imputed data sets + pool using Rubin's rules
estimating_spec_per_study <- function( study_index,
                                       input_data){
  
  # take study subset
  single_study_data <- input_data[ which( input_data$study == study_index), ]
  
  # true negatives per study across all imputations
  true_negatives <-  c( by( single_study_data,
                            # subset by imputed data set
                            INDICES = as.factor( single_study_data$.imp),
                            # log-binomial model
                            function(y) { sum( y[ , "PEdiagn"] == 0 & y[ , "gestalt"] == 0)
                            }))
  
  # false positives per study across all imputations
  false_positives <-  c( by( single_study_data,
                             # subset by imputed data set
                             INDICES = as.factor( single_study_data$.imp),
                             # log-binomial model
                             function(y) { sum( y[ , "PEdiagn"] == 0 & y[ , "gestalt"] == 1)
                             }))
  
  spec <- true_negatives / ( true_negatives + false_positives)
  
  se_spec <- sqrt( (spec * ( 1- spec)) / ( true_negatives + false_positives))
  
  pooled_spec <- mean( spec)
  
  pooled_var_spec <- rubins_rules_var( estimates = spec,
                                       ses = se_spec)
  
  return( list( pooled_spec = pooled_spec,
                pooled_var_spec = pooled_var_spec,
                TN = round( mean( true_negatives)),
                FP = round( mean( false_positives))))
}


# estimating_sensspec_input_per_study ----
# compute true positives, false positives, true negative, and
# false negatives per study averaged over imputations as input for bivariate 
# model
estimating_sensspec_input_per_study <- function( study_index,
                                                 input_data){
  # take study subset
  single_study_data <- input_data[ which( input_data$study == study_index), ]
  
  # true positives per study across all imputations (averaged)
  true_positives <-  mean( c( by( single_study_data,
                                  # subset by imputed data set
                                  INDICES = as.factor( single_study_data$.imp),
                                  # log-binomial model
                                  function(y) { sum( y[ , "PEdiagn"] == 1 & y[ , "gestalt"] == 1)
                                  })))
  
  # false positives per study across all imputations (averaged)
  false_positives <-  mean( c( by( single_study_data,
                                   # subset by imputed data set
                                   INDICES = as.factor( single_study_data$.imp),
                                   # log-binomial model
                                   function(y) { sum( y[ , "PEdiagn"] == 0 & y[ , "gestalt"] == 1)
                                   })))
  
  # true negatives per study across all imputations (averaged)
  true_negatives <-  mean( c( by( single_study_data,
                                  # subset by imputed data set
                                  INDICES = as.factor( single_study_data$.imp),
                                  # log-binomial model
                                  function(y) { sum( y[ , "PEdiagn"] == 0 & y[ , "gestalt"] == 0)
                                  })))
  
  # false negatives per study across all imputations (averaged)
  false_negatives <-  mean( c( by( single_study_data,
                                   # subset by imputed data set
                                   INDICES = as.factor( single_study_data$.imp),
                                   # log-binomial model
                                   function(y) { sum( y[ , "PEdiagn"] == 1 & y[ , "gestalt"] == 0)
                                   })))
  
  
  study_data <- data.frame( "SID" = rep( study_index, each = 2),
                            "TP" = rep( true_positives, each = 2),
                            "FN" = rep( false_negatives, each = 2),
                            "FP" = rep( false_positives, each = 2),
                            "TN" = rep( true_negatives, each = 2))
  
  bivariate_data <- study_data
  bivariate_data$dis <- c( 1, 0)
  bivariate_data$nondis <- c( 0, 1)
  bivariate_data$pos <- bivariate_data$dis * bivariate_data$TP + bivariate_data$nondis * bivariate_data$TN
  bivariate_data$n <- bivariate_data$dis * (bivariate_data$TP + bivariate_data$FN) + 
    bivariate_data$nondis * (bivariate_data$TN + bivariate_data$FP)
  
  return( bivariate_data)
}


# extract_subset_input_sensspec ----
# estimating sensitivity and specificity in subsets
# define function to subset data per study and per subset
# compute dis / nondis input averaged over imputed data sets
# run bivariate model
extract_subset_input_sensspec <- function( subset_data_level){
  # estimate log RR per study pooled across imputed data sets
  ### use function X above
  bivariate_input_per_study <- lapply( unique( subset_data_level$study),
                                       function(x) estimating_sensspec_input_per_study( study_index = x,
                                                                                        input_data = subset_data_level))
  
  bivariate_input <- do.call( rbind, bivariate_input_per_study)
  bivariate_input <- cbind( "study" = bivariate_input$SID,
                            round( bivariate_input[ ,-1])) # studyID is a factor
  
  model_biv <- glmer( cbind( pos, n - pos) ~ -1 + dis + nondis + (-1 + dis + nondis | study),
                      family = binomial,
                      data = bivariate_input)   
  
  # sensitivity
  ## point estimate
  mean_logit_sens <- fixef( model_biv)["dis"]
  mean_sens <- plogis( mean_logit_sens)
  
  ## standard error
  se_logit_sens <- sqrt( diag( vcov( model_biv)))["dis"]
  
  ## confidence interval
  lci_sens <- plogis( mean_logit_sens - 1.96 * se_logit_sens)
  uci_sens <- plogis( mean_logit_sens + 1.96 * se_logit_sens)
  
  # specificity
  ## point estimate
  mean_logit_spec <- fixef( model_biv)["nondis"]
  mean_spec <- plogis( mean_logit_spec)
  
  ## standard error
  se_logit_spec <- sqrt( diag( vcov( model_biv)))["nondis"]
  
  ## confidence interval
  lci_spec <- plogis( mean_logit_spec - 1.96 * se_logit_spec)
  uci_spec <- plogis( mean_logit_spec + 1.96 * se_logit_spec)
  
  sens <- data.frame( "subgroup_sens" = "placeholder",
                      "sensitivity (95% CI)" = paste0( round( mean_sens, digits = 2),
                                                       " (",
                                                       round( lci_sens, digits = 2),
                                                       " to ",
                                                       round( uci_sens, digits = 2),
                                                       ")"))
  
  spec <- data.frame( "subgroup_spec" = "placeholder",
                      "specificity (95% CI)" = paste0( round( mean_spec, digits = 2),
                                                       " (",
                                                       round( lci_spec, digits = 2),
                                                       " to ",
                                                       round( uci_spec, digits = 2),
                                                       ")"))
  
  out_table <- cbind( sens, spec)
  
  return( out_table)
  
}


# extract_subset_input_sensspec_bayesian ----
# estimating sensitivity and specificity in subsets
# define function to subset data per study and per subset
# compute dis / nondis input averaged over imputed data sets
# run bivariate model using bayesian approach
extract_subset_input_sensspec_bayesian <- function( subset_data_level){
  # estimate log RR per study pooled across imputed data sets
  ### use function X above
  bivariate_input_per_study <- lapply( unique( subset_data_level$study),
                                       function(x) estimating_sensspec_input_per_study( study_index = x,
                                                                                        input_data = subset_data_level))
  
  bivariate_input <- do.call( rbind, bivariate_input_per_study)
  bivariate_input <- cbind( "study" = bivariate_input$SID,
                            round( bivariate_input[ ,-1])) # studyID is a factor
  
  model_biv <- stan_glmer( cbind( pos, n - pos) ~ -1 + dis + nondis + (-1 + dis + nondis | study),
                           family = binomial,
                           data = bivariate_input)   
  
  # sensitivity
  ## point estimate
  mean_logit_sens <- fixef( model_biv)["dis"]
  mean_sens <- plogis( mean_logit_sens)
  
  ## standard error
  se_logit_sens <- sqrt( diag( vcov( model_biv)))["dis"]
  
  ## confidence interval
  lci_sens <- plogis( mean_logit_sens - 1.96 * se_logit_sens)
  uci_sens <- plogis( mean_logit_sens + 1.96 * se_logit_sens)
  
  # specificity
  ## point estimate
  mean_logit_spec <- fixef( model_biv)["nondis"]
  mean_spec <- plogis( mean_logit_spec)
  
  ## standard error
  se_logit_spec <- sqrt( diag( vcov( model_biv)))["nondis"]
  
  ## confidence interval
  lci_spec <- plogis( mean_logit_spec - 1.96 * se_logit_spec)
  uci_spec <- plogis( mean_logit_spec + 1.96 * se_logit_spec)
  
  sens <- data.frame( "subgroup_sens" = "placeholder",
                      "sensitivity (95% CI)" = paste0( round( mean_sens, digits = 2),
                                                       " (",
                                                       round( lci_sens, digits = 2),
                                                       " to ",
                                                       round( uci_sens, digits = 2),
                                                       ")"))
  
  spec <- data.frame( "subgroup_spec" = "placeholder",
                      "specificity (95% CI)" = paste0( round( mean_spec, digits = 2),
                                                       " (",
                                                       round( lci_spec, digits = 2),
                                                       " to ",
                                                       round( uci_spec, digits = 2),
                                                       ")"))
  
  out_table <- cbind( sens, spec)
  
  return( out_table)
  
}


# extract_subset_input_sensspec_fixed ----
# estimating sensitivity and specificity in subsets
# define function to subset data per study and per subset
# compute dis / nondis input averaged over imputed data sets
# run bivariate model with fixed effects only
extract_subset_input_sensspec_fixed <- function( subset_data_level){
  # estimate log RR per study pooled across imputed data sets
  ### use function X above
  bivariate_input_per_study <- lapply( unique( subset_data_level$study),
                                       function(x) estimating_sensspec_input_per_study( study_index = x,
                                                                                        input_data = subset_data_level))
  
  bivariate_input <- do.call( rbind, bivariate_input_per_study)
  bivariate_input <- cbind( "study" = as.numeric( bivariate_input$SID),
                            round( bivariate_input[ ,-1])) # change studyID to numeric
  
  model_biv <- glm( cbind( pos, n - pos) ~ -1 + dis + nondis + (-1 + dis + nondis | study),
                    family = binomial,
                    data = bivariate_input)   
  
  # sensitivity
  ## point estimate
  mean_logit_sens <- model_biv$coefficients["dis"]
  mean_sens <- plogis( mean_logit_sens)
  
  ## standard error
  se_logit_sens <- sqrt( diag( vcov( model_biv)))["dis"]
  
  ## confidence interval
  lci_sens <- plogis( mean_logit_sens - 1.96 * se_logit_sens)
  uci_sens <- plogis( mean_logit_sens + 1.96 * se_logit_sens)
  
  # specificity
  ## point estimate
  mean_logit_spec <- model_biv$coefficients["nondis"]
  mean_spec <- plogis( mean_logit_spec)
  
  ## standard error
  se_logit_spec <- sqrt( diag( vcov( model_biv)))["nondis"]
  
  ## confidence interval
  lci_spec <- plogis( mean_logit_spec - 1.96 * se_logit_spec)
  uci_spec <- plogis( mean_logit_spec + 1.96 * se_logit_spec)
  
  sens <- data.frame( "subgroup_sens" = "placeholder",
                      "sensitivity (95% CI)" = paste0( round( mean_sens, digits = 2),
                                                       " (",
                                                       round( lci_sens, digits = 2),
                                                       " to ",
                                                       round( uci_sens, digits = 2),
                                                       ")"))
  
  spec <- data.frame( "subgroup_spec" = "placeholder",
                      "specificity (95% CI)" = paste0( round( mean_spec, digits = 2),
                                                       " (",
                                                       round( lci_spec, digits = 2),
                                                       " to ",
                                                       round( uci_spec, digits = 2),
                                                       ")"))
  
  out_table <- cbind( sens, spec)
  
  return( out_table)
  
}

# prop_miss_table ----
# function for missings written by Toshihiko Takada
prop_miss_table <- function( dataframe) {
  m <- sapply( dataframe, function(x) {
    data.frame(
      nmiss = sum( is.na( x)), 
      n = length( x), 
      propmiss = sum( is.na( x))/ length( x)
    )
  })
  d <- data.frame( t( m))
  d <- sapply( d, unlist)
  d <- as.data.frame( d)
  d$variable <- row.names( d)
  row.names(d) <- NULL
  d <- cbind( d[ ncol( d)], d[ -ncol( d)])
  return(d)
}
