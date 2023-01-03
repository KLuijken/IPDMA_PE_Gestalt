# ---------------------------------------------------------------------------- #
# Accuracy of the physicians’ intuitive risk estimation in the diagnostic 
# management of pulmonary embolism
# Script author: Kim Luijken
# 
# The current script generates simulated data to illustrate use of the code
# ---------------------------------------------------------------------------- #

set.seed( 24112022)

# specify sample size and number of variables involved
n <- 10000
n_vars <- 15

# introduce a dependency between variables
sig_mat <- matrix( 0.4, nrow = n_vars, ncol = n_vars)
diag( sig_mat) <- 1
correlation <- mvrnorm( n, mu = rep( 0, n_vars), Sigma = sig_mat)

# generate binary variables
PEdiagn         <- rbinom( n, 1, plogis( -1 + correlation[ , 1]))
gestalt         <- rbinom( n, 1, plogis( -1.6 + correlation[ , 2]))
pub_before_2010 <- rbinom( n, 1, plogis( -0.8 + correlation[ , 3]))
prev_VTE        <- rbinom( n, 1, plogis( -1 + correlation[ , 4]))
HR_100          <- rbinom( n, 1, plogis( -1.7 + correlation[ , 5]))
hemoptysis      <- rbinom( n, 1, plogis( -1.2 + correlation[ , 6]))
cancer          <- rbinom( n, 1, plogis( -0.9 + correlation[ , 7]))
symp_dvt        <- rbinom( n, 1, plogis( -1 + correlation[ , 8]))
immosurg        <- rbinom( n, 1, plogis( -2 + correlation[ , 9]))
sex             <- rbinom( n, 1, plogis( -0.4 + correlation[ , 10]))
chf             <- rbinom( n, 1, plogis( -0.9 + correlation[ , 11]))
cld             <- rbinom( n, 1, plogis( -2.3 + correlation[ , 12]))

# generate categorical variables
setting         <- rbinom( n, 3, plogis( -2 + correlation[ , 13]))
study           <- rbinom( n, 6, plogis( correlation[ , 14]))

# generate continuous variables
age <- correlation[ , 15] + 65


# store in data frame
single_simulated_data <- data.frame( "PEdiagn" = PEdiagn,
                              "gestalt" = gestalt,
                              "pub_before_2010" = pub_before_2010,
                              "prev_VTE" = prev_VTE,
                              "HR_100" = HR_100,
                              "hemoptysis" = hemoptysis,
                              "cancer" = cancer,
                              "symp_dvt" = symp_dvt,
                              "immosurg" = immosurg,
                              "sex" = sex,
                              "chf" = chf,
                              "cld" = cld,
                              "setting"  = setting,
                              "study" = study,
                              "age" = age
                              )

# introduce missing values
missing_simulated_data <- ampute( single_simulated_data,
                                  prop = 0.3,   # 30% missing
                                  mech = "MAR", # missing at random
                                  pattern = c( rep( 0, times = 13),
                                               1, 0) # no missings in "study"
                                  )$amp

# make data right type for imputation (ignore warning)
missing_simulated_data <- data.frame( "PEdiagn" = factor( missing_simulated_data$PEdiagn),
                                     "gestalt" = factor( missing_simulated_data$gestalt),
                                     "pub_before_2010" = factor( missing_simulated_data$pub_before_2010),
                                     "prev_VTE" = factor( missing_simulated_data$prev_VTE),
                                     "HR_100" = factor( missing_simulated_data$HR_100),
                                     "hemoptysis" = factor( missing_simulated_data$hemoptysis),
                                     "cancer" = factor( missing_simulated_data$cancer),
                                     "symp_dvt" = factor( missing_simulated_data$symp_dvt),
                                     "immosurg" = factor( missing_simulated_data$immosurg),
                                     "sex" = factor( missing_simulated_data$sex),
                                     "chf" = factor( missing_simulated_data$chf),
                                     "cld" = factor( missing_simulated_data$cld),
                                     "setting"  = factor( missing_simulated_data$setting),
                                     "study" = factor( missing_simulated_data$study),
                                     "age" = age
                                     )
# save incompete data
saveRDS( missing_simulated_data, "./data/simulated_data_incomplete.rds")

# impute missing data
imputed_simulated_data <- mice( missing_simulated_data, m = 10)

# store data in a way similar to provided dataset of current IPDMA study
simulated_data <- complete( imputed_simulated_data, "long")

simulated_data <- cbind( data.frame( apply( simulated_data[ , 1:2], 
                                     MARGIN = 2,
                                     as.numeric)),
                         data.frame( "study" = simulated_data$study,
                                     "age" = simulated_data$age),
                         data.frame( apply( simulated_data[ , 3:13], 
                                     MARGIN = 2,
                                     as.numeric))
                              )
saveRDS( simulated_data, "./data/simulated_data.rds")