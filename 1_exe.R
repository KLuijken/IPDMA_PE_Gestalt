# ---------------------------------------------------------------------------- #
# Accuracy of the physicians’ intuitive risk estimation in the diagnostic 
# management of pulmonary embolism
# Script author: Kim Luijken
#
# The current script executes the analysis
# ---------------------------------------------------------------------------- #


# read in data ----
simulated_data <- readRDS( "./data/simulated_data.rds")
simulated_data_incomplete <- readRDS( "./data/simulated_data_incomplete.rds")

# read in user-defined functions
source( "helpers.R")

# perform descriptive analyses ----
# outputs two tables under the /tables folder with missingness and sample size 
# descriptions; as well as three figures with the RR, sensitivity, and 
# specificity of gestalt in individual studies
source( "analysis_descriptive.R")


# perform main analyses ----
# outputs six tables under the /tables folder with the RR, sensitivity, and 
# specificity overall and for a single subgroup
source( "analysis_main.R")
