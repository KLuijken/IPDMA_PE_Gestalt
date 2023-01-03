# ---------------------------------------------------------------------------- #
# Accuracy of the physicians’ intuitive risk estimation in the diagnostic 
# management of pulmonary embolism
# Script author: Kim Luijken
#
# The current script executes the analysis
# ---------------------------------------------------------------------------- #

# load packages ----
source( "2_packages.R")

# read in data ----
# alternatively, code can be simulated using 3_simulate_data.R
# note that this script can be replaced with a data cleaning script in a new
# empirical analysis
simulated_data <- readRDS( "./data/simulated_data.rds")
simulated_data_incomplete <- readRDS( "./data/simulated_data_incomplete.rds")

# read in user-defined functions
source( "4_helpers.R")

# perform descriptive analyses ----
# outputs two tables under the /tables folder with missingness and sample size 
# descriptions; as well as three figures with the RR, sensitivity, and 
# specificity of gestalt in individual studies
source( "5_analysis_descriptive.R")


# perform main analyses ----
# outputs six tables under the /tables folder with the RR, sensitivity, and 
# specificity overall and for a single subgroup
source( "6_analysis_main.R")
