# ---------------------------------------------------------------------------- #
# Accuracy of the physicians’ intuitive risk estimation in the diagnostic 
# management of pulmonary embolism
# Script author: Kim Luijken
#
# The current script describes dependencies
# ---------------------------------------------------------------------------- #

library( MASS)    # for data simulation
library( mice)    # for missing data in simulated data
library( metafor) # for forest function
library( dplyr)   # for data wrangling
library( lme4)    # for main analysis
library( rstanarm)# for non-converging analysis of sens and spec
library( gt)      # for presentation of tables in rmd
