library(tidyverse) 
library(odin) 

# Setup
source('code/utils.R')                   # CHECKED 
source('code/global_parameters.R')       # CHECKED
source('code/parameters.R')              # CHECKED

# Uncontrolled epidemics
source('code/episims.R')                 # CHECKED
source('code/overdispersion.R')          # LOGICKED
source('code/survival.R')                # CHECKED 

# Controlled epidemics 
source('code/isolation.R')               # LOGICKED
source('code/gatheringsize.R')           # LOGICKED

# Inference
source('code/growthrate.R')              # LOGICKED
source('code/g_identifiability_2.R')     # LOGICKED 

source('code/psi_empirical_ligti.R')     # LOGICKED 