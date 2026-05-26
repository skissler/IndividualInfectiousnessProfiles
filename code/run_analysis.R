library(tidyverse) 
library(odin) 

# Setup
source('code/utils.R') 
source('code/global_parameters.R')
source('code/parameters.R')

# Uncontrolled epidemics
source('code/episims.R')                 # CHECKED
source('code/survival.R')                # 
source('code/growthrate.R')
source('code/g_identifiability_2.R')
source('code/overdispersion.R')

# Controlled epidemics 
# source('code/isolation.R')
source('code/gatheringsize.R')

# Inference
# source('code/identifiability.R')
source('code/psi_identifiability.R')

