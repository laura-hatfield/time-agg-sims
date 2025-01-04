library(tidyverse)
library(doParallel)

source("sim_funs.R")
# Other dependencies:
# MASS for mvrnorm()
# car linearHypothesis() tests
# estimatr for lm_robust()
# did  for Callaway and Sant'Anna estimators
# Matrix for nearPD()

#### Set parameters ####
myparams <- list(
  n.year = 6,
  n.units = 50,
  sigma_quarter = 1, # Variance of the quarter effects
  phi = c(.1,0), 
  const = 0,   # Overall (non-stationary) trend of month effect
  # sigma_month is the white noise
  sigma_month = 1,
  
  # Cons Treatment effect 
  tau = 3)


# Simulation in parallel
cl <- parallel::makeCluster(12)
doParallel::registerDoParallel(cl)
nsims <- 2000

## This one is already done
if (F){
  #### Parametric, common adoption ####
  com.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
    single.iter.param(myparams,staggered=F) %>% mutate(iter = i)
  }
  
  saveRDS(com.results, file="common_sim_results.rds")
}

#### Parametric, staggered adoption ####
stag.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  single.iter.param(myparams,staggered=T) %>% mutate(iter = i)
}

saveRDS(stag.results, file="stag_sim_results.rds")
rm(stag.results); gc()

#### Resampling, common adoption ####
load("cleaned_force_data.RData")

resamp.com.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  single.iter.resamp(myparams,force.dat,starts=3:4,staggered=F) %>% mutate(iter=i)
}

saveRDS(resamp.com.results,file="resamp_com_sim_results.rds")
rm(resamp.com.results);gc()

#### Resampling, staggered adoption ####
# Supply the real distribution of start years (but require 2 years pre and post)
last.start <- max(force.dat$year)-2
starts <- (force.dat %>% group_by(unitID) %>% slice(1) %>% 
     ungroup() %>% filter(start.year %in% 3:last.start))$start.year

resamp.stag.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
 single.iter.resamp(myparams,force.dat,starts=starts,staggered=T) %>% mutate(iter=i)
}

saveRDS(resamp.stag.results,file="resamp_stag_sim_results.rds")
rm(resamp.stag.results); gc()

parallel::stopCluster(cl)
