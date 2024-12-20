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
  n.month = 12, # each year has 12 months
  n.quarter = 4, # each year has 4 quarters
  n.year = 6,
  n.units = 50,
  sigma_week = 5,
  sigma_quarter = 1,
  phi = c(.1,0), 
  # Overall (non-stationary) trend 
  const = 0,
  # sigma_month is the white noise
  sigma_month = 1,
  
  # Cons Treatment effect 
  tau = 3)


# Simulation in parallel infrastructure
nsims <- 2
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)

#### Run simulation at scale: common adoption ####
com.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  single_iteration(myparams,staggered=F) %>% mutate(iter = i)
}

saveRDS(com.results, file=paste0(Sys.Date(),"_common_sim_results.rds"))


#### Run simulation at scale: staggered adoption ####
stag.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  single_iteration(myparams,staggered=T) %>% mutate(iter = i)
}

saveRDS(stag.results, file=paste0(Sys.Date(),"_stag_sim_results.rds"))

#### Run resampling simulation at scale: staggered adoption ####
load("cleaned_force_data.RData")

resample.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  resamp.dat <- resample(force.dat,n.units=100)
  # Inject treatment effects
  inject.analyze(resamp.dat,myparams,staggered=T) %>% mutate(iter=i)
}
parallel::stopCluster(cl)
