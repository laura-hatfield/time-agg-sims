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
## Set this 
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
  resamp.dat <- resample(force.dat,n.units=50)
  # Inject treatment effects
  inject.analyze(resamp.dat,myparams,staggered=T) %>% mutate(iter=i)
}

saveRDS(resample.results,file=paste0(Sys.Date(),"_resamp_stag_sim_results.rds"))
parallel::stopCluster(cl)
