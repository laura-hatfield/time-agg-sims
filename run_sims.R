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
cl <- parallel::makeCluster(10)
doParallel::registerDoParallel(cl)
nsims <- 5000

#### Parametric, common adoption ####
com.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  single.iter(params=myparams,staggered=F,resampling=F) %>% mutate(iter = i)
}
saveRDS(com.results, file="common_sim_results.rds")
rm(com.results); gc()

#### Parametric, staggered adoption ####
stag.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  single.iter(params=myparams,staggered=T,resampling=F) %>% mutate(iter = i)
}
saveRDS(stag.results, file="stag_sim_results.rds")
rm(stag.results); gc()

#### Resampling, staggered adoption ####
load("cleaned_force_data.RData")

# Need a much smaller treatment effect for this outcome scale:
myparams$tau <- 0.08
## Need more treatment units with a rare outcome
myparams$n.units <- 100 

# Supply the real distribution of start years (but require 2 years pre and post)
# So the only possible start years are 3 or 4
last.start <- max(force.dat$year)-2
mystarts <- (force.dat %>% group_by(unitID) %>% slice(1) %>% 
             ungroup() %>% filter(start.year %in% 3:last.start))$start.year

resamp.stag.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine = 'rbind') %dopar% {
  single.iter(params=myparams,staggered=T,resampling=T,data=force.dat,starts=mystarts) %>% mutate(iter=i)
}
saveRDS(resamp.stag.results,file="resamp_stag_sim_results.rds")
rm(resamp.stag.results); gc()

#### Resampling, common adoption #### # .combine = 'rbind'
resamp.com.results <- foreach::foreach(i = 1:nsims,.packages=c('tidyverse','estimatr','car','did'), .combine='rbind') %dopar% {
  single.iter(params=myparams,staggered=F,resampling=T,data=force.dat,starts=mystarts) %>% mutate(iter=i)
}
saveRDS(resamp.com.results,file="resamp_com_sim_results.rds")
rm(resamp.com.results);gc()

parallel::stopCluster(cl)
