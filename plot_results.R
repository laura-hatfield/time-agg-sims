library(tidyverse)
theme_set(theme_minimal()+theme(legend.position="bottom"))
source("plot_funs.R")
source("sim_funs.R")

## this step is expensive on large simulation runs, so just do it once
if (F){
  com.results <- readRDS("time-agg-sims-results/common_sim_results.rds")
  com.summaries <- process_results_com(com.results)
  save(com.summaries,file="common_sim_summaries.RData")
  rm(com.results);gc()
  
  stag.results <- readRDS("time-agg-sims-results/stag_sim_results.rds")
  stag.summaries <- process_results_stag(stag.results)
  save(stag.summaries,file="stag_sim_summaries.RData")
  rm(stag.results);gc()  
  
  resamp.com.results <- readRDS("time-agg-sims-results/resamp_com_sim_results.rds")
  resamp.com.summaries <- process_results_com(resamp.com.results)
  save(resamp.com.summaries,file="resamp_com_sim_summaries.RData")
  rm(resamp.com.results);gc()  
  
  resamp.stag.results <- readRDS("time-agg-sims-results/resamp_stag_sim_results.rds")
  resamp.stag.summaries <- process_results_stag(resamp.stag.results)
  save(resamp.stag.summaries,file="resamp_stag_sim_summaries.RData")
  rm(resamp.stag.results);gc()  
}

#### Common, parametric ####
load("common_sim_summaries.RData")

# Performance of different models
com_mod_perf <- score_mod_perf(com.summaries)
plot_mod_perf(com_mod_perf,save.prefix="figures/com_param")

# Performance of different time aggregations
com_agg_perf <- score_agg_perf(com.summaries,staggered=F)
plot_agg_perf(com_agg_perf,save.prefix="figures/com_param")

#### Staggered, parametric ####
load("stag_sim_summaries.RData")

stag_agg_perf <- score_agg_perf(stag.summaries,staggered=T)
plot_agg_perf(stag_agg_perf,save.prefix="figures/stag_param")

#### Common, resampling ####
load("resamp_com_sim_summaries.RData")

resamp_com_mod_perf <- score_mod_perf(resamp.com.summaries)
plot_mod_perf(resamp_com_mod_perf,save.prefix="figures/com_resamp")

resamp_com_agg_perf <- score_agg_perf(resamp.com.summaries,staggered=F)
plot_agg_perf(resamp_com_agg_perf,save.prefix="figures/com_resamp")

#### Staggered, resampling ####
load("resamp_stag_sim_summaries.RData")

resamp_stag_agg_perf <- score_agg_perf(resamp.stag.summaries,staggered=T)
plot_agg_perf(resamp_stag_agg_perf,save.prefix="figures/stag_resamp")

## Plot the true treatment effects:
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

# Time-varying truth in parametric, common adoption scenarios:
com.param.truth <- plot_truth(myparams,starts=rep(4,2),resample=F,staggered=F) %>% mutate(time=as.numeric(time))
ggplot(com.param.truth,aes(x=time,group=interaction(grp,panel))) + geom_line(aes(y=truth,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("truth_param_common.png",width=5,height=5)

# Time-varying truth in parametric, staggered adoption scenarios:
stag.param.truth <- plot_truth(myparams,starts=3:4,resample=F,staggered=T) %>% filter(start.year>0)
ggplot(stag.param.truth,aes(x=time,group=interaction(grp,start.year,panel))) + geom_line(aes(y=truth,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("truth_param_staggered.png",width=5,height=5)

load("cleaned_force_data.RData")
# Need a much smaller treatment effect for this outcome scale:
myparams$tau <- 0.04
## Need more treatment units with a rare outcome
myparams$n.units <- 100 

# Time-varying truth in resampling, common adoption scenarios:
com.resamp.truth <- plot_truth(myparams,starts=rep(4,2),resample=T,staggered=F,data=force.dat) %>% mutate(time=as.numeric(time))
ggplot(com.resamp.truth,aes(x=time,group=interaction(grp,panel))) + geom_line(aes(y=truth,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("truth_resamp_common.png",width=5,height=5)

# Time-varying truth in resampling, staggered adoption scenarios
stag.resamp.truth <- plot_truth(myparams,starts=3:4,resample=T,staggered=T,data=force.dat) %>% filter(start.year>0)
ggplot(stag.resamp.truth,aes(x=time,group=interaction(grp,start.year,panel))) + geom_line(aes(y=truth,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("truth_resamp_staggered.png",width=5,height=5)
