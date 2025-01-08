library(tidyverse)
theme_set(theme_minimal()+theme(legend.position="bottom"))
source("plot_funs.R")

## this step is expensive on large simulation runs, so just do it once
if (F){
  com.results <- readRDS("common_sim_results.rds")
  com.summaries <- process_results_com(com.results)
  save(com.summaries,file="common_sim_summaries.RData")
  rm(com.results);gc()
  
  stag.results <- readRDS("stag_sim_results.rds")
  stag.summaries <- process_results_stag(stag.results)
  save(stag.summaries,file="stag_sim_summaries.RData")
  rm(stag.results);gc()  
  
  resamp.com.results <- readRDS("resamp_com_sim_results.rds")
  resamp.com.summaries <- process_results_com(resamp.com.results)
  save(resamp.com.summaries,file="resamp_com_sim_summaries.RData")
  rm(resamp.com.results);gc()  
  
  resamp.stag.results <- readRDS("resamp_stag_sim_results.rds")
  resamp.stag.summaries <- process_results_stag(resamp.stag.results)
  save(resamp.stag.summaries,file="resamp_stag_sim_summaries.RData")
  rm(resamp.stag.results);gc()  
}

#### Common, parametric ####
load("common_sim_summaries.RData")

# Performance of different models
com_mod_perf <- score_mod_perf(com.summaries)
plot_mod_perf(com_mod_perf,save.prefix="com_param")

# Performance of different time aggregations
com_agg_perf <- score_agg_perf(com.summaries,staggered=F)
plot_agg_perf(com_agg_perf,save.prefix="com_param")

#### Staggered, parametric ####
load("stag_sim_summaries.RData")

stag_agg_perf <- score_agg_perf(stag.summaries,staggered=T)
plot_agg_perf(stag_agg_perf,save.prefix="stag_param")

#### Common, resampling ####
load("resamp_com_sim_summaries.RData")

resamp_com_mod_perf <- score_mod_perf(resamp.com.summaries)
plot_mod_perf(resamp_com_mod_perf,save.prefix="com_resamp")

resamp_com_agg_perf <- score_agg_perf(resamp.com.summaries,staggered=F)
plot_agg_perf(resamp_com_agg_perf,save.prefix="comp_resamp")


#### Staggered, resampling ####
load("resamp_stag_sim_summaries.RData")

resamp_stag_agg_perf <- score_agg_perf(resamp.stag.summaries,staggered=T)
plot_agg_perf(resamp_stag_agg_perf)
