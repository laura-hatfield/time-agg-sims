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

## Just load the summaries of the simulation results
load("common_sim_summaries.RData")

com_winners <- score_model_winners(com.summaries)

ggplot(filter(com_winners,panel=="balanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=pct.decrement)) + geom_text(aes(label=signif(value,2),col=pct.decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none") +
  scale_fill_gradient(low='#762a83',high='white') +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_winners,panel=="unbalanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=pct.decrement)) + geom_text(aes(label=signif(value,2),col=pct.decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none") +
  scale_fill_gradient(low='#762a83',high='white') +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_winners,panel=="balanced",estimand=="Time-varying",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=pct.decrement)) + geom_text(aes(label=signif(value,2),col=pct.decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none") +
  scale_fill_gradient(low='#762a83',high='white') +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_winners,panel=="unbalanced",estimand=="Time-varying",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=pct.decrement)) + geom_text(aes(label=signif(value,2),col=pct.decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none") +
  scale_fill_gradient(low='#762a83',high='white') +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
