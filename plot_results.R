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

com_mod_perf <- score_model_perf(com.summaries)

ggplot(filter(com_mod_perf,panel=="balanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_mod_perf,panel=="unbalanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_mod_perf,panel=="balanced",estimand=="Time-varying",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_mod_perf,panel=="unbalanced",estimand=="Time-varying",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))


com_agg_perf <- score_agg_perf(com.summaries,staggered=F)

ggplot(filter(com_agg_perf,estimand=="Overall",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_agg_perf,estimand=="Time-varying",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

load("stag_sim_summaries.RData")

stag_agg_perf <- score_agg_perf(stag.summaries,staggered=T)

ggplot(filter(stag_agg_perf,estimand=="Overall",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(stag_agg_perf,estimand=="Time-varying",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))


load("resamp_com_sim_summaries.RData")

resamp_com_agg_perf <- score_agg_perf(resamp.com.summaries,staggered=F)

ggplot(filter(resamp_com_agg_perf,estimand=="Overall",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(resamp_com_agg_perf,estimand=="Time-varying",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement)) + 
  scale_color_gradient(low='white',high='black',guide="none",limits=c(0,.25),oob=scales::squish) +
  scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish) +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
