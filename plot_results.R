library(tidyverse)
theme_set(theme_minimal())
source("plot_funs.R")

## this step is expensive on large simulation runs, so just do it once
if (F){
  com.results <- readRDS("common_sim_results.rds")
  com.summaries <- summarize_results(com.results)
  save(com.summaries,file="com_results_summaries.RData")
  rm(com.results);gc()
}

## Just load the summaries of the simulation results
load("com_results_summaries.RData")

com_winners <- score_winners(com.summaries)

ggplot(filter(com_winners,panel=="balanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=round(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(mod_win_plot_dat,panel=="unbalanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=round(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
