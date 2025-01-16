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
plot_agg_perf(resamp_com_agg_perf,save.prefix="com_resamp")


#### Staggered, resampling ####
load("resamp_stag_sim_summaries.RData")

resamp_stag_agg_perf <- score_agg_perf(resamp.stag.summaries,staggered=T)
plot_agg_perf(resamp_stag_agg_perf,save.prefix="stag_resamp")

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

bal.data   <- generate.data(myparams,month.byunit=T,quarter.byunit=T,unbalanced=F)
bal.trt.dat <- assign.treatment(bal.data,starts=4,staggered=F)
these.bal <- make.data(bal.data,bal.trt.dat)

cons.data <- add.trt.effects(these.bal,tau=myparams$tau,grp.var = F, time.var = F)
cons.agg <- agg.data.common(cons.data,is.tv=T,is.gv=F)

tv.data   <- add.trt.effects(these.bal,tau=myparams$tau,grp.var = F, time.var = T)
tv.agg <- agg.data.common(tv.data,is.tv=T,is.gv=F)

gv.data   <- add.trt.effects(these.bal,tau=myparams$tau,grp.var = T, time.var = F)
gv.agg <- agg.data.common(gv.data,is.tv=T,is.gv=T)

gtv.data  <- add.trt.effects(these.bal,tau=myparams$tau,grp.var = T, time.var = T)
gtv.agg <- agg.data.common(gtv.data,is.tv=T,is.gv=T)

bal.truth <- bind_rows(cons.agg$month_true %>% mutate(agg="month"),cons.agg$quarter_true %>% mutate(agg="quarter"),cons.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="constant",grp=1) %>%
  bind_rows(bind_rows(tv.agg$month_true %>% mutate(agg="month"),tv.agg$quarter_true %>% mutate(agg="quarter"),tv.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="time-varying",grp=1)) %>%
  bind_rows(bind_rows(gv.agg$month_true %>% mutate(agg="month"),gv.agg$quarter_true %>% mutate(agg="quarter"),gv.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="group-varying")) %>%
  bind_rows(bind_rows(gtv.agg$month_true %>% mutate(agg="month"),gtv.agg$quarter_true %>% mutate(agg="quarter"),gtv.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="group- and time-varying")) %>%
  mutate(panel="balanced")
rm(cons.data,tv.data,tv.agg,gv.agg,gtv.agg)

unbal.data <- generate.data(myparams,month.byunit=T,quarter.byunit=T,unbalanced=T)
unbal.trt.dat <- assign.treatment(unbal.data,starts=4,staggered=F)
these.unbal <- make.data(unbal.data,unbal.trt.dat)

cons.data <- add.trt.effects(these.unbal,tau=myparams$tau,grp.var = F, time.var = F)
cons.agg <- agg.data.common(cons.data,is.tv=T,is.gv=F)

tv.data   <- add.trt.effects(these.unbal,tau=myparams$tau,grp.var = F, time.var = T)
tv.agg <- agg.data.common(tv.data,is.tv=T,is.gv=F)

gv.data   <- add.trt.effects(these.unbal,tau=myparams$tau,grp.var = T, time.var = F)
gv.agg <- agg.data.common(gv.data,is.tv=T,is.gv=T)

gtv.data  <- add.trt.effects(these.unbal,tau=myparams$tau,grp.var = T, time.var = T)
gtv.agg <- agg.data.common(gtv.data,is.tv=T,is.gv=T)

unbal.truth <- bind_rows(cons.agg$month_true %>% mutate(agg="month"),cons.agg$quarter_true %>% mutate(agg="quarter"),cons.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="constant",grp=1) %>%
  bind_rows(bind_rows(tv.agg$month_true %>% mutate(agg="month"),tv.agg$quarter_true %>% mutate(agg="quarter"),tv.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="time-varying",grp=1)) %>%
  bind_rows(bind_rows(gv.agg$month_true %>% mutate(agg="month"),gv.agg$quarter_true %>% mutate(agg="quarter"),gv.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="group-varying")) %>%
  bind_rows(bind_rows(gtv.agg$month_true %>% mutate(agg="month"),gtv.agg$quarter_true %>% mutate(agg="quarter"),gtv.agg$year_true %>% mutate(agg="year")) %>% mutate(trteff="group- and time-varying")) %>%
  mutate(panel="unbalanced")
rm(cons.data,tv.data,tv.agg,gv.agg,gtv.agg)

truth <- bind_rows(bal.truth,unbal.truth) %>%
  mutate(grp=factor(grp),
         trteff=factor(trteff,levels=c('constant','time-varying','group-varying','group- and time-varying')))

# Time-varying estimates in common adoption case:
ggplot(truth,aes(x=time,group=interaction(grp,panel))) + geom_line(aes(y=truth,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("truth_param_common.png",width=5,height=5)

these <- single.iter.param(myparams,staggered=T) %>%
  mutate(trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying')))
# Time-varying estimates in common adoption case:
ggplot(filter(these,group=="monthly_0",grepl("monthtrt",name),trteff!="null"),aes(x=time,group=panel)) + geom_line(aes(y=truth,col=panel)) +
  facet_wrap(~trteff,ncol=2)
