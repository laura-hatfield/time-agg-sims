library(tidyverse)
theme_set(theme_minimal()+theme(legend.position="bottom"))
source("plot_funs.R")
source("sim_funs.R")

## this step is expensive on large simulation runs, so just do it once
if (F){
  com.results <- readRDS("results/common_sim_results.rds")
  com.summaries <- process_results_com(com.results)
  save(com.summaries,file="common_sim_summaries.RData")
  rm(com.results);gc()
  
  stag.results <- readRDS("results/stag_sim_results.rds")
  stag.summaries <- process_results_stag(stag.results)
  save(stag.summaries,file="stag_sim_summaries.RData")
  rm(stag.results);gc()  
  
  resamp.com.results <- readRDS("results/resamp_com_sim_results.rds")
  resamp.com.results <- resamp.com.results %>%
    ## Multiply the estimates (and SEs) by 100 for readability
    mutate(across(.cols=c(est,se,truth,lb,ub),.fns=~.x*100))
  resamp.com.summaries <- process_results_com(resamp.com.results)
  save(resamp.com.summaries,file="resamp_com_sim_summaries.RData")
  rm(resamp.com.results);gc()  
  
  resamp.stag.results <- readRDS("results/resamp_stag_sim_results.rds")
  resamp.stag.results <- resamp.stag.results %>%
    ## Multiply the estimates (and SEs) by 100 for readability
    mutate(across(.cols=c(est,truth,se),.fns=~.x*100))    
  resamp.stag.summaries <- process_results_stag(resamp.stag.results)
  save(resamp.stag.summaries,file="resamp_stag_sim_summaries.RData")
  rm(resamp.stag.results);gc()  
}

#### Heat plots of summary statistics ####
load("common_sim_summaries.RData")
load("stag_sim_summaries.RData")
load("resamp_com_sim_summaries.RData")
load("resamp_stag_sim_summaries.RData")

# Common, parametric, models
com_mod_perf <- score_mod_perf(com.summaries)
plot_mod_perf(com_mod_perf,save.prefix="figures/com_param")

# Common, parametric, aggregations
com_agg_perf <- score_agg_perf(com.summaries,staggered=F)
plot_agg_perf(com_agg_perf,save.prefix="figures/com_param")
table_winners(com_agg_perf,bias=.01,reject=.01,rmse=.1,prefix="figures/com_param")

# Staggered, parametric, aggregations
stag_agg_perf <- score_agg_perf(stag.summaries,staggered=T)
plot_agg_perf(stag_agg_perf,save.prefix="figures/stag_param")
table_winners(stag_agg_perf,bias=.01,reject=.01,rmse=.1,prefix="figures/stag_param")

# Common, resampling, models
resamp_com_mod_perf <- score_mod_perf(resamp.com.summaries)
plot_mod_perf(resamp_com_mod_perf,save.prefix="figures/com_resamp")

# Common, resampling, aggregations
resamp_com_agg_perf <- score_agg_perf(resamp.com.summaries,staggered=F)
plot_agg_perf(resamp_com_agg_perf,save.prefix="figures/com_resamp")
table_winners(resamp_com_agg_perf,bias=.01,reject=.01,rmse=.1,prefix="figures/com_resamp")

# Staggered, resampling, aggregations
resamp_stag_agg_perf <- score_agg_perf(resamp.stag.summaries,staggered=T)
plot_agg_perf(resamp_stag_agg_perf,save.prefix="figures/stag_resamp")
table_winners(resamp_stag_agg_perf,bias=.01,reject=.01,rmse=.1,prefix="figures/stag_resamp")

#### True treatment effects by time, group ####
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
com.param.truth <- plot_truth(myparams,starts=rep(4,2),resample=F,staggered=F) %>% mutate(time=as.numeric(time)) %>%
  mutate(trteff=factor(trteff,levels=c("constant","time-varying","group-varying","group- and time-varying"),
                       labels=c('constant','time','group','group-time')))
ggplot(com.param.truth,aes(x=time,group=interaction(grp,panel))) + geom_line(aes(y=truth,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("figures/truth_param_common.png",width=4,height=4)

# Time-varying truth in parametric, staggered adoption scenarios:
stag.param.truth <- plot_truth(myparams,starts=3:4,resample=F,staggered=T) %>% filter(start.year>0) %>%
  mutate(trteff=factor(trteff,levels=c("constant","time-varying","group-varying","group- and time-varying"),
                       labels=c('constant','time','group','group-time')))
ggplot(stag.param.truth,aes(x=time,group=interaction(grp,start.year,panel))) + geom_line(aes(y=truth,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("figures/truth_param_staggered.png",width=4,height=4)

load("cleaned_force_data.RData")
# Need a much smaller treatment effect for this outcome scale:
myparams$tau <- 0.08
## Need more treatment units with a rare outcome
myparams$n.units <- 100 

# Time-varying truth in resampling, common adoption scenarios:
com.resamp.truth <- plot_truth(myparams,starts=rep(4,2),resample=T,staggered=F,data=force.dat) %>% mutate(time=as.numeric(time))
ggplot(com.resamp.truth,aes(x=time,group=interaction(grp,panel))) + geom_line(aes(y=truth*100,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("figures/truth_resamp_common.png",width=5,height=5)

# Time-varying truth in resampling, staggered adoption scenarios
stag.resamp.truth <- plot_truth(myparams,starts=3:4,resample=T,staggered=T,data=force.dat) %>% filter(start.year>0)
ggplot(stag.resamp.truth,aes(x=time,group=interaction(grp,start.year,panel))) + geom_line(aes(y=truth*100,col=grp,lty=panel)) +
  facet_grid(trteff~agg,scale="free_x") + theme(panel.grid.minor=element_blank()) +
  scale_color_brewer(palette="Dark2",guide="none")
ggsave("figures/truth_resamp_staggered.png",width=5,height=5)

#### Distribution of truth and estimates ####
com.results <- readRDS("results/common_sim_results.rds")
plot_mc.dist(com.results,staggered=F,filename="figures/com_param_est_dist.png")
rm(com.results); gc()

stag.results <- readRDS("results/stag_sim_results.rds")
plot_mc.dist(stag.results,staggered=T,filename="figures/stag_param_est_dist.png")
rm(stag.results); gc()

resamp.com.results <- readRDS("results/resamp_com_sim_results.rds")
resamp.com.results <- resamp.com.results %>%
  ## Multiply the estimates (and SEs) by 100 for readability
  mutate(across(.cols=c(est,se,truth,lb,ub),.fns=~.x*100))
plot_mc.dist(resamp.com.results,staggered=F,filename="figures/com_resamp_est_dist.png")
rm(resamp.com.results); gc()

resamp.stag.results <- readRDS("results/resamp_stag_sim_results.rds")
resamp.stag.results <- resamp.stag.results %>%
  ## Multiply the estimates (and SEs) by 100 for readability
  mutate(across(.cols=c(est,truth,se),.fns=~.x*100))    
plot_mc.dist(resamp.stag.results,staggered=T,filename="figures/stag_resamp_est_dist.png")
rm(resamp.stag.results); gc()

#### Details in a single simulation ####
load("cleaned_force_data.RData")
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

set.seed(2283)
## make one draw from parametric and resampling and compare their variances
resamp.data <- resample(force.dat,n.units=100,starts=3:4,unbalanced=T) %>%
  mutate(Y=Y*100)
resamp.trt <- assign.treatment(resamp.data,starts=3:4,staggered=T)
resamp.trt.dat <- make.data(resamp.data,resamp.trt)
resamp.null.data <- add.trt.effects(resamp.trt.dat,tau=0, grp.var = F, time.var = F)

ggplot(resamp.null.data,aes(x=m,y=Y)) + geom_line(aes(group=unitID,col=factor(treatment))) + 
  geom_smooth(aes(group=factor(treatment),col=factor(treatment)),se=F)
# These are mostly 0s at the unit level

param.data <- generate.data(myparams,month.byunit=T,quarter.byunit=T,unbalanced=F)
param.trt <- assign.treatment(param.data,starts=3:4,staggered=T)
param.trt.dat <- make.data(param.data,param.trt)
param.null.data <- add.trt.effects(param.trt.dat,tau=0, grp.var=F,time.var=F)

ggplot(param.null.data,aes(x=m,y=Y)) + geom_line(aes(group=unitID,col=factor(treatment))) + 
  geom_smooth(aes(group=factor(treatment),col=factor(treatment)),se=F)
# because of the "wandering away" autocorrelation structure, the variance increases over time

# Plot mean and SD of the outcome over time at different time aggregations
resamp.stats <- resamp.null.data %>% 
  group_by(start.year,year) %>% mutate(y.mean=mean(Y),y.se=sqrt(var(Y))) %>% 
  group_by(start.year,q) %>% mutate(q.mean=mean(Y),q.se=sqrt(var(Y))) %>% 
  group_by(start.year,m) %>% mutate(m.mean=mean(Y),m.se=sqrt(var(Y))) %>%
  group_by(start.year,year,q,m) %>% slice(1) %>% ungroup() %>% 
  select(start.year,year,m,q,y.mean:m.se) %>% pivot_longer(y.mean:m.se) %>% mutate(outcomes="resampling") %>%
  separate(name,into=c('agg','stat')) %>% pivot_wider(names_from=stat,values_from=value)

param.stats <- param.null.data %>% 
  group_by(start.year,year) %>% mutate(y.mean=mean(Y),y.se=sqrt(var(Y))) %>% 
  group_by(start.year,q) %>% mutate(q.mean=mean(Y),q.se=sqrt(var(Y))) %>% 
  group_by(start.year,m) %>% mutate(m.mean=mean(Y),m.se=sqrt(var(Y))) %>%
  group_by(start.year,year,q,m) %>% slice(1) %>% ungroup() %>% 
  select(start.year,year,m,q,y.mean:m.se) %>% pivot_longer(y.mean:m.se) %>% mutate(outcomes="parametric") %>%
  separate(name,into=c('agg','stat')) %>% pivot_wider(names_from=stat,values_from=value)

sum.stats <- rbind(resamp.stats,param.stats) %>% 
  mutate(lb=mean-qnorm(.975)*se,
         ub=mean+qnorm(.975)*se,
         grp=ifelse(start.year==4,'late',ifelse(start.year==3,"early","never-treated")),
         time=ifelse(agg=="y",year,ifelse(agg=="q",q,m)))

# Show mean plus or minus SD
ggplot(sum.stats,aes(x=time)) + geom_hline(yintercept=0,col='darkgrey') + 
  geom_point((aes(y=mean,col=grp)),position=position_dodge(width=0.5)) + 
  geom_segment(aes(y=lb,yend=ub,col=grp),position=position_dodge(width=0.5)) + 
  facet_grid(outcomes~agg,scale="free") +
  scale_y_continuous("Mean outcome")
ggsave("figures/mean_se_of_outcome.png",width=8,height=6)

## Show how this translates into estimates:
resamp.ests <- analyze.data.CS(resamp.null.data)
param.ests <- analyze.data.CS(param.null.data)
  
sum.ests <- rbind(resamp.ests %>% mutate(outcomes="resampling"),
                  param.ests %>% mutate(outcomes="parametric")) %>%
  filter(time!="post") %>%
  mutate(time=as.numeric(time),
         lb=est-qnorm(.975)*se,
         ub=est+qnorm(.975)*se) %>%
  filter(estimand=="Time-varying")

# Show estimate plus or minus SE
ggplot(sum.ests,aes(x=time)) + geom_hline(yintercept=0,col='darkgrey') +
  geom_point(aes(y=est)) + geom_segment(aes(y=lb,yend=ub)) +facet_grid(outcomes~agg,scale="free") 
ggsave("figures/mean_se_of_estimate.png",width=8,height=6)

## Plot the same summary stats, but for the *difference* between treated and comparison
resamp.stats <- resamp.null.data %>% 
  group_by(start.year,year) %>% mutate(y.mean=mean(Y),y.var=var(Y),y.n=n()) %>%
  group_by(start.year,q) %>% mutate(q.mean=mean(Y),q.var=var(Y),q.n=n()) %>% 
  group_by(start.year,m) %>% mutate(m.mean=mean(Y),m.var=var(Y),m.n=n()) %>%
  group_by(start.year,year,q,m) %>% slice(1) %>% ungroup() %>%
  select(m,q,year,start.year,y.mean:m.n) %>%
  pivot_wider(names_from=start.year,values_from=y.mean:m.n) %>%
  mutate(y.mean.late=y.mean_4-y.mean_0,q.mean.late=q.mean_4-q.mean_0,m.mean.late=m.mean_4-m.mean_0,
         y.mean.early=y.mean_3-y.mean_0,q.mean.early=q.mean_3-q.mean_0,m.mean.early=m.mean_3-m.mean_0,
         y.se.late=sqrt(y.var_4/y.n_4 + y.var_0/y.n_0),y.se.early=sqrt(y.var_3/y.n_3 + y.var_0/y.n_0),
         q.se.late=sqrt(q.var_4/q.n_4 + q.var_0/q.n_0),q.se.early=sqrt(q.var_3/q.n_3 + q.var_0/q.n_0),
         m.se.late=sqrt(m.var_4/m.n_4 + m.var_0/m.n_0),m.se.early=sqrt(m.var_3/m.n_3 + m.var_0/m.n_0)) %>%
  select(year,m,q,y.mean.late:m.se.early) %>% pivot_longer(y.mean.late:m.se.early) %>% mutate(outcomes="resampling") %>%
  separate(name,into=c('agg','stat','grp')) %>% pivot_wider(names_from=stat,values_from=value)

param.stats <- param.null.data %>% 
  group_by(start.year,year) %>% mutate(y.mean=mean(Y),y.var=var(Y),y.n=n()) %>%
  group_by(start.year,q) %>% mutate(q.mean=mean(Y),q.var=var(Y),q.n=n()) %>% 
  group_by(start.year,m) %>% mutate(m.mean=mean(Y),m.var=var(Y),m.n=n()) %>%
  group_by(start.year,year,q,m) %>% slice(1) %>% ungroup() %>%
  select(m,q,year,start.year,y.mean:m.n) %>%
  pivot_wider(names_from=start.year,values_from=y.mean:m.n) %>%
  mutate(y.mean.late=y.mean_4-y.mean_0,q.mean.late=q.mean_4-q.mean_0,m.mean.late=m.mean_4-m.mean_0,
         y.mean.early=y.mean_3-y.mean_0,q.mean.early=q.mean_3-q.mean_0,m.mean.early=m.mean_3-m.mean_0,
         y.se.late=sqrt(y.var_4/y.n_4 + y.var_0/y.n_0),y.se.early=sqrt(y.var_3/y.n_3 + y.var_0/y.n_0),
         q.se.late=sqrt(q.var_4/q.n_4 + q.var_0/q.n_0),q.se.early=sqrt(q.var_3/q.n_3 + q.var_0/q.n_0),
         m.se.late=sqrt(m.var_4/m.n_4 + m.var_0/m.n_0),m.se.early=sqrt(m.var_3/m.n_3 + m.var_0/m.n_0)) %>%
  select(year,m,q,y.mean.late:m.se.early) %>% pivot_longer(y.mean.late:m.se.early) %>% mutate(outcomes="parametric") %>%
  separate(name,into=c('agg','stat','grp')) %>% pivot_wider(names_from=stat,values_from=value)

diff.stats <- rbind(resamp.stats,param.stats) %>% 
  mutate(lb=mean-qnorm(.975)*se,
         ub=mean+qnorm(.975)*se,
         time=ifelse(agg=="y",year,ifelse(agg=="q",q,m)))

# Show the differences 
ggplot(diff.stats,aes(x=time)) + geom_hline(yintercept=0,col='darkgrey') + 
  geom_point((aes(y=mean,col=grp)),position=position_dodge(width=0.4)) + geom_segment(aes(y=lb,yend=ub,col=grp),position=position_dodge(width=0.5)) + 
  facet_grid(outcomes~agg,scale="free") +
  scale_y_continuous("Difference from never-treated")
ggsave("figures/mean_se_of_tx-ctrl_diff.png",width=8,height=6)

