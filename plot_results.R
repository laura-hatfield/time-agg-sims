com.results <- readRDS("common_sim_results.rds")

com.est.summaries <- filter(com.results,!name%in%c('Joint F','Trend')) %>% 
  # Compute the metrics relevant to individual effect estimates
  mutate(estimand=ifelse(time=="post","Overall","Time-varying"),
         covers=ifelse((lb<truth)&(ub>truth),1,0),
         abs.bias=abs(est-truth),
         stat.sig=pval<0.05,
         len.ci=ub-lb,
         sq.err=(est-truth)^2) %>%
  separate(group,into=c('agg','model'),sep="_") %>%
  group_by(panel,assignment,trteff,estimand,agg,model) %>% 
  summarize(mc.err=sqrt(var(est)),
            across(c(est,truth,covers:sq.err),\(x) mean(x,,na.rm=T)),.groups="keep") %>% 
  ungroup() %>% 
  mutate(rmse=sqrt(sq.err),
         agg=factor(agg,levels=c('weekly','monthly','quarterly','yearly'),
                    labels=c('week','month','quarter','year')),
         model=factor(model,levels=c('2x2','0','uy','uq','sq','um','sm','uw','sw')),
         trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                       labels=c('null','constant','time-varying','group-varying','group-time')),
         label=paste(estimand,panel,sep="\n"))
  
com.test.summaries <- filter(com.results,name%in%c('Joint F','Trend')) %>% 
  # Compute the metrics relevant to individual effect estimates
  mutate(stat.sig=pval<.05,
         estimand="Time-varying") %>%
  separate(group,into=c('agg','model'),sep="_") %>%
  group_by(name,panel,assignment,trteff,estimand,agg,model) %>% 
  summarize(stat.sig=mean(stat.sig),.groups="keep") %>% ungroup() %>%
  mutate(agg=factor(agg,levels=c('weekly','monthly','quarterly','yearly'),
                    labels=c('week','month','quarter','year')),
         model=factor(model,levels=c('2x2','0','uy','uq','sq','um','sm','uw','sw')),
         trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                       labels=c('null','constant','time-varying','group-varying','group-time'))) %>%
  # Bind on the individual time-varying effects averaged over post-period times
  bind_rows(com.est.summaries %>% mutate(name="Individual") %>%
              select(name,panel,assignment,trteff,estimand,agg,model,stat.sig)) %>%
  mutate(power=ifelse(name=="Trend",
                      ifelse(grepl("time",trteff),stat.sig,NA),
                      ifelse(trteff!="null",stat.sig,NA)),
         t1e=ifelse(name=="Trend",
                    ifelse(!grepl("time",trteff),stat.sig,NA),
                    ifelse(trteff=="null",stat.sig,NA)),
         label=paste(name,panel,sep="\n"))




# Results from the joint tests and single overall tests
com.test.summaries <- bind_rows(test.prelim) %>%
  

save(com.est.summaries,com.test.summaries,file="2024-08-16_common_sim_summaries.RData")
```

```{r load_stag_sim_results, eval=FALSE, message = FALSE, warning=FALSE, include=FALSE}
stag_effects <- readRDS("2024-08-16_stag_sim_results.rds")

# The sims object is really big, so do it in a loop over the four balance x assignment conditions:
these <- expand.grid(panel = c("balanced","unbalanced"),
                     assignment=c("random","trended"))
smaller <- list()
for (i in 1:dim(these)[1]){
  smaller[[i]] <- stag_effects %>% filter(panel==these[i,'panel'],assignment==these[i,'assignment'])
}
rm(stag_effects); gc()
est.prelim <- list()

if (F){
  for (i in 1:dim(these)[1]){
    # Exploring variability in the CS'A estimator for week aggregation
    ggplot(smaller[[4]] %>% filter(estimand=="Time-varying",trteff=="null") %>%
             mutate(agg=factor(agg,levels=c('week','month','quarter','year'))),aes(x=time,y=est,group=time)) + geom_boxplot() + facet_wrap(~agg,scale="free_x") +
      ggtitle(paste(unlist(these[i,c('panel','assignment')]),collapse=", "))
    if (save.figs) ggsave(file.path(figs.dir,paste(c(paste(unlist(these[i,c('panel','assignment')]),collapse="_"),"csa_estimate_boxplots.png"),collapse="_")),width=11,height=8.5)
    
  }
}


for (i in 1:dim(these)[1]){
  est.prelim[[i]] <- smaller[[i]] %>%
    # Compute the metrics relevant to individual effect estimates
    mutate(lb=est-qnorm(.975)*se,
           ub=est+qnorm(.975)*se,
           covers=ifelse((lb<truth)&(ub>truth),1,0),
           abs.bias=abs(est-truth),
           stat.sig=((lb<0)&(ub<0))|((lb>0)&(ub>0)),
           len.ci=ub-lb,
           sq.err=(est-truth)^2,
           model="CSA") %>%
    group_by(panel,assignment,trteff,estimand,agg,model) %>% 
    summarize(mc.err=sqrt(var(est)),across(c(est,truth,covers:sq.err),\(x) mean(x,,na.rm=T)),.groups="keep") %>% ungroup()
}
rm(smaller); gc()

# Bind over the four combinations, reconcile factors 
stag.est.summaries <- bind_rows(est.prelim) %>%
  mutate(agg=factor(agg,levels=c('week','month','quarter','year')),
         trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                       labels=c('null','constant','time-varying','group-varying','group-time')),
         power=ifelse(trteff!="null",stat.sig,NA),
         t1e=ifelse(trteff=="null",stat.sig,NA),
         rmse=sqrt(sq.err))

save(stag.est.summaries,file="2024-08-16_stag_sim_summaries.RData")
```

# Common timing heatplots
```{r heat_common_data}
load("2024-06-23_common_sim_summaries.RData")
```

```{r mod_winners}
mod_test_perf <- com.test.summaries %>% 
  select(panel,assignment,trteff,estimand,agg,model,name,power,t1e) %>%
  mutate(name=ifelse(name=="Joint F","JointF",name)) %>%
  pivot_wider(names_from=name,values_from=c(power,t1e),names_sep=".") %>%
  pivot_longer(cols=c(power.Trend,power.Individual,power.JointF,
                      t1e.Trend,t1e.Individual,t1e.JointF),names_to="metric")
mod_est_perf <- com.est.summaries %>% 
  select(panel,assignment,trteff,estimand,agg,model,covers,abs.bias,len.ci,rmse,mc.err) %>%
  pivot_longer(cols=c(covers,abs.bias,len.ci,rmse,mc.err),names_to="metric")

mod_performance <-  bind_rows(mod_test_perf,mod_est_perf) %>%
  arrange(panel,assignment,trteff,estimand,metric,agg,value) %>%
  group_by(panel,assignment,trteff,estimand,metric,agg) %>% 
  mutate(best=ifelse(metric%in%c('covers','power.Individual','power.JointF','power.Trend'),
                     last(value,na_rm=TRUE),first(value,na_rm=TRUE))) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(value-best)/abs(best),abs(value-best)))

mod_winners_prelim <- mod_performance %>%
  # Change the definition of winner so that for metrics with nominal levels, winning = achieving nominal
  # and for all others, we use the 5% difference from the best definition
  mutate(winner=ifelse(metric%in%c("rmse","abs.bias"),pct.decrement<.05, # smaller is better
                       ifelse(grepl("t1e",metric),signif(value,2)<=.05, # nominal wins
                              ifelse(metric=="covers",signif(value,2)>=.95,FALSE)))) # nominal wins


mod_ci_winners <- mod_winners_prelim %>% filter(metric%in%c('covers','len.ci')) %>%
  select(panel,assignment,trteff,estimand,agg,model,metric,value) %>% 
  pivot_wider(names_from=metric,values_from=value) %>% 
  filter(signif(covers,2)>=0.95) %>%
  arrange(panel,assignment,trteff,estimand,agg,len.ci) %>%
  group_by(panel,assignment,trteff,estimand,agg) %>%
  # CI length only wins among levels with nominal coverage
  mutate(best=first(len.ci,na_rm=TRUE)) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(len.ci-best)/abs(best),abs(len.ci-best)),
         ci.len.winner=pct.decrement<.05) %>% 
  mutate(metric="len.ci") %>% rename(value=len.ci) %>%
  select(panel,assignment,trteff,estimand,agg,model,metric,value,ci.len.winner)


mod_power_winners <- left_join(mod_winners_prelim %>% filter(trteff!="null",grepl('power',metric),!is.na(value)) %>%
                                 separate(metric,into=c("metric","specification")) %>%
                                 select(panel,assignment,trteff,estimand,agg,model,specification,metric,value),
                               mod_winners_prelim %>% filter(trteff=="null",grepl('t1e',metric),!is.na(value)) %>% 
                                 separate(metric,into=c("metric","specification")) %>%
                                 select(panel,assignment,estimand,agg,model,specification,value) %>%
                                 rename(t1e=value),by=c('panel','assignment','estimand','model','agg','specification')) %>%
  filter(signif(t1e,2)<=.2) %>%
  arrange(panel,assignment,trteff,estimand,agg,specification,value) %>%
  group_by(panel,assignment,trteff,estimand,agg,specification) %>%
  mutate(best=last(value,na_rm=T)) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(value-best)/abs(best),abs(value-best)),
         power.winner=pct.decrement<.05) %>%
  unite("metric",metric,specification,sep=".") %>%
  select(panel,assignment,trteff,estimand,agg,model,metric,value,power.winner)


# Join them back on:
mod_winners <- left_join(mod_winners_prelim,mod_ci_winners) %>%
  mutate(winner=ifelse(metric=="len.ci",ifelse(is.na(ci.len.winner),FALSE,ci.len.winner),winner)) %>%
  select(-ci.len.winner) %>%
  left_join(mod_power_winners) %>%
  mutate(winner=ifelse(grepl("power",metric),ifelse(is.na(power.winner),FALSE,power.winner),winner))

rm(mod_winners_prelim,mod_ci_winners,mod_est_perf,mod_test_perf,mod_performance)
```


```{r plot_model_winners}
# Put the value in the cell, 
mod_win_plot_dat <- filter(mod_winners,assignment=="random",!metric%in%c("mc.err","abs.bias")) %>%
  mutate(label=factor(metric,levels=c('rmse','len.ci','covers',
                                      'power.Individual','power.JointF','power.Trend',
                                      't1e.Individual','t1e.JointF','t1e.Trend'),
                      labels=c('RMSE','CI Length','Coverage',
                               'Power (indiv)','Power (joint)','Power (trend)',
                               'Alpha (indiv)','Alpha (joint)','Alpha (trend)')))

ggplot(filter(mod_win_plot_dat,panel=="balanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=round(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_model_winners_overall_balanced.png"),width=14,height=8.5)

ggplot(filter(mod_win_plot_dat,panel=="unbalanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=round(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_model_winners_overall_unbalanced.png"),width=14,height=8.5)

ggplot(filter(mod_win_plot_dat,panel=="balanced",estimand=="Time-varying",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=round(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_model_winners_time-varying_balanced.png"),width=14,height=8.5)

ggplot(filter(mod_win_plot_dat,panel=="unbalanced",estimand=="Time-varying",!is.na(value)),aes(x=model,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=round(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_model_winners_time-varying_unbalanced.png"),width=14,height=8.5)

```


# Common timing simplified
```{r plot_common_sq_error_simple}
ggplot(filter(com.est.summaries,model=="0",assignment=="random"),aes(x=agg,y=trteff)) + geom_tile(aes(fill=rmse)) +
  facet_grid(estimand~panel) +
  scale_fill_gradient("Root Mean Squared Error",low='#762a83',high='white') + xlab("") + ylab("") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_sq_error_simple.png"),width=6,height=6)
```

```{r plot_common_coverage_simple}
ggplot(filter(com.est.summaries,model=="0",assignment=="random"),aes(x=agg,y=trteff)) + geom_tile(aes(fill=covers)) +
  facet_grid(estimand~panel) +
  scale_fill_gradient("Coverage",high='#762a83',low='white') + xlab("") + ylab("") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_coverage_simple.png"),width=6,height=6)
```

```{r plot_common_ci_len_simple}
ggplot(filter(com.est.summaries,model=="0",assignment=="random"),aes(x=agg,y=trteff)) + geom_tile(aes(fill=len.ci)) +
  facet_grid(estimand~panel) +
  scale_fill_gradient("CI Length",low='#762a83',high='white') + xlab("") + ylab("") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_ci_len_simple.png"),width=6,height=6)
```

```{r plot_common_power_simple}
ggplot(filter(com.test.summaries,model=="0",assignment=="random",!is.na(power)),aes(x=agg,y=trteff)) + geom_tile(aes(fill=power)) +
  scale_fill_gradient("Power",high='#762a83',low='white') + xlab("") + ylab("") +
  facet_grid(name~panel,scale="free_y") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_power_simple.png"),width=6,height=6)
```

```{r plot_common_typeIerror_simple}
ggplot(filter(com.test.summaries,model=="0",assignment=="random",!is.na(t1e)),aes(x=agg,y=trteff)) + geom_tile(aes(fill=t1e)) +
  scale_fill_gradient("Type I error",low='#762a83',high='white') + xlab("") + ylab("") +
  facet_grid(name~panel,scale="free_y") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_typeIerror_simple.png"),width=6,height=6)
```

# Common winners plots
```{r com_winners}
## Rearrange the data into the table format:
com_test_perf <- com.test.summaries %>% select(panel,assignment,trteff,estimand,agg,model,name,power,t1e) %>%
  mutate(name=ifelse(name=="Joint F","JointF",name)) %>%
  filter(model=="0") %>% select(-model) %>%
  pivot_wider(names_from=name,values_from=c(power,t1e),names_sep=".") %>%
  pivot_longer(cols=c(power.Trend,power.Individual,power.JointF,
                      t1e.Trend,t1e.Individual,t1e.JointF),names_to="metric")
com_est_perf <- com.est.summaries %>% 
  select(panel,assignment,trteff,estimand,agg,model,covers,abs.bias,len.ci,rmse) %>%
  filter(model=="0") %>% select(-model) %>%
  pivot_longer(cols=c(covers,abs.bias,len.ci,rmse),names_to="metric")

com_performance <-  bind_rows(com_test_perf,com_est_perf) %>%
  arrange(panel,assignment,trteff,estimand,metric,value) %>%
  group_by(panel,assignment,trteff,estimand,metric) %>% 
  mutate(best=ifelse(grepl('covers|power',metric),last(value,na_rm=TRUE), #bigger is better
                     first(value,na_rm=TRUE))) %>% #smaller is better
  ungroup() %>% 
  # difference from best
  mutate(pct.decrement=ifelse(best!=0,abs(value-best)/abs(best),abs(value-best))) 

com_winners_prelim <- com_performance %>%
  mutate(winner=ifelse(metric%in%c("rmse","abs.bias"),pct.decrement<.05, # smaller is better
                       ifelse(grepl("t1e",metric),signif(value,2)<=.05, # nominal wins
                              ifelse(metric=="covers",signif(value,2)>=.95,FALSE)))) # nominal wins

com_ci_winners <- com_winners_prelim %>% filter(metric%in%c('covers','len.ci')) %>%
  select(panel,assignment,trteff,estimand,agg,metric,value) %>% 
  pivot_wider(names_from=metric,values_from=value) %>% 
  filter(signif(covers,2)>=0.95) %>%
  arrange(panel,assignment,trteff,estimand,len.ci) %>%
  group_by(panel,assignment,trteff,estimand) %>%
  # CI length only wins among levels with nominal coverage
  mutate(best=first(len.ci,na_rm=TRUE)) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(len.ci-best)/abs(best),abs(len.ci-best)),
         ci.len.winner=pct.decrement<.05) %>% 
  mutate(metric="len.ci") %>% rename(value=len.ci) %>%
  select(panel,assignment,trteff,estimand,agg,metric,value,ci.len.winner)

# Power only wins among levels with type I error <20%
com_power_winners <- left_join(com_winners_prelim %>% filter(trteff!="null",grepl('power',metric),!is.na(value)) %>%
                                 separate(metric,into=c("metric","specification")) %>%
                                 select(panel,assignment,trteff,estimand,agg,specification,metric,value),
                               com_winners_prelim %>% filter(trteff=="null",grepl('t1e',metric),!is.na(value)) %>% 
                                 separate(metric,into=c("metric","specification")) %>%
                                 select(panel,assignment,estimand,agg,specification,value) %>%
                                 rename(t1e=value),by=c('panel','assignment','estimand','agg','specification')) %>%
  filter(signif(t1e,2)<=.2) %>%
  arrange(panel,assignment,trteff,estimand,specification,value) %>%
  group_by(panel,assignment,trteff,estimand,specification) %>%
  mutate(best=last(value,na_rm=T)) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(value-best)/abs(best),abs(value-best)),
         power.winner=pct.decrement<.05) %>%
  unite("metric",metric,specification,sep=".") %>%
  select(panel,assignment,trteff,estimand,agg,metric,value,power.winner)

# Join them back on:
com_winners <- left_join(com_winners_prelim,com_ci_winners) %>%
  mutate(winner=ifelse(metric=="len.ci",ifelse(is.na(ci.len.winner),FALSE,ci.len.winner),winner)) %>%
  select(-ci.len.winner) %>%
  left_join(com_power_winners) %>%
  mutate(winner=ifelse(grepl("power",metric),ifelse(is.na(power.winner),FALSE,power.winner),winner))
rm(com_winners_prelim,com_ci_winners,com_est_perf,com_test_perf,com_performance)
```

```{r plot_common_winners}
# Put the value in the cell, 
com_win_plot_dat <- filter(com_winners,assignment=="random",metric!="abs.bias") %>%
  mutate(label=factor(metric,levels=c('rmse','len.ci','covers',
                                      'power.Individual','power.JointF','power.Trend',
                                      't1e.Individual','t1e.JointF','t1e.Trend'),
                      labels=c('RMSE','CI Length','Coverage',
                               'Power (indiv)','Power (joint)','Power (trend)',
                               'Alpha (indiv)','Alpha (joint)','Alpha (trend)')))

ggplot(filter(com_win_plot_dat,estimand=="Overall",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_winners_overall.png"),width=5.5,height=5)

ggplot(filter(com_win_plot_dat,estimand=="Time-varying",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_winners_time-varying.png"),width=5.5,height=6)

## Try a version that separates testing from estimation:
ggplot(filter(com_win_plot_dat,estimand=="Time-varying",!is.na(value),metric%in%c('covers','len.ci','rmse')),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(filter(com_win_plot_dat,estimand=="Time-varying",!is.na(value),!metric%in%c('covers','len.ci','rmse')),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))

```

# Staggered timing heatplots

```{r heat_staggered_data}
load("2024-06-23_stag_sim_summaries.RData")
```

```{r plot_stag_sq_error}
ggplot(filter(stag.est.summaries,assignment=="random"),aes(x=agg,y=trteff)) + geom_tile(aes(fill=rmse)) +
  facet_grid(estimand~panel) +
  scale_fill_gradient("Root Mean Squared Error",low='#762a83',high='white') + xlab("") + ylab("") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_sq_error.png"),width=6,height=6)
```

```{r plot_stag_coverage}
ggplot(filter(stag.est.summaries,assignment=="random"),aes(x=agg,y=trteff)) + geom_tile(aes(fill=covers)) +
  facet_grid(estimand~panel) +
  scale_fill_gradient("Coverage",high='#762a83',low='white') + xlab("") + ylab("") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_coverage.png"),width=6,height=6)
```

```{r plot_stag_ci_len}
ggplot(filter(stag.est.summaries,assignment=="random"),aes(x=agg,y=trteff)) + geom_tile(aes(fill=len.ci)) +
  facet_grid(estimand~panel) +
  scale_fill_gradient("CI Length",low='#762a83',high='white') + xlab("") + ylab("") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_ci_len.png"),width=6,height=6)
```

```{r plot_stag_power}
ggplot(filter(stag.est.summaries,assignment=="random",!is.na(power)),aes(x=agg,y=trteff)) + geom_tile(aes(fill=power)) +
  scale_fill_gradient("Power",high='#762a83',low='white') + xlab("") + ylab("") +
  facet_grid(estimand~panel,scale="free_y") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_power.png"),width=6,height=6)
```

```{r plot_stag_typeIerror}
ggplot(filter(stag.est.summaries,assignment=="random",!is.na(t1e)),aes(x=agg,y=trteff)) + geom_tile(aes(fill=t1e)) +
  scale_fill_gradient("Type I error",low='#762a83',high='white') + xlab("") + ylab("") +
  facet_grid(estimand~panel,scale="free_y") +
  theme(legend.position="bottom",axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_typeIerror.png"),width=6,height=6)
```

# Staggered winners plots
```{r stag_winners}
stag_performance <-  stag.est.summaries %>% 
  select(panel,assignment,trteff,estimand,agg,model,covers,abs.bias,len.ci,rmse,mc.err,t1e,power) %>%
  pivot_longer(cols=c(covers,abs.bias,len.ci,rmse,mc.err,t1e,power),names_to="metric") %>%
  arrange(panel,assignment,trteff,estimand,metric,value) %>%
  group_by(panel,assignment,trteff,estimand,metric) %>% 
  mutate(best=ifelse(metric%in%c('covers','power'),
                     last(value,na_rm=TRUE),first(value,na_rm=TRUE))) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(value-best)/abs(best),abs(value-best)))

stag_winners_prelim <- stag_performance %>%
  # Change the definition of winner so that for metrics with nominal levels, winning = achieving nominal
  # and for all others, we use the 5% difference from the best definition
  mutate(winner=ifelse(metric%in%c("rmse","abs.bias","mc.err"),pct.decrement<.05,
                       ifelse(grepl("power",metric),ifelse(best>=.8,signif(value,2)>=.8,pct.decrement<.05),
                              ifelse(grepl("t1e",metric),signif(value,2)<=.05,
                                     ifelse(metric=="covers",signif(value,2)>=.95,FALSE)))))

stag_ci_winners <- stag_winners_prelim %>% filter(metric%in%c('covers','len.ci')) %>%
  select(panel,assignment,trteff,estimand,agg,metric,value) %>% 
  pivot_wider(names_from=metric,values_from=value) %>% 
  filter(signif(covers,2)>=0.95) %>%
  arrange(panel,assignment,trteff,estimand,len.ci) %>%
  group_by(panel,assignment,trteff,estimand) %>%
  # CI length only wins among levels with nominal coverage
  mutate(best=first(len.ci,na_rm=TRUE)) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(len.ci-best)/abs(best),abs(len.ci-best)),
         ci.len.winner=pct.decrement<.05) %>% 
  mutate(metric="len.ci") %>% rename(value=len.ci) %>%
  select(panel,assignment,trteff,estimand,agg,metric,value,ci.len.winner)

# Join them back on:
stag_winners <- left_join(stag_winners_prelim,stag_ci_winners) %>%
  mutate(winner=ifelse(metric=="len.ci",ifelse(is.na(ci.len.winner),FALSE,ci.len.winner),winner)) %>%
  select(-ci.len.winner)
```

```{r plot_stag_winners}
# Put the value in the cell, 
stag_win_plot_dat <- filter(stag_winners,assignment=="random",!metric%in%c("abs.bias","mc.err")) %>%
  mutate(label=factor(metric,levels=c('rmse','len.ci','covers','power','t1e'),
                      labels=c('RMSE','CI Length','Coverage','Power','Alpha')))

ggplot(filter(stag_win_plot_dat,estimand=="Overall",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_winners_overall.png"),width=5.5,height=5)

ggplot(filter(stag_win_plot_dat,estimand=="Time-varying",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_winners_time-varying.png"),width=5.5,height=5)

```



# Common winners (biased assignment)
```{r plot_common_winners_biased}
# Put the value in the cell, 
com_win_plot_dat <- filter(com_winners,assignment=="trended",metric!="abs.bias") %>%
  mutate(label=factor(metric,levels=c('rmse','len.ci','covers',
                                      'power.Individual','power.JointF','power.Trend',
                                      't1e.Individual','t1e.JointF','t1e.Trend'),
                      labels=c('RMSE','CI Length','Coverage',
                               'Power (indiv)','Power (joint)','Power (trend)',
                               'Alpha (indiv)','Alpha (joint)','Alpha (trend)')))

ggplot(filter(com_win_plot_dat,estimand=="Overall",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_winners_overall_trended.png"),width=5.5,height=5)

ggplot(filter(com_win_plot_dat,estimand=="Time-varying",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"common_winners_time-varying_trended.png"),width=5.5,height=6)

```

# Staggered winners (biased assignment)
```{r plot_stag_winners_biased}
# Put the value in the cell, 
stag_win_plot_dat <- filter(stag_winners,assignment=="trended",!metric%in%c("abs.bias",'mc.err')) %>%
  mutate(label=factor(metric,levels=c('rmse','len.ci','covers','power','t1e'),
                      labels=c('RMSE','CI Length','Coverage','Power','Alpha')))

ggplot(filter(stag_win_plot_dat,estimand=="Overall",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_winners_overall_trended.png"),width=5.5,height=5)

ggplot(filter(stag_win_plot_dat,estimand=="Time-varying",!is.na(value)),aes(x=agg,y=label)) +
  geom_tile(aes(fill=winner)) + geom_text(aes(label=signif(value,2),col=winner)) + 
  scale_color_manual(values=c('TRUE'='white','FALSE'='darkgrey'),guide="none") +
  scale_fill_manual(values=c('TRUE'='#762a83','FALSE'='white'),guide="none") +
  facet_grid(trteff~panel,scale="free") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
if (save.figs) ggsave(file.path(figs.dir,"stag_winners_time-varying_trended.png"),width=5.5,height=5)

```

