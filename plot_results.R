library(tidyverse)
theme_set(theme_minimal())
com.results <- readRDS("common_sim_results.rds")

com.est.summaries <- filter(com.results,!name%in%c('Joint F','Trend')) %>% 
  # Compute the metrics relevant to individual effect estimates
  mutate(estimand=ifelse(time=="post","Overall","Time-varying"),
         covers=ifelse((lb<truth)&(ub>truth),1,0),
         bias=est-truth,
         stat.sig=pval<0.05,
         len.ci=ub-lb,
         sq.err=(est-truth)^2) %>%
  separate(group,into=c('agg','model'),sep="_") %>%
  group_by(panel,assignment,trteff,estimand,agg,model) %>% 
  summarize(mc.sd=sqrt(var(est)),
            across(c(est,truth,covers:sq.err),\(x) mean(x,na.rm=T)),.groups="keep") %>% 
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

mod_test_perf <- com.test.summaries %>% 
  select(panel,assignment,trteff,estimand,agg,model,name,power,t1e) %>%
  mutate(name=ifelse(name=="Joint F","JointF",name)) %>%
  pivot_wider(names_from=name,values_from=c(power,t1e),names_sep=".") %>%
  pivot_longer(cols=c(power.Trend,power.Individual,power.JointF,
                      t1e.Trend,t1e.Individual,t1e.JointF),names_to="metric")
mod_est_perf <- com.est.summaries %>% 
  select(panel,assignment,trteff,estimand,agg,model,covers,bias,len.ci,rmse,mc.sd) %>%
  pivot_longer(cols=c(covers,bias,len.ci,rmse,mc.sd),names_to="metric")

mod_performance <-  bind_rows(mod_test_perf,mod_est_perf) %>%
  arrange(panel,assignment,trteff,estimand,metric,agg,value) %>%
  group_by(panel,assignment,trteff,estimand,metric,agg) %>% 
  mutate(best=ifelse(metric%in%c('covers','power.Individual','power.JointF','power.Trend'),
                     last(value,na_rm=TRUE),first(abs(value),na_rm=TRUE))) %>% ungroup() %>%
  mutate(pct.decrement=ifelse(best!=0,abs(abs(value)-abs(best))/abs(best),abs(abs(value)-abs(best))))

mod_winners_prelim <- mod_performance %>%
  # Change the definition of winner so that for metrics with nominal levels, winning = achieving nominal
  # and for all others, we use the 5% difference from the best definition
  mutate(winner=ifelse(metric%in%c("rmse","bias"),pct.decrement<.05, # smaller is better
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

mod_win_plot_dat <- filter(mod_winners,assignment=="random") %>%
  mutate(label=factor(metric,levels=c('mc.sd','bias','rmse','len.ci','covers',
                                      'power.Individual','power.JointF','power.Trend',
                                      't1e.Individual','t1e.JointF','t1e.Trend'),
                      labels=c('MC SD','Bias','RMSE','CI Length','Coverage',
                               'Power (indiv)','Power (joint)','Power (trend)',
                               'Alpha (indiv)','Alpha (joint)','Alpha (trend)')))

ggplot(filter(mod_win_plot_dat,panel=="balanced",estimand=="Overall",!is.na(value)),aes(x=model,y=label)) +
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
