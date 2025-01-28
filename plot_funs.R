process_results_com <- function(results){
  est.prelim <- filter(results,!name%in%c('Joint F','Trend')) %>% 
      mutate(estimand=ifelse(time=="post","Overall","Time-varying"),
             covers=ifelse((lb<truth)&(ub>truth),1,0),
             bias=est-truth,
             stat.sig=pval<0.05,
             sq.err=(est-truth)^2) %>%
      separate(group,into=c('agg','model'),sep="_") %>%
      group_by(panel,assignment,trteff,estimand,agg,model) %>% 
      summarize(mc.sd=sqrt(var(est)),
                across(c(est,se,truth,covers:sq.err),\(x) mean(x,na.rm=T)),.groups="keep") %>% 
      ungroup() %>% mutate(rmse=sqrt(sq.err))
    
    test.prelim <- filter(results,name%in%c('Joint F','Trend')) %>% 
      # Compute the metrics relevant to individual effect estimates
      mutate(stat.sig=pval<.05,
             estimand="Time-varying") %>%
      separate(group,into=c('agg','model'),sep="_") %>%
      group_by(name,panel,assignment,trteff,estimand,agg,model) %>% 
      summarize(stat.sig=mean(stat.sig),.groups="keep") %>% ungroup()
 
  com.est.summaries <- est.prelim %>%
    mutate(agg=factor(agg,levels=c('weekly','monthly','quarterly','yearly'),
                      labels=c('week','month','quarter','year')),
           model=factor(model,levels=c('2x2','0','uy')),
           trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                         labels=c('null','constant','time-varying','group-varying','group-time')),
           label=paste(estimand,panel,sep="\n"))
  
  # Results from the joint tests and single overall tests
  com.test.summaries <- test.prelim %>%
    mutate(agg=factor(agg,levels=c('weekly','monthly','quarterly','yearly'),
                      labels=c('week','month','quarter','year')),
           model=factor(model,levels=c('2x2','0','uy')),
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
  return(list("est"=com.est.summaries,"test"=com.test.summaries))
}

process_results_stag <- function(results){
  stag.est.summaries <- results %>%
    # For comparability with the common scenarios, just look at post-period effects
    filter(time >= 0) %>%
    # Compute the metrics relevant to individual effect estimates
    mutate(lb=est-qnorm(.975)*se,
           ub=est+qnorm(.975)*se,
           covers=ifelse((lb<truth)&(ub>truth),1,0),
           bias=est-truth,
           stat.sig=((lb<0)&(ub<0))|((lb>0)&(ub>0)),
           sq.err=(est-truth)^2,
           model="CSA") %>%
    # Average over iterations, but not time (yet)
    group_by(panel,assignment,trteff,estimand,agg,model) %>% 
    summarize(mc.sd=sqrt(var(est)),
              across(c(est,se,truth,covers:sq.err),\(x) mean(x,na.rm=T)),.groups="keep") %>% ungroup() %>%
    mutate(agg=factor(agg,levels=c('month','quarter','year')),
           trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                         labels=c('null','constant','time-varying','group-varying','group-time')),
           power=ifelse(trteff!="null",stat.sig,NA),
           t1e=ifelse(trteff=="null",stat.sig,NA)) %>%
    mutate(rmse=sqrt(sq.err))
  return(list('est'=stag.est.summaries))
}

score_mod_perf <- function(summaries){
  test_perf <- summaries$test %>% 
    select(panel,assignment,trteff,estimand,agg,model,name,power,t1e) %>%
    mutate(name=ifelse(name=="Joint F","JointF",name),
           truth=NA) %>%
    pivot_wider(names_from=name,values_from=c(power,t1e),names_sep=".") %>%
    pivot_longer(cols=c(power.Trend,power.Individual,power.JointF,
                        t1e.Trend,t1e.Individual,t1e.JointF),names_to="metric")
  
  est_perf <- summaries$est %>% 
    select(panel,assignment,trteff,estimand,agg,model,truth,covers,bias,se,rmse,mc.sd) %>%
    pivot_longer(cols=c(covers,bias,se,rmse,mc.sd),names_to="metric")
  
  perf <-  bind_rows(test_perf,est_perf)
  
  # this scores models at each aggregation level
  scored_perf <- perf %>%
    arrange(panel,assignment,trteff,estimand,metric,agg,abs(value)) %>%
    group_by(panel,assignment,trteff,estimand,metric,agg) %>%
    # Find the best-performing (or nominal) level: 
    mutate(best=ifelse(metric=="covers",.95,
                       ifelse(grepl("t1e",metric),.05,
                              ifelse(grepl('power',metric),last(value,na_rm=TRUE),
                                     # For the rest, first (smallest) is best
                                     ifelse(metric%in%c('bias','mc.sd','rmse','se'),first(value,na_rm=TRUE),NA))))) %>% ungroup() %>%
    mutate(decrement=abs(abs(value)-abs(best)),
           # Percent decrement if on SD scale
           pct.decrement=ifelse(metric%in%c('mc.sd','rmse','se'),decrement/abs(best),decrement),
           label=factor(metric,levels=c('bias','mc.sd','rmse','se','covers',
                                        'power.Individual','power.JointF','power.Trend',
                                        't1e.Individual','t1e.JointF','t1e.Trend'),
                        labels=c('Bias','MC SD','RMSE','SE','Coverage',
                                 'Power (indiv)','Power (joint)','Power (trend)',
                                 'Alpha (indiv)','Alpha (joint)','Alpha (trend)')))
  return(scored_perf)
}

score_agg_perf <- function(summaries,staggered){
  est_perf <- summaries$est %>% 
    select(panel,assignment,trteff,estimand,agg,model,truth,covers,bias,se,rmse,mc.sd) %>%
    pivot_longer(cols=c(covers,bias,se,rmse,mc.sd),names_to="metric")
  
  if (!staggered){
    # Only look at the models with unit FE
    est_perf <- filter(est_perf,model=="0")
    
    test_perf <- summaries$test %>% filter(model=="0") %>%
      select(panel,assignment,trteff,estimand,agg,model,name,power,t1e) %>%
      mutate(name=ifelse(name=="Joint F","JointF",name),
             truth=NA) %>%
      pivot_wider(names_from=name,values_from=c(power,t1e),names_sep=".") %>%
      pivot_longer(cols=c(power.Trend,power.Individual,power.JointF,
                          t1e.Trend,t1e.Individual,t1e.JointF),names_to="metric")
    
    perf <-  bind_rows(test_perf,est_perf)
  } else {
    test_perf <- summaries$est %>% 
      select(panel,assignment,trteff,estimand,agg,model,power,t1e) %>%
      rename(power.Individual=power,t1e.Individual=t1e) %>%
      pivot_longer(cols=c(power.Individual,t1e.Individual),names_to="metric")
    perf <- bind_rows(test_perf,est_perf)
  }
  scored_perf <- perf %>%
    arrange(panel,assignment,trteff,estimand,metric,abs(value)) %>%
    group_by(panel,assignment,trteff,estimand,metric) %>%
    # Find the best-performing (or nominal) level: 
    mutate(best=ifelse(metric=="covers",.95,
                       ifelse(grepl("t1e",metric),.05,
                              ifelse(grepl('power',metric),last(value,na_rm=TRUE),
                                     # For the rest, first (smallest) is best
                                     ifelse(metric%in%c('bias','mc.sd','rmse','se'),first(value,na_rm=TRUE),NA))))) %>% ungroup() %>%
    # Don't give the joint test credit for power when its t1e is so bad
    mutate(decrement=abs(abs(value)-abs(best)),
           # Percent decrement if on SD scale
           pct.decrement=ifelse(metric%in%c('mc.sd','rmse','se'),decrement/abs(best),decrement),
           label=factor(metric,levels=c('bias','mc.sd','rmse','se',
                                        'covers',
                                        'power.Individual','power.JointF','power.Trend',
                                        't1e.Individual','t1e.JointF','t1e.Trend'),
                        labels=c('Bias','MC SD','RMSE','SE','Coverage',
                                 'Power (indiv)','Power (joint)','Power (trend)',
                                 'Alpha (indiv)','Alpha (joint)','Alpha (trend)')))
  return(scored_perf)
}

plot_mod_perf <- function(perf_dat,save.prefix){
  panel <- c('balanced','unbalanced')
  estimand <- c('Overall','Time-varying')
  combos <- expand.grid('panel'=panel,'estimand'=estimand)  
  for (i in 1:dim(combos)[1]){
    p1 <- ggplot(filter(perf_dat,panel==combos[i,'panel'],estimand==combos[i,'estimand'],
                  !is.na(value)),aes(x=model,y=label)) +
      geom_tile(aes(fill=pct.decrement)) + geom_text(aes(label=signif(value,2))) + 
      scale_color_manual(values=c('FALSE'='black','TRUE'='white'),guide="none") +
      scale_fill_gradient(low='white',high='darkgrey',guide="none",limits=c(0,.2),oob=scales::squish) +
      facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle=45,hjust=1))
    ggsave(p1,file=paste0(c(paste0(c(save.prefix,'mod',as.character(combos[i,'panel']),
                            as.character(combos[i,'estimand'])),collapse="_"),".png"),collapse=""),
           width=6.5,height=9)
    print(p1)
  }
}

plot_agg_perf <- function(perf_dat,save.prefix){
  this.est <- c('Overall','Time-varying')
  for (i in 1:2){
    p1 <- ggplot(filter(perf_dat,estimand==this.est[i],!is.na(value),metric!="power.JointF"),aes(x=agg,y=label)) +
      geom_tile(aes(fill=pct.decrement)) + geom_text(aes(label=signif(value,2))) + 
      scale_color_manual(values=c('FALSE'='black','TRUE'='white'),guide="none") +
      scale_fill_gradient(low='white',high='darkgrey',guide="none",limits=c(0,.2),oob=scales::squish) +
      facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle=45,hjust=1))  
    ggsave(p1,file=paste0(c(paste0(c(save.prefix,'agg',as.character(this.est[i])),collapse="_"),".png"),collapse=""),
           width=6.5,height=9)
    print(p1)
  }
}

agg_truth <- function(data,params,staggered){
  cons.data <- add.trt.effects(data,tau=params$tau,grp.var = F, time.var = F)
  tv.data   <- add.trt.effects(data,tau=params$tau,grp.var = F, time.var = T)
  gv.data   <- add.trt.effects(data,tau=params$tau,grp.var = T, time.var = F)
  gtv.data  <- add.trt.effects(data,tau=params$tau,grp.var = T, time.var = T)
  
  if (staggered){
    cons.agg <- agg.data.staggered(cons.data,is.gv=F)
    tv.agg <- agg.data.staggered(tv.data,is.gv=F)
    gv.agg <- agg.data.staggered(gv.data,is.gv=T)
    gtv.agg <- agg.data.staggered(gtv.data,is.gv=T)
    truth <- bind_rows(cons.agg$month_true %>% mutate(start.year=ceiling(start.month/12)),cons.agg$quarter_true %>% mutate(start.year=ceiling(start.quarter/4)),cons.agg$year_true) %>% 
      mutate(trteff="constant",grp=1) %>% select(-start.month,-start.quarter) %>%
      bind_rows(bind_rows(tv.agg$month_true %>% mutate(start.year=ceiling(start.month/12)),tv.agg$quarter_true %>% mutate(start.year=ceiling(start.quarter/4)),tv.agg$year_true) %>% 
                  mutate(trteff="time-varying",grp=1) %>% select(-start.month,-start.quarter)) %>%
      bind_rows(bind_rows(gv.agg$month_true %>% mutate(start.year=ceiling(start.month/12)),gv.agg$quarter_true %>% mutate(start.year=ceiling(start.quarter/4)),gv.agg$year_true) %>% 
                  mutate(trteff="group-varying") %>% select(-start.month,-start.quarter)) %>%
      bind_rows(bind_rows(gtv.agg$month_true %>% mutate(start.year=ceiling(start.month/12)),gtv.agg$quarter_true %>% mutate(start.year=ceiling(start.quarter/4)),gtv.agg$year_true) %>% 
                  mutate(trteff="group- and time-varying")  %>% select(-start.month,-start.quarter))
  } else {
    cons.agg <- agg.data.common(cons.data,is.tv=T,is.gv=F)
    tv.agg <- agg.data.common(tv.data,is.tv=T,is.gv=F)
    gv.agg <- agg.data.common(gv.data,is.tv=T,is.gv=T)
    gtv.agg <- agg.data.common(gtv.data,is.tv=T,is.gv=T)
    truth <- bind_rows(cons.agg$month_true,cons.agg$quarter_true,cons.agg$year_true) %>% mutate(trteff="constant",grp=1) %>%
      bind_rows(bind_rows(tv.agg$month_true,tv.agg$quarter_true,tv.agg$year_true) %>% mutate(trteff="time-varying",grp=1)) %>%
      bind_rows(bind_rows(gv.agg$month_true,gv.agg$quarter_true,gv.agg$year_true) %>% mutate(trteff="group-varying")) %>%
      bind_rows(bind_rows(gtv.agg$month_true,gtv.agg$quarter_true,gtv.agg$year_true) %>% mutate(trteff="group- and time-varying"))
  }
  return(truth)
}

plot_truth <- function(params,data,starts,resample=F,staggered=F){
  if (resample){
    bal.data <- resample(data,n.units=myparams$n.units,starts,unbalanced=F) 
    unbal.data <- resample(data,n.units=params$n.units,starts,unbalanced=T)
  } else {
    bal.data   <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=F)
    unbal.data <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=T)
  }
  bal.trt.dat <- assign.treatment(bal.data,starts=starts,staggered=staggered)
  bal.truth <- agg_truth(make.data(bal.data,bal.trt.dat),params,staggered) %>% mutate(panel="balanced")

  unbal.trt.dat <- assign.treatment(unbal.data,starts=starts,staggered)
  unbal.truth <- agg_truth(make.data(unbal.data,unbal.trt.dat),params,staggered) %>% mutate(panel="unbalanced")
  
  truth <- bind_rows(bal.truth,unbal.truth) %>%
    mutate(grp=factor(grp),
           trteff=factor(trteff,levels=c('constant','time-varying','group-varying','group- and time-varying')))
  return(truth)
}

plot_mc.dist <- function(results,staggered,filename){
  if (staggered){
    these.res <- results %>% filter(estimand=="Time-varying") %>% select(agg,trteff,adoption,panel,assignment,iter,time,est,truth) %>%
      mutate(trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying')))
  } else {
    these.res <- results %>% filter(time!="post",grepl("0",group)) %>% mutate(agg=ifelse(group=="monthly_0","month",ifelse(group=="quarterly_0","quarter","year"))) %>% 
      select(agg,trteff,adoption,panel,assignment,iter,time,est,truth) %>%
      mutate(trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying')))
  }
  
  tv.truth <- these.res %>% group_by(agg,trteff,panel,time) %>% summarize(truth=mean(truth),.groups="drop") %>% mutate(time=as.numeric(time))
  tv.est.quantiles <- these.res %>% group_by(agg,trteff,panel,time) %>% 
    summarize(q025=quantile(est,.025),q25=quantile(est,.25),q50=quantile(est,.50),q75=quantile(est,.75),q975=quantile(est,.975),.groups="drop") %>% mutate(time=as.numeric(time))
  
  thisplot <- ggplot(filter(tv.est.quantiles,panel=="balanced"),aes(x=time,group=panel)) + 
    geom_hline(yintercept = 0,col='darkgrey') + 
    geom_line(aes(y=q50)) + geom_line(aes(y=q25),lty=2) + geom_line(aes(y=q75),lty=2) +
    geom_line(aes(y=q025),lty=3) + geom_line(aes(y=q975),lty=3) +
    geom_line(data=filter(tv.truth,panel=="balanced"),aes(x=time,y=truth),col='red',lty=2) + 
    facet_grid(trteff~agg,scale="free_x") + scale_y_continuous("Treatment effect")
  ggsave(filename,thisplot,width=6.5,height=8)
}

table_winners <- function(summaries,bias,rmse,reject,prefix){
  winners <- summaries %>% filter(metric %in%c('power.Individual','bias','rmse','t1e.Individual')) %>% 
    # Apply separate thresholds for each metric
    mutate(win=ifelse(metric%in%c("power.Individual","t1e.Individual"),decrement < reject,
                      ifelse(metric=="bias",decrement < bias,
                             ifelse(metric=="rmse",decrement < bias,FALSE)))) %>%
    filter(win) %>%
    select(panel,estimand,trteff,metric,agg) %>% mutate(agg=factor(agg,levels=c('month','quarter','year'),labels=c('m','q','y'))) %>%
    arrange(panel,estimand,trteff,metric,agg) %>%
    pivot_wider(names_from=metric,values_from = agg,values_fn=~paste0(.x,collapse="")) %>% 
    mutate(reject=ifelse(trteff=="null",t1e.Individual,power.Individual)) %>%
    select(panel,estimand,trteff,bias,reject,rmse)
  print(winners)
  write_csv(winners,file=paste0(prefix,"_agg_winners.csv"))
}



