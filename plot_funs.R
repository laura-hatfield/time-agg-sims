process_results_com <- function(results){
  est.prelim <- filter(results,!name%in%c('Joint F','Trend')) %>% 
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
           model=factor(model,levels=c('2x2','0','uy','uq','sq','um','sm','uw','sw')),
           trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                         labels=c('null','constant','time-varying','group-varying','group-time')),
           label=paste(estimand,panel,sep="\n"))
  
  # Results from the joint tests and single overall tests
  com.test.summaries <- test.prelim %>%
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
  
  return(list("est"=com.est.summaries,"test"=com.test.summaries))
}

process_results_stag <- function(results){
  stag.est.summaries <- results %>%
    # Compute the metrics relevant to individual effect estimates
    mutate(lb=est-qnorm(.975)*se,
           ub=est+qnorm(.975)*se,
           covers=ifelse((lb<truth)&(ub>truth),1,0),
           bias=est-truth,
           stat.sig=((lb<0)&(ub<0))|((lb>0)&(ub>0)),
           len.ci=ub-lb,
           sq.err=(est-truth)^2,
           model="CSA") %>%
    group_by(panel,assignment,trteff,estimand,agg,model) %>% 
    summarize(mc.sd=sqrt(var(est)),
              across(c(est,truth,covers:sq.err),\(x) mean(x,na.rm=T)),.groups="keep") %>% ungroup() %>%
    mutate(agg=factor(agg,levels=c('week','month','quarter','year')),
           trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                         labels=c('null','constant','time-varying','group-varying','group-time')),
           power=ifelse(trteff!="null",stat.sig,NA),
           t1e=ifelse(trteff=="null",stat.sig,NA),
           rmse=sqrt(sq.err))
  return(list('est'=stag.est.summaries))
}

score_model_winners <- function(summaries){
    test_perf <- summaries$test %>% 
      select(panel,assignment,trteff,estimand,agg,model,name,power,t1e) %>%
      mutate(name=ifelse(name=="Joint F","JointF",name)) %>%
      pivot_wider(names_from=name,values_from=c(power,t1e),names_sep=".") %>%
      pivot_longer(cols=c(power.Trend,power.Individual,power.JointF,
                          t1e.Trend,t1e.Individual,t1e.JointF),names_to="metric")


  est_perf <- summaries$est %>% 
    select(panel,assignment,trteff,estimand,agg,model,covers,bias,len.ci,rmse,mc.sd) %>%
    pivot_longer(cols=c(covers,bias,len.ci,rmse,mc.sd),names_to="metric")
  
  perf <-  bind_rows(test_perf,est_perf)
  
  # this scores models, not time aggregation levels 
  scored_perf <- perf %>%
    arrange(panel,assignment,trteff,estimand,metric,agg,abs(value)) %>%
    group_by(panel,assignment,trteff,estimand,metric,agg) %>%
    # Find the best-performing level: for these, last is best
    mutate(best=ifelse(metric%in%c('covers','power.Individual','power.JointF','power.Trend'),last(value,na_rm=TRUE),
                       # For the rest, first is best
                       first(value,na_rm=TRUE))) %>% ungroup() %>%
    # Use pct.decrement instead of winingness for colors in plots
    mutate(pct.decrement=ifelse(best!=0,abs(abs(value)-abs(best))/abs(best),abs(abs(value)-abs(best))),
           label=factor(metric,levels=c('mc.sd','bias','rmse','len.ci','covers',
                                        'power.Individual','power.JointF','power.Trend',
                                        't1e.Individual','t1e.JointF','t1e.Trend'),
                        labels=c('MC SD','Bias','RMSE','CI Length','Coverage',
                                 'Power (indiv)','Power (joint)','Power (trend)',
                                 'Alpha (indiv)','Alpha (joint)','Alpha (trend)')))
  return(scored_perf)
}


