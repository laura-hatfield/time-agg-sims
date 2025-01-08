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
           sq.err=(est-truth)^2,
           model="CSA") %>%
    group_by(panel,assignment,trteff,estimand,agg,model) %>% 
    summarize(mc.sd=sqrt(var(est)),
              across(c(est,se,truth,covers:sq.err),\(x) mean(x,na.rm=T)),.groups="keep") %>% ungroup() %>%
    mutate(agg=factor(agg,levels=c('week','month','quarter','year')),
           trteff=factor(trteff,levels=c('null','constant','time-varying','group-varying','group- and time-varying'),
                         labels=c('null','constant','time-varying','group-varying','group-time')),
           power=ifelse(trteff!="null",stat.sig,NA),
           t1e=ifelse(trteff=="null",stat.sig,NA),
           rmse=sqrt(sq.err))
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
    arrange(panel,assignment,trteff,estimand,metric,agg,truth,abs(value)) %>%
    group_by(panel,assignment,trteff,estimand,metric,agg,truth) %>%
    # Find the best-performing level: for these, last is best
    mutate(best=ifelse(metric%in%c('covers','power.Individual','power.JointF','power.Trend'),last(value,na_rm=TRUE),
                       # For the rest, first is best
                       first(value,na_rm=TRUE))) %>% ungroup() %>%
    # Scale the difference in units of the treatment effect if not already bounded
    mutate(decrement=ifelse((metric%in%c('bias','rmse','se','mc.sd'))&(trteff!="null"),abs(abs(value)-abs(best))/abs(truth),abs(abs(value)-abs(best))),
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
  } else perf <- est_perf
  
  scored_perf <- perf %>%
    arrange(panel,assignment,trteff,estimand,metric,truth,abs(value)) %>%
    group_by(panel,assignment,trteff,estimand,metric,truth) %>%
    # Find the best-performing level: for these, last is best
    mutate(best=ifelse(metric%in%c('covers','power.Individual','power.JointF','power.Trend'),last(value,na_rm=TRUE),
                       # For the rest, first is best
                       first(value,na_rm=TRUE))) %>% ungroup() %>%
    mutate(scale=ifelse(metric%in%c('bias','rmse','se','mc.sd'),"Trt eff","Prop"),
           # Scale the difference in units of the treatment effect
           decrement=ifelse((scale=="Trt eff")&(trteff!="null"),abs(abs(value)-abs(best))/abs(truth),
                            # Or as absolute difference if truth is 0 or on scale of proportions
                            abs(abs(value)-abs(best))),
           label=factor(metric,levels=c('bias','mc.sd','rmse','se','covers',
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
      geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement>0.25)) + 
      scale_color_manual(values=c('FALSE'='white','TRUE'='black'),guide="none") +
      scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish,guide="none") +
      facet_grid(trteff~agg,scale="free") + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle=45,hjust=1))
    ggsave(p1,file=paste0(c(paste0(c('mod',save.prefix,as.character(combos[i,'panel']),
                            as.character(combos[i,'estimand'])),collapse="_"),".png"),collapse=""),
           width=6.5,height=9)
  }
}

plot_agg_perf <- function(perf_dat,save.prefix){
  this.est <- c('Overall','Time-varying')
  for (i in 1:2){
    p1 <- ggplot(filter(perf_dat,estimand==this.est[i],!is.na(value)),aes(x=agg,y=label)) +
      geom_tile(aes(fill=decrement)) + geom_text(aes(label=signif(value,2),col=decrement>0.25)) + 
      scale_color_manual(values=c('FALSE'='white','TRUE'='black'),guide="none") +
      scale_fill_gradient(low='#762a83',high='white',limits=c(0,1),oob=scales::squish,guide="none") +
      facet_grid(trteff~panel,scale="free_y") + xlab("") + ylab("") +
      theme(axis.text.x = element_text(angle=45,hjust=1))  
    ggsave(p1,file=paste0(c(paste0(c('agg',save.prefix,as.character(this.est[i])),collapse="_"),".png"),collapse=""),
           width=6.5,height=9)
  }
}


