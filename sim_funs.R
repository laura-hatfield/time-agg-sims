#' Generate monthly shocks 
#' 
#' Generate monthly shocks from iid normal
#' 
#' @description
#' `draw.months` generates monthly shocks from an ARIMA process
#'
#' @param params list with named elements: total_num_month,sigma_month,const,phi,total_num_quarter,n.year
#' @returns tibble with the following variables: Y.month, m, q
#' @details The parameter `phi` is a 2-vector containing the 2 autoregressive parameters; 
#' `phi[1]` controls likelihood of wandering away while `phi[2]` controls "acceleration" (how fast)
#' The parameter `sigma_month` is a white noise error SD
#' The parameter `const` is an overall non-stationary trend

draw.months <- function(params){
  # Create empty object of the correct dimension
  y <- Y <- rep(NA,params$total_num_month)
  # Initialize the time series
  y[1] <- params$const + rnorm(1,0,params$sigma_month)
  Y[1] <- y[1]
  y[2] <- params$const + rnorm(1,0,params$sigma_month)
  Y[2] <- y[2] + 2*Y[1]
  # Iterate over months
  for (m in 3:params$total_num_month){
    y[m] <- params$const + params$phi[1]*y[m-1] + 
      params$phi[2]*y[m-2] + rnorm(1,0,params$sigma_month)
    Y[m] <- y[m] + Y[m-1] 
  }
  months.per.quarter <- params$total_num_month/params$total_num_quarter
  months <- tibble(Y.month=Y,m=1:params$total_num_month,
                   q=rep(1:params$total_num_quarter,each=months.per.quarter))
  return(months)
}


#' Generate quarterly shocks from iid normal
#' 
#' @description
#' `draw.quarters` generates iid normal quarter shocks but sorted so that they are always increasing
#'
#' @param params list with named elements: n.quarter,sigma_quarter,n.year,
#' @returns tibble with the following variables: Y.quarter,q,year

draw.quarters <- function(params){
  ## Draw a random vector of quarter effects
  temp <- rnorm(4,0,params$sigma_quarter)
  # Sort so that the seasonal effects is always increasing over the year
  quarter.effects <- temp[order(temp)]
  quarters <- tibble(Y.quarter=rep(quarter.effects,params$n.year),
                     q=rep(1:params$total_num_quarter),
                     year=rep(1:params$n.year,each=4))
  return(quarters)
}

#' Generate quarterly shocks from iid normal
#' 
#' @description
#' `generate.data` generates autocorrelated data for n.units over n.year
#'
#' @details Generates quarterly shocks
#' @param params list with named elements: n.units, n.year, sigma_quarter, sigma_month
#' @param month.byunit logical indicating whether each unit has its own monthly shocks
#' @param quarter.byunit logical indicating whether each unit has its own quarterly shocks
#' @param unbalanced logical indicating whether the data should be artificially unbalanced
#' @returns tibble with the following variables: Y, m, q, year

generate.data <- function(params,
                          month.byunit=F,
                          quarter.byunit=F,
                          unbalanced=F){
  # Convenience computations for use in function calls
  params$total_num_month = 12 * params$n.year
  params$total_num_quarter = 4 * params$n.year
  
  # Draw month vectors
  if (month.byunit){
    monthdat <- replicate(params$n.units,draw.months(params),simplify=FALSE)
    monthdat <- bind_rows(monthdat, .id="unitID")
  } else {
    months <- draw.months(params)
    monthdat <- tibble(unitID=as.character(1:params$n.units)) %>%
      cross_join(months) %>% arrange(unitID,m)
  }
  data <- monthdat
  
  # Draw quarter vectors
  if (quarter.byunit){
    quarterdat <- replicate(params$n.units,draw.quarters(params),simplify=FALSE)
    quarterdat <- bind_rows(quarterdat,.id="unitID")
  } else {
    quarters <- draw.quarters(params)
    quarterdat <- tibble(unitID=as.character(1:params$n.units)) %>% 
      cross_join(quarters) %>% arrange(unitID,q)
  }
  # Bind to the previous data
  data <- left_join(data,quarterdat,by=c('unitID','q'))
  
  # Add up the effects to get the outcome data
  data <- data %>% mutate(Y=Y.month+Y.quarter)
  
  # Unbalance the panel
  if (unbalanced){
    # For each unit, draw a random first month (with prob = 0.5) and last month (with prob = 0.5)
    data <- data %>% group_by(unitID) %>% 
      mutate(trunc.begin=rbinom(1,1,.5),
             first.month=ifelse(trunc.begin,sample(1:((params$total_num_month/2)-1),1),1),
             trunc.end=rbinom(1,1,.5),
             last.month=ifelse(trunc.end,sample(((params$total_num_month/2)+1):(params$total_num_month-1),1),params$total_num_month)) %>%
      ungroup() %>% filter((m>=first.month) & (m<=last.month))
  }
  
  return(data %>% select(unitID, Y, m, q, year))
}

resample <- function(data,n.units,starts,unbalanced){
  # Make long data for randomly sampled (with replacement) units
  units <- unique(data$unitID)
  sampled <- tibble('unitID'= sample(units,n.units,replace=T),
                    'newid'=1:n.units) %>% 
    left_join(data %>% select(unitID,Y,m,q,year),by="unitID",
              relationship="many-to-many") %>% 
    select(-unitID) %>% rename(unitID=newid)
  if (unbalanced){
    # For each unit, draw a random first month (with prob = 0.5) 
    # and last month (with prob = 0.5)
    first.months <- 1:((min(starts)*12)-1)
    end.month <- max(sampled$year)*12
    last.months <- ((max(starts)*12)+1):end.month
    sampled <- sampled %>% group_by(unitID) %>% 
      mutate(trunc.begin=rbinom(1,1,.5),
             first.month=ifelse(trunc.begin,sample(first.months,1),1),
             trunc.end=rbinom(1,1,.5),
             last.month=ifelse(trunc.end,sample(last.months,1),end.month)) %>%
      ungroup() %>% filter((m>=first.month) & (m<=last.month)) %>%
      select(unitID,Y,m,q,year)
  }
  return(sampled)
}

add.trt.effects <- function(data, tau, grp.var, time.var){
  if (grp.var){
    data <- data %>% group_by(unitID) %>%
      mutate(grp=rbinom(1,1,.5)) %>% ungroup() %>%
      mutate(tau.g = ifelse(grp==1,tau-(.5*tau),tau+(.5*tau)))
  } else {
    data <- data %>% mutate(tau.g = tau)
  }
  if (time.var){
    data <- data %>% mutate(truth = posttrt*(tau.g-tau.g*exp(-.05*(m-start.month))))
  } else {
    data <- data %>% mutate(truth = posttrt*tau.g)
  }
  data <- data %>%
    rename(Y.unt = Y) %>%
    mutate(Y = Y.unt + truth)
}

joint.tests <- function(fitted.model){
  results <- tibble(name=c("Joint F","Trend"),
                    test.stat=c(NA,NA),
                    pval=c(NA,NA))
  # Get the clustered var-cov matrix
  # with a little cheat to make sure it's invert-able
  vcovCL <- Matrix::nearPD(fitted.model$vcov)$mat 
  param.names <- rownames(vcovCL)
  # Just the treatment effects
  which.params <- grepl("trt",param.names)
  
  # If there are multiple coefficients, do a joint F test and trend test:
  if ((sum(which.params)>1)&(sum(which.params)<100)){
    joint.F <- car::linearHypothesis(fitted.model,param.names[which.params],vcov.=vcovCL,test="F",singular.ok = TRUE)
    # Construct the hypothesis matrix for linear trend test:
    contr <- rep(0,length(param.names))
    contr[which.params] <- contr.poly(sum(which.params))[,".L"]
    trend <- car::linearHypothesis(fitted.model,contr,vcov.=vcovCL,test="F",singular.ok = TRUE)
    results$test.stat <- c(joint.F[2,'F'],trend[2,'F'])
    results$pval <- c(joint.F[2,'Pr(>F)'],trend[2,'Pr(>F)'])
  }
  return(results)
}

safe.joint.tests <- safely(joint.tests,
                           otherwise=tibble(name=c("Joint F","Trend"),
                                            test.stat=c(NA,NA),
                                            pval=c(NA,NA)))

# Post-processing the fitted model object
filter.fun <- function(fitted.model,group) {
  data <- data.frame(summary(fitted.model)$coefficients)
  nice_results <- tibble(name=rownames(data),data) %>% 
    filter(grepl("trt", name)) %>% # just the treatment effects
    mutate(time.counter = str_extract(name, "\\d+"),# extract time counters
           time=ifelse(is.na(time.counter),"post",time.counter),
           group=group) %>%
    rename(est=Estimate,se=Std..Error,test.stat=t.value,pval=Pr...t..,lb=CI.Lower,ub=CI.Upper) %>%
    select(-time.counter,-DF)
  return(nice_results)
}

agg.data.common <- function(data,is.tv,is.gv){
  monthly_data <- data %>% 
    mutate(unitID=factor(unitID),
           m=factor(m),
           q=factor(q),
           year=factor(year),
           post=as.numeric(post),
           postmonth = factor(ifelse(post>0,m,0)),
           postmonthtrt = factor(ifelse(posttrt==1,m,0)))
  
   if (is.gv){
     quarterly_data = monthly_data %>% group_by(unitID,grp,q,year)
     yearly_data = monthly_data %>% group_by(unitID,grp,year)
   } else {
     quarterly_data = monthly_data %>% group_by(unitID,q,year)
     yearly_data = monthly_data %>% group_by(unitID,year)
  }
 
  ## Aggregate to quarter
  quarterly_data = quarterly_data %>%
    summarise(across(c(Y,treatment,post,posttrt,truth),mean),.groups="keep") %>% ungroup() %>%
    mutate(postquarter = factor(ifelse(post>0,q,0)), 
           postquartertrt = factor(ifelse(posttrt==1,q,0)))
  
  # Aggregate to year
  yearly_data = yearly_data %>%
    summarise(across(c(Y,treatment,post,posttrt,truth),mean),.groups="keep") %>% ungroup() %>%
    mutate(postyear = factor(ifelse(post>0,year,0)),
           postyeartrt = factor(ifelse(posttrt==1,year,0)))
  
  if (is.tv){
    if (is.gv){
      month_true <- monthly_data %>% filter(treatment==1) %>% group_by(m,grp) 
      quarter_true <- quarterly_data %>% filter(treatment==1) %>% group_by(q,grp) 
      year_true <- yearly_data %>% filter(treatment==1) %>% group_by(year,grp) 
    } else {
      month_true <- monthly_data %>% filter(treatment==1) %>% group_by(m) 
      quarter_true <- quarterly_data %>% filter(treatment==1) %>% group_by(q)
      year_true <- yearly_data %>% filter(treatment==1) %>% group_by(year) 
    }
    month_true <- month_true %>% summarize(truth=mean(truth),.groups="drop") %>% rename(time=m) 
    quarter_true <- quarter_true  %>% summarize(truth=mean(truth),.groups="drop") %>% rename(time=q)
    year_true <- year_true %>% summarize(truth=mean(truth),.groups="drop") %>% rename(time=year)
  } else {
    if (is.gv){
      month_true <- monthly_data %>% filter(treatment==1) %>% group_by(post,grp) %>%
        summarize(truth=mean(truth),.groups="drop") %>% mutate(time=ifelse(post==1,"post","pre")) %>% select(-post)
      quarter_true <- quarterly_data %>% filter(treatment==1) %>% group_by(post,grp) %>%
        summarize(truth=mean(truth),.groups="drop") %>% mutate(time=ifelse(post==1,"post","pre"))  %>% select(-post)
      year_true <- yearly_data %>% filter(treatment==1) %>% group_by(post,grp) %>%
        summarize(truth=mean(truth),.groups="drop") %>% mutate(time=ifelse(post==1,"post","pre"))  %>% select(-post)
    } else {
      month_true <- monthly_data %>% filter(treatment==1) %>% group_by(post) %>%
        summarize(truth=mean(truth),.groups="drop") %>% mutate(time=ifelse(post==1,"post","pre")) %>% select(-post)
      quarter_true <- quarterly_data %>% filter(treatment==1) %>% group_by(post) %>%
        summarize(truth=mean(truth),.groups="drop") %>% mutate(time=ifelse(post==1,"post","pre"))  %>% select(-post)
      year_true <- yearly_data %>% filter(treatment==1) %>% group_by(post) %>%
        summarize(truth=mean(truth),.groups="drop") %>% mutate(time=ifelse(post==1,"post","pre")) %>% select(-post)
    }
  }
  return(list('monthly'=monthly_data,'quarterly'=quarterly_data,'yearly'=yearly_data,
         'month_true'=month_true %>% mutate(agg="month"),
         'quarter_true'=quarter_true %>% mutate(agg="quarter"),
         'year_true'=year_true %>% mutate(agg="year")))
}

agg.data.staggered <- function(data,is.gv=F){
  monthly_data <- data %>% ungroup() %>%
    mutate(unitID=as.numeric(unitID)) 
  
  ## Aggregate to quarter
  quarterly_data = monthly_data %>% group_by(unitID,q,start.quarter) %>%
    summarise(Y=mean(Y),.groups="drop")
  
  # Aggregate to year
  yearly_data = monthly_data %>% group_by(unitID,year,start.year) %>%
    summarise(Y=mean(Y),.groups="drop") 
  
  # Truth by treatment timing group at each post-treatment time point
  # For use in plotting functions, can also group by heterogeneous treatment groups
  if (is.gv){
    month_true <- data %>% group_by(start.month,grp,m)
    quarter_true <- data %>% group_by(start.quarter,grp,q) 
    year_true <- data %>% group_by(start.year,grp,year)
  } else {
    month_true <- data %>% group_by(start.month,m)
    quarter_true <- data %>% group_by(start.quarter,q)
    year_true <- data %>% group_by(start.year,year)
  }
  month_true <- month_true %>% summarize(truth=mean(truth),count=length(unique(unitID)),.groups="drop") %>% rename(time=m) %>% mutate(agg="month")
  quarter_true <- quarter_true %>% summarize(truth=mean(truth),count=length(unique(unitID)),.groups="drop") %>% rename(time=q) %>% mutate(agg="quarter")
  year_true <-  year_true %>%  summarize(truth=mean(truth),count=length(unique(unitID)),.groups="drop") %>% rename(time=year) %>% mutate(agg="year")
    
  # Using relative time, compute the true treatment effect by weighting the truth by the proportion of units contributing at each time
  dynamic_truth <-  bind_rows(month_true %>% filter(start.month!=0) %>% mutate(time=time-start.month) %>% group_by(time) %>%
                                mutate(total=sum(count)) %>% ungroup() %>% mutate(wt=count/total) %>% group_by(time,agg) %>%
                                summarize(truth=weighted.mean(truth,w=wt),.groups="drop"),
                              quarter_true %>% filter(start.quarter!=0) %>% mutate(time=time-start.quarter) %>% group_by(time) %>%
                                mutate(total=sum(count)) %>% ungroup() %>% mutate(wt=count/total) %>% group_by(time,agg) %>%
                                summarize(truth=weighted.mean(truth,w=wt),.groups="drop"),
                              year_true %>% filter(start.year!=0) %>% mutate(time=time-start.year) %>% group_by(time) %>% 
                                mutate(total=sum(count)) %>% ungroup() %>% mutate(wt=count/total) %>% group_by(time,agg) %>%
                                summarize(truth=weighted.mean(truth,w=wt),.groups="drop"))
  
  static_truth <- bind_rows(month_true %>% filter(time>=start.month,start.month!=0) %>% 
                              group_by(start.month,agg) %>% summarize(truth=mean(truth),.groups="drop") %>% rename(time=start.month),
                            quarter_true %>% filter(time>=start.quarter,start.quarter!=0) %>% 
                              group_by(start.quarter,agg)  %>% summarize(truth=mean(truth),.groups="drop") %>% rename(time=start.quarter),
                            year_true %>% filter(time>=start.year,start.year!=0) %>% 
                              group_by(start.year,agg) %>% summarize(truth=mean(truth),.groups="drop") %>% rename(time=start.year))
  
  return(list('monthly'=monthly_data,'quarterly'=quarterly_data,'yearly'=yearly_data,
              'month_true'=month_true,'quarter_true'=quarter_true,'year_true'=year_true,
              'dynamic_truth'=dynamic_truth,'static_truth'=static_truth))
}

analyze.data.common <- function(data,is.tv){
  agg.dat <- agg.data.common(data,is.tv,is.gv=F)
  
  if (is.tv){
    monthly_DID_2x2= estimatr::lm_robust(Y ~ treatment + postmonthtrt, fixed_effects = ~postmonth, data = agg.dat$monthly, clusters=unitID, se_type="stata")
    monthly_DID_0  = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth, data = agg.dat$monthly, clusters=unitID,se_type="stata")
    monthly_DID_uy = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + year + postmonth,data = agg.dat$monthly, clusters=unitID,se_type="stata")
    
    quarterly_DID_2x2= estimatr::lm_robust(Y ~ treatment + postquartertrt, fixed_effects = ~postquarter, data = agg.dat$quarterly, clusters=unitID,se_type="stata")
    quarterly_DID_0  = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter, data = agg.dat$quarterly, clusters=unitID,se_type="stata")
    quarterly_DID_uy = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + year + postquarter, data = agg.dat$quarterly, clusters=unitID,se_type="stata")
    
    yearly_DID_2x2 = estimatr::lm_robust(Y ~ treatment + postyeartrt, fixed_effects = ~postyear, data = agg.dat$yearly, clusters=unitID,se_type="stata")
    yearly_DID_0   = estimatr::lm_robust(Y ~ postyeartrt, fixed_effects = ~unitID + postyear, data = agg.dat$yearly, clusters=unitID,se_type="stata")
    yearly_DID_uy  = estimatr::lm_robust(Y ~ postyeartrt, fixed_effects = ~unitID + year + postyear, data = agg.dat$yearly, clusters=unitID,se_type="stata")
    
  } else {
    monthly_DID_2x2 =estimatr::lm_robust(Y ~ post + treatment + posttrt, data = agg.dat$monthly, clusters=unitID, se_type="stata")
    monthly_DID_0 =  estimatr::lm_robust(Y ~ post + posttrt, fixed_effects = ~unitID, data = agg.dat$monthly, clusters=unitID,se_type="stata")
    monthly_DID_uy = estimatr::lm_robust(Y ~ posttrt, fixed_effects = ~unitID + year,data = agg.dat$monthly, clusters=unitID,se_type="stata")
    
    quarterly_DID_2x2 =estimatr::lm_robust(Y ~ post + treatment + posttrt, data = agg.dat$quarterly, clusters=unitID,se_type="stata")
    quarterly_DID_0  = estimatr::lm_robust(Y ~ post + posttrt, fixed_effects = ~unitID, data = agg.dat$quarterly, clusters=unitID,se_type="stata")
    quarterly_DID_uy = estimatr::lm_robust(Y ~ posttrt, fixed_effects = ~unitID + year, data = agg.dat$quarterly, clusters=unitID,se_type="stata")
    
    yearly_DID_2x2 =estimatr::lm_robust(Y ~ post + treatment + posttrt, data = agg.dat$yearly, clusters=unitID, se_type="stata")
    yearly_DID_0 =  estimatr::lm_robust(Y ~ post + posttrt, fixed_effects = ~unitID, data = agg.dat$yearly, clusters=unitID,se_type="stata")
    yearly_DID_uy = estimatr::lm_robust(Y ~ posttrt, fixed_effects = ~unitID + year, data = agg.dat$yearly, clusters=unitID,se_type="stata")
  }
  # Joint model estimates with truth
  monthly_results <- left_join(rbind(filter.fun(monthly_DID_2x2, "monthly_2x2"),
                                     filter.fun(monthly_DID_0,"monthly_0"),
                                     filter.fun(monthly_DID_uy, "monthly_uy")),agg.dat$month_true,by="time")
  quarterly_results <- left_join(rbind(filter.fun(quarterly_DID_2x2, "quarterly_2x2"),
                                       filter.fun(quarterly_DID_0, "quarterly_0"),
                                       filter.fun(quarterly_DID_uy, "quarterly_uy")),agg.dat$quarter_true,by="time")
  yearly_results <- left_join(rbind(filter.fun(yearly_DID_2x2, "yearly_2x2"),
                                    filter.fun(yearly_DID_0, "yearly_0"),
                                    filter.fun(yearly_DID_uy, "yearly_uy")),agg.dat$year_true,by="time")
  
  results <- bind_rows(monthly_results,quarterly_results,yearly_results)
  
  # Do the joint tests
  if (is.tv){
    test.results <- bind_rows(safe.joint.tests(monthly_DID_2x2)$result %>% mutate(group= "monthly_2x2",time="post"),
                              safe.joint.tests(monthly_DID_0)$result %>% mutate(group= "monthly_0",time="post"),
                              safe.joint.tests(monthly_DID_uy)$ result %>% mutate(group= "monthly_uy",time="post"),
                              safe.joint.tests(quarterly_DID_2x2)$ result %>% mutate(group= "quarterly_2x2",time="post"),
                              safe.joint.tests(quarterly_DID_0)$ result %>% mutate(group= "quarterly_0",time="post"),
                              safe.joint.tests(quarterly_DID_uy)$result %>% mutate(group= "quarterly_uy",time="post"),
                              safe.joint.tests(yearly_DID_2x2)$result %>% mutate(group= "yearly_2x2",time="post"),
                              safe.joint.tests(yearly_DID_0)$result %>% mutate(group= "yearly_0",time="post"),
                              safe.joint.tests(yearly_DID_uy)$result %>% mutate(group= "yearly_uy",time="post"))
    results <- bind_rows(results,test.results)
  }
  return(results)
}

analyze.data.CS <- function(data){
  agg.dat <- agg.data.staggered(data)
  monthly_DID <- did::att_gt(yname = "Y",
                             tname = "m",
                             idname = "unitID",
                             gname = "start.month",
                             xformla = ~1,
                             data = agg.dat$monthly,
                             allow_unbalanced_panel = TRUE)
  
  quarterly_DID <- did::att_gt(yname = "Y",
                               tname = "q",
                               idname = "unitID",
                               gname = "start.quarter",
                               xformla = ~1,
                               data = agg.dat$quarterly,
                               allow_unbalanced_panel = TRUE)
  
  yearly_DID <-  did::att_gt(yname = "Y",
                             tname = "year",
                             idname = "unitID",
                             gname = "start.year",
                             xformla = ~1,
                             data = agg.dat$yearly,
                             allow_unbalanced_panel = TRUE)
  
  month.agg.es <- did::aggte(monthly_DID, type = "dynamic", na.rm=TRUE)
  month.agg.gs <- did::aggte(monthly_DID, type = "group", na.rm=TRUE)
  quarter.agg.es <- did::aggte(quarterly_DID, type = "dynamic", na.rm=TRUE)
  quarter.agg.gs <- did::aggte(quarterly_DID, type = "group", na.rm=TRUE) 
  year.agg.es <- did::aggte(yearly_DID, type = "dynamic", na.rm=TRUE)
  year.agg.gs <- did::aggte(yearly_DID, type = "group", na.rm=TRUE) 
  
  dynamic <- bind_rows(tibble('est'= month.agg.es$att.egt, 'se'= month.agg.es$se.egt, 'time'=month.agg.es$egt, 'agg'="month"),
                       tibble('est'= quarter.agg.es$att.egt, 'se'= quarter.agg.es$se.egt, 'time'=quarter.agg.es$egt, 'agg'="quarter"),
                       tibble('est'= year.agg.es$att.egt, 'se'= year.agg.es$se.egt, 'time'=year.agg.es$egt, 'agg'="year"))
  static <- bind_rows(tibble('est'= month.agg.gs$att.egt, 'se'= month.agg.gs$se.egt, 'time'=month.agg.gs$egt, 'agg'="month"),
                      tibble('est'= quarter.agg.gs$att.egt, 'se'= quarter.agg.gs$se.egt, 'time'=quarter.agg.gs$egt, 'agg'="quarter"),
                      tibble('est'= year.agg.gs$att.egt, 'se'= year.agg.gs$se.egt, 'time'=year.agg.gs$egt, 'agg'="year"))
  
  dynamic_results <- left_join(dynamic,agg.dat$dynamic_truth,by=c('time','agg')) %>% mutate(estimand="Time-varying")
  static_results <- left_join(static,agg.dat$static_truth,by=c("time","agg")) %>% mutate(estimand="Overall")
  
  results <- bind_rows(dynamic_results,static_results)
  return(results)
}

analyze.and.combine <- function(data){
  overall <- analyze.data.common(data,is.tv=F)
  timevar <- analyze.data.common(data,is.tv=T)
  return(rbind(overall,timevar))
}

inject.analyze <- function(data,params,staggered){
  null.data <- add.trt.effects(data,tau=0,         grp.var = F, time.var = F)
  cons.data <- add.trt.effects(data,tau=params$tau,grp.var = F, time.var = F)
  tv.data   <- add.trt.effects(data,tau=params$tau,grp.var = F, time.var = T)
  gv.data   <- add.trt.effects(data,tau=params$tau,grp.var = T, time.var = F)
  gtv.data  <- add.trt.effects(data,tau=params$tau,grp.var = T, time.var = T)
  
  if (staggered) {
    null.res <- analyze.data.CS(null.data)
    cons.res <- analyze.data.CS(cons.data)
    tv.res   <- analyze.data.CS(tv.data)
    gv.res   <- analyze.data.CS(gv.data)
    gtv.res  <- analyze.data.CS(gtv.data)
  } else {
    null.res <- analyze.and.combine(null.data)
    cons.res <- analyze.and.combine(cons.data)
    tv.res   <- analyze.and.combine(tv.data)
    gv.res   <- analyze.and.combine(gv.data)
    gtv.res  <- analyze.and.combine(gtv.data)
  }
  results <- bind_rows(null.res %>% mutate(trteff='null'),
                       cons.res %>% mutate(trteff='constant'),
                       tv.res %>% mutate(trteff='time-varying'),
                       gv.res %>% mutate(trteff='group-varying'),
                       gtv.res %>% mutate(trteff='group- and time-varying')) %>%
    mutate(adoption=ifelse(staggered,'staggered','common'))
  return(results)
}

# returns a unit-level data frame with binary indicators for treatment
# and start group
assign.treatment <- function(data,starts,staggered){
  # Randomly assign treatment to each unit
  trt.dat <- data %>% group_by(unitID) %>% slice(1) %>% mutate(treatment = rbinom(1,1,0.5))
  if (staggered){
    trt.dat <- trt.dat %>%
      # Choose a random start year;
      # Untreated units have start times of 0
      mutate(start.year=treatment*sample(starts,1),
             start.month=treatment*(12*(start.year-1)+1),
             start.quarter=treatment*(4*(start.year-1)+1)) 
  } else {
    # Choose a single random start year
    single.start <- sample(starts,1)
    # apply it to all units
    trt.dat <- trt.dat %>%
      mutate(start.year=single.start,
             start.month=12*(start.year-1)+1,
             start.quarter=4*(start.year-1)+1)
  }
  trt.dat <- trt.dat %>%
    select(unitID,treatment,start.year,start.month,start.quarter)
  return(trt.dat)
}

make.data <- function(long.dat,trt.dat){
  left_join(long.dat,trt.dat,by='unitID') %>% 
    mutate(post = (m >= start.month) & (start.month!=0),
           posttrt=post*treatment)
}


single.iter <- function(params,staggered,resampling=F,starts,data){
  if (resampling){
    bal.data <- resample(data,n.units=myparams$n.units,starts,unbalanced=F) 
    unbal.data <- resample(data,n.units=myparams$n.units,starts,unbalanced=T) 
  } else {
    bal.data   <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=F)
    unbal.data <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=T)
    starts <- 3:(params$n.year-2)
  }

  bal.trt.dat <- assign.treatment(bal.data,starts=starts,staggered=staggered)
  bal.rand.results <- inject.analyze(make.data(bal.data,bal.trt.dat),params,staggered) %>%
    mutate(panel='balanced',assignment='random')  

  unbal.trt.dat <- assign.treatment(unbal.data,starts=starts,staggered=staggered)
  unbal.rand.results <- inject.analyze(make.data(unbal.data,unbal.trt.dat),params,staggered) %>%
    mutate(panel='unbalanced',assignment='random')
  
  return(rbind(bal.rand.results,unbal.rand.results))
}
