#' Generate monthly shocks 
#' 
#' Generate monthly shocks from iid normal
#' 
#' @description
#' `draw.months` generates monthly shocks from an ARIMA process
#'
#' @param params list with named elements: total_num_month,sigma_month,const,phi,total_num_quarter,n.year
#' @returns tibble with the following variables: Y.month, monthofyear, monthofquarter, m, q
#' @details The parameter `phi` is a 2-vector containing the 2 autoregressive parameters; 
#' `phi[1]` controls likelihood of wandering away while `phi[2]` controls "acceleration" (how fast)
#' The parameter `sigma_month` is a white noise error SD
#' The parameter `const` is an overall non-stationary trend

draw.months <- function(params){
  # Create empty object of the correct dimension
  y <- Y <- rep(NA,params$total_num_month)
  # Start the timeseries
  y[1] <- params$const + rnorm(1,0,params$sigma_month)
  Y[1] <- y[1]
  y[2] <- params$const + rnorm(1,0,params$sigma_month)
  Y[2] <- y[2] + 2*Y[1]
  # Iterate over the month index to generate the rest of the autoregressive series
  for (m in 3:params$total_num_month){
    y[m] <- params$const + params$phi[1]*y[m-1] + 
      params$phi[2]*y[m-2] + rnorm(1,0,params$sigma_month)
    Y[m] <- y[m] + Y[m-1] 
  }
  months.per.quarter <- params$total_num_month/params$total_num_quarter
  months <- tibble(Y.month=Y,
                   # Seasonal month indices (within year)
                   monthofyear = rep(1:12,params$n.year),
                   # Seasonal month indices (within quarter)
                   monthofquarter = rep(1:months.per.quarter,params$total_num_quarter), 
                   # Overall month index
                   m=1:params$total_num_month,
                   # Overall quarter index
                   # (generated here so we can match to quarter effects)
                   q=rep(1:params$total_num_quarter,each=months.per.quarter))
  return(months)
}


#' Generate quarterly shocks from iid normal
#' 
#' @description
#' `draw.quarters` generates iid normal quarter shocks but sorted so that they are always increasing
#'
#' @param params list with named elements: n.quarter,sigma_quarter,n.year,
#' @returns tibble with the following variables: Y.quarter,quarterofyear,q,year

draw.quarters <- function(params){
  ## Draw a random vector of quarter effects
  temp <- rnorm(4,0,params$sigma_quarter)
  # Sort so that the seasonal effects is always increasing over the year
  quarter.effects <- temp[order(temp)]
  quarters <- tibble(Y.quarter=rep(quarter.effects,params$n.year),
                     # Seasonal quarter indices (within year)
                     quarterofyear = rep(1:4,params$n.year),
                     # Overall quarter indices
                     q=rep(1:params$total_num_quarter),
                     # Year indices
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
#' @returns tibble with the following variables: Y, m, q, year, monthofyear, quarterofyear

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
  
  return(data %>% select(unitID, Y, m, q, year, monthofyear, quarterofyear))
}

#### Resample data ####
resample <- function(data,n.units){
  # Randomly select (with replacement) individual units
  units <- unique(data$unitID)
  these.units <- sample(units,n.units,replace=T)
  # Randomly select (with replacement) start months from the unit-level distribution of start months
  starts <- (data %>% group_by(unitID) %>% slice(1) %>% ungroup())$start.month
  these.starts <- sample(starts,n.units,replace=T)
  
  sampled <- tibble('unitID'=these.units,
                    'start.month'=these.starts,
                    'newid'=1:n.units) %>% 
    left_join(data %>% select(-c(start.month,start.quarter,start.year)),
              by="unitID",relationship="many-to-many") %>%
    # Recreate the start.quarter and start.year variables:
    mutate(treatment=ifelse(start.month<4*12,1,0), # never-treated are units with start months after year 4
           start.month=treatment*start.month,
           start.quarter = treatment*ceiling(start.month/3),
           start.year = treatment*ceiling(start.month/12),
           unitID=newid,
           post = m >= start.month & start.month!=0,
           posttrt=post*treatment,
           start.grp=factor(start.month)) %>%
    select(unitID,Y,m,q,year,start.month,start.quarter,start.year,treatment,post,posttrt,start.grp)
}

#### Inject treatment effects ####
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

#### Analyze data functions ####
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
    contr[which.params] <- contr.poly(sum(which.param))[,".L"]
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
filter_fun <- function(fitted.model,group) {
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

analyze.data.common <- function(data,is.tv){
  monthly_data <- data %>% 
    mutate(unitID=factor(unitID),
           monthofyear=factor(monthofyear),
           quarterofyear=factor(quarterofyear),
           m=factor(m),
           q=factor(q),
           year=factor(year),
           postmonth = factor(ifelse(post>0,m,0)),
           postmonthtrt = factor(ifelse(posttrt==1,m,0)))
  
  ## Aggregate to quarter
  quarterly_data = monthly_data %>%
    group_by(unitID,q,quarterofyear,year) %>%
    summarise(across(c(Y,treatment,post,posttrt,truth),mean),.groups="keep") %>% ungroup() %>%
    mutate(postquarter = factor(ifelse(post>0,q,0)), 
           postquartertrt = factor(ifelse(posttrt==1,q,0)))
  
  # Aggregate to year
  yearly_data = monthly_data %>%
    group_by(unitID,year) %>%
    summarise(across(c(Y,treatment,post,posttrt,truth),mean),.groups="keep") %>% ungroup() %>%
    mutate(postyear = factor(ifelse(post>0,year,0)),
           postyeartrt = factor(ifelse(posttrt==1,year,0)))
  
  if (is.tv){
    monthly_DID_2x2= estimatr::lm_robust(Y ~ post + treatment + postmonthtrt, data = monthly_data, clusters=unitID, se_type="stata")
    monthly_DID_0  = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth, data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_sm = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + monthofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_sq = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + quarterofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_um = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + m,data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_uq = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + q,data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_uy = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + year,data = monthly_data, clusters=unitID,se_type="stata")
    
    quarterly_DID_2x2= estimatr::lm_robust(Y ~ post + treatment + postquartertrt, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_0  = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter, data = quarterly_data, clusters=unitID,se_type="stata")
    #quarterly_DID_sq = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter + quarterofyear + year, data = quarterly_data, clusters=unitID,se_type="stata")
    #quarterly_DID_uq = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter + q, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_uy = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter + year, data = quarterly_data, clusters=unitID,se_type="stata")
    
    yearly_DID_2x2 = estimatr::lm_robust(Y ~ post + treatment + postyeartrt, data = yearly_data, clusters=unitID,se_type="stata")
    yearly_DID_0   = estimatr::lm_robust(Y ~ postyeartrt, fixed_effects = ~unitID + postyear, data = yearly_data, clusters=unitID,se_type="stata")
    yearly_DID_uy  = estimatr::lm_robust(Y ~ postyeartrt, fixed_effects = ~unitID + postyear + year, data = yearly_data, clusters=unitID,se_type="stata")
    
    # Average the truth over counties at each post-treatment time point
    # for later merging on 
    month_true <- monthly_data %>% filter(posttrt==1) %>% group_by(postmonth) %>%
      summarize(truth=mean(truth)) %>% rename(time=postmonth)
    quarter_true <- quarterly_data %>% filter(posttrt==1) %>% group_by(postquarter) %>%
      summarize(truth=mean(truth)) %>% rename(time=postquarter)
    year_true <- yearly_data %>% filter(posttrt==1) %>% group_by(postyear) %>%
      summarize(truth=mean(truth)) %>% rename(time=postyear)
  } else {
    monthly_DID_2x2 =estimatr::lm_robust(Y~post+treatment+posttrt, data = monthly_data, clusters=unitID, se_type="stata")
    monthly_DID_0 =  estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID, data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_sm = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + monthofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_sq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + quarterofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_um = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + m,data = monthly_data, clusters=unitID,se_type="stata")
    #monthly_DID_uq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + q,data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_uy = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + year,data = monthly_data, clusters=unitID,se_type="stata")
    
    quarterly_DID_2x2 =estimatr::lm_robust(Y~post+treatment+posttrt, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_0  = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID, data = quarterly_data, clusters=unitID,se_type="stata")
    #quarterly_DID_sq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + quarterofyear + year, data = quarterly_data, clusters=unitID,se_type="stata")
    #quarterly_DID_uq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + q, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_uy = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + year, data = quarterly_data, clusters=unitID,se_type="stata")
    
    yearly_DID_2x2 =estimatr::lm_robust(Y~post+treatment+posttrt, data = yearly_data, clusters=unitID,se_type="stata")
    yearly_DID_0 =  estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID, data = yearly_data, clusters=unitID,se_type="stata")
    yearly_DID_uy = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + year, data = yearly_data, clusters=unitID,se_type="stata")
    
    month_true <- monthly_data %>% filter(posttrt==1) %>% group_by(post) %>%
      summarize(truth=mean(truth)) %>% mutate(time="post") %>% select(-post)
    quarter_true <- quarterly_data %>% filter(posttrt==1) %>% group_by(post) %>%
      summarize(truth=mean(truth)) %>% mutate(time="post") %>% select(-post)
    year_true <- yearly_data %>% filter(posttrt==1) %>% group_by(post) %>%
      summarize(truth=mean(truth)) %>% mutate(time="post") %>% select(-post)
  }
  monthly_results <- left_join(rbind(filter_fun(monthly_DID_2x2, "monthly_2x2"),
                                     filter_fun(monthly_DID_0,"monthly_0"),
                                     #filter_fun(monthly_DID_sm, "monthly_sm"),
                                     #filter_fun(monthly_DID_sq, "monthly_sq"),
                                     #filter_fun(monthly_DID_um, "monthly_um"),
                                     #filter_fun(monthly_DID_uq, "monthly_uq"),
                                     filter_fun(monthly_DID_uy, "monthly_uy")),month_true,by="time")
  quarterly_results <- left_join(rbind(filter_fun(quarterly_DID_2x2, "quarterly_2x2"),
                                       filter_fun(quarterly_DID_0, "quarterly_0"),
                                       #filter_fun(quarterly_DID_sq, "quarterly_sq"),
                                       #filter_fun(quarterly_DID_uq, "quarterly_uq"),
                                       filter_fun(quarterly_DID_uy, "quarterly_uy")),quarter_true,by="time")
  yearly_results <- left_join(rbind(filter_fun(yearly_DID_2x2, "yearly_2x2"),
                                    filter_fun(yearly_DID_0, "yearly_0"),
                                    filter_fun(yearly_DID_uy, "yearly_uy")),year_true,by="time")
  
  results <- bind_rows(monthly_results,quarterly_results,yearly_results)
  
  # Do the joint tests
  if (is.tv){
    test.results <- bind_rows(safe.joint.tests(monthly_DID_2x2)$result %>% mutate(group= "monthly_2x2",time="post"),
                              safe.joint.tests(monthly_DID_0)$result %>% mutate(group= "monthly_0",time="post"),
                              #safe.joint.tests(monthly_DID_sm)$result %>% mutate(group= "monthly_sm",time="post"),
                              #safe.joint.tests(monthly_DID_sq)$result %>% mutate(group= "monthly_sq",time="post"),
                              #safe.joint.tests(monthly_DID_um)$result %>% mutate(group= "monthly_um",time="post"),
                              #safe.joint.tests(monthly_DID_uq)$ result %>% mutate(group= "monthly_uq",time="post"),
                              safe.joint.tests(monthly_DID_uy)$ result %>% mutate(group= "monthly_uy",time="post"),
                              safe.joint.tests(quarterly_DID_2x2)$ result %>% mutate(group= "quarterly_2x2",time="post"),
                              safe.joint.tests(quarterly_DID_0)$ result %>% mutate(group= "quarterly_0",time="post"),
                              #safe.joint.tests(quarterly_DID_sq)$result %>% mutate(group= "quarterly_sq",time="post"),
                              #safe.joint.tests(quarterly_DID_uq)$result %>% mutate(group= "quarterly_uq",time="post"),
                              safe.joint.tests(quarterly_DID_uy)$result %>% mutate(group= "quarterly_uy",time="post"),
                              safe.joint.tests(yearly_DID_2x2)$result %>% mutate(group= "yearly_2x2",time="post"),
                              safe.joint.tests(yearly_DID_0)$result %>% mutate(group= "yearly_0",time="post"),
                              safe.joint.tests(yearly_DID_uy)$result %>% mutate(group= "yearly_uy",time="post"))
    results <- bind_rows(results,test.results)
  }
  return(results)
}

analyze.data.CS <- function(data){
  monthly_data <- data %>% ungroup() %>%
    mutate(unitID=as.numeric(unitID)) # Requires a numeric unit ID variable
  
  ## Aggregate to quarter
  quarterly_data = monthly_data %>%
    group_by(unitID,q,start.quarter) %>%
    summarise(Y=mean(Y),.groups="keep") %>% ungroup() 
  
  # Aggregate to year
  yearly_data = monthly_data %>%
    group_by(unitID,year,start.year) %>%
    summarise(Y=mean(Y),.groups="keep") %>% ungroup() 
  
  monthly_DID <- did::att_gt(yname = "Y",
                             tname = "m",
                             idname = "unitID",
                             gname = "start.month",
                             xformla = ~1,
                             data = monthly_data,
                             control_group="notyettreated",
                             allow_unbalanced_panel = TRUE)
  
  quarterly_DID <- did::att_gt(yname = "Y",
                               tname = "q",
                               idname = "unitID",
                               gname = "start.quarter",
                               xformla = ~1,
                               data = quarterly_data,
                               control_group="notyettreated",
                               allow_unbalanced_panel = TRUE)
  
  yearly_DID <-  did::att_gt(yname = "Y",
                             tname = "year",
                             idname = "unitID",
                             gname = "start.year",
                             xformla = ~1,
                             data = yearly_data,
                             control_group="notyettreated",
                             allow_unbalanced_panel = TRUE)
  
  # Truth by treatment timing group at each post-treatment time point
  # These are the group-time estimates
  month_true <- data %>% group_by(start.month,m) %>%
    summarize(truth=mean(truth),count=length(unique(unitID)),.groups="drop") %>% rename(time=m)
  quarter_true <- data %>% group_by(start.quarter,q) %>%
    summarize(truth=mean(truth),count=length(unique(unitID)),.groups="drop") %>% rename(time=q)
  year_true <- data %>% group_by(start.year,year) %>%
    summarize(truth=mean(truth),count=length(unique(unitID)),.groups="drop") %>% rename(time=year)

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
  dynamic_truth <-  bind_rows(month_true %>% filter(start.month!=0) %>% mutate(time=time-start.month) %>% group_by(time) %>%
                                mutate(total=sum(count)) %>% ungroup() %>% mutate(wt=count/total) %>% group_by(time) %>%
                                summarize(truth=weighted.mean(truth,w=wt)) %>% mutate(agg="month"),
                              quarter_true %>% filter(start.quarter!=0) %>% mutate(time=time-start.quarter) %>% group_by(time) %>%
                                mutate(total=sum(count)) %>% ungroup() %>% mutate(wt=count/total) %>% group_by(time) %>%
                                summarize(truth=weighted.mean(truth,w=wt)) %>% mutate(agg="quarter"),
                              year_true %>% filter(start.year!=0) %>% mutate(time=time-start.year) %>% group_by(time) %>% 
                                mutate(total=sum(count)) %>% ungroup() %>% mutate(wt=count/total) %>% group_by(time) %>%
                                summarize(truth=weighted.mean(truth,w=wt)) %>% mutate(agg="year"))
  static_truth <- bind_rows(month_true %>% group_by(start.month) %>% filter(time>=start.month) %>% summarize(truth=mean(truth)) %>% 
                              mutate(agg="month") %>% rename(time=start.month),
                            quarter_true %>% group_by(start.quarter) %>% filter(time>=start.quarter) %>% summarize(truth=mean(truth)) %>%
                              mutate(agg="quarter") %>% rename(time=start.quarter),
                            year_true %>% group_by(start.year) %>% filter(time>=start.year) %>% summarize(truth=mean(truth)) %>%
                              mutate(agg="year") %>% rename(time=start.year))
  dynamic_results <- left_join(dynamic,dynamic_truth,by=c('time','agg')) %>% mutate(estimand="Time-varying")
  static_results <- left_join(static,static_truth,by=c("time","agg")) %>% mutate(estimand="Overall")
  
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
assign.treatment <- function(data,params,by.trend,staggered){
  if (by.trend) {
    trt.dat <- data %>% group_by(unitID) %>%
      # fit linear model to each unit's data
      group_modify(~broom::tidy(lm(Y~year,data=.x))) %>% filter(term=="year") %>%
      # center the slopes (to net out overall trend)
      mutate(mean.slope=mean(estimate),
             ctr.slope=estimate-mean.slope,
             # Units with more positive slopes have higher pr(treat)
             treatment = rbinom(1,1,plogis(ctr.slope)))
  } else {
    # Randomly assign treatment to each unit
    trt.dat <- data %>% group_by(unitID) %>% slice(1) %>%
      mutate(treatment = rbinom(1,1,0.5))
  }
  if (staggered){
    trt.dat <- trt.dat %>%
      # Choose a random year between 3 and n.year-2 for each unit
      # so there are at least 2 pre and 2 post for every timing group
      mutate(start.year=sample(3:(params$n.year-2),1))
  } else {
    trt.dat <- trt.dat %>%
      # Every units starts at the middle year (rounded up)
      mutate(start.year=ceiling(params$n.year/2))
  }
  trt.dat <- trt.dat %>% 
    mutate(start.month=treatment*(12*(start.year-1)+1),
           start.quarter=treatment*(4*(start.year-1)+1)) %>%
    select(unitID,treatment,start.year,start.month,start.quarter)
  return(trt.dat)
}

trt.assignments <- function(data,params,staggered){
  rand.trt <- assign.treatment(data,params,by.trend=F,staggered=staggered)
  trend.trt <-  assign.treatment(data,params,by.trend=T,staggered=staggered)
  
  return(list('rand.trt'=rand.trt,
              'trend.trt'=trend.trt))
}

make.data <- function(long.dat,trt.dat){
  left_join(long.dat,trt.dat,by='unitID') %>%
    mutate(post = (m >= start.month) & (start.month!=0),
           posttrt=post*treatment)
}


single_iteration <- function(params,staggered){
  bal.data   <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=F)
  unbal.data <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=T)
  
  bal.trt.dat <- trt.assignments(bal.data,params,staggered=staggered)
  unbal.trt.dat <- trt.assignments(unbal.data,params,staggered=staggered)
  
  bal.rand.results <- inject.analyze(make.data(bal.data,bal.trt.dat$rand.trt),params,staggered) %>%
    mutate(panel='balanced',assignment='random')
  
  #bal.trend.results <- inject.analyze(make.data(bal.data,bal.trt.dat$trend.trt),params,staggered) %>%
  #  mutate(panel='balanced',assignment='trended')
  
  unbal.rand.results <- inject.analyze(make.data(unbal.data,unbal.trt.dat$rand.trt),params,staggered) %>%
    mutate(panel='unbalanced',assignment='random')
  
  #unbal.trend.results <- inject.analyze(make.data(unbal.data,bal.trt.dat$trend.trt),params,staggered) %>%
  #  mutate(panel='unbalanced',assignment='trended')
  
  return(rbind(bal.rand.results,unbal.rand.results))
}
