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
                   monthofyear = rep(1:params$n.month,params$n.year),
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
  temp <- rnorm(params$n.quarter,0,params$sigma_quarter)
  # Sort so that the seasonal effects is always increasing over the year
  quarter.effects <- temp[order(temp)]
  quarters <- tibble(Y.quarter=rep(quarter.effects,params$n.year),
                     # Seasonal quarter indices (within year)
                     quarterofyear = rep(1:params$n.quarter,params$n.year),
                     # Overall quarter indices
                     q=rep(1:params$total_num_quarter),
                     # Year indices
                     year=rep(1:params$n.year,each=params$n.quarter))
  return(quarters)
}

#' Generate quarterly shocks from iid normal
#' 
#' @description
#' `generate.data` generates autocorrelated data for n.units over n.years
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
  params$total_num_month = params$n.month * params$n.year
  params$total_num_quarter = params$n.quarter * params$n.year

  # Draw month vectors
  if (month.byunit){
    monthdat <- replicate(params$n.units,draw.months(params),simplify=FALSE)
    monthdat <- bind_rows(monthdat, .id="unitID")
  } else {
    months <- draw.months(params)
    monthdat <- cross_join(months,tibble(unitID=as.character(1:params$n.units))) %>% 
      arrange(unitID,m)
  }
  data <- monthdat
  
  # Draw quarter vectors
  if (quarter.byunit){
    quarterdat <- replicate(params$n.units,draw.quarters(params),simplify=FALSE)
    quarterdat <- bind_rows(quarterdat,.id="unitID")
  } else {
    quarters <- draw.quarters(params)
    quarterdat <- cross_join(quarters,tibble(unitID=as.character(1:params$n.units))) %>% 
      arrange(unitID,q)
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
cons.trt.effects <- function(data, tau, eq.by.grp){
  # For treatment effects that are equal in all units
  if (eq.by.grp){
    data <- data %>%
      mutate(truth = posttrt*tau) %>%
      rename(Y.unt = Y) %>% 
      mutate(Y = Y.unt + truth)
    return(data)
  } else {
    # Randomly assign units to be in a "high" or "low" treatment response group
    data <- data %>% group_by(unitID) %>%
      mutate(grp = rbinom(1,1,.5)) %>% ungroup() %>%
      mutate(tau.g = ifelse(grp==1,tau-(.5*tau),tau+(.5*tau)),
             truth = posttrt*tau.g) %>%
      rename(Y.unt = Y) %>%
      mutate(Y = Y.unt + truth)
    return(data)
  }
}

cons.trt.effects.stag <- function(data, tau, eq.by.grp){
  # For treatment effects that are the same in all treatment timing groups
  if (eq.by.grp){
    data <- data %>%
      mutate(truth = posttrt*tau) %>%
      mutate(Y.unt = Y,
             Y = Y.unt + truth)
    return(data)
  } else {
    # Each treatment timing group has a random treatment response 
    # Centered around tau (plus or minus 1/2 tau SD on either side)
    data <- data %>%
      group_by(start.grp) %>%
      mutate(tau.g = rnorm(1,tau,.5*tau)) %>%
      ungroup() %>%
      mutate(truth = posttrt*tau.g,
             Y.unt = Y,
             Y = Y.unt + truth)
    return(data)
  }
}

time.trt.effects <- function(data, tau, eq.by.grp){
  # Find the start of the post-treatment period
  t0 <- (filter(data,posttrt==1) %>% mutate(start=first(m)) %>% slice(1))$m
  
  # For treatment effects that are equal in all units
  if (eq.by.grp){
    data <- data %>%
      mutate(truth=posttrt*(tau-tau*exp(-.05*(m-t0)))) %>%
      rename(Y.unt = Y) %>%
      mutate(Y = Y.unt + truth)
    return(data)
  } else {
    # Randomly assign units to be in a "high" or "low" treatment response group
    data <- data %>% 
      group_by(unitID) %>%
      mutate(grp=rbinom(1,1,.5)) %>% ungroup() %>%
      mutate(tau.g = ifelse(grp==1,tau-(.5*tau),tau+(.5*tau)),
             truth = posttrt*(tau.g-tau.g*exp(-.05*(m-t0)))) %>%
      rename(Y.unt = Y) %>%
      mutate(Y = Y.unt + truth)
    return(data)
  }
}

time.trt.effects.stag <- function(data, tau, eq.by.grp){
  # For treatment effects that are the same in all treatment timing groups
  if (eq.by.grp){
    data <- data %>%
      mutate(truth=posttrt*(tau-tau*exp(-.05*(m-start.month)))) %>%
      rename(Y.unt = Y) %>%
      mutate(Y = Y.unt + truth)
    return(data)
  } else {
    data <- data %>% 
      # Treatment timing group has either "low" (early) or "high" (late) treatment response
      mutate(tau.g = ifelse(start.grp==1,tau-(.5*tau),tau+(.5*tau)),
             truth = posttrt*(tau.g-tau.g*exp(-.05*(m-start.month)))) %>%
      rename(Y.unt = Y) %>%
      mutate(Y = Y.unt + truth)
    return(data)
  }
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
  data <- data %>% 
    mutate(unitID=factor(unitID),
           monthofyear=factor(monthofyear),
           quarterofyear=factor(quarterofyear),
           m=factor(m),
           q=factor(q),
           year=factor(year))
  
  ## Aggregate to month
  monthly_data = data %>%
    group_by(unitID,m,monthofyear,q,quarterofyear,year) %>%
    summarise(across(c(Y,treatment,post,posttrt,truth),mean),.groups="keep") %>% ungroup() %>%
    mutate(postmonth = factor(ifelse(post>0,m,0)), 
           postmonthtrt = factor(ifelse(posttrt==1,m,0)))
  
  ## Aggregate to quarter
  quarterly_data = data %>%
    group_by(unitID,q,quarterofyear,year) %>%
    summarise(across(c(Y,treatment,post,posttrt,truth),mean),.groups="keep") %>% ungroup() %>%
    mutate(postquarter = factor(ifelse(post>0,q,0)), 
           postquartertrt = factor(ifelse(posttrt==1,q,0)))
  
  # Aggregate to year
  yearly_data = data %>%
    group_by(unitID,year) %>%
    summarise(across(c(Y,treatment,post,posttrt,truth),mean),.groups="keep") %>% ungroup() %>%
    mutate(postyear = factor(ifelse(post>0,year,0)),
           postyeartrt = factor(ifelse(posttrt==1,year,0)))
  
  if (is.tv){
    monthly_DID_2x2= estimatr::lm_robust(Y ~ post + treatment + postmonthtrt, data = monthly_data, clusters=unitID, se_type="stata")
    monthly_DID_0  = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth, data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_sm = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + monthofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_sq = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + quarterofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_um = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + m,data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_uq = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + q,data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_uy = estimatr::lm_robust(Y ~ postmonthtrt, fixed_effects = ~unitID + postmonth + year,data = monthly_data, clusters=unitID,se_type="stata")
    
    quarterly_DID_2x2= estimatr::lm_robust(Y ~ post + treatment + postquartertrt, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_0  = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_sq = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter + quarterofyear + year, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_uq = estimatr::lm_robust(Y ~ postquartertrt, fixed_effects = ~unitID + postquarter + q, data = quarterly_data, clusters=unitID,se_type="stata")
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
    monthly_DID_sm = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + monthofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_sq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + quarterofyear + year, data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_um = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + m,data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_uq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + q,data = monthly_data, clusters=unitID,se_type="stata")
    monthly_DID_uy = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + year,data = monthly_data, clusters=unitID,se_type="stata")
    
    quarterly_DID_2x2 =estimatr::lm_robust(Y~post+treatment+posttrt, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_0  = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_sq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + quarterofyear + year, data = quarterly_data, clusters=unitID,se_type="stata")
    quarterly_DID_uq = estimatr::lm_robust(Y~post+treatment+posttrt, fixed_effects = ~unitID + q, data = quarterly_data, clusters=unitID,se_type="stata")
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
                                     filter_fun(monthly_DID_sm, "monthly_sm"),
                                     filter_fun(monthly_DID_sq, "monthly_sq"),
                                     filter_fun(monthly_DID_um, "monthly_um"),
                                     filter_fun(monthly_DID_uq, "monthly_uq"),
                                     filter_fun(monthly_DID_uy, "monthly_uy")),month_true,by="time")
  quarterly_results <- left_join(rbind(filter_fun(quarterly_DID_2x2, "quarterly_2x2"),
                                       filter_fun(quarterly_DID_0, "quarterly_0"),
                                       filter_fun(quarterly_DID_sq, "quarterly_sq"),
                                       filter_fun(quarterly_DID_uq, "quarterly_uq"),
                                       filter_fun(quarterly_DID_uy, "quarterly_uy")),quarter_true,by="time")
  
  yearly_results <- left_join(rbind(filter_fun(yearly_DID_2x2, "yearly_2x2"),
                                    filter_fun(yearly_DID_0, "yearly_0"),
                                    filter_fun(yearly_DID_uy, "yearly_uy")),year_true,by="time")
  
  results <- bind_rows(monthly_results,quarterly_results,yearly_results)
  
  # Do the joint tests
  if (is.tv){
    test.results <- bind_rows(safe.joint.tests(monthly_DID_2x2)$result %>% mutate(group= "monthly_2x2",time="post"),
                              safe.joint.tests(monthly_DID_0)$result %>% mutate(group= "monthly_0",time="post"),
                              safe.joint.tests(monthly_DID_sm)$result %>% mutate(group= "monthly_sm",time="post"),
                              safe.joint.tests(monthly_DID_sq)$result %>% mutate(group= "monthly_sq",time="post"),
                              safe.joint.tests(monthly_DID_um)$result %>% mutate(group= "monthly_um",time="post"),
                              safe.joint.tests(monthly_DID_uq)$ result %>% mutate(group= "monthly_uq",time="post"),
                              safe.joint.tests(monthly_DID_uy)$ result %>% mutate(group= "monthly_uy",time="post"),
                              safe.joint.tests(quarterly_DID_2x2)$ result %>% mutate(group= "quarterly_2x2",time="post"),
                              safe.joint.tests(quarterly_DID_0)$ result %>% mutate(group= "quarterly_0",time="post"),
                              safe.joint.tests(quarterly_DID_sq)$result %>% mutate(group= "quarterly_sq",time="post"),
                              safe.joint.tests(quarterly_DID_uq)$result %>% mutate(group= "quarterly_uq",time="post"),
                              safe.joint.tests(quarterly_DID_uy)$result %>% mutate(group= "quarterly_uy",time="post"),
                              safe.joint.tests(yearly_DID_2x2)$result %>% mutate(group= "yearly_2x2",time="post"),
                              safe.joint.tests(yearly_DID_0)$result %>% mutate(group= "yearly_0",time="post"),
                              safe.joint.tests(yearly_DID_uy)$result %>% mutate(group= "yearly_uy",time="post"))
    results <- bind_rows(results,test.results)
  }
  return(results)
}

analyze.data.CS <- function(data){
  data <- data %>% 
    mutate(unitID=as.numeric(unitID)) # Requires a numeric unit ID variable
  
  ## Aggregate to month
  monthly_data = data %>%
    group_by(unitID,m,start.month) %>%
    summarise(Y=mean(Y),.groups="keep") %>% ungroup() 
  
  ## Aggregate to quarter
  quarterly_data = data %>%
    group_by(unitID,q,start.quarter) %>%
    summarise(Y=mean(Y),.groups="keep") %>% ungroup() 
  
  # Aggregate to year
  yearly_data = data %>%
    group_by(unitID,year,start.year) %>%
    summarise(Y=mean(Y),.groups="keep") %>% ungroup() 
  
  monthly_DID <- did::att_gt(yname = "Y",
                             tname = "m",
                             idname = "unitID",
                             gname = "start.month",
                             xformla = ~1,
                             data = monthly_data,
                             allow_unbalanced_panel = TRUE)
  
  quarterly_DID <- did::att_gt(yname = "Y",
                               tname = "q",
                               idname = "unitID",
                               gname = "start.quarter",
                               xformla = ~1,
                               data = quarterly_data,
                               allow_unbalanced_panel = TRUE)
  
  yearly_DID <-  did::att_gt(yname = "Y",
                             tname = "year",
                             idname = "unitID",
                             gname = "start.year",
                             xformla = ~1,
                             data = yearly_data,
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
  
  if (staggered) {
    null.data <- cons.trt.effects.stag(data,tau=0,         eq.by.grp = T)
    cons.data <- cons.trt.effects.stag(data,tau=params$tau,eq.by.grp = T)
    tv.data   <- time.trt.effects.stag(data,tau=params$tau,eq.by.grp = T)
    gv.data   <- cons.trt.effects.stag(data,tau=params$tau,eq.by.grp = F)
    gtv.data  <- time.trt.effects.stag(data,tau=params$tau,eq.by.grp = F)
    
    null.res <- analyze.data.CS(null.data)
    cons.res <- analyze.data.CS(cons.data)
    tv.res   <- analyze.data.CS(tv.data)
    gv.res   <- analyze.data.CS(gv.data)
    gtv.res  <- analyze.data.CS(gtv.data)
    
    results <- bind_rows(null.res %>% mutate(trteff='null'),
                         cons.res %>% mutate(trteff='constant'),
                         tv.res %>% mutate(trteff='time-varying'),
                         gv.res %>% mutate(trteff='group-varying'),
                         gtv.res %>% mutate(trteff='group- and time-varying')) %>%
      mutate(adoption='staggered')
  } else {
    null.data <- cons.trt.effects(data,tau=0, eq.by.grp = T)
    cons.data <- cons.trt.effects(data, tau=params$tau, eq.by.grp = T)
    tv.data   <- time.trt.effects(data, tau=params$tau, eq.by.grp = T)
    gv.data   <- cons.trt.effects(data, tau=params$tau, eq.by.grp = F)
    gtv.data  <- time.trt.effects(data, tau=params$tau, eq.by.grp = F)
    
    null.res <- analyze.and.combine(null.data)
    cons.res <- analyze.and.combine(cons.data)
    tv.res   <- analyze.and.combine(tv.data)
    gv.res   <- analyze.and.combine(gv.data)
    gtv.res  <- analyze.and.combine(gtv.data)
    
    results <- bind_rows(null.res %>% mutate(trteff='null'),
                         cons.res %>% mutate(trteff='constant'),
                         tv.res %>% mutate(trteff='time-varying'),
                         gv.res %>% mutate(trteff='group-varying'),
                         gtv.res %>% mutate(trteff='group- and time-varying')) %>%
      mutate(adoption='common')
  }
  return(results)
}

make.data <- function(params,staggered){
  bal.data   <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=F)
  unbal.data <- generate.data(params,month.byunit=T,quarter.byunit=T,unbalanced=T)
  
  if (staggered){
    ## Randomly assign treatment to units in balanced data 
    bal.rand.trt.dat <- bal.data %>% group_by(unitID) %>% 
      mutate(treatment=rbinom(1,1,.5),
             # randomly assign to one of two start groups: early (0) and late (1)
             start.grp=rbinom(1,1,.5)) %>%
      ungroup() %>%
      # Early group starts at year 3, late group starts at year 4
      mutate(start.year=treatment*ifelse(start.grp==1,3,4),
             start.month=treatment*(params$n.month*(start.year-1)+1),
             start.quarter=treatment*(params$n.quarter*(start.year-1)+1),
             post = m >= start.month & start.month!=0,
             posttrt=post*treatment)
    
    ## Assign treatment on the basis of pre-period trends in balanced data
    bal.trend.trt <- bal.data %>% group_by(unitID) %>%
      group_modify(~broom::tidy(lm(Y~year,data=.x))) %>% filter(term=="year") %>%
      # center the slopes (i.e., ignore any overall time trend)
      ungroup() %>% mutate(mean.slope=mean(estimate),
                           ctr.slope=estimate-mean.slope) %>% rowwise() %>%
      # Draw treatment assignment such that units with more positive slopes
      # have higher probability of being treated
      mutate(treatment=rbinom(1,1,plogis(ctr.slope)),
             start.grp=rbinom(1,1,.5))
    
    bal.trend.trt.dat <- left_join(bal.data,bal.trend.trt %>%
                                     select(unitID,treatment,start.grp),by="unitID") %>%
      ungroup() %>%
      mutate(start.year=treatment*ifelse(start.grp==1,3,4),
             start.month=treatment*(params$n.month*(start.year-1)+1),
             start.quarter=treatment*(params$n.quarter*(start.year-1)+1),
             post = m >= start.month & start.month != 0,
             posttrt=post*treatment)
    
    # Randomly assign treatment in unbalanced data
    unbal.rand.trt.dat <- unbal.data %>% group_by(unitID) %>% 
      mutate(treatment=rbinom(1,1,.5),
             start.grp=rbinom(1,1,.5)) %>%
      ungroup() %>%
      mutate(start.year=treatment*ifelse(start.grp==1,3,4),
             start.month=treatment*(params$n.month*(start.year-1)+1),
             start.quarter=treatment*(params$n.quarter*(start.year-1)+1),
             post = m >= start.month & start.month != 0,
             posttrt=post*treatment)
    
    ## Assign treatment on the basis of pre-period trends in unbalanced data
    unbal.trend.trt <- unbal.data %>% group_by(unitID) %>%
      group_modify(~broom::tidy(lm(Y~year,data=.x))) %>% filter(term=="year") %>%
      # center the slopes (i.e., ignore any overall time trend)
      ungroup() %>% mutate(mean.slope=mean(estimate,na.rm=T),
                           ctr.slope=estimate-mean.slope) %>% rowwise() %>%
      # Draw treatment assignment such that units with more positive slopes
      # have higher probability of being treated
      mutate(pr.trt=ifelse(!is.na(ctr.slope),plogis(ctr.slope),.5),
             treatment=rbinom(1,1,pr.trt),
             start.grp=rbinom(1,1,.5))
    
    unbal.trend.trt.dat <- left_join(unbal.data,unbal.trend.trt %>%
                                       select(unitID,treatment,start.grp),by="unitID") %>%
      ungroup() %>%
      mutate(start.year=treatment*ifelse(start.grp==1,3,4),
             start.month=treatment*(params$n.month*(start.year-1)+1),
             start.quarter=treatment*(params$n.quarter*(start.year-1)+1),
             post = m >= start.month & start.month != 0,
             posttrt=post*treatment)
  } else {
    ## Randomly assign treatment to units in balanced data 
    bal.rand.trt.dat <- bal.data %>% group_by(unitID) %>% 
      mutate(treatment=rbinom(1,1,.5)) %>%
      ungroup() %>%
      mutate(post=ifelse(year>params$n.year/2,1,0),
             posttrt=post*treatment)
    
    ## Assign treatment on the basis of pre-period trends in balanced data
    bal.trend.trt <- bal.data %>% group_by(unitID) %>%
      group_modify(~broom::tidy(lm(Y~year,data=.x))) %>% filter(term=="year") %>%
      # center the slopes (i.e., ignore any overall time trend)
      ungroup() %>% mutate(mean.slope=mean(estimate),
                           ctr.slope=estimate-mean.slope) %>% rowwise() %>%
      # Draw treatment assignment such that units with more positive slopes
      # have higher probability of being treated
      mutate(treatment=rbinom(1,1,plogis(ctr.slope)))
    
    bal.trend.trt.dat <- left_join(bal.data,bal.trend.trt %>%
                                     select(unitID,treatment),by="unitID") %>%
      mutate(post=ifelse(year>params$n.year/2,1,0),
             posttrt=post*treatment)
    
    # Randomly assign treatment in unbalanced data
    unbal.rand.trt.dat <- unbal.data %>% group_by(unitID) %>% 
      mutate(treatment=rbinom(1,1,.5)) %>%
      ungroup() %>%
      # Treatment starts at the midpoint of the timeseries
      mutate(post=ifelse(year>params$n.year/2,1,0),
             posttrt=post*treatment)
    
    ## Assign treatment on the basis of pre-period trends in unbalanced data
    unbal.trend.trt <- unbal.data %>% group_by(unitID) %>%
      group_modify(~broom::tidy(lm(Y~year,data=.x))) %>% filter(term=="year") %>%
      # center the slopes (i.e., ignore any overall time trend)
      ungroup() %>% mutate(mean.slope=mean(estimate,na.rm=T),
                           ctr.slope=estimate-mean.slope) %>% rowwise() %>%
      # Draw treatment assignment such that units with more positive slopes
      # have higher probability of being treated
      mutate(pr.trt=ifelse(!is.na(ctr.slope),plogis(ctr.slope),.5),
             treatment=rbinom(1,1,pr.trt))
    
    unbal.trend.trt.dat <- left_join(unbal.data,unbal.trend.trt %>%
                                       select(unitID,treatment),by="unitID") %>%
      # Treatment starts at the midpoint of the timeseries
      mutate(post=ifelse(year>params$n.year/2,1,0),
             posttrt=post*treatment)
  }
  
  return(list('bal.rand.trt.dat'=bal.rand.trt.dat,
              'bal.trend.trt.dat'=bal.trend.trt.dat,
              'unbal.rand.trt.dat'=unbal.rand.trt.dat,
              'unbal.trend.trt.dat'=unbal.trend.trt.dat))
}

single_iteration <- function(params,staggered){
  all_data <- make.data(params,staggered)
  
  bal.rand.results <-  inject.analyze(all_data$bal.rand.trt.dat,params,staggered) %>%
    mutate(panel='balanced',assignment='random')
  
  bal.trend.results <- inject.analyze(all_data$bal.trend.trt.dat,params,staggered) %>%
    mutate(panel='balanced',assignment='trended')
  
  unbal.rand.results <- inject.analyze(all_data$unbal.rand.trt.dat,params,staggered) %>%
    mutate(panel='unbalanced',assignment='random')
  
  unbal.trend.results <- inject.analyze(all_data$unbal.trend.trt.dat,params,staggered) %>%
    mutate(panel='unbalanced',assignment='trended')
  
  return(rbind(bal.rand.results,bal.trend.results,unbal.rand.results,unbal.trend.results))
}