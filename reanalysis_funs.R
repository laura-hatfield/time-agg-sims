analyze.CS.applied <- function(data,ant.pers=0){
  agg.dat <- agg.data.staggered(data)
  
  monthly_DID <- did::att_gt(yname = "Y",
                             tname = "m",
                             idname = "unitID",
                             gname = "start.month",
                             control_group = "notyettreated",
                             anticipation=ant.pers,
                             cband=FALSE,
                             xformla = ~1,
                             data = agg.dat$monthly,
                             allow_unbalanced_panel = TRUE)
  
  
  quarterly_DID <- did::att_gt(yname = "Y",
                               tname = "q",
                               idname = "unitID",
                               gname = "start.quarter",
                               control_group = "notyettreated",
                               cband=FALSE,
                               anticipation=ant.pers,
                               xformla = ~1,
                               data = agg.dat$quarterly,
                               allow_unbalanced_panel = TRUE)
  
  yearly_DID <-  did::att_gt(yname = "Y",
                             tname = "year",
                             idname = "unitID",
                             gname = "start.year",
                             control_group = "notyettreated",
                             cband=FALSE,
                             anticipation=ant.pers,
                             xformla = ~1,
                             data = agg.dat$yearly,
                             allow_unbalanced_panel = TRUE)
  
  month.agg.es <- did::aggte(monthly_DID, type = "dynamic", na.rm=TRUE)
  quarter.agg.es <- did::aggte(quarterly_DID, type = "dynamic", na.rm=TRUE)
  year.agg.es <- did::aggte(yearly_DID, type = "dynamic", na.rm=TRUE)
  
  dynamic <- bind_rows(tibble('est'= month.agg.es$att.egt, 'se'= month.agg.es$se.egt, 'time'=month.agg.es$egt, 'agg'="month"),
                       tibble('est'= quarter.agg.es$att.egt, 'se'= quarter.agg.es$se.egt, 'time'=quarter.agg.es$egt, 'agg'="quarter"),
                       tibble('est'= year.agg.es$att.egt, 'se'= year.agg.es$se.egt, 'time'=year.agg.es$egt, 'agg'="year")) %>%
    mutate(est=est*100,
           lb=est-qnorm(.975)*se*100,
           ub=est+qnorm(.975)*se*100)
  
  # Aggregating to get overall effects
  # 2 years
  month.agg.gs.2y <- did::aggte(monthly_DID, type = "group", na.rm=TRUE, max_e=24)
  quarter.agg.gs.2y <- did::aggte(quarterly_DID, type = "group", na.rm=TRUE, max_e=8) 
  year.agg.gs.2y <- did::aggte(yearly_DID, type = "group", na.rm=TRUE, max_e=2) 
  # 3 years 
  month.agg.gs.3y <- did::aggte(monthly_DID, type = "group", na.rm=TRUE, max_e=36)
  quarter.agg.gs.3y <- did::aggte(quarterly_DID, type = "group", na.rm=TRUE, max_e=12) 
  year.agg.gs.3y <- did::aggte(yearly_DID, type = "group", na.rm=TRUE, max_e=3) 
  # All periods
  month.agg.gs.all <- did::aggte(monthly_DID, type = "group", na.rm=TRUE)
  quarter.agg.gs.all <- did::aggte(quarterly_DID, type = "group", na.rm=TRUE) 
  year.agg.gs.all <- did::aggte(yearly_DID, type = "group", na.rm=TRUE) 
  
  static <- bind_rows(tibble('est'= month.agg.gs.2y$overall.att, 'se'= month.agg.gs.2y$overall.se, 'agg'="month",period="2 years"),
                      tibble('est'= quarter.agg.gs.2y$overall.att, 'se'= quarter.agg.gs.2y$overall.se, 'agg'="quarter", period="2 years"),
                      tibble('est'= year.agg.gs.2y$overall.att, 'se'= year.agg.gs.2y$overall.se, 'agg'="year", period="2 years"),
                      tibble('est'= month.agg.gs.3y$overall.att, 'se'= month.agg.gs.3y$overall.se, 'agg'="month",period="3 years"),
                      tibble('est'= quarter.agg.gs.3y$overall.att, 'se'= quarter.agg.gs.3y$overall.se, 'agg'="quarter", period="3 years"),
                      tibble('est'= year.agg.gs.3y$overall.att, 'se'= year.agg.gs.3y$overall.se, 'agg'="year", period="3 years"),
                      tibble('est'= month.agg.gs.all$overall.att, 'se'= month.agg.gs.all$overall.se, 'agg'="month",period="All"),
                      tibble('est'= quarter.agg.gs.all$overall.att, 'se'= quarter.agg.gs.all$overall.se, 'agg'="quarter", period="All"),
                      tibble('est'= year.agg.gs.all$overall.att, 'se'= year.agg.gs.all$overall.se, 'agg'="year", period="All")) %>%
    mutate(est=est*100,
           lb=est-qnorm(.975)*se*100,
           ub=est+qnorm(.975)*se*100)
  
  return(list('dynamic'=dynamic,'static'=static))
}
