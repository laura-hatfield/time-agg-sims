library(did)
library(tidyverse)
theme_set(theme_minimal()+theme(legend.position="bottom"))
source("sim_funs.R")
source("reanalysis_funs.R")
load("cleaned_force_data.RData")

## Aggregate data at each level
## Truth is not used; just to make data compatible with the agg.data.staggered function from the simulations
force.results <- analyze.CS.applied(force.dat %>% mutate(truth=0))

dynamic <- force.results$dynamic %>%
  mutate(inscope=ifelse(agg=="month",(time>=-36)&(time<=36),ifelse(agg=="quarter",(time>=-12)&(time<=12),ifelse(agg=="year",(time>=-3)&(time<=3),FALSE)))) %>%
  filter(inscope)

ggplot(dynamic,aes(x=time)) +  geom_point(aes(y=est)) +
  geom_segment(aes(y = lb, yend = ub),col = grey(0.2)) +
  facet_wrap(~agg,scale="free_x") +
  geom_hline(yintercept = 0, color = grey(0.4), linetype = "dashed") +
  geom_vline(xintercept = 0, color = grey(0.4), linetype = "dashed") +
  scale_x_continuous(name= "Time relative to training") +
  scale_y_continuous("Estimated ATT") +
  theme(panel.grid.major = element_blank())
ggsave("figures/reanalysis_dynamic.png",width=6.5,height=4)


ggplot(force.results$static,aes(x=agg,group=period)) +   
  geom_segment(aes(y=lb, yend=ub,col=period),position=position_dodge(width=0.75)) +
  geom_label(aes(y=est, label = round(est,2),col=period),position=position_dodge(width=0.75)) +
  scale_color_grey(start=0,end=0.6,guide="none") +
  geom_hline(yintercept = 0, color = grey(0.4), linetype = "dashed") +
  scale_x_discrete("Aggregation", labels = function(x) str_to_title(x)) +
  scale_y_continuous("Estimated ATT")
ggsave("figures/reanalysis_static.png",width=6.5,height=3)

## Re-fitting with a different reference period
force.results.1antper <- analyze.CS.applied(force.dat %>% mutate(truth=0),ant.pers=1)
dynamic <- force.results.1antper$dynamic %>%
  mutate(inscope=ifelse(agg=="month",(time>=-36)&(time<=36),ifelse(agg=="quarter",(time>=-12)&(time<=12),ifelse(agg=="year",(time>=-3)&(time<=3),FALSE)))) %>%
  filter(inscope)

ggplot(dynamic,aes(x=time)) +  geom_point(aes(y=est)) +
  geom_segment(aes(y = lb, yend = ub),col = grey(0.2)) +
  facet_wrap(~agg,scale="free_x") +
  geom_hline(yintercept = 0, color = grey(0.4), linetype = "dashed") +
  geom_vline(xintercept = 0, color = grey(0.4), linetype = "dashed") +
  scale_x_continuous(name= "Time relative to training") +
  scale_y_continuous("Estimated ATT") +
  theme(panel.grid.major = element_blank())
ggsave("figures/reanalysis_dynamic_1antper.png",width=6.5,height=4)


ggplot(force.results.1antper$static,aes(x=agg,group=period)) +   
  geom_segment(aes(y=lb, yend=ub,col=period),position=position_dodge(width=0.75)) +
  geom_label(aes(y=est, label = round(est,2),col=period),position=position_dodge(width=0.75)) +
  scale_color_grey(start=0,end=0.6,guide="none") +
  geom_hline(yintercept = 0, color = grey(0.4), linetype = "dashed") +
  scale_x_discrete("Aggregation", labels = function(x) str_to_title(x)) +
  scale_y_continuous("Estimated ATT")
ggsave("figures/reanalysis_static_1antper4.png",width=6.5,height=3)