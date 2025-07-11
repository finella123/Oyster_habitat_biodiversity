# This is my Master's thesis project- exploring the effects of oyster trophic interactions on the benthic community
#This script is to anlyze the data
source("scripts2/01_install_packages.R")

com1 <- readRDS("wdata2/wd_community_data_m^2.rds")

reef1 <- readRDS("wdata2/wd_reef_characteristics_m^2.rds")


reef2 <- reef1%>%
  select(site, sample, season, oys.density)%>%
  rename(Site=site,Sample=sample)%>%
  filter(Site!=1)

com2 <- com1%>%
  group_by(Site, Sample, season) %>%
  summarize(
    spr= n_distinct(TaxaID),
    abm2=sum(Abundance.m2))
  
rac <- reef2 %>%
  left_join(com2, by = c("Site", "Sample", "season"))






##### ABUNDANCE ----
### sample level LOV
lod.ab.samp <- glmmTMB(log(t.abund) ~ mlov*season,data = acon%>% mutate(season = relevel(season,ref="w")))
summary(mlovt.abund2)
glmm.resids(mlovt.abund2)
Anova(mlovt.abund2, type="III")