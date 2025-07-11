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
    sprm2= n_distinct(TaxaID),
    abm2=sum(Abundance.m2))
  
rac <- reef2 %>%
  left_join(com2, by = c("Site", "Sample", "season"))%>%
  drop_na() # I have 16 sites with no com data because I only processed 6/7 per site
any(is.na(rac))# final check for issues- looks good!


##### ABUNDANCE ----
### sample level LOD
lod.ab.samp <- glmmTMB(abm2 ~ oys.density*season,data = rac)


summary(lod.ab.samp)
glmm.resids(lod.ab.samp)
Anova(lod.ab.samp, type="III")