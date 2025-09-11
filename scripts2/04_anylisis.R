# This is my Master's thesis project- exploring the effects of oyster trophic interactions on the benthic community
#This script is to anlyze the data
source("scripts2/01_install_packages.R")

#loading data and getting into correct format for analysis-----
com1 <- readRDS("wdata2/wd_community_data_m^2.rds")

reef1 <- readRDS("wdata2/wd_reef_characteristics_m^2.rds")

ffg3 <- readRDS("wdata2/functional_feeding_groups_data.rds")

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
  drop_na()%>% # I have 16 sites with no com data because I only processed 6/7 per site
  group_by(Site, season) %>%
  mutate(moys.density = mean(oys.density)) %>%
  ungroup()

rac$season <- as.factor(rac$season)
any(is.na(rac))# final check for issues- looks good!

comcomp <- com1%>%
  select(Site, Sample, season, scientific, Abundance.m2)%>%
  ungroup()%>%
  select(-TaxaID)%>%
  group_by(Site,Sample,season)%>%
  left_join(reef2)%>%
  group_by(Site, season) %>%
  mutate(moys.density = mean(oys.density)) %>%
  ungroup()
comcomp2 <- comcomp[, c("Site", "Sample", "season", "oys.density", "moys.density", "scientific", "Abundance.m2")]%>%
  pivot_wider(names_from = scientific, values_from = Abundance.m2, values_fill = 0)


#data set to anylyze functional feeding groups

ffg4 <- com1%>%
  left_join(ffg3, by = "TaxaID")
#function for residuals
glmm.resids<-function(model){
  t1 <- simulateResiduals(model)
  print(testDispersion(t1))
  plot(t1)
}
##### ABUNDANCE ----
###! microhabitat abundance-----
lod.ab.samp <- glmmTMB(log(abm2) ~ oys.density*season+ (1|Site) ,data = rac%>% mutate(season = relevel(season,ref="w")))
summary(lod.ab.samp)
glmm.resids(lod.ab.samp)
Anova(lod.ab.samp, type="III")


# Extract residuals
res <- residuals(lod.ab.samp)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(lod.ab.samp, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Benthic Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

###! reef level abundance-----
mlod.ab.reef <- glmmTMB(log(abm2) ~ moys.density*season+ (1|Site) ,data = rac%>% mutate(season = relevel(season,ref="w")))
#using this random effect to account for the shared structure of multiple samples having the sample mlod- so the model doesnt think each mlod is a unique value
summary(mlod.ab.reef)
glmm.resids(mlod.ab.reef)
Anova(mlod.ab.reef, type="III")

# Extract residuals
res <- residuals(mlod.ab.reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(mlod.ab.reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = " Reef Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Benthic Abundance by reef Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
### site level abundance-----
#### SPECIES RICHNESS-----
###! microhabitat richness-----
lod.sprm2.samp <- glmmTMB(sprm2 ~ oys.density*season+ (1|Site) ,data = rac%>% mutate(season = relevel(season,ref="w")))
summary(lod.sprm2.samp)
glmm.resids(lod.sprm2.samp)
Anova(lod.sprm2.samp, type="III")

# Extract residuals
res <- residuals(lod.sprm2.samp)
shapiro.test(res)

# plot sprm2 vs oyster density, color by season
preds <- ggpredict(lod.sprm2.samp, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Benthic Richness by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


###! reef level richness-----
mlod.spr.reef <- glmmTMB(sprm2 ~ moys.density*season+ (1|Site) ,data = rac%>% mutate(season = relevel(season,ref="w")))
#using this random effect to account for the shared structure of multiple samples having the sample mlod- so the model doesnt think each mlod is a unique value
summary(mlod.spr.reef)
glmm.resids(mlod.spr.reef)
Anova(mlod.spr.reef, type="III")

# Extract residuals
res <- residuals(mlod.spr.reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(mlod.spr.reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = " Reef Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Benthic Abundance by reef Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
### site level richness-----
##### COMMUNITY COMPOSITION-----
###! microhabitat comp-----
#nmds for abundance between seasons and live oyster density

lod.ab.nmds<-metaMDS(comcomp2[,-1:-5], distance="bray", k=2, trymax = 1000)

w.lod.ab.nmds<-data.frame(scores(lod.ab.nmds,choices=c(1,2),display="sites"),LOD=comcomp2$oys.density, sample=comcomp2$Sample,season=comcomp2$season)

lod.ab.nmdsplot <- ggplot(w.lod.ab.nmds, aes(y=NMDS2,x=NMDS1,color=LOD, shape=season))+
  geom_point(size=2)+
   theme(plot.title = element_text(hjust=0.5)) +
   ggtitle("Abundance")
lod.ab.nmdsplot

ggplot(w.lod.ab.nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = LOD, shape = season), size = 3, alpha = 0.9) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "NMDS of Benthic Community Composition",
    x = "NMDS1",
    y = "NMDS2",
    color = "Live Oyster Density",
    shape = "Season"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")+ stat_ellipse(aes(group = season), linetype = "dashed", color = "gray40")+ scale_color_viridis_c(option = "plasma")+ facet_wrap(~season)
###! reef level comp------
#nmds for abundance between seasons and live oyster density

mlod.ab.nmds<-metaMDS(comcomp2[,-1:-5], distance="bray", k=2, trymax = 1000)

w.mlod.ab.nmds<-data.frame(scores(mlod.ab.nmds,choices=c(1,2),display="sites"),mLOD=comcomp2$moys.density, sample=comcomp2$Sample,season=comcomp2$season, site=comcomp2$Site)

ggplot(w.mlod.ab.nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = mLOD, shape = season), size = 3, alpha = 0.9) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "NMDS of Benthic Community Composition",
    x = "NMDS1",
    y = "NMDS2",
    color = "mLOD",
    shape = "Season"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")+ 
  stat_ellipse(aes(group = season), linetype = "dashed", color = "gray40")+
  scale_color_viridis_c(option = "plasma")+ 
    facet_wrap(~site)
### site level comp-----




##### FFG ABUNDANCE ----
### microhabitat FFG abundance-----
lod.ffg.samp <- glmmTMB(log(abm2) ~ oys.density*season+ (1|Site) ,data = rac%>% mutate(season = relevel(season,ref="w")))
summary(lod.ffg.samp)
glmm.resids(lod.ffg.samp)
Anova(lod.ffg.samp, type="III")

# Extract residuals
res <- residuals(lod.ffg.samp)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(lod.ab.samp, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Benthic Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### reef level FFG abundance-----
### site level FFG abndance-----
#### FFG SPECIES RICHNESS-----
### microhabitat FFG richness-----
lod.sprm2.samp <- glmmTMB(sprm2 ~ oys.density*season+ (1|Site) ,data = rac)
summary(lod.sprm2.samp)
glmm.resids(lod.sprm2.samp)
Anova(lod.sprm2.samp, type="III")

# Extract residuals
res <- residuals(lod.sprm2.samp)
shapiro.test(res)

# plot sprm2 vs oyster density, color by season
preds <- ggpredict(lod.sprm2.samp, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Benthic Richness by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


### reef level FFG richness-----
### site level FFG richness-----
##### FFG COMMUNITY COMPOSITION-----
### microhabitat FFG comp-----
### reef level FFG comp-----
### site level FFG comp-----