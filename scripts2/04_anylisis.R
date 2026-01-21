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
  group_by(Site,season) %>%
  mutate(moys.density = mean(oys.density))

rac$season <- as.factor(rac$season)
rac$Site <- as.factor(rac$Site)
rac$Sample <- as.factor(rac$Sample)
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
  
ffg5 <- ffg4%>%
  group_by(Site, Sample, season,ffg) %>%
  summarize(
    ffabm2=sum(Abundance.m2))%>%
  left_join(reef2)%>%
  group_by(Site) %>%
  mutate(moys.density = mean(oys.density))

ffgcomcomp <- ffg5[, c("Site", "Sample", "season", "oys.density", "moys.density", "ffg", "ffabm2")]%>%
  pivot_wider(names_from = ffg, values_from = ffabm2, values_fill = 0)

ffg_long <- ffgcomcomp %>%
  pivot_longer(
    cols = c(c,d,o,p,s,g),   
    names_to = "FFG",
    values_to = "ffabm2"
  )
ffg_long$FFG <- factor(ffg_long$FFG)
ffg_long$season <- factor(ffg_long$season)
ffg_long$Site <- factor(ffg_long$Site)
ffg_long$Sample <- factor(ffg_long$Sample)

#oyster predator dataset
oys.pred <- ffg4%>%
  group_by(Site, Sample, season,oyster_predator) %>%
  summarize(
    opabm2=sum(Abundance.m2))%>%
  left_join(reef2)%>%
  group_by(Site) %>%
  mutate(moys.density = mean(oys.density))%>%
  filter(oyster_predator=='y')%>%
  select(-oyster_predator)
oys.pred$Site <- factor(oys.pred$Site)
oys.pred$Sample <- factor(oys.pred$Sample)
oys.pred$season <- factor(oys.pred$season)

#function for residuals
glmm.resids<-function(model){
  t1 <- simulateResiduals(model)
  print(testDispersion(t1))
  plot(t1)
}
#live oyster density rabbit hole----- 
#CANNOT GROUP SAMPLES ACROSS SEASON BECAUSE THEY ARE UNIQUE SAMPPLES
# 
# site_season_mod <- glmmTMB(
#   oys.density ~ season + (1 | Site),
#   data = rac
# )
# Anova(site_season_mod, type = "III")
# 
# 
# sample_season_mod <- glmmTMB(
#   oys.density ~ season + (1 | Site/Sample),
#   data = rac
# )
# Anova(sample_season_mod, type = "III")
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
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA)+
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

w.lod.ab.nmds<-data.frame(vegan::scores(lod.ab.nmds,choices=c(1,2),display="sites"),LOD=comcomp2$oys.density, sample=comcomp2$Sample,season=comcomp2$season)


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
  theme(legend.position = "right")+ stat_ellipse(aes(group = season), linetype = "dashed", color = "gray40")+ scale_color_viridis_c(option = "plasma")

# 
# 
# lod.ab.nmdsplot <- ggplot(w.lod.ab.nmds, aes(y=NMDS2,x=NMDS1,color=LOD, shape=season))+
#   geom_point(size=2)+
#   theme(plot.title = element_text(hjust=0.5)) +
#   ggtitle("Abundance")
# lod.ab.nmdsplot
# ggplot(w.lod.ab.nmds, aes(NMDS1, NMDS2, color = LOD)) +
#   geom_point(size = 3) +
#   scale_color_viridis_c() +
#   facet_wrap(~ season) +
#   theme_classic() +
#   labs(
#     title = "NMDS by Season",
#     color = "Live Oyster Density (m⁻²)"
#   )
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
###! microhabitat FFG abundance----
### c= carnivore-----
carnivore_mod <- glmmTMB(
  log(ffabm2 + 1) ~ oys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "c")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(carnivore_mod)
summary(carnivore_mod)
Anova(carnivore_mod, type = "III")

# Extract residuals
res <- residuals(carnivore_mod)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(carnivore_mod, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Carnivore Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### d= deposit feeders-----
deposit_mod <- glmmTMB(
  log(ffabm2 + 1) ~ oys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "d")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(deposit_mod)
summary(deposit_mod)
Anova(deposit_mod, type = "III")

# Extract residuals
res <- residuals(deposit_mod)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(deposit_mod, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Deposit feeders Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### o=omnivore-----
omnivore_mod <- glmmTMB(
  log(ffabm2 + 1) ~ oys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "o")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(omnivore_mod)
summary(omnivore_mod)
Anova(omnivore_mod, type = "III")

# Extract residuals
res <- residuals(omnivore_mod)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(omnivore_mod, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted omnivore Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### p=parasite-----
#tried lots but low abundance 
parasite_mod <- glmmTMB(
  log(ffabm2+1) ~ oys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "p")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(parasite_mod)
summary(parasite_mod)
Anova(parasite_mod, type = "III")

# Extract residuals
res <- residuals(parasite_mod)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(parasite_mod, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted parasite Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### s = suspension feeder-----
suspension_mod <- glmmTMB(
  log(ffabm2+1)~ oys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "s")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(suspension_mod)
summary(suspension_mod)
Anova(suspension_mod, type = "III")

# Extract residuals
res <- residuals(suspension_mod)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(suspension_mod, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted suspension feeder Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### g=grazers-----
grazer_mod <- glmmTMB(
  log(ffabm2+1) ~ oys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "g")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(grazer_mod)
summary(grazer_mod)
Anova(grazer_mod, type = "III")

# Extract residuals
res <- residuals(grazer_mod)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(grazer_mod, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted grazer Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


### oyster consumers ---------
oys_pred_mod <- glmmTMB(
  opabm2 ~ oys.density * season + (1 | Site),
  data = oys.pred %>% mutate(season = relevel(season,ref="w")))
glmm.resids(oys_pred_mod)
summary(oys_pred_mod)
Anova(oys_pred_mod, type = "III")

# Extract residuals
res <- residuals(oys_pred_mod)
shapiro.test(res)

# plot oys pred abm2 vs oyster density, color by season
preds <- ggpredict(oys_pred_mod, terms = c("oys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Oyster Predator Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
### reef level FFG abundance-----
### c= carnivore-----
carnivore_mod_reef <- glmmTMB(
  log(ffabm2 + 1) ~ moys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "c")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(carnivore_mod_reef)
summary(carnivore_mod_reef)
Anova(carnivore_mod_reef, type = "III")

# Extract residuals
res <- residuals(carnivore_mod_reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(carnivore_mod_reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Carnivore Abundance by Mean Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### d= deposit feeders-----
deposit_mod_reef <- glmmTMB(
  log(ffabm2 + 1) ~ moys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "d")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(deposit_mod_reef)
summary(deposit_mod_reef)
Anova(deposit_mod_reef, type = "III")

# Extract residuals
res <- residuals(deposit_mod_reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(deposit_mod_reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Deposit feeders Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### o=omnivore-----
omnivore_mod_reef <- glmmTMB(
  log(ffabm2 + 1) ~ moys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "o")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(omnivore_mod_reef)
summary(omnivore_mod_reef)
Anova(omnivore_mod_reef, type = "III")

# Extract residuals
res <- residuals(omnivore_mod_reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(omnivore_mod_reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted omnivore Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### p=parasite-----
#tried lots but low abundance 
parasite_mod_reef <- glmmTMB(
  log(ffabm2+1) ~ moys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "p")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(parasite_mod_reef)
summary(parasite_mod_reef)
Anova(parasite_mod_reef, type = "III")

# Extract residuals
res <- residuals(parasite_mod_reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(parasite_mod_reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted parasite Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### s = suspension feeder-----
suspension_mod_reef <- glmmTMB(
  log(ffabm2+1)~ moys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "s")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(suspension_mod_reef)
summary(suspension_mod_reef)
Anova(suspension_mod_reef, type = "III")

# Extract residuals
res <- residuals(suspension_mod_reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(suspension_mod_reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted suspension feeder Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

### g=grazers-----
grazer_mod_reef <- glmmTMB(
  log(ffabm2+1) ~ moys.density * season + (1 | Site),
  data = ffg_long %>% filter(FFG == "g")%>% mutate(season = relevel(season,ref="w")))
glmm.resids(grazer_mod_reef)
summary(grazer_mod_reef)
Anova(grazer_mod_reef, type = "III")

# Extract residuals
res <- residuals(grazer_mod_reef)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
preds <- ggpredict(grazer_mod_reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted grazer Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


### oyster consumers ---------
oys_pred_mod_reef <- glmmTMB(
  log(opabm2) ~ moys.density * season + (1 | Site),
  data = oys.pred %>% mutate(season = relevel(season,ref="w")))
glmm.resids(oys_pred_mod_reef)
summary(oys_pred_mod_reef)
Anova(oys_pred_mod_reef, type = "III")

# plot oys pred abm2 vs oyster density, color by season
preds <- ggpredict(oys_pred_mod_reef, terms = c("moys.density", "season"))

ggplot(preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "Live Oyster Density",
    y = expression(Predicted~abundance~(m^2)),
    color = "Season",
    fill = "Season",
    title = "Predicted Oyster Predator Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
### site level FFG abndance-----
##### FFG COMMUNITY COMPOSITION-----




####simper mess around-------
comm_mat <- ffgcomcomp %>%
  mutate(
    Site = as.character(Site),
    Sample = as.character(Sample),
    season = as.character(season),
    SampleID = paste(Site, Sample, season, sep = "_")
  ) %>%
  select(-Site, -Sample, -season) %>%
  column_to_rownames("SampleID")%>%
  ungroup()%>%
  select(-Site, -oys.density, -moys.density)
any(duplicated(rownames(comm_mat)))
str(comm_mat)

comm_tr <- decostand(comm_mat, method = "hellinger")

season_vec <- sub(".*_", "", rownames(comm_tr))  # takes everything after the last "_"
table(season_vec)  # sanity check counts per season

simper_out <- simper(
  comm_tr,
  group = season_vec,
  permutations = 999
)

summary(simper_out)
simper_summary <- summary(simper_out)

# Loop over each pairwise comparison
top_contributors <- lapply(names(simper_summary), function(pair) {
  df <- as.data.frame(simper_summary[[pair]])
  df$FFG <- rownames(df)
  
  # Find the FFG with the highest average contribution
  top <- df[which.max(df$average), c("FFG", "average", "cumsum")]
  top$comparison <- pair
  top
})

# Combine into one data frame
top_contributors_df <- do.call(rbind, top_contributors)
row.names(top_contributors_df) <- NULL

top_contributors_df







####NMDS mess around -------
mlod.ab.nmds<-metaMDS(comcomp2[,-1:-5], distance="bray", k=2, trymax = 1000)

# ---------------------------------------------------------
# 1. Separate metadata and species matrix
# ---------------------------------------------------------

meta <- comcomp2 %>%
  select(Site, Sample, season, oys.density, moys.density)

species <- comcomp2 %>%
  select(-(Site: moys.density))   # drops first 5 columns, keeps species only
# 1. Hellinger transform 
species <- decostand(species, method = "hellinger")
# ---------------------------------------------------------
# 2. Run NMDS (Bray-Curtis)
# ---------------------------------------------------------

set.seed(123)
nmds <- metaMDS(species, distance = "bray", k = 2, trymax = 200)

# Extract NMDS scores for plotting
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores <- bind_cols(meta, nmds_scores)

# ---------------------------------------------------------
# 3. Fit species vectors (envfit)
# ---------------------------------------------------------

fit <- envfit(nmds, species, permutations = 999)

# Extract significant vectors (p < 0.05)
vec <- as.data.frame(scores(fit, display = "vectors"))
vec$species <- rownames(vec)
#vec_sig <- vec[fit$vectors$pvals < 0.003, ]
# Calculate vector lengths 
vec$length <- sqrt(vec$NMDS1^2 + vec$NMDS2^2) 
# Select top 5 longest vectors 
vec_top5 <- vec %>% arrange(desc(length)) %>% slice(1:5)
# ---------------------------------------------------------
# 4. Plot NMDS with ggplot2
# ---------------------------------------------------------

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = season, color = oys.density), size = 3) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  geom_segment(data = vec_sig,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "red", linewidth = 0.8) +
  geom_text(data = vec_sig,
            aes(x = NMDS1, y = NMDS2, label = species),
            color = "red", vjust = -0.5, size = 3) +
  theme_minimal(base_size = 14) +
  labs(title = "NMDS of Community Composition",
       color = "Oyster Density",
       shape = "Season")


ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(shape = season, color = oys.density), size = 3) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  
  # Ellipses by season
  stat_ellipse(aes(group = season), type = "t", linewidth = 1) +
  
  # Species vectors
  geom_segment(data = vec_sig,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "red", linewidth = 0.8) +
  geom_text(data = vec_sig,
            aes(label = species),
            color = "red", vjust = -0.5, size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(title = "NMDS (Hellinger-transformed community)",
       color = "Oyster Density",
       shape = "Season")

ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(shape = season, color = oys.density), size = 3, alpha = 1) +
  #scale_color_gradient(low = "lightblue", high = "darkblue") +
  stat_ellipse(aes(group = season),  linetype = "dashed", color = "gray40")+ 
  scale_color_viridis_c(option = "plasma") +
  
  # Top 5 species vectors
  geom_segment(data = vec_top5,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "red", linewidth = 0.8) +
  geom_text(data = vec_top5,
            aes(label = species),
            color = "red", vjust = -0.5, size = 3) +
  
  theme_minimal(base_size = 14) +
  labs(title = "NMDS (Hellinger-transformed community)",
       color = "Live Oyster Density",
       shape = "Season")




