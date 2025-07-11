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

#function for residuals
glmm.resids<-function(model){
  t1 <- simulateResiduals(model)
  print(testDispersion(t1))
  plot(t1)
}
##### ABUNDANCE ----

### sample level LOD
lod.ab.samp <- glmmTMB(log(abm2) ~ oys.density*season,data = rac)

summary(lod.ab.samp)
glmm.resids(lod.ab.samp)
Anova(lod.ab.samp, type="III")

# Extract residuals
res <- residuals(lod.ab.samp)
shapiro.test(res)

# plot log(abm2) vs oyster density, color by season
ggplot(rac, aes(x = oys.density, y = abm2, color = season)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  facet_wrap(~ season) +
  labs(
    x = "Oyster Density (per m²)",
    y = "Abundance per m²",
    title = "Relationship Between Abundance and Oyster Density by Season"
  ) +
  theme_minimal()

# Create new data frame for prediction
newdata <- expand.grid(
  oys.density = seq(min(rac$oys.density), max(rac$oys.density), length.out = 100),
  season = unique(rac$season)
)

# Predict log(abm2)
newdata$log_abm2 <- predict(lod.ab.samp, newdata = newdata)

# Back-transform
newdata$abm2_pred <- exp(newdata$log_abm2)

ggplot(newdata, aes(x = oys.density, y = abm2_pred, color = season)) +
  geom_line(size = 1.2) +
  labs(
    x = "Oyster Density (per m²)",
    y = "Predicted Abundance per m²",
    title = "Predicted Abundance vs Oyster Density by Season"
  ) +
  theme_minimal()+ geom_point(data = rac, aes(x = oys.density, y = abm2), alpha = 0.3, inherit.aes = FALSE)

#look into the high outlier 180000\
# Sort and view top 5
head(rac[order(-rac$abm2), ], 5)




#SPECIES RICHNESS
### sample level LOD
lod.spr.samp <- glmmTMB(sprm2 ~ oys.density*season,data = rac)

summary(lod.spr.samp)
glmm.resids(lod.spr.samp)
Anova(lod.spr.samp, type="III")

# Extract residuals
res <- residuals(lod.spr.samp)
shapiro.test(res)

ggplot(data = rac, aes(x = oys.density, y = sprm2, color = season)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(
    title = "Effect of Oyster Density on Species Richness by Season",
    x = "Oyster Density",
    y = "Species Richness (per m²)"
  ) +
  theme_minimal()
