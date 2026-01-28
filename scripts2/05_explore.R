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

print(comcomp2)
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
library(ggplot2)

ggplot(rac, aes(x = Site, y = oys.density)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Oyster Density by Site",
       x = "Site",
       y = "Oyster Density") +
  theme_minimal()

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
ggplot(rac, aes(x = oys.density,
                y = abm2,
                color = season)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = "Live Oyster Density",
    y = expression(Abundance~(m^2)),
    color = "Season",
    title = "Raw Benthic Abundance by Oyster Density and Season"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


ggplot(rac, aes(x = oys.density,
                y = log(abm2),
                color = Site)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    color = "Season",
    title = "Raw Log-Abundance by Oyster Density and site"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggplot(rac, aes(x = oys.density)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    x = "Live Oyster Density",
    y = "Number of Samples",
    title = "Distribution of Samples Across Live Oyster Density"
  ) +
  theme_minimal()

ggplot(rac, aes(x = oys.density, y = abm2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(y = abm2), method = "loess")+
  labs(
    x = "Live Oyster Density",
    y = expression(Abundance~(m^2)),
    color = "Season",
    title = "Smoothed Relationship Between Abundance and Oyster Density"
  ) +
  theme_minimal()


ggplot(rac, aes(x = oys.density, y = log(abm2), color = season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    color = "Season",
    title = "Linear Trend Between Log-Abundance and Oyster Density"
  ) +
  theme_minimal()

ggplot(rac, aes(x = oys.density)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    x = "Live Oyster Density",
    y = "Number of Samples",
    title = "Distribution of Live Oyster Density Samples"
  ) +
  theme_minimal()+
  facet_wrap(~season)

###! reef level abundance-----
mlod.ab.reef <- glmmTMB(log(abm2) ~ moys.density*season,data = rac%>% mutate(season = relevel(season,ref="w")))
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
mlod.spr.reef <- glmmTMB(sprm2 ~ moys.density*season ,data = rac%>% mutate(season = relevel(season,ref="w")))
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
# 1. Separate metadata and species matrix
meta <- comcomp2 %>%
  select(Site, Sample, season, oys.density, moys.density)

species <- comcomp2 %>%
  select(-(Site: moys.density))   # drops first 5 columns, keeps species only
# 1. Hellinger transform 
species <- decostand(species, method = "hellinger")

# 2. Run NMDS (Bray-Curtis)

set.seed(123)
nmds <- metaMDS(species, distance = "bray", k = 2, trymax = 200)

# Extract NMDS scores for plotting
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores <- bind_cols(meta, nmds_scores)

# 3. Fit species vectors (envfit)

fit <- envfit(nmds, species, permutations = 999)

# Extract significant vectors (p < 0.05)
vec <- as.data.frame(scores(fit, display = "vectors"))
vec$species <- rownames(vec)
# Calculate vector lengths 
vec$length <- sqrt(vec$NMDS1^2 + vec$NMDS2^2) 
# Select top 5 longest vectors 
vec_top5 <- vec %>% arrange(desc(length)) %>% slice(1:5)

# 4. Plot NMDS with ggplot2

ggplot(nmds_scores, aes(NMDS1, NMDS2)) + 
  geom_point(aes(shape = season, color = oys.density), size = 5, alpha = 1) + 
  scale_shape_manual(values = c( "f" = 16, "s" = 17, "u" = 15, "w" = 10 ))+
  stat_ellipse(aes(group = season), linetype = "dashed", color = "gray40") + 
  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 0.8) + 
  # Top 5 species vectors (green arrows) 
  geom_segment( 
    data = vec_top5, 
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
    arrow = arrow(length = unit(0.3, "cm")), 
    color = "forestgreen", linewidth = 0.9 ) + 
  # Species labels: larger, green, moved inward 
  geom_text( 
    data = vec_top5, 
    aes( 
      x = NMDS1 * 0.85+ (NMDS2 * 0.23), # pulls labels 15% closer to origin 
      y = NMDS2 * 0.85- (NMDS1 * 0.05),
      label = species ), 
    color = "forestgreen", 
    size = 5, # larger text 
    fontface = "bold" )+
  theme_minimal(base_size = 14) + 
  labs( title = "NMDS (Hellinger-transformed community)", 
        color = "Live Oyster Density", shape = "Season" )+
  facet_wrap(~season)


#Sample × Species × Abundance × oys.density
species_long <- comcomp2 %>%
  pivot_longer(
    cols = 6:47,              # all species columns
    names_to = "species",
    values_to = "abundance"
  )

ggplot(species_long, aes(x = abundance)) +
  geom_histogram(fill = "steelblue", color = "white") +
  facet_wrap(~ species, scales = "free") +
  theme_minimal() +
  labs(
    x = "Abundance",
    y = "Frequency",
    title = "Distribution of Abundance for Each Species"
  )

ggplot(species_long, aes(x = abundance, fill = oys.density)) +
  geom_histogram(color = "white", bins = 20) +
  facet_wrap(~ species, scales = "free") +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    x = "Abundance",
    y = "Frequency",
    fill = "Live Oyster Density",
    title = "Species Abundance Distributions by Live Oyster Density"
  )


species_long <- species_long %>%
  mutate(oys.bin = cut(oys.density,
                       breaks = c(-Inf, 0, 10, 50, Inf),
                       labels = c("0", "1–10", "11–75", "75+")))

ggplot(species_long, aes(x = abundance, fill = oys.bin)) +
  geom_histogram(color = "white", bins = 20, alpha = 0.7) +
  facet_wrap(~ species, scales = "free") +
  theme_minimal() +
  labs(
    x = "Abundance",
    y = "Frequency",
    fill = "Oyster Density Bin",
    title = "Species Abundance Distributions by Oyster Density Category"
  )

ggplot(species_long, aes(x = oys.density, y = abundance)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ species, scales = "free") +
  theme_minimal() +
  labs(
    x = "Live Oyster Density",
    y = "Abundance",
    title = "Abundance vs Live Oyster Density for Each Species"
  )

species_long2 <- species_long %>%
  mutate(oys.bin = cut(
    oys.density,
    breaks = c(-Inf, 0, 10, 50, Inf),
    labels = c("0", "1–10", "11–50", "50+")
  ))
ggplot(species_long2, aes(x = oys.bin, y = abundance)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ species, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Live Oyster Density (binned)",
    y = "Abundance",
    title = "How Abundance Changes with Live Oyster Density (Boxplots)"
  )

species_long <- species_long %>%
  mutate(oys.bin = cut(
    oys.density,
    breaks = pretty(oys.density, n = 2),
    include.lowest = TRUE
  ))




ggplot(species_long, aes(x = oys.bin, y = abundance)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ species, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Live Oyster Density (binned by distribution)",
    y = "Abundance",
    title = "How Abundance Changes with Live Oyster Density"
  )



species_long <- species_long %>%
  mutate(oys.bin = ntile(oys.density, 3))


# Compute quantile cutpoints
qs <- quantile(species_long$oys.density, probs = seq(0, 1, length.out = 4), na.rm = TRUE)

# Create pretty labels like "0–12", "12–40", "40–80"
pretty_labels <- paste0(
  round(qs[-length(qs)], 1), "–", round(qs[-1], 1)
)

# Apply labels
species_long <- species_long %>%
  mutate(oys.bin = factor(oys.bin, labels = pretty_labels))


ggplot(species_long, aes(x = oys.bin, y = abundance)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ species, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Live Oyster Density (3 equal‑sized bins)",
    y = "Abundance",
    title = "How Abundance Changes with Live Oyster Density"
  )


ggplot(species_long, aes(x = abundance)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_grid(species ~ oys.bin, scales = "free") +
  theme_minimal() +
  labs(
    x = "Abundance",
    y = "Frequency",
    title = "Species Abundance Distributions Across Oyster Density Bins",
    subtitle = "Oyster density bins are equal‑sized (quantile‑based) with pretty numeric labels"
  )

ggplot(species_long, aes(x = abundance, fill = oys.bin)) +
  geom_histogram(color = "white", bins = 20, alpha = 0.7) +
  facet_wrap(~ species, scales = "free") +
  theme_minimal() +
  labs(
    x = "Abundance",
    y = "Frequency",
    fill = "Oyster Density Bin",
    title = "Species Abundance Distributions by Oyster Density Category"
  )







###! reef level comp------
#nmds for abundance between seasons and mean live oyster density
# 4. Plot NMDS with ggplot2

ggplot(nmds_scores, aes(NMDS1, NMDS2)) + 
  geom_point(aes(shape = season, color = moys.density), size = 5, alpha = 1) + 
  scale_shape_manual(values = c( "f" = 16, "s" = 17, "u" = 15, "w" = 10 ))+
  stat_ellipse(aes(group = season), linetype = "dashed", color = "gray40") + 
  scale_color_viridis_c(option = "plasma", begin = 0.2, end = 0.8) + 
  # Top 5 species vectors (green arrows) 
  geom_segment( 
    data = vec_top5, 
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
    arrow = arrow(length = unit(0.3, "cm")), 
    color = "forestgreen", linewidth = 0.9 ) + 
  # Species labels: larger, green, moved inward 
  geom_text( 
    data = vec_top5, 
    aes( 
      x = NMDS1 * 0.85+ (NMDS2 * 0.23), # pulls labels 15% closer to origin 
      y = NMDS2 * 0.85- (NMDS1 * 0.05),
      label = species ), 
    color = "forestgreen", 
    size = 3, # larger text 
    fontface = "bold" )+
  theme_minimal(base_size = 14) + 
  labs( title = "NMDS (Hellinger-transformed community)", 
        color = "Mean Live Oyster Density", shape = "Season" )+
  facet_wrap(~season)



### site level comp-----













##### FFG ABUNDANCE ----
###! microhabitat FFG abundance----
###! c= carnivore-----
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

###! d= deposit feeders-----
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

###! o=omnivore-----
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
###! reef level FFG abundance-----
###! c= carnivore-----
carnivore_mod_reef <- glmmTMB(
  log(ffabm2 + 1) ~ moys.density * season,
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

###! d= deposit feeders-----
deposit_mod_reef <- glmmTMB(
  log(ffabm2 + 1) ~ moys.density * season,
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

###! o=omnivore-----
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
