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

rac<- reef2 %>%
  left_join(com2, by = c("Site", "Sample", "season"))%>%
  drop_na()%>% # I have 16 sites with no com data because I only processed 6/7 per site
  group_by(Site,season) %>%
  mutate(moys.density = mean(oys.density))

rac$season <- as.factor(rac$season)
rac$Site <- as.factor(rac$Site)
rac$Sample <- as.factor(rac$Sample)
any(is.na(rac))# final check for issues- looks good!






#com_noamps= commmunity metrics without amphipods and isopods :'(
com_noamps <- com1 %>%
  # Remove the unwanted TaxaIDs
  filter(!TaxaID %in% c("amp-1", "amp-11", "amp-4", "amp-7", "iso-1", "iso-4", "iso-8")) %>%
  group_by(Site, Sample, season) %>%
  summarize(
    sprm2 = n_distinct(TaxaID),
    abm2  = sum(Abundance.m2))

#rac_noamps= commmunity metrics and oyster data without amphipods and isopods :'(
rac_noamps <- reef2 %>%
  left_join(com_noamps, by = c("Site", "Sample", "season"))%>%
  drop_na()%>% # I have 16 sites with no com data because I only processed 6/7 per site
  group_by(Site,season) %>%
  mutate(moys.density = mean(oys.density))

rac_noamps$season <- as.factor(rac_noamps$season)
rac_noamps$Site <- as.factor(rac_noamps$Site)
rac_noamps$Sample <- as.factor(rac_noamps$Sample)
any(is.na(rac_noamps))# final check for issues- looks good!




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



#community metric data set with samples above 150 LOD removed
rac_filtered <- rac %>%
  filter(oys.density <= 150)

#function for residuals
glmm.resids<-function(model){
  t1 <- simulateResiduals(model)
  print(testDispersion(t1))
  plot(t1)
}


ggplot(rac, aes(x = Site, y = oys.density)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Oyster Density by Site",
       x = "Site",
       y = "Oyster Density") +
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


#oysters removed 150
ggplot(rac_filtered, aes(x = Site, y = oys.density)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Oyster Density by Site",
       x = "Site",
       y = "Oyster Density") +
  theme_minimal()

ggplot(rac_filtered, aes(x = oys.density)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    x = "Live Oyster Density",
    y = "Number of Samples",
    title = "Distribution of Live Oyster Density Samples"
  ) +
  theme_minimal()+
  facet_wrap(~season)



ggplot(rac_filtered, aes(x = oys.density)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white", alpha = 0.8) +
  labs(
    x = "Live Oyster Density",
    y = "Number of Samples",
    title = "Distribution of Samples Across Live Oyster Density filtered"
  ) +
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
                color = season)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    color = "Season",
    title = "Raw Log-Abundance by Oyster Density and season"
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

#explore data with removal of samples with LOD above 150---------

#distribution of abundance 
ggplot(rac_filtered, aes(x = abm2)) +
  geom_histogram(bins = 300, color = "black", fill = "grey70") +
  labs(
    x = "Abundance (m²)",
    y = "Count",
    title = "Distribution of Abundance"
  ) +
  theme_minimal()
#log abundance because heavy right tail skew
rac_filtered <- rac_filtered %>%
  mutate(log_abm2 = log(abm2))

ggplot(rac_filtered, aes(x = log_abm2)) +
  geom_histogram(bins = 300, color = "black", fill = "grey70") +
  labs(
    x = "LOGAbundance (m²)",
    y = "Count",
    title = "Distribution of Abundance"
  ) +
  theme_minimal()

ggplot(rac_filtered, aes(x = log_abm2)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 300,
    color = "black",
    fill = "grey70"
  ) +
  stat_function(
    fun = dnorm,
    args = list(
      mean = mean(rac_filtered$log_abm2, na.rm = TRUE),
      sd   = sd(rac_filtered$log_abm2, na.rm = TRUE)
    ),
    linewidth = 1
  ) +
  labs(
    x = "Log Abundance (m²)",
    y = "Density",
    title = "Distribution of Log Abundance with Normal Fit"
  )

    
#using rac_filtered

#just data points logged
ggplot(rac_filtered, aes(x = oys.density,
                y = log(abm2),
                color = season)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    color = "Season",
    title = "Raw Log-Abundance by Oyster Density and season  filtered"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


#smooth relationship
ggplot(rac_filtered, aes(x = oys.density, y = abm2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(y = abm2), method = "loess")+
  labs(
    x = "Live Oyster Density",
    y = expression(Abundance~(m^2)),
    color = "Season",
    title = "Smoothed Relationship Between Abundance and Oyster Density filtered"
  ) +
  theme_minimal()
#smooth relationship logged
ggplot(rac_filtered, aes(x = oys.density, y = log(abm2))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    title = "Smoothed Relationship Between Log Abundance and Oyster Density (filtered)"
  ) +
  theme_minimal()



#smooth relationship w season
ggplot(rac_filtered, aes(x = oys.density, y = abm2, color=season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(y = abm2), method = "loess")+
  labs(
    x = "Live Oyster Density",
    y = expression(Abundance~(m^2)),
    color = "Season",
    title = "Smoothed Relationship Between Abundance and Oyster Density  by season filtered"
  ) +
  theme_minimal()+
  facet_wrap(~season, scales="free_y")
#smooth relationship logged w season
ggplot(rac_filtered, aes(x = oys.density, y = log(abm2), color=season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    title = "Smoothed Relationship Between Log Abundance and Oyster Density (filtered)"
  ) +
  theme_minimal()+
  facet_wrap(~season, scales="free_y")


#linear logged
ggplot(rac_filtered, aes(x = oys.density, y = log(abm2), color = season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    color = "Season",
    title = "Linear Trend Between Log-Abundance and Oyster Density filtered"
  ) +
  theme_minimal()+
  facet_wrap(~season, scales="free_y")


#linear no log

ggplot(rac_filtered, aes(x = oys.density, y = abm2, color = season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Live Oyster Density",
    y = expression((Abundance~(m^2))),
    color = "Season",
    title = "Linear Trend Between Abundance and Oyster Density filtered"
  ) +
  theme_minimal()+
  facet_wrap(~season, scales="free_y")


#explore data with removal of amphipods and isopods-------

#distribution of abundance 
ggplot(rac_noamps, aes(x = abm2)) +
  geom_histogram(bins = 300, color = "black", fill = "grey70") +
  labs(
    x = "Abundance (m²)",
    y = "Count",
    title = "Distribution of Abundance"
  ) +
  theme_minimal()


ggplot(rac_noamps, aes(x = abm2)) +
  # Histogram scaled to density
  geom_histogram(aes(y = after_stat(density)),
                 bins = 300,
                 color = "black",
                 fill = "grey70") +
  # Normal distribution curve
  stat_function(
    fun = dnorm,
    args = list(
      mean = mean(rac_noamps$abm2, na.rm = TRUE),
      sd   = sd(rac_noamps$abm2, na.rm = TRUE)
    ),
    color = "red",
    linewidth = 1
  ) +
  labs(
    x = "Abundance (m²)",
    y = "Density",
    title = "Distribution of Abundance with Normal Fit"
  ) +
  theme_minimal()

#logged
ggplot(rac_noamps, aes(x = log(abm2))) +
  geom_histogram(aes(y = after_stat(density)), bins = 300, color = "black", fill = "grey70") +
  stat_function(fun = dnorm,
                args = list(mean = mean(log(rac_noamps$abm2), na.rm = TRUE),
                            sd   = sd(log(rac_noamps$abm2), na.rm = TRUE)),
                color = "red", linewidth = 1) +
  labs(x = "Log Abundance (m²)",
       y = "Density",
       title = "Distribution of Log Abundance with Normal Fit") +
  theme_minimal()



#using rac_noamps

#just data points logged
ggplot(rac_noamps, aes(x = oys.density,
                         y = log(abm2),
                         color = season)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    color = "Season",
    title = "Raw Log-Abundance by Oyster Density and season  filtered no amphipods/isopods"
  ) +
  theme_minimal() +
  theme(legend.position = "top")


#smooth relationship
ggplot(rac_noamps, aes(x = oys.density, y = abm2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(y = abm2), method = "loess")+
  labs(
    x = "Live Oyster Density",
    y = expression(Abundance~(m^2)),
    color = "Season",
    title = "Smoothed Relationship Between Abundance and Oyster Density filtered"
  ) +
  theme_minimal()
#smooth relationship logged
ggplot(rac_noamps, aes(x = oys.density, y = log(abm2))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    title = "Smoothed Relationship Between Log Abundance and Oyster Density (filtered)"
  ) +
  theme_minimal()



#smooth relationship w season
ggplot(rac_noamps, aes(x = oys.density, y = abm2, color=season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(y = abm2), method = "loess")+
  labs(
    x = "Live Oyster Density",
    y = expression(Abundance~(m^2)),
    color = "Season",
    title = "Smoothed Relationship Between Abundance and Oyster Density  by season filtered"
  ) +
  theme_minimal()
# +
#   facet_wrap(~season, scales="free_y")
#smooth relationship logged w season
ggplot(rac_noamps, aes(x = oys.density, y = log(abm2), color=season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    title = "Smoothed Relationship Between Log Abundance and Oyster Density (filtered)"
  ) +
  theme_minimal()+
  facet_wrap(~season, scales="free_y")


#linear logged
ggplot(rac_noamps, aes(x = oys.density, y = log(abm2), color = season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Live Oyster Density",
    y = expression(log(Abundance~(m^2))),
    color = "Season",
    title = "Linear Trend Between Log-Abundance and Oyster Density filtered"
  ) +
  theme_minimal()+
  facet_wrap(~season, scales="free_y")


#linear no log

ggplot(rac_noamps, aes(x = oys.density, y = abm2, color = season)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Live Oyster Density",
    y = expression((Abundance~(m^2))),
    color = "Season",
    title = "Linear Trend Between Abundance and Oyster Density filtered"
  ) +
  theme_minimal()+
  facet_wrap(~season, scales="free_y")
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









