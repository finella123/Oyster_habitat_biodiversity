# This is my Master's thesis project- exploring the effects of oyster trophic interactions on the benthic community
#This script is to reorganize my data 
source("scripts2/01_installpackagesM.R")
source("scripts2/02_import_data-EX.R")

#####Calculate sample abundance and dry weight from the sub sample amphipod and isopod data-----
#all data from subsampling came from tray methods
#there were no vouchers taken 
#remove unneeded columns, get true dry weight, and get rid of negative values
allaiwd <- allaisub %>%
  select(-Project, -Date,-Mesh.Size,-Method)%>%
  filter(Site!="M")%>%
  mutate(dw=Dryw-WeighBoat, dw= ifelse(dw<0.0001, 0.001, dw))%>%
  select(-Dryw,-WeighBoat, -Wetw)%>%
  rename(Sample=Bag.ID)

#Isolate the aibulks from the subsamples-
allaibulk <- allaiwd%>%
  filter(Taxa.ID %in% "ai-bulk")%>%
  mutate(aib.dw=dw)%>%
  select(-Taxa.ID,-Abundance,-dw, -Voucher)

#Isolate amphipod and isopod subsample of 50
aisub <- allaiwd%>%
  filter(!Taxa.ID %in% "ai-bulk")

#Attain biomass for individuals and attain proportions by subsample
#add a new column in aisub that takes the dry weight of each taxa and divides it by the abundance to get an individual weight
aisub$Abundance <- as.numeric(aisub$Abundance)
aisub2 <- aisub%>%
  mutate(ind.dw=dw/Abundance)%>%
  group_by(Site, Sample, season) %>%
  mutate(total_abundance = sum(Abundance, na.rm = TRUE))%>%
  mutate(proportion=Abundance/total_abundance)

#Use the proportions and individual dry weights to back calculate the abundance and biomass of the taxa for the bulk subsamples
subcalc <- aisub2%>%
  left_join(allaibulk)%>%
  mutate(pxb=proportion*aib.dw)%>% # proportion of biomass in bulk sample by taxa
  mutate(b.abund=pxb/ind.dw)%>% # abundance of each taxa in the bulk sample
  mutate(Total.Abundance=round(Abundance+b.abund))%>% #total abundance of each taxa by sample
  mutate(Total.Biomass=dw+pxb)%>% #total biomass of each taxa by sample
  select(-Abundance, -dw, -ind.dw, -total_abundance, -proportion, -aib.dw, -pxb, -b.abund,-Voucher)%>%
  rename(Abundance=Total.Abundance)%>%
  rename(true.dw=Total.Biomass)%>%
  rename(TaxaID=Taxa.ID)

#####Create working data for all community abundance and biomass data-----
#true dry weight of the samples and remove negative values. use wet weights when dry weights are not available. Add 1 to abundance and 1 average dry weight per taxa for each row that the voucher specimen in Y. Remove site 1 and suction samples
biomass$Abundance<- as.numeric(biomass$Abundance)
wd1 <- biomass%>%
  mutate(ww=Wetw-WeighBoat, ww= ifelse(ww<0.0001, 0.001, ww))%>%
  mutate(dw=Dryw-WeighBoat, dw= ifelse(dw<0.0001, 0.001, dw))%>%
  filter(Site!=1)%>%
  filter(Method!='s')%>%
  select(-Project, -Date, -WeighBoat,-Wetw, -Dryw,-Method,-mesh.size)%>%
  rename(Sample=BagID)%>%
  mutate(Abundance=case_when(
    Voucher=="Y"~Abundance+1,
    Voucher=="N"~Abundance))

#calculate average dry weight for each taxa to add to biomass after 
adw<-wd1%>%
  group_by(TaxaID)%>%
  mutate(dryweight=sum(dw,na.rm = TRUE),
         abund=sum(Abundance,na.rm = TRUE))%>%
  mutate(adw=dryweight/abund, adw=ifelse(adw==0,ww,adw))%>%
  select(TaxaID,adw)%>%
  distinct()

#combine average dw with dw for voucher specimens because they did not get their dry weights recorded
wd2<-wd1 %>%
  group_by(season,Site, Sample,TaxaID)%>%
  left_join(adw)%>%
  mutate(true.dw=case_when(
    Voucher=="Y"~dw+adw,
    Voucher=="N"~dw,
    TRUE~adw), true.dw=ifelse(is.na(true.dw),adw,true.dw))%>%
  select(-Voucher,-ww,-dw,-adw)%>%
  bind_rows(subcalc)%>%
  summarize(Abundance=sum(Abundance),
            Biomass=sum(true.dw))%>%
  mutate(Abundance.m2=Abundance/0.2304)%>% #standardized to 1m^2
  mutate(Biomass.m2=Biomass/0.2304) #standardized to 1m^2

#Change taxa ID's to the most updated and confirmed taxa ID
upid <- taxaid%>%
  select(-TaxaID2,-notes,-scientific,-common)

wd3 <- wd2%>%
  left_join(upid)%>%
  select(-TaxaID)%>%
  rename(TaxaID=upID)%>%
  group_by(Site, Sample, TaxaID,season)%>%
  summarise(Abundance=sum(Abundance),
            Biomass=sum(Biomass),
            Abundance.m2=sum(Abundance.m2),
            Biomass.m2=sum(Biomass.m2))
wd3 <- wd3[-1163, ] # remove WTF-7 since it will not be used in this dataset
any(is.na(wd3))# final check for outliers


#this is where I need to get all my taxa ids and make a data sheet with the updated id's 
unique(wd3$TaxaID)
# Get unique values from a column (example: "column_name" in data frame "df")
unique_values <- unique(wd3$TaxaID)

# Turn into a data frame
unique_df <- data.frame(unique_values)

# Save to CSV
write.csv(unique_df, "unique_values.csv", row.names = FALSE)
#LEFT OFF HERE



write_rds(wd3,"wdata2/wd_community_data_m^2.rds")

##### Create working data for reef characteristics-----
# 1 shell height= ">115" because it maxed out the caliper.make it 116 
shellheight[133, 8] = 116
#Remove spaces in data and make shell height numeric
shellheight$shell.height <- as.numeric(gsub(" ", "", shellheight$shell.height))
write_rds(shellheight,"wdata2/wd_reef_shell_height.rds") 
#Calculate average cluster vol, cluster count, total substrate, oyster count, mussel count, and barnacle count
reef1 <- reefquads%>%
  mutate(clust.v=cluster.volume-starting.volume,
         t.sub.v=total.reef.substrate.volume-starting.volume)%>%
  select(-month,-day,-year,-project,-starting.volume,-cluster.volume,-total.reef.substrate.volume)%>%
  rename(sample=quadrat)%>%
  group_by(site,sample, season)%>%
  mutate(clust.density=(cluster.count)/.0625, #standardize m^-2
         oys.density=(oyster.count)/.0625,#standardize m^-2
         mus.density=(mussel.count)/.0625,#standardize m^-2
         barn.density=(barnacle.count)/.0625,#standardize m^-2
         clust.vol=(clust.v)/.0625,#standardize m^-2
         reef.vol=(t.sub.v)/.0625)%>%#standardize m^-2
  select(-cluster.count, -oyster.count, -mussel.count,-barnacle.count,-clust.v,-t.sub.v)
write_rds(reef1,"wdata2/wd_reef_characteristics_m^2.rds") 

##### Sandra paired ind shell height and volume oyster-----
#join all data sheets and pull out needed columns 
pair.sh.v <- paired.sh.v%>%
  rename(shell.height=Height..mm., volume= Volume.Oyster..mL.)%>% select(-Whole.Oys.weight..g.)
#Remove 1 outlier due to human error during data entry and no way to confirm.
pair.sh.v <- pair.sh.v[-c(537), ]
write_rds(pair.sh.v,"wdata2/wd_sh_ov_predictor_data.rds") 
