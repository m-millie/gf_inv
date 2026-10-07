###################################################################
###
### 0_dataCleaning.R : Importing and cleaning data for upload to EDI.
###
### Authors: Millie Ortiz, Kimberly Komatsu
###
###################################################################

# Packages and Set-Up ----------------------------------------------------------------

library(readxl)
library(tidyverse)


# Read Data ---------------------------------------------------------------

### abundance data

abundance2014 <- read.csv("inv_data/ghost_fire_invert_community_2014.csv") %>% 
  select(year, month, watershed, block, plot, burn_trt, order, family, arthropod_ID, stage, collected, count)

abundance2019 <- read.csv("inv_data/ghost_fire_invert_community_2019.csv") %>% 
  select(year, month, watershed, block, plot, burn_trt, order, family, arthropod_ID, stage, collected, count)

abundance2024 <- read.csv("inv_data/ghost_fire_invert_community_2024.csv") %>% 
  select(year, month, watershed, block, plot, burn_trt, order, family, arthropod_ID, stage, collected, count)

abundance <- rbind(abundance2014, abundance2019, abundance2024) %>% 
  mutate(order=str_to_lower(order),
         family=str_to_lower(family),
         collected=str_to_lower(collected), 
         burn_trt=str_to_sentence(burn_trt)) %>% 
  mutate(order=ifelse(order %in% c('hemiptera:heteroptera', 'hemiptera:auchenorrhyncha', 'hemiptera:sternorrhyncha'), 'hemiptera', 
               ifelse(family %in% c('scathophagidae', 'culicidae', 'psilidae', 'chloropidae', 'pipunculidae'), 'diptera', 
               ifelse(family=='chrysomelidae', 'coleoptera', 
               ifelse(family %in% c('encrytidae', 'formicidae'), 'hymenoptera', 
               ifelse(family %in% c('blissidae', 'lygaeidae', 'tingidae'), 'hemiptera', order))))),
         family=ifelse(family=='tettigoniidae', 'tettigonidae', family)) %>% 
  select(-arthropod_ID)


### biomass data (only 2014 and 2024; missing samples from 2019 prevent analysis)

biomass2014 <- read.csv("inv_data/ghost_fire_invert_biomass_2014.csv") %>% 
  mutate(tube_mass=1000*tube_mass,
         combined_mass=1000*combined_mass,
         biomass=combined_mass-tube_mass) %>% 
  select(year, month, watershed, block, plot, burn_trt, contents, biomass) %>% 
  mutate(contents=ifelse(contents %in% c('Acrididae', 'Tettigonidae'), 'orthoptera', 
                  ifelse(contents=='Heteronemiidae', 'other', contents))) #change to more accurately reflect contents

biomass2024 <- read.csv("inv_data/ghost_fire_invert_biomass_2024.csv") %>% 
  mutate(tube_mass=ifelse(tube_mass==113.87, 1013.87, tube_mass),
         biomass=as.numeric(combined_mass)-tube_mass) %>% 
  select(year, month, watershed, block, plot, burn_trt, contents, biomass) %>% 
  mutate(contents=ifelse(contents=='oher', 'other', contents)) #fix spelling error

biomass <- rbind(biomass2014, biomass2024) %>% 
  mutate(burn_trt=str_to_sentence(burn_trt))

ggplot(biomass, aes(x=biomass, fill=as.factor(year))) + geom_histogram()

# Orthoptera only, generating estimates of biomass of observed but not collected individuals (hopped away too fast)

orthopteraCollected <- abundance %>% 
  filter(order=='orthoptera', collected=='collected') %>% 
  select(year, watershed, block, plot, order, family, count) %>% 
  group_by(year, watershed, block, plot, order) %>% 
  summarise(count=sum(count), .groups='drop')

orthopteraBiomass <- biomass %>% 
  filter(contents=="orthoptera") %>% 
  left_join(orthopteraCollected) %>%
  mutate(count=ifelse(is.na(count), 1, count)) %>% 
  group_by(year, burn_trt) %>% 
  summarise(orthoptera_biomass = mean(biomass/count)) %>% 
  ungroup()

observedBiomass <- abundance %>% 
  filter(collected=='observed',
         order=='orthoptera',
         year!=2019) %>% 
  left_join(orthopteraBiomass) %>% 
  mutate(biomass=count*orthoptera_biomass) %>% 
  select(year, month, watershed, block, plot, burn_trt, biomass)

biomassAll <- biomass %>% 
  select(year, month, watershed, block, plot, burn_trt, biomass) %>% 
  rbind(observedBiomass) %>% 
  group_by(year, month, watershed, block, plot, burn_trt) %>% 
  summarise(invertebrate_biomass=sum(biomass), .groups='drop')


### plant %C and %N data

CN2019 <- read_xlsx('inv_data/2019_ghost_fire_CN.xlsx') %>% 
  mutate(year=2019)

CN2024 <- read_xlsx('inv_data/2024_ghost_fire_CN.xlsx') %>% 
  mutate(year=2024)

CN <- rbind(CN2019, CN2024) %>% 
  rename(watershed=Watershed,
         block=Block,
         plot=Plot,
         growth_form=Type) %>% 
  mutate(burn_trt=ifelse(watershed %in% c('1D', 'SpB'), 'Annual', 'Unburned'))


# Write Data for EDI Project ---------------------------------------------------------------

write.csv(abundance, 'inv_data/GF_invertAbundance.csv', row.names=F)
write.csv(biomassAll, 'inv_data/GF_invertBiomass.csv', row.names=F)
write.csv(CN, 'inv_data/GF_plantCN.csv', row.names=F)