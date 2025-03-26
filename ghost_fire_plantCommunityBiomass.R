################################################################################
##  ghost_fire_plantCommunityBiomass.R: Getting plant community diversity and biomass for each plot.
##
##  Authors: Kim Komatsu
##  Date created: March 26, 2025
################################################################################

library(codyn)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\konza projects\\Ghost Fire\\DATA')


##### import and bind plant community data #####

comm2014 <- read.csv('GhostFire2014_Data\\Species Comp\\GhostFire_SpComp_2014.csv') %>% 
  select(-Watershed, -Species)
comm2019 <- read.csv('GhostFire2019_Data\\SpeciesComp\\GhostFire_SpComp_2019.csv') %>% 
  select(-Comments)
comm2024 <- read.csv('GhostFire2024_Data\\SpeciesComp\\GhostFire_SpComp_2024.csv') %>% 
  select(-Watershed, -Comments)

commAll <- rbind(comm2014, comm2019, comm2024) %>% 
  pivot_longer(cols=c(June, August), names_to='month', values_to='cover') %>% 
  group_by(Year, Experiment, Site, Burn.Trt, Block, Plot, spnum) %>% 
  summarise(max_cover=max(cover)) %>% 
  ungroup()


##### calculate plant species richness #####

richnessAll <- commAll %>% 
  group_by(Year, Experiment, Site, Burn.Trt, Block, Plot) %>% 
  summarise(plant_richness=length(spnum)) %>% 
  ungroup() %>% 
  rename(year=Year,
         experiment=Experiment,
         site=Site,
         burn_trt=Burn.Trt,
         block=Block,
         plot=Plot)


##### import plant biomass data #####

bio2014 <- read.csv('GhostFire2014_Data\\Biomass\\GhostFire_Biomass_2014.csv') %>% 
  mutate(Year=2014) %>% 
  rename(burn_trt=BurnFreq) %>% 
  select(Year, burn_trt, Watershed, Block, Plot, Replicate, Grass, Forb, Woody, P.Dead) 
bio2019 <- read.csv('GhostFire2019_Data\\Biomass\\GhostFire_Biomass_DataEntry2019.csv') %>% 
  rename(Watershed=Wateshed) %>% 
  mutate(Year=2019,
         burn_trt=ifelse(BurnFreq==20, 'Annual', 'Unburned')) %>% 
  select(Year, burn_trt, Watershed, Block, Plot, Replicate, Grass, Forb, Woody, P.Dead) 
bio2024 <- read.csv('GhostFire2024_Data\\Biomass\\GhostFire_ANPP_2024.csv') %>% 
  rename(Replicate=Rep,
         P.Dead=Pdead) %>% 
  mutate(Year=2024,
         burn_trt=ifelse(BurnFreq==20, 'Annual', 'Unburned')) %>% 
  select(Year, burn_trt, Watershed, Block, Plot, Replicate, Grass, Forb, Woody, P.Dead) 


##### calculate average plot live and litter biomass #####
bioAll <- rbind(bio2014, bio2019, bio2024) %>% 
  mutate_at(c('Grass', 'Forb', 'Woody'), ~replace(., is.na(.), 0)) %>% 
  mutate(live_biomass=(Grass+Forb+Woody)) %>% 
  rename(litter_biomass=P.Dead) %>% 
  group_by(Year, burn_trt, Watershed, Block, Plot) %>% 
  summarise(live_biomass=mean(live_biomass)*100,
            litter_biomass=mean(litter_biomass)*100) %>% 
  ungroup() %>% 
  # select(-Watershed) %>% 
  rename(year=Year, 
         block=Block, 
         plot=Plot)


##### merge plant diversity and biomass data #####

plantData <- full_join(bioAll, richnessAll)

# write.csv(plantData, 'C:\\Users\\kjkomatsu\\Desktop\\R files\\gf_inv\\inv_data\\ghost_fire_plant_data.csv')