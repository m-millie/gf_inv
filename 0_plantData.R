################################################################################
##  0_plantData.R: Getting plant community diversity and biomass for each plot.
##
##  Authors: Kim Komatsu
##  Date created: March 26, 2025
################################################################################

library(codyn)
library(tidyverse)


##### import plant community data and calculate plant species richness #####

richnessAll <- read.csv('https://pasta.lternet.edu/package/data/eml/knb-lter-knz/101/4/17ec1f61e5234d4931cd4c395f4bd643') %>% 
  filter(RecYear %in% c(2014, 2019, 2024)) %>% 
  mutate(Watershed=ifelse(Watershed=='SPB', 'SpB', Watershed)) %>% 
  group_by(RecYear, BurnTrt, Watershed, Block, Plot) %>% 
  summarise(plant_richness=length(Spnum)) %>% 
  ungroup() %>% 
  rename(year=RecYear,
         burn_trt=BurnTrt,
         watershed=Watershed,
         block=Block,
         plot=Plot)


##### import plant biomass data and calculate average plot live and litter biomass #####
bioAll <- read.csv('https://pasta.lternet.edu/package/data/eml/knb-lter-knz/101/4/053e369d68886f36f9f15c175749c59f') %>% 
  filter(RecYear %in% c(2014, 2019, 2024)) %>% 
  mutate(burn_trt=ifelse(BurnFreq==1, 'Annual', 'Unburned')) %>% 
  mutate_at(c('Grass', 'Forb', 'Woody'), ~replace(., is.na(.), 0)) %>% 
  mutate(live_biomass=(Grass+Forb+Woody)) %>% 
  rename(litter_biomass=PreDead) %>% 
  group_by(RecYear, burn_trt, Watershed, Block, Plot) %>% 
  summarise(live_biomass=mean(live_biomass)*100,
            litter_biomass=mean(litter_biomass)*100) %>% 
  ungroup() %>% 
  rename(year=RecYear, 
         block=Block, 
         plot=Plot, 
         watershed=Watershed)


##### merge plant diversity and biomass data #####

plantData <- full_join(bioAll, richnessAll)

# saveRDS(plantData, 'plantData.RDS')