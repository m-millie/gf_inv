#' ---
#' 
#' 
#' ---


# Packages ----------------------------------------------------------------

library(tidyverse)
library(gridExtra)
library(lme4)
library(vegan)
library(lmerTest)
library(piecewiseSEM)


# Read Data ---------------------------------------------------------------

acrididae_14 <- read.csv("inv_data/ghost_fire_invert_biomass_2014.csv") %>% 
  filter(contents == "Acrididae") %>% 
  group_by(year, burn_trt) %>% 
  summarise(acrididae_biomass = mean(biomass)) %>% 
  ungroup()

acrididae_obs_14 <- read.csv("inv_data/ghost_fire_invert_community_2014.csv") %>% 
  filter(collected == "observed", order == "Orthoptera") %>% 
  left_join(acrididae_14) %>% 
  mutate(biomass_2 = count * acrididae_biomass) %>% 
  select(year, month, watershed, block, plot, burn_trt, biomass_2)

biomass_14 <- read.csv("inv_data/ghost_fire_invert_biomass_2014.csv") %>% 
  group_by(year, month, watershed, block, plot, burn_trt) %>% 
  summarise(biomass_2 = sum(biomass)) %>% 
  ungroup() %>% 
  rbind(acrididae_obs_14) %>% 
  group_by(year, month, watershed, block, plot, burn_trt) %>% 
  summarise(total_biomass = sum(biomass_2)) %>% 
  ungroup()


acrididae_24 <- read.csv("inv_data/ghost_fire_invert_biomass_2024.csv") %>% 
  filter(contents == "orthoptera") %>% 
  group_by(year, burn_trt, plot_trt, litter_trt) %>% 
  summarise(acrididae_biomass = mean(biomass)) %>% 
  ungroup()

acrididae_obs_24 <- read.csv("inv_data/ghost_fire_invert_community_2024.csv") %>% 
  filter(collected == "Observed", order == "Orthoptera") %>% 
  left_join(acrididae_24) %>% 
  mutate(biomass_2 = count * acrididae_biomass) %>% 
  select(year, month, watershed, block, plot, burn_trt, biomass_2, plot_trt, litter_trt)

biomass_24 <- read.csv("inv_data/ghost_fire_invert_biomass_2024.csv") %>% 
  group_by(year, month, watershed, block, plot, burn_trt, plot_trt, litter_trt) %>% 
  summarise(biomass_2 = sum(biomass)) %>% 
  ungroup() %>% 
  rbind(acrididae_obs_24) %>% 
  group_by(year, month, watershed, block, plot, burn_trt, plot_trt, litter_trt) %>% 
  summarise(total_biomass = sum(biomass_2)) %>% 
  ungroup()


comm_14 <- read.csv("inv_data/ghost_fire_invert_community_2014.csv")
comm_19 <- read.csv("inv_data/ghost_fire_invert_community_2019.csv")
comm_24 <- read.csv("inv_data/ghost_fire_invert_community_2024.csv")
funct <- read.csv("inv_data/gf_funct_groups.csv")
plant_data <- read.csv("inv_data/ghost_fire_plant_data.csv") %>% 
  select(-burn_trt)
soil <- read.csv("inv_data/ghost_fire_soil_moisture_2024.csv") %>% 
  group_by(watershed, block, plot) %>% 
  summarise(mean_soil = mean(soil_moisture)) %>% 
  ungroup()


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=20, vjust=2),
             strip.text.x = element_text(size=20), 
             strip.text.y = element_text(size=20),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

# Organize WS -------------------------------------------------------------

comm_14$watershed <- factor(comm_14$watershed, levels = c("1D", "SpB", "20C", "20B"))
biomass_14$watershed <- factor(biomass_14$watershed, levels = c("1D", "SpB", "20C", "20B"))

comm_19$watershed <- factor(comm_19$watershed, levels = c("1D", "SpB", "20C", "20B"))

comm_24$watershed <- factor(comm_24$watershed, levels = c("1D", "SpB", "20C", "20B"))
biomass_24$watershed <- factor(biomass_24$watershed, levels = c("1D", "SpB", "20C", "20B"))



# 2014 Analysis -----------------------------------------------------------

abun_14 <- comm_14 %>%
  group_by(burn_trt, plot, watershed, block) %>%
  summarise(total_abun = sum(count), .groups = 'drop')

## mixed model for abundance
abun_mod_14 <- lmer(total_abun ~ burn_trt + (1 | watershed), data = abun_14)


summary(abun_mod_14)

anova(abun_mod_14)

count_14<- ggplot(abun_14, aes(x = burn_trt, y = total_abun)) +
  geom_boxplot() +
  xlab("") +
  ylab("Total Abundance")
  

## mixed model for biomass
bio_mod_14 <- lmer(total_biomass ~ burn_trt + (1 | watershed), data = biomass_14)

summary(bio_mod_14)

anova(bio_mod_14)

bio_14<- ggplot(biomass_14, aes(x = burn_trt, y = total_biomass)) +
  geom_boxplot() +
  xlab("") +
  ylab("Total Biomass (g)")

## model for total family richness
# calculate total family richness by plot
family_richness_14 <- comm_14 %>%
  group_by(plot, burn_trt, watershed, block) %>%
  summarise(total_family_richness = n_distinct(family)) %>% 
  ungroup()

# merge family richness back to the original df
#comm_14 <- comm_14 %>%
 # left_join(family_richness_14, by = "plot")

# Fit a mixed-effects model with burn as a fixed effect and watershed as a random effect
richness_mod_14 <- lmer(total_family_richness ~ as.factor(burn_trt) + (1 | watershed), data = family_richness_14)

summary(richness_mod_14)

anova(richness_mod_14) #figure out df for richness

## model for family evenness

family_evenness_14 <- comm_14 %>% 
  group_by(plot, burn_trt, watershed, block) %>% 
  summarise(family_evenness = diversity(count, index = "shannon") / log(specnumber(count))) %>% 
  ungroup()

# Merge family evenness back to the original dataframe 
#comm_14 <- comm_14 %>% 
# left_join(family_evenness_14, by = "plot")

# Fit a mixed-effects model for family evenness 
evenness_mod_14 <- lmer(family_evenness ~ burn_trt + (1 | watershed), data = family_evenness_14)

# Display the model summary for evenness 
summary(evenness_mod_14)

anova(evenness_mod_14)


rich_14 <- ggplot(family_richness_14, aes(x = burn_trt, y = total_family_richness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Family Richness")

even_14 <- ggplot(family_evenness_14, aes(x = burn_trt, y = family_evenness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Family Evenness")

grid.arrange(count_14, bio_14, rich_14, even_14, nrow = 2)


## merged functional groups
merged_fun_14 <- merge(comm_14, funct, by.x = c("family", "order"), by.y = c("family", "order")) %>% 
  group_by(year, month, watershed, block, plot, burn_trt, functional_group) %>% 
  summarise(total_count = sum(count)) %>% 
  ungroup()

herb_mod_14 <- lmer(total_count ~ burn_trt + (1|watershed),
                   data = merged_fun_14[merged_fun_14$functional_group == "Herbivore", ])

anova(herb_mod_14)

pred_mod_14 <- lmer(total_count ~ burn_trt + (1|watershed),
                    data = merged_fun_14[merged_fun_14$functional_group == "Predator", ])

anova(pred_mod_14)

para_mod_14 <- lmer(total_count ~ burn_trt + (1|watershed),
                    data = merged_fun_14[merged_fun_14$functional_group == "Parasitoid", ])

anova(para_mod_14)

omni_mod_14 <- lmer(total_count ~ burn_trt + (1|watershed),
                    data = merged_fun_14[merged_fun_14$functional_group == "Omnivore", ])

anova(omni_mod_14)


# 2019 Analysis -----------------------------------------------------------

abun_19 <- comm_19 %>% #figure out why 20C-B-4 have an extra space after nutrient
  #mutate(plot_trt2 = ifelse(plot_trt == "C ", "C", plot_trt)) %>% 
  #separate(plot_trt, into = c("plot_trt", "drop"),sep = " ")
  #select(-plot_trt) %>% 
  #rename(plot_trt = plot_trt2) %>% 
  group_by(year, month, burn_trt, watershed, block, plot, plot_trt, litter_trt) %>% 
  summarise(total_abun = sum(count)) %>% 
  ungroup()

## abundance mixed model
abun_mod_19 <- lmer(total_abun ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = abun_19)

summary(abun_mod_19)

anova(abun_mod_19)

## family richness model
df_richness_19 <- comm_19 %>%
  ungroup() %>% 
  #mutate(block_plot = paste(block, plot, sep = "::")) %>% 
  group_by(watershed, block, plot, burn_trt, plot_trt, litter_trt, year) %>%
  summarise(family_richness = n_distinct(family))

richness_mod_19 <- lmer(family_richness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_richness_19)

summary(richness_mod_19)

anova(richness_mod_19)

## family evenness model

df_evenness_19 <- comm_19 %>%
  group_by(block, plot, burn_trt, plot_trt, litter_trt, watershed) %>%
  summarise(family_evenness = diversity(count, index = "shannon") / log(specnumber(count)))


# Fit the mixed-effects model
evenness_mod_19 <- lmer(family_evenness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_evenness_19)

summary(evenness_mod_19)

anova(evenness_mod_19)

count_19 <- ggplot(abun_19, aes(x = plot_trt, y = total_abun, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Total Abundance") +
  scale_fill_manual(values = c("#337539", "#dccd7d")) +
  facet_wrap(~burn_trt)

rich_19 <- ggplot(df_richness_19, aes(x = plot_trt, y = family_richness, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Family Richness") +
  scale_fill_manual(values = c("#337539", "#dccd7d")) +
  facet_wrap(~burn_trt)

even_19 <- ggplot(df_evenness_19, aes(x = plot_trt, y = family_evenness, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Family Evenness") +
  scale_fill_manual(values = c("#337539", "#dccd7d")) +
  facet_wrap(~burn_trt)

grid.arrange(arrangeGrob(count_19),
             arrangeGrob(rich_19, even_19, ncol = 1),
             ncol = 2, widths = c(2,2))



merged_funct_19 <- merge(comm_19, funct, by.x = c("family", "order"), by.y = c("family", "order"))

merged_fun_19 <- merge(comm_19, funct, by.x = c("family", "order"), by.y = c("family", "order")) %>% 
  group_by(year, month, watershed, block, plot, burn_trt, plot_trt, litter_trt, functional_group) %>% 
  summarise(total_count = sum(count)) %>% 
  ungroup()

herb_mod_19 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                    data = merged_fun_19[merged_fun_19$functional_group == "Herbivore", ])

anova(herb_mod_19)


pred_mod_19 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                    data = merged_fun_19[merged_fun_19$functional_group == "Predator", ])

anova(pred_mod_19)

para_mod_19 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                    data = merged_fun_19[merged_fun_19$functional_group == "Parasitoid", ])

anova(para_mod_19)

det_mod_19 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                    data = merged_fun_19[merged_fun_19$functional_group == "Detritivore", ])

anova(det_mod_19)

par_mod_19 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                   data = merged_fun_19[merged_fun_19$functional_group == "Parasite", ])

anova(par_mod_19)

pol_mod_19 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                   data = merged_fun_19[merged_fun_19$functional_group == "Pollinator", ])

anova(pol_mod_19)


ggplot(merged_fun_19 %>% filter(functional_group == "Predator"), 
       aes(x = burn_trt, y = total_count)) +
  geom_boxplot() +
  xlab("") +
  ylab("Predator Abundance")
  #scale_fill_manual(values = c("#337539", "#dccd7d")) +
  #facet_wrap(~burn_trt)

herb_nut <- ggplot(merged_fun_19 %>% filter(functional_group == "Herbivore"), 
       aes(x = plot_trt, y = total_count)) +
  geom_boxplot() +
  xlab("") +
  ylab("Herbivore Abundance")

herb_lit <- ggplot(merged_fun_19 %>% filter(functional_group == "Herbivore"), 
       aes(x = litter_trt, y = total_count)) +
  geom_boxplot() +
  xlab("") +
  ylab("Herbivore Abundance")


grid.arrange(herb_nut, herb_lit, nrow = 1)

# 2024 analysis ------------------------------------------------------
abun_24 <- comm_24 %>% 
  group_by(year, month, burn_trt, watershed, block, plot, plot_trt, litter_trt) %>% 
  summarise(total_abun = sum(count)) %>% 
  ungroup()

## abundance mixed model
abun_mod_24 <- lmer(total_abun ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = abun_24)

summary(abun_mod_24)

anova(abun_mod_24)

## abundance mixed model
bio_mod_24 <- lmer(total_biomass ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = biomass_24)

summary(bio_mod_24)

anova(bio_mod_24)

## family richness model
df_richness_24 <- comm_24 %>%
  group_by(plot, burn_trt, plot_trt, litter_trt, watershed, block, year) %>%
  summarise(family_richness = n_distinct(family))

richness_mod_24 <- lmer(family_richness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_richness_24)

summary(richness_mod_24)

anova(richness_mod_24)

## family evenness model

df_evenness_24 <- comm_24 %>%
  group_by(plot, burn_trt, plot_trt, litter_trt, watershed, block) %>%
  summarise(family_evenness = diversity(count, index = "shannon") / log(specnumber(count)))

# Fit the mixed-effects model
evenness_mod_24 <- lmer(family_evenness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_evenness_24)

summary(evenness_mod_24)

anova(evenness_mod_24)

count_24 <- ggplot(abun_24, aes(x = plot_trt, y = total_abun, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Total Abundance") +
  scale_fill_manual(values = c("#337539", "#dccd7d")) +
  facet_wrap(~burn_trt)

bio_24 <- ggplot(biomass_24, aes(x = plot_trt, y = total_biomass, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Total Biomass (g)") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_fill_manual(values = c("#337539", "#dccd7d")) +
  facet_wrap(~burn_trt)

rich_24 <- ggplot(df_richness_24, aes(x = plot_trt, y = family_richness, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Family Richness") +
  scale_fill_manual(values = c("#337539", "#dccd7d")) +
  facet_wrap(~burn_trt)

even_24 <- ggplot(df_evenness_24, aes(x = plot_trt, y = family_evenness, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Family Evenness") +
  scale_fill_manual(values = c)("#337539", "#dccd7d") +
  facet_wrap(~burn_trt)

grid.arrange(count_24, bio_24, rich_24, even_24, nrow = 2)
 

merged_fun_24 <- merge(comm_24, funct, by.x = c("family", "order"), by.y = c("family", "order"))

merged_fun_24 <- merge(comm_24, funct, by.x = c("family", "order"), by.y = c("family", "order")) %>% 
  group_by(year, month, watershed, block, plot, burn_trt, plot_trt, litter_trt, functional_group) %>% 
  summarise(total_count = sum(count)) %>% 
  ungroup()

herb_mod_24 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                    data = merged_fun_24[merged_fun_24$functional_group == "Herbivore", ])

anova(herb_mod_24)

pred_mod_24 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                    data = merged_fun_24[merged_fun_24$functional_group == "Predator", ])

anova(pred_mod_24)

omni_mod_24 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                    data = merged_fun_24[merged_fun_24$functional_group == "Omnivore", ])

anova(omni_mod_24)

para_mod_24 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                   data = merged_fun_24[merged_fun_24$functional_group == "Parasitoid", ])

anova(para_mod_24)

pol_mod_24 <- lmer(total_count ~ burn_trt * plot_trt * litter_trt + (1 | watershed),
                   data = merged_fun_24[merged_fun_24$functional_group == "Polyphagous", ])

anova(pol_mod_24)

ggplot(merged_fun_24 %>% filter(functional_group == "Omnivore"),
       aes(x = plot_trt, y = total_count, fill = litter_trt)) +
  geom_boxplot() +
  xlab("") +
  ylab("Omnivore Abundance") +
  scale_fill_manual(values = c("#337539", "#dccd7d")) +
  facet_wrap(~burn_trt)


# regressions -----------------------------------------------
## treatment plant regressions
abun_trt <- rbind(abun_19, abun_24) %>% 
  full_join(plant_data) %>% 
  na.omit()

summary(lm(total_abun ~ plant_richness, data = subset(abun_trt, year %in% c(2019, 2024))))


rich_trt <- ggplot(abun_trt %>% filter(year != 2014),
       aes(x = plant_richness, y = total_abun)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Plant Richness") +
  ylab("Arthropod Abundance")

summary(lm(total_abun ~ live_biomass, data = subset(abun_trt, year %in% c(2019, 2024))))


live_trt <- ggplot(abun_trt %>% filter(year != 2014),
       aes(x = live_biomass, y = total_abun)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Live Biomass") +
  ylab("Arthropod Abundance")


summary(lm(total_abun ~ litter_biomass, data = subset(abun_trt, year %in% c(2019, 2024))))


litter_trt <- ggplot(abun_trt %>% filter(year != 2014), 
       aes(x = litter_biomass, y = total_abun)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +  
  geom_smooth(method = "lm",  se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Litter Biomass") +
  ylab("Arthropod Abundance")


## pretreatment regressions
abun_pre <- abun_14 %>% 
  full_join(plant_data)


summary(lm(total_abun ~ plant_richness, data = subset(abun_pre, year %in% c(2014))))


rich_pre <- ggplot(abun_pre %>% filter(!year %in% c(2019, 2024)), 
       aes(x = plant_richness, y = total_abun)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  xlab("Plant Richness") +
  ylab("Arthropod Abundance")


summary(lm(total_abun ~ live_biomass, data = subset(abun_pre, year %in% c(2014))))


live_pre <- ggplot(abun_pre %>% filter(!year %in% c(2019, 2024)), 
       aes(x = live_biomass, y = total_abun)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  xlab("Live Biomass") +
  ylab("Arthropod Abundance")


summary(lm(total_abun ~ litter_biomass, data = subset(abun_pre, year %in% c(2014))))


litter_pre <- ggplot(abun_pre %>% filter(!year %in% c(2019, 2024)), 
       aes(x = litter_biomass, y = total_abun)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  xlab("Litter Biomass") +
  ylab("Arthropod Abundance")

grid.arrange(arrangeGrob(rich_pre, live_pre, litter_pre, ncol = 1),
             arrangeGrob(rich_trt, live_trt, litter_trt, ncol = 1),
             ncol = 2, widths = c(2,2))

## soil moisture regression

abun_soil <- abun_24 %>% 
  full_join(soil)

summary(lm(total_abun ~ mean_soil, data = abun_soil))


sm_trt <- ggplot(abun_soil,
       aes(x = mean_soil, y = total_abun)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  xlab("Soil Moisture") +
  ylab("Arthropod Abundance")


biomass_soil <- biomass_24 %>% 
  left_join(soil)

summary(lm(total_biomass ~ mean_soil, data = biomass_soil))

ggplot(biomass_soil,
       aes(x = mean_soil, y = total_biomass)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Soil Moisture") +
  ylab("Arthropod Biomass")


r_soil <- df_richness_24 %>% 
  left_join(soil)

summary(lm(family_richness ~ mean_soil, data = r_soil))

ggplot(r_soil,
       aes(x = mean_soil, y = family_richness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Soil Moisture") +
  ylab("Arthropod Family Richness")


e_soil <- df_evenness_24 %>% 
  left_join(soil)

summary(lm(family_evenness ~ mean_soil, data = e_soil))

ggplot(e_soil,
       aes(x = mean_soil, y = family_evenness)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Soil Moisture") +
  ylab("Arthropod Family Evenness")


sem_data <- abun_trt %>% 
  mutate(plot_trt_num = ifelse(plot_trt == "Carbon", -1, ifelse(plot_trt == "Control", 0, 1)), 
         litter_trt_num = ifelse(litter_trt == "Absent", 0, 1),
         burn_trt_num = ifelse(burn_trt == "Annual", 0, 1))

## richness regressions
summary(lm(total_abun ~ plant_richness, data = subset(abun_trt, year %in% c(2019, 2024))))


rich_trt <- ggplot(abun_trt %>% filter(year != 2014),
                   aes(x = plant_richness, y = total_abun)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Plant Richness") +
  ylab("Arthropod Abundance")

summary(lm(total_abun ~ live_biomass, data = subset(abun_trt, year %in% c(2019, 2024))))


live_trt <- ggplot(abun_trt %>% filter(year != 2014),
                   aes(x = live_biomass, y = total_abun)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Live Biomass") +
  ylab("Arthropod Abundance")


summary(lm(total_abun ~ litter_biomass, data = subset(abun_trt, year %in% c(2019, 2024))))


litter_trt <- ggplot(abun_trt %>% filter(year != 2014), 
                     aes(x = litter_biomass, y = total_abun)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +  
  geom_smooth(method = "lm",  se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Litter Biomass") +
  ylab("Arthropod Abundance")


###all data, all years-----------

###difference metrics not through composition model
#all years

richness_trt <- rbind(df_richness_19, df_richness_24) %>% 
  left_join(plant_data) %>% 
  mutate(plot_trt_num = ifelse(plot_trt == "Carbon", -1, ifelse(plot_trt == "Control", 0, 1)), 
         litter_trt_num = ifelse(litter_trt == "Absent", 0, 1),
         burn_trt_num = ifelse(burn_trt == "Annual", 0, 1)) %>% 
  ungroup()

summary(div_sem <- psem(
  lm(family_richness ~ burn_trt_num + litter_trt_num + litter_biomass + live_biomass + plant_richness, data = richness_trt),
  lm(litter_biomass ~ burn_trt_num + litter_trt_num + plot_trt_num, data = richness_trt),
  lm(live_biomass ~ burn_trt_num + litter_trt_num + plot_trt_num, data = richness_trt),
  lm(plant_richness ~ burn_trt_num + litter_trt_num + plot_trt_num, data = richness_trt),
  plant_richness %~~% live_biomass,
  plant_richness %~~% litter_biomass,
  live_biomass %~~% litter_biomass,
  data = richness_trt 
))

sem_data <- abun_trt %>% 
  mutate(plot_trt_num = ifelse(plot_trt == "Carbon", -1, ifelse(plot_trt == "Control", 0, 1)), 
         litter_trt_num = ifelse(litter_trt == "Absent", 0, 1),
         burn_trt_num = ifelse(burn_trt == "Annual", 0, 1))

summary(abundance_sem <- psem(
  lm(total_abun ~ litter_biomass + live_biomass + plant_richness + litter_trt_num + burn_trt_num, data = sem_data),
  lm(litter_biomass ~ burn_trt_num + litter_trt_num + plot_trt_num, data = sem_data),
  lm(live_biomass ~ burn_trt_num + litter_trt_num + plot_trt_num, data = sem_data),
  lm(plant_richness ~ burn_trt_num + litter_trt_num + plot_trt_num, data = sem_data),
  plant_richness %~~% live_biomass,
  plant_richness %~~% litter_biomass,
  live_biomass %~~% litter_biomass,
  data = sem_data 
))

coefs1 <- coefs(div_sem, standardize = "scale", standardize.type = "latent.linear", intercepts = FALSE)


summary(lm(family_richness ~ plant_richness, data = subset(richness_trt, year %in% c(2019, 2024))))


rich_rich <- ggplot(richness_trt %>% filter(year != 2014),
                   aes(x = plant_richness, y = family_richness)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +
  #geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Plant Richness") +
  ylab("Arthropod Richness")

summary(lm(family_richness ~ live_biomass, data = subset(richness_trt, year %in% c(2019, 2024))))


live_rich <- ggplot(richness_trt %>% filter(year != 2014),
                   aes(x = live_biomass, y = family_richness)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Live Biomass") +
  ylab("Arthropod Richness")


summary(lm(family_richness ~ litter_biomass, data = subset(richness_trt, year %in% c(2019, 2024))))


litter_rich <- ggplot(richness_trt %>% filter(year != 2014), 
                     aes(x = litter_biomass, y = family_richness)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +  
  #geom_smooth(method = "lm",  se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Litter Biomass") +
  ylab("Arthropod Richness")



richness_pre <- family_richness_14 %>% 
  left_join(plant_data)


summary(lm(total_family_richness ~ plant_richness, data = subset(richness_pre, year %in% c(2019, 2024))))


rich_rich_pre <- ggplot(richness_pre %>% filter(!year %in% c(2019, 2024)),
                    aes(x = plant_richness, y = total_family_richness)) +
  geom_point() +
  #geom_smooth(method = "lm", se = F, color = "black") +
  xlab("Plant Richness") +
  ylab("Arthropod Richness")

summary(lm(total_family_richness ~ live_biomass, data = subset(richness_pre, year %in% c(2019, 2024))))


live_rich_pre <- ggplot(richness_pre %>% filter(!year %in% c(2019, 2024)),
                        aes(x = live_biomass, y = total_family_richness)) +
  geom_point() +
  #geom_smooth(method = "lm", se = F, color = "black") +
  xlab("Live Biomass") +
  ylab("Arthropod Richness")


summary(lm(total_family_richness ~ litter_biomass, data = subset(richness_pre, year %in% c(2019, 2024))))


litter_rich_pre <- ggplot(richness_pre %>% filter(!year %in% c(2019, 2024)),
                          aes(x = litter_biomass, y = total_family_richness)) +
  geom_point() +  
  #geom_smooth(method = "lm",  se = F, color = "black") +
  xlab("Litter Biomass") +
  ylab("Arthropod Richness")


grid.arrange(arrangeGrob(rich_rich_pre, live_rich_pre, litter_rich_pre, ncol = 1),
             arrangeGrob(rich_rich, live_rich, litter_rich, ncol = 1),
             ncol = 2, widths = c(2,2))



# permanova 2024---------------------------------------------------------------

fam_abun <- comm_24 %>% 
  group_by(year, watershed, block, plot, burn_trt, litter_trt, plot_trt, arthropod_ID) %>% 
  summarise(total_count = sum(count)) %>% 
  ungroup() %>% 
  mutate(trt = paste(burn_trt, litter_trt, plot_trt, sep = "_")) %>% #you have to make your dataframe wide form for this
  select(year, watershed, block, plot, burn_trt, litter_trt, plot_trt, arthropod_ID, total_count, trt) %>% #you want some replicate variable, treatment variable, and your taxonomic identifier and count columns
  pivot_wider(names_from='arthropod_ID', values_from = 'total_count', values_fill = 0)  
  #pivot_wider so that species are the column names and the counts are filled in, with 0's put in if a species wasn't found in a plot

permanova <- adonis(formula = fam_abun[,9:112] ~ litter_trt * plot_trt * burn_trt, data=fam_abun, permutations=999, method="bray") #this runs the PERMANOVA test on the relCover2021 data with only the columns related to the species as the response variable, the trt as the dependent variable, 999 permutations of the test using bray curtis dissimilarity as your distance metric

print(permanova) #print the permanova output

results_table <- as.data.frame(permanova$aov.tab)


#all the code below is for plotting the NMDS (a non-metric dimensional scaling plot) that shows differences between treatments in terms of their community composition
sppBC <- metaMDS(fam_abun[,9:112])

plotData <- fam_abun[,1:8]

#Use the vegan ellipse function to make ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group= fam_abun$trt)
BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)               

ord3 <- data.frame(plotData,scores(sppBC,display="sites"))%>%
  group_by(trt)

BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)
BC_Ellipses <- data.frame() #Make a new empty data frame called BC_Ellipses  
for(g in unique(BC_NMDS$group)){
  BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                             veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
                                          ,group=g))
} #Generate ellipses points
BC_Ellipses2 <- BC_Ellipses %>% 
  separate(col = group, into= c("burn_trt", "litter_trt", "plot_trt"), sep = "_", remove = F)

nmds1 <- ggplot(subset(BC_NMDS_Graph, burn_trt  = "Annual"), aes(x=MDS1, y=MDS2, color=plot_trt,linetype = litter_trt)) +
  geom_point(size=6)+ 
  geom_path(data = filter(BC_Ellipses2, group%in%c("Annual_Absent_Carbon","Annual_Absent_Control", "Annual_Absent_Nitrogen", "Annual_Present_Carbon", "Annual_Present_Control", "Annual_Present_Nitrogen")), aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("#CC79A7", "#D55E00", "#009E73"), name = "") +
  #scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 22))
  

nmds2 <- ggplot(subset(BC_NMDS_Graph, burn_trt  = "Unburned"), aes(x=MDS1, y=MDS2, color=plot_trt,linetype = litter_trt)) +
  geom_point(size=6)+ 
  geom_path(data = filter(BC_Ellipses2, group%in%c("Unburned_Absent_Carbon","Unburned_Absent_Control", "Unburned_Absent_Nitrogen", "Unburned_Present_Carbon", "Unburned_Present_Control", "Unburned_Present_Nitrogen")), aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("#CC79A7", "#D55E00", "#009E73"), name = "") +
  #scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 22))


grid.arrange(nmds1, nmds2, nrow = 1)



# permanova 2019 ----------------------------------------------------------

fam_abun_19 <- comm_19 %>% 
  group_by(year, watershed, block, plot, burn_trt, litter_trt, plot_trt, arthropod_ID) %>% 
  summarise(total_count = sum(count)) %>% 
  ungroup() %>% 
  mutate(trt = paste(burn_trt, litter_trt, plot_trt, sep = "_")) %>% #you have to make your dataframe wide form for this
  select(year, watershed, block, plot, burn_trt, litter_trt, plot_trt, arthropod_ID, total_count, trt) %>% #you want some replicate variable, treatment variable, and your taxonomic identifier and count columns
  pivot_wider(names_from='arthropod_ID', values_from = 'total_count', values_fill = 0)  
#pivot_wider so that species are the column names and the counts are filled in, with 0's put in if a species wasn't found in a plot

permanova <- adonis(formula = fam_abun_19[,9:133] ~ litter_trt * plot_trt * burn_trt, data=fam_abun, permutations=999, method="bray") #this runs the PERMANOVA test on the relCover2021 data with only the columns related to the species as the response variable, the trt as the dependent variable, 999 permutations of the test using bray curtis dissimilarity as your distance metric

print(permanova) #print the permanova output

results_table <- as.data.frame(permanova$aov.tab)


#all the code below is for plotting the NMDS (a non-metric dimensional scaling plot) that shows differences between treatments in terms of their community composition
sppBC <- metaMDS(fam_abun_19[,9:133])

plotData <- fam_abun_19[,1:8]

#Use the vegan ellipse function to make ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group= fam_abun_19$trt)
BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)               

ord3 <- data.frame(plotData,scores(sppBC,display="sites"))%>%
  group_by(trt)

BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)
BC_Ellipses <- data.frame() #Make a new empty data frame called BC_Ellipses  
for(g in unique(BC_NMDS$group)){
  BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                             veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
                                          ,group=g))
} #Generate ellipses points
BC_Ellipses2 <- BC_Ellipses %>% 
  separate(col = group, into= c("burn_trt", "litter_trt", "plot_trt"), sep = "_", remove = F)

nmds1 <- ggplot(subset(BC_NMDS_Graph, burn_trt  = "Annual"), aes(x=MDS1, y=MDS2, color=plot_trt,linetype = litter_trt)) +
  geom_point(size=6)+ 
  geom_path(data = filter(BC_Ellipses2, group%in%c("Annual_Absent_Carbon","Annual_Absent_Control", "Annual_Absent_Nitrogen", "Annual_Present_Carbon", "Annual_Present_Control", "Annual_Present_Nitrogen")), aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("#CC79A7", "#D55E00", "#009E73"), name = "") +
  #scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 22))


nmds2 <- ggplot(subset(BC_NMDS_Graph, burn_trt  = "Unburned"), aes(x=MDS1, y=MDS2, color=plot_trt,linetype = litter_trt)) +
  geom_point(size=6)+ 
  geom_path(data = filter(BC_Ellipses2, group%in%c("Unburned_Absent_Carbon","Unburned_Absent_Control", "Unburned_Absent_Nitrogen", "Unburned_Present_Carbon", "Unburned_Present_Control", "Unburned_Present_Nitrogen")), aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("#CC79A7", "#D55E00", "#009E73"), name = "") +
  #scale_linetype_manual(values = c("twodash", "solid", "twodash", "solid", "twodash", "solid"), name = "") +
  xlab("NMDS1")+ 
  ylab("NMDS2")+ 
  theme(axis.text.x=element_text(size=24, color = "black"), axis.text.y = element_text(size = 24, color = "black"), legend.text = element_text(size = 22))


grid.arrange(nmds1, nmds2, nrow = 1)



# permanova 2014 ----------------------------------------------------------

fam_abun_14 <- comm_14 %>% 
  group_by(year, watershed, block, plot, burn_trt, arthropod_ID) %>% 
  summarise(total_count = sum(count)) %>% 
  ungroup() %>% 
  #mutate(trt = paste(burn_trt, litter_trt, plot_trt, sep = "_")) %>% #you have to make your dataframe wide form for this
  select(year, watershed, block, plot, burn_trt, arthropod_ID, total_count) %>% #you want some replicate variable, treatment variable, and your taxonomic identifier and count columns
  pivot_wider(names_from='arthropod_ID', values_from = 'total_count', values_fill = 0)  
#pivot_wider so that species are the column names and the counts are filled in, with 0's put in if a species wasn't found in a plot

permanova <- adonis(formula = fam_abun_14[,6:47] ~ burn_trt, data=fam_abun, permutations=999, method="bray") #this runs the PERMANOVA test on the relCover2021 data with only the columns related to the species as the response variable, the trt as the dependent variable, 999 permutations of the test using bray curtis dissimilarity as your distance metric

print(permanova) #print the permanova output

results_table <- as.data.frame(permanova$aov.tab)


#all the code below is for plotting the NMDS (a non-metric dimensional scaling plot) that shows differences between treatments in terms of their community composition
sppBC <- metaMDS(fam_abun_14[,6:47])

plotData <- fam_abun_14[,1:5]

#Use the vegan ellipse function to make ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

BC_NMDS = data.frame(MDS1 = sppBC$points[,1], MDS2 = sppBC$points[,2],group= fam_abun_14$burn_trt)
BC_NMDS_Graph <- cbind(plotData,BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$burn_trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)               

ord3 <- data.frame(plotData,scores(sppBC,display="sites"))%>%
  group_by(burn_trt)

BC_Ord_Ellipses<-ordiellipse(sppBC, plotData$burn_trt, display = "sites",
                             kind = "se", conf = 0.95, label = T)
BC_Ellipses <- data.frame() #Make a new empty data frame called BC_Ellipses  
for(g in unique(BC_NMDS$group)){
  BC_Ellipses <- rbind(BC_Ellipses, cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                             veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,BC_Ord_Ellipses[[g]]$center,BC_Ord_Ellipses[[g]]$scale)))
                                          ,group=g))
} #Generate ellipses points

ggplot(subset(BC_NMDS_Graph), aes(x=MDS1, y=MDS2)) +
  geom_point(size=6, aes(color=burn_trt)) +  # Color points by burn_trt
  geom_path(data = filter(BC_Ellipses), 
            aes(x = NMDS1, y = NMDS2, color = group),  # Color ellipses by burn_trt
            size = 3) +
  labs(color="Burn Treatment", linetype = "", shape = "") +
  scale_color_manual(values=c("#de1a24", "#056517")) +  # Custom colors for treatments
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme(axis.text.x = element_text(size=24, color = "black"),
        axis.text.y = element_text(size = 24, color = "black"),
        legend.text = element_text(size = 22))


