#' ---
#' 
#' 
#' ---


# Packages ----------------------------------------------------------------

library(tidyverse)
library(lme4)
library(vegan)


# Read Data ---------------------------------------------------------------

biomass_14 <- read.csv("inv_data/ghost_fire_invert_biomass_2014.csv")
comm_14 <- read.csv("inv_data/ghost_fire_invert_community_2014.csv")
comm_19 <- read.csv("inv_data/ghost_fire_invert_community_2019.csv")
comm_24 <- read.csv("inv_data/ghost_fire_invert_community_2024.csv")
biomass_24 <- read.csv("inv_data/ghost_fire_invert_biomass_2024.csv")
funct <- read.csv("inv_data/gf_funct_groups.csv")
anpp_24 <- read.csv("inv_data/GhostFire_ANPP_2024.csv")

# Organize WS -------------------------------------------------------------

comm_14$watershed <- factor(comm_14$watershed, levels = c("1D", "SpB", "20C", "20B"))
biomass_14$watershed <- factor(biomass_14$watershed, levels = c("1D", "SpB", "20C", "20B"))

comm_19$watershed <- factor(comm_19$watershed, levels = c("1D", "SpB", "20C", "20B"))

comm_24$watershed <- factor(comm_24$watershed, levels = c("1D", "SpB", "20C", "20B"))
biomass_24$watershed <- factor(biomass_24$watershed, levels = c("1D", "SpB", "20C", "20B"))



# ANOVAS and other analyses ------------------------------------------------------------------

## Count by burn_trt 
m5 <- aov(count ~ burn_trt, data = comm_24) 
summary(m5)

##Count by watershed 
m6 <- aov(count ~ watershed, data = comm_24)
summary(m6)

## Biomass by burn_trt 
m7 <- aov(biomass ~ burn_trt, data = biomass_24)
summary(m7)

## Biomass by watershed 
m8 <- aov(biomass ~ watershed, data = biomass_24)
summary(m8)

## Count by nutrients
m9 <- aov(count ~ nutrient, data = comm_24)
summary(m9)

## Count by litter
m10 <- aov(count ~ litter, data = comm_24)
summary(m10)


# Summary statistics for count by litter
litter_sum <- comm_24 %>%
  group_by(litter) %>%
  summarise(mean_count = mean(count), sd_count = sd(count), n = n())

# Summary statistics for count by nutrient
nutrient_sum <- comm_24 %>%
  group_by(nutrient) %>%
  summarise(mean_count = mean(count), sd_count = sd(count), n = n())


## Calculate richness 
 richness_24 <- comm_24 %>%
  group_by(burn_trt, watershed, plot) %>%
  summarise(richness = length(unique(arthropod_ID)))

# Calculate evenness
evenness_24 <- comm_24 %>%
  group_by(burn_trt, watershed, plot) %>%
  summarise(shannon = diversity(count, index = "shannon"))

# Combine richness and evenness into a single data frame
re_24 <- full_join(richness_24, evenness_24, by = c("burn_trt", "watershed"))

re_24$watershed <- factor(re_24$watershed, levels = c('1D', 'SpB', '20C', '20B'))


# Boxplots for richness by burn_trt and watershed
richness_burn_trt_plot <- ggplot(re_24, aes(x = burn_trt, y = richness, fill = burn_trt)) +
  geom_boxplot() +
  theme_bw() +
  theme_classic() +
  labs(x = "Burn Treatment", y = "Richness")

richness_watershed_plot <- ggplot(re_24, aes(x = watershed, y = richness, fill = watershed)) +
  geom_boxplot() +
  theme_bw() +
  theme_classic() +
  labs(x = "Watershed", y = "Richness")

# Boxplots for evenness by burn_trt and watershed
evenness_burn_trt_plot <- ggplot(re_24, aes(x = burn_trt, y = shannon, fill = burn_trt)) +
  geom_boxplot() +
  theme_bw() +
  theme_classic() +
  labs(x = "Burn Treatment", y = "Evenness")

evenness_watershed_plot <- ggplot(re_24, aes(x = watershed, y = shannon, fill = watershed)) +
  geom_boxplot() +
  theme_bw() +
  theme_classic() +
  labs(x = "Watershed", y = "Evenness")

# Combine plots into one
gridExtra::grid.arrange(
  richness_burn_trt_plot, richness_watershed_plot, 
  evenness_burn_trt_plot, evenness_watershed_plot, nrow = 2)

#  ANOVA for richness by burn
anova_richness_burn_trt <- aov(richness ~ burn_trt, data = re_24)
summary(anova_richness_burn_trt)

# ANOVA for richness by watershed
anova_richness_watershed <- aov(richness ~ watershed, data = re_24)
summary(anova_richness_watershed)

# ANOVA for evenness by burn_trt
anova_evenness_burn_trt <- aov(shannon ~ burn_trt, data = re_24)
summary(anova_evenness_burn_trt)

#  ANOVA for evenness by watershed
anova_evenness_watershed <- aov(shannon ~ watershed, data = re_24)
summary(anova_evenness_watershed)



# Figures? ----------------------------------------------------------------

