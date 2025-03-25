#' ---
#' 
#' 
#' ---


# Packages ----------------------------------------------------------------

library(tidyverse)
library(lme4)
library(vegan)


# Read Data ---------------------------------------------------------------

biomass_14 <- read.csv("inv_data/ghost_fire_invert_biomass_2014.csv") %>% 
  group_by(year, month, watershed, block, plot, burn_trt) %>% 
  summarise(total_biomass = sum(biomass)) %>% 
  ungroup()
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



# 2014 Analysis -----------------------------------------------------------
gml_14 <- lmer(count ~ burn_trt + (1 | watershed), data = comm_14)
glm_14_results <- Anova(model, type = "III")


df_richness_evenness_14 <- comm_14 %>%
  group_by(plot, burn_trt, watershed) %>%
  summarise(
    richness = n_distinct(family),
    evenness = diversity(count) / log(richness)
  )

# Fit the mixed effects model for richness
model_richness_14 <- lmer(richness ~ burn_trt + (1 | watershed), data = df_richness_evenness_14)

# Perform ANOVA for richness
anova_richness_14 <- Anova(model_richness_14, type = "III")

# Fit the mixed effects model for evenness
model_evenness_14 <- lmer(evenness ~ burn_trt + (1 | watershed), data = df_richness_evenness_14)

# Perform ANOVA for evenness
anova_evenness_14 <- Anova(model_evenness_14, type = "III")

# Display the results
list(richness_anova_14 = anova_richness_14, evenness_anova_14 = anova_evenness_14)


df_richness_evenness <- comm_14 %>% 
  group_by(plot, burn_trt, watershed) %>% 
  summarise( richness = n_distinct(family), 
             evenness = diversity(count) / log(n_distinct(family)) )

plot_data <- df_richness_evenness_14 %>%
  pivot_longer(cols = c(richness, evenness), names_to = "metric", values_to = "value")

ggplot(plot_data, aes(x = burn_trt, y = value, fill = burn_trt)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Family Richness and Evenness by Burn Treatment",
       x = "Burn Treatment",
       y = "Value") +
  theme(legend.position = "none")


abun_14 <- comm_14 %>%
  group_by(burn_trt, plot, watershed, block) %>%
  summarise(total_count = sum(count), .groups = 'drop')

## mixed model for abundance
abun_mod_14 <- lmer(total_count ~ burn_trt + (1 | watershed), data = abun_14)


summary(abun_mod_14)

anova(abun_mod_14)



ggplot(abun_14, aes(x = burn_trt, y = total_count, fill = burn_trt)) +
  geom_boxplot() +
  labs(title = "Total Invertebrate Count by Burn Treatment",
       x = "Burn Treatment",
       y = "Number of Individuals") +
  theme_minimal() +
  theme(legend.position = "none")

## mixed model for biomass
bio_mod_14 <- lmer(total_biomass ~ burn_trt + (1 | watershed), data = biomass_14)

summary(bio_mod_14)

anova(bio_mod_14)

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
richness_mod_14 <- lmer(total_family_richness ~ burn_trt + (1 | watershed), data = family_richness_14)

summary(richness_mod_14)

anova(richness_mod_14) #figure out df for richness

## model for family evenness
family_evenness_14 <- comm_14 %>% 
  group_by(plot, burn_trt, watershed, block) %>% 
  summarise(family_evenness = sum(count^2) / (sum(count)^2))

# Merge family evenness back to the original dataframe 
#comm_14 <- comm_14 %>% 
# left_join(family_evenness_14, by = "plot")

# Fit a mixed-effects model for family evenness 
evenness_mod_14 <- lmer(family_evenness ~ burn_trt + (1 | watershed), data = family_evenness_14)

# Display the model summary for evenness 
summary(evenness_mod_14)

anova(evenness_mod_14)


ggplot(family_richness_14, aes(x = burn_trt, y = total_family_richness)) +
  geom_boxplot() +
  labs(title = "Boxplot of Family Richness by Burn Treatment") +
  theme_minimal()

ggplot(family_evenness_14, aes(x = burn_trt, y = family_evenness)) +
  geom_boxplot() +
  labs(title = "Boxplot of Family Evenness by Burn Treatment") +
  theme_minimal()


# 2019 Analysis -----------------------------------------------------------

abun_19 <- comm_19 %>% 
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
  group_by(watershed, as.factor(block), plot, burn_trt, plot_trt, litter_trt) %>%
  summarise(family_richness = n_distinct(family))

richness_mod_19 <- lmer(family_richness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_richness_19)

# Display the summary of the model
summary(richness_mod_19)

## family evenness model
df_evenness_19 <- comm_19 %>%
  group_by(plot, burn_trt, plot_trt, litter_trt, watershed) %>%
  summarise(family_evenness = diversity(count, index = "shannon") / log(n_distinct(family)))

# Fit the mixed-effects model
evenness_mod_19 <- lmer(family_evenness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_evenness_19)

# Display the summary of the model
summary(evenness_mod_19)



# new 2024 analysis? ------------------------------------------------------

## abundance mixed model
abun_mod_24 <- lmer(count ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = comm_24)

summary(abun_mod_24)

## abundance mixed model
bio_mod_24 <- lmer(biomass ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = biomass_24)

summary(bio_mod_24)


## family richness model
df_richness_24 <- comm_24 %>%
  group_by(plot, burn_trt, plot_trt, litter_trt, watershed) %>%
  summarise(family_richness = n_distinct(family))

richness_mod_24 <- lmer(family_richness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_richness_24)

# Display the summary of the model
summary(richness_mod_24)

## family evenness model
df_evenness_24 <- comm_24 %>%
  group_by(plot, burn_trt, plot_trt, litter_trt, watershed) %>%
  summarise(family_evenness = diversity(count, index = "shannon") / log(n_distinct(family)))

# Fit the mixed-effects model
evenness_mod_24 <- lmer(family_evenness ~ burn_trt * plot_trt * litter_trt + (1 | watershed), data = df_evenness_24)

# Display the summary of the model
summary(evenness_mod_24)









# 2024 Analysis -----------------------------------------------------------


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


richness_data <- comm_24 %>%
  group_by(litter, plot) %>%
  summarise(family_richness = n_distinct(family)) %>%
  ungroup()

# Perform ANOVA to analyze family richness based on litter
anova_result <- aov(family_richness ~ litter, data = richness_data)

# Display ANOVA results
summary(anova_result)

richness_dat <- comm_24 %>%
  group_by(nutrient, plot) %>%
  summarise(family_richness = n_distinct(family)) %>%
  ungroup()

# Perform ANOVA to analyze family richness based on nutrient
anova_res <- aov(family_richness ~ nutrient, data = richness_dat)

# Display ANOVA results
summary(anova_res)

# Create a boxplot to visualize family richness by nutrient
ggplot(richness_dat, aes(x = nutrient, y = family_richness, fill = nutrient)) +
  geom_boxplot() +
  labs(title = "Family Richness by Nutrient",
       x = "Nutrient",
       y = "Family Richness") +
  theme_minimal()


evenness_24 <- comm_24 %>%
  group_by(plot, litter, nutrient) %>%
  summarise(family_count = n_distinct(family),
            total_count = sum(count)) %>%
  mutate(evenness = family_count / total_count)

# Perform ANOVA to analyze family evenness based on litter
anova_litter <- aov(evenness ~ litter, data = evenness_24)

# Perform ANOVA to analyze family evenness based on nutrient
anova_nutrient <- aov(evenness ~ nutrient, data = evenness_24)

# Display ANOVA summaries
summary(anova_litter)
summary(anova_nutrient)




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


richness_data <- comm_24 %>%
  group_by(litter, plot) %>%
  summarise(family_richness = n_distinct(family)) %>%
  ungroup()

# Perform ANOVA to analyze family richness based on litter
anova_result <- aov(family_richness ~ litter, data = richness_data)

# Display ANOVA results
summary(anova_result)

richness_dat <- comm_24 %>%
  group_by(nutrient, plot) %>%
  summarise(family_richness = n_distinct(family)) %>%
  ungroup()

# Perform ANOVA to analyze family richness based on nutrient
anova_res <- aov(family_richness ~ nutrient, data = richness_dat)

# Display ANOVA results
summary(anova_res)

# Create a boxplot to visualize family richness by nutrient
ggplot(richness_dat, aes(x = nutrient, y = family_richness, fill = nutrient)) +
  geom_boxplot() +
  labs(title = "Family Richness by Nutrient",
       x = "Nutrient",
       y = "Family Richness") +
  theme_minimal()


evenness_24 <- comm_24 %>%
  group_by(plot, litter, nutrient) %>%
  summarise(family_count = n_distinct(family),
            total_count = sum(count)) %>%
  mutate(evenness = family_count / total_count)

# Perform ANOVA to analyze family evenness based on litter
anova_litter <- aov(evenness ~ litter, data = evenness_24)

# Perform ANOVA to analyze family evenness based on nutrient
anova_nutrient <- aov(evenness ~ nutrient, data = evenness_24)

# Display ANOVA summaries
summary(anova_litter)
summary(anova_nutrient)




## 2014 data mixed model
gml_14 <- lmer(count ~ burn_trt + (1 | watershed), data = comm_14)
glm_14_results <- Anova(model, type = "III")


df_richness_evenness_14 <- comm_14 %>%
  group_by(plot, burn_trt, watershed) %>%
  summarise(
    richness = n_distinct(family),
    evenness = diversity(count) / log(richness)
  )

# Fit the mixed effects model for richness
model_richness_14 <- lmer(richness ~ burn_trt + (1 | watershed), data = df_richness_evenness_14)

# Perform ANOVA for richness
anova_richness_14 <- Anova(model_richness_14, type = "III")

# Fit the mixed effects model for evenness
model_evenness_14 <- lmer(evenness ~ burn_trt + (1 | watershed), data = df_richness_evenness_14)

# Perform ANOVA for evenness
anova_evenness_14 <- Anova(model_evenness_14, type = "III")

# Display the results
list(richness_anova_14 = anova_richness_14, evenness_anova_14 = anova_evenness_14)


df_richness_evenness <- comm_14 %>% 
  group_by(plot, burn_trt, watershed) %>% 
  summarise( richness = n_distinct(family), 
             evenness = diversity(count) / log(n_distinct(family)) )

plot_data <- df_richness_evenness_14 %>%
  pivot_longer(cols = c(richness, evenness), names_to = "metric", values_to = "value")

ggplot(plot_data, aes(x = burn_trt, y = value, fill = burn_trt)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Family Richness and Evenness by Burn Treatment",
       x = "Burn Treatment",
       y = "Value") +
  theme(legend.position = "none")


## Same thing but to get f values? 
# Fit the mixed effects model for richness 
model_richness <- lmer(richness ~ burn_trt + (1 | watershed), data = df_richness_evenness)

# Perform ANOVA for richness 
anova_richness <- anova(model_richness)

# Fit the mixed effects model for evenness 
model_evenness <- lmer(evenness ~ burn_trt + (1 | watershed), data = df_richness_evenness)

# Perform ANOVA for evenness 
anova_evenness <- anova(model_evenness)

# Extract F values 
f_values_richness <- anova_richness$`F value` 
f_values_evenness <- anova_evenness$`F value`

# Display the results with F values 
list( richness_anova = anova_richness, richness_f_values = f_values_richness, evenness_anova = anova_evenness, evenness_f_values = f_values_evenness )




# Figures? ----------------------------------------------------------------
