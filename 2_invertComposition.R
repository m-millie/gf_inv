###################################################################
###
### 2_invertComposition.R : Importing and cleaning data for analysis.
###
### Authors: Millie Ortiz, Kimberly Komatsu
###
###################################################################

# Packages and Set-Up ----------------------------------------------------------------

library(lme4)
library(lmerTest)
library(codyn)
library(vegan)
library(tidyverse)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=20, vjust=2),
             strip.text.x = element_text(size=20), 
             strip.text.y = element_text(size=20),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  


# Read Data ---------------------------------------------------------------

abundance <- readRDS('abundance.RDS') %>% # invertebrate counts
  mutate(replicate=paste(burn_trt, watershed, block, plot, sep='::'))
biomass <- readRDS('biomassAll.RDS') # invertebrate biomass
CN <- readRDS('CN.RDS') # plant %C and %N
plant <- read.csv('inv_data/ghost_fire_plant_data.csv') # plant biomass and richness
functionalGroups <- read.csv('inv_data/gf_funct_groups.csv')

# Community Metrics ---------------------------------------------------------------

totalAbundance <- abundance %>%
  group_by(year, watershed, block, plot, burn_trt) %>%
  summarise(total_count = sum(count), .groups = 'drop')

communityStructure <- community_structure(abundance, time.var='year', abundance.var='count', replicate.var='replicate', metric='Evar') %>%
  separate(col=replicate, into=c('burn_trt','watershed','block','plot'), sep='::') %>% 
  mutate(plot=as.integer(plot)) %>% 
  left_join(totalAbundance)

functionalStructure <- abundance %>% 
  left_join(functionalGroups) %>% 
  group_by(year, watershed, block, plot, burn_trt, eco_functional_group) %>% 
  summarise(funct_count=sum(count), .groups='drop')


# Pre-Treatment (2014) Analysis -----------------------------------------------------------

# Composition

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

# Richness
richness2014 <- lmer(richness ~ burn_trt + (1 | watershed), data = subset(communityStructure, year==2014))
summary(richness2014)
anova(richness2014)

richnessFig2014 <- ggplot(communityStructure, aes(x = burn_trt, y = richness)) +
  geom_boxplot() +
  xlab("") +
  ylab("Invertebrate Richness")

# Evenness
evenness2014 <- lmer(Evar ~ burn_trt + (1 | watershed), data = subset(communityStructure, year==2014))
summary(evenness2014)
anova(evenness2014)

evennessFig2014 <- ggplot(communityStructure, aes(x = burn_trt, y = Evar)) +
  geom_boxplot() +
  xlab("") +
  ylab("Invertebrate Evenness")

# Abundance
count2014 <- lmer(total_count ~ burn_trt + (1 | watershed), data = subset(communityStructure, year==2014))
summary(count2014)
anova(count2014)

countFig2014 <- ggplot(communityStructure, aes(x = burn_trt, y = total_count)) +
  geom_boxplot() +
  xlab("") +
  ylab("Invertebrate Abundance")


# Biomass
biomass2014 <- lmer(invertebrate_biomass ~ burn_trt + (1 | watershed), data = subset(biomass, year==2014))
summary(biomass2014)
anova(biomass2014)

biomassFig2014 <- ggplot(subset(biomass, year==2014), aes(x = burn_trt, y = invertebrate_biomass)) +
  geom_boxplot() +
  xlab("") +
  ylab("Invertebrate Biomass (mg)")


# Functional Groups
functionalGroup2014 <- lmer(funct_count ~ burn_trt + (1|watershed), 
                            data = subset(functionalStructure, year==2014 & eco_functional_group=='herbivore'))
summary(functionalGroup2014)
anova(functionalGroup2014)

functionalGroup2014 <- lmer(funct_count ~ burn_trt + (1|watershed), 
                            data = subset(functionalStructure, year==2014 & eco_functional_group=='predator'))
summary(functionalGroup2014)
anova(functionalGroup2014)

functionalGroup2014 <- lmer(funct_count ~ burn_trt + (1|watershed), 
                            data = subset(functionalStructure, year==2014 & eco_functional_group=='parasitoid'))
summary(functionalGroup2014)
anova(functionalGroup2014)

functionalGroup2014 <- lmer(funct_count ~ burn_trt + (1|watershed), 
                            data = subset(functionalStructure, year==2014 & eco_functional_group=='omnivore'))
summary(functionalGroup2014)
anova(functionalGroup2014)






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
