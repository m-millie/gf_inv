###################################################################
###
### 2_invertComposition.R : Data analysis and figure generation 
###                         related to invertebrate responses to 
###                         burning, soil N, and plant litter.
###
### Authors: Millie Ortiz, Kimberly Komatsu
###
###################################################################

# Packages and Set-Up ----------------------------------------------------------------

library(EDIutils)
library(lme4)
library(lmerTest)
library(emmeans)
library(codyn)
library(vegan)
library(cowplot)
library(piecewiseSEM)
library(tidyverse)

myEDIAccessKey="PASTE-KEY-HERE"

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

trt <- read.csv('GF_PlotList.csv') %>% 
  mutate(plot=as.integer(str_extract(Plot, "\\d+"))) %>% 
  select(-plot_id, -Plot, -Burn.Trt2) %>% 
  rename(watershed=Watershed,
         burn_trt=Burn.Trt,
         block=Block,
         litter=Litter,
         nutrient=Nutrient) %>% 
  mutate(nutrient= recode(nutrient,
         'C' = 'Control',
         'S' = 'Carbon',
         'U' = 'Nitrogen')) %>% 
  mutate(litter=recode(litter,
         'P' = 'Present',
         'A' = 'Absent'))

abundance <- read.csv(paste0("https://pasta.lternet.edu/package/data/eml/knb-lter-knz/101/4/12f22afc603c70b1174455e8636115bb",
                             "?key=", myEDIAccessKey)) %>% # invertebrate counts
  rename(year=RecYear, watershed=Watershed, block=Block, plot=Plot, burn_trt=BurnTrt, 
         order=invert_order, family=invert_family, collected=Collected, count=Count) %>% 
  mutate(block=str_squish(block)) %>% 
  group_by(year, watershed, block, plot, burn_trt, order, family) %>% 
  summarise(count = sum(count), .groups='drop') %>%  # combine collected and observed counts
  left_join(trt) %>%
  mutate(replicate=paste(burn_trt, watershed, block, plot, litter, nutrient, sep='::'))

biomass <- read.csv(paste0("https://pasta.lternet.edu/package/data/eml/knb-lter-knz/101/4/c773c9eb9480a03851ff9e818150d91f",
                           "?key=", myEDIAccessKey)) %>%  # invertebrate biomass
  rename(year=RecYear, watershed=Watershed, block=Block, plot=Plot, burn_trt=BurnTrt) %>% 
  select(-DataCode, -RecType, -RecMonth, -Comments) %>% 
  left_join(trt) %>% 
  filter(invertebrate_biomass>0) # drop single missing datapoint

CN <- read.csv(paste0("https://pasta.lternet.edu/package/data/eml/knb-lter-knz/101/4/4665458445755e751deeac090ddb6fb0",
                      "?key=", myEDIAccessKey)) %>% # plant %C and %N
  select(-DataCode, -RecType, -Comments) %>% 
  mutate(BurnTrt=str_to_sentence(BurnTrt)) %>% 
  rename(year=RecYear, watershed=Watershed, block=Block, plot=Plot, burn_trt=BurnTrt) %>% 
  mutate(CN=as.numeric(Total_C)/as.numeric(Total_N)) %>% 
  select(-Total_N, -Total_C) %>% 
  pivot_wider(names_from=Growth_form, values_from=CN) %>% 
  mutate(avg_CN=rowMeans(across(c(forb, grass)), na.rm=T)) %>% 
  select(-grass, -forb, -woody)

plant <- readRDS('plantData.RDS') # plant biomass and richness

functionalGroups <- read.csv('gf_funct_groups.csv')

# Community Metrics ---------------------------------------------------------------

totalAbundance <- abundance %>%
  group_by(year, watershed, block, plot, burn_trt, litter, nutrient) %>%
  summarise(total_count = sum(count), .groups = 'drop')

communityStructure <- community_structure(abundance, time.var='year', abundance.var='count', replicate.var='replicate', metric='Evar') %>%
  separate(col=replicate, into=c('burn_trt','watershed','block','plot', 'litter', 'nutrient'), sep='::') %>% 
  mutate(plot=as.integer(plot)) %>% 
  left_join(totalAbundance)

functionalStructure <- abundance %>% 
  left_join(functionalGroups) %>% 
  group_by(year, watershed, block, plot, burn_trt, litter, nutrient, eco_functional_group) %>% 
  summarise(funct_count=sum(count), .groups='drop')


# Pre-Treatment (2014) Analysis -----------------------------------------------------------

# Composition - PERMANOVA

abundance2014 <- abundance %>% 
  mutate(taxa=paste(order, family, sep='_')) %>% 
  select(year, watershed, block, plot, burn_trt, taxa, count) %>% 
  filter(year==2014) %>% 
  pivot_wider(names_from='taxa', values_from = 'count', values_fill = 0)  
  
permanova2014 <- adonis2(formula = abundance2014[,6:46] ~ burn_trt, data=abundance2014[,1:5], permutations=999, method='bray')
print(permanova2014)


# SIMPER and bar graph

summary(simper2014 <- simper(abundance2014[, 6:46], abundance2014$burn_trt, permutations = 999))

simper2014Table <- summary(simper2014)$Annual_Unburned %>%
  as.data.frame() %>%
  arrange(desc(average)) %>%
  rownames_to_column(var='species')

ggplot(simper2014Table, aes(x=reorder(species, average), y=average)) +
  geom_col() +
  coord_flip() +
  labs(x='Species', y='Contribution to Dissimilarity')

sppBC <- metaMDS(abundance2014[,6:46])

nmds_df <- data.frame(
  scores(sppBC, display = "sites"),
  burn_trt = abundance2014$burn_trt
)

sp_scores <- as.data.frame(scores(sppBC, display='species')) %>% 
  rownames_to_column(var='species')

nmdsSpecies <- left_join(simper2014Table, sp_scores) %>% 
  mutate(species=str_to_title(str_replace(species, "_", " "))) %>% 
  separate(species, into=c('order', 'family'), remove=F)

# function to create ellipse coordinates
veganCovEllipse <- function(cov, center = c(0,0), scale = 1, npoints = 100) {
  theta <- seq(0, 2 * pi, length.out = npoints)
  circle <- cbind(cos(theta), sin(theta))
  ellipse <- t(center + scale * t(circle %*% chol(cov)))
  as.data.frame(ellipse)
}

ord_ell <- ordiellipse(sppBC, groups = abundance2014$burn_trt, display = "sites", kind = "se", conf = 0.95, draw = "none")

ellipse_df <- bind_rows(lapply(names(ord_ell), 
                               function(g) {
                                 df <- veganCovEllipse(ord_ell[[g]]$cov, ord_ell[[g]]$center, ord_ell[[g]]$scale)
                                 colnames(df) <- c("NMDS1", "NMDS2")
                                 df$burn_trt <- g
                                 df
                                 }))

nmdsFig2014 <- ggplot() +
  geom_point(data=nmds_df, aes(x = NMDS1, y = NMDS2, color = burn_trt), size = 6) + # plot loadings
  geom_path(data = ellipse_df, aes(x = NMDS1, y = NMDS2, color = burn_trt), linewidth = 1.5) +
  scale_color_manual(values = c("#de1a24", "#056517")) +
  labs(color = "Burn Treatment") +
  # geom_point(data=nmdsSpecies, aes(x=NMDS1, y=NMDS2, size=average), color='darkgrey') + 
  geom_text(data=subset(nmdsSpecies, cumsum<0.8), aes(x=NMDS1, y=NMDS2, label=family), color='black', size=5) + # species loadings
  annotate("text", x=-Inf, y=Inf, label='(e)', hjust=-0.2, vjust=1.2, size=6) +
  theme(axis.text = element_text(size = 24, color = "black"),
        legend.text = element_text(size = 18),
        legend.position='none')



# Richness
richness2014 <- lmer(richness ~ burn_trt + (1 | watershed), data = subset(communityStructure, year==2014))
summary(richness2014)
anova(richness2014)

richnessFig2014 <- ggplot(communityStructure, aes(x = burn_trt, y = richness, color=burn_trt)) +
  geom_boxplot() +
  scale_color_manual(values = c("#de1a24", "#056517")) +
  xlab("") +
  ylab("Family Richness") +
  annotate("text", x=-Inf, y=Inf, label='(c)', hjust=-0.2, vjust=1.2, size=6) +
  theme(legend.position='none',
        plot.margin = margin(5, 20, 5, 5))


# Evenness
evenness2014 <- lmer(Evar ~ burn_trt + (1 | watershed), data = subset(communityStructure, year==2014))
summary(evenness2014)
anova(evenness2014)

evennessFig2014 <- ggplot(communityStructure, aes(x = burn_trt, y = Evar, color=burn_trt)) +
  geom_boxplot() +
  scale_color_manual(values = c("#de1a24", "#056517")) +
  xlab("") +
  ylab("Family Evenness") +
  annotate("text", x=-Inf, y=Inf, label='(d)', hjust=-0.2, vjust=1.2, size=6) +
  theme(legend.position='none',
        plot.margin = margin(5, 5, 5, 20))


# Abundance
count2014 <- lmer(total_count ~ burn_trt + (1 | watershed), data = subset(communityStructure, year==2014))
summary(count2014)
anova(count2014)

countFig2014 <- ggplot(communityStructure, aes(x = burn_trt, y = total_count, color=burn_trt)) +
  geom_boxplot() +
  scale_color_manual(values = c("#de1a24", "#056517")) +
  xlab("") +
  ylab("Total Abundance") +
  annotate("text", x=-Inf, y=Inf, label='(a)', hjust=-0.2, vjust=1.2, size=6) +
  theme(legend.position='none',
        plot.margin = margin(5, 20, 5, 5))


# Biomass
biomass2014 <- lmer(invertebrate_biomass ~ burn_trt + (1 | watershed), data = subset(biomass, year==2014))
summary(biomass2014)
anova(biomass2014)

biomassFig2014 <- ggplot(subset(biomass, year==2014), aes(x = burn_trt, y = invertebrate_biomass, color=burn_trt)) +
  geom_boxplot() +
  scale_color_manual(values = c("#de1a24", "#056517")) +
  xlab("") +
  ylab("Total Biomass (mg)") +
  annotate("text", x=-Inf, y=Inf, label='(b)', hjust=-0.2, vjust=1.2, size=6) +
  theme(legend.position='none',
        plot.margin = margin(5, 5, 5, 20))


# Functional Groups
functionalGroup2014 <- lmer(funct_count ~ burn_trt + (1|watershed),
                            data = subset(functionalStructure, year==2014 & eco_functional_group=='herbivore'))
summary(functionalGroup2014)
anova(functionalGroup2014)

functionalGroup2014 <- lmer(funct_count ~ burn_trt + (1|watershed),
                            data = subset(functionalStructure, year==2014 & eco_functional_group=='predator'))
summary(functionalGroup2014)
anova(functionalGroup2014)


### Combined pre-treatment figure ###

top <- plot_grid(
  countFig2014, biomassFig2014,
  richnessFig2014, evennessFig2014,
  ncol = 2,
  rel_spacing = 0.1
)

plot_grid(
  top,
  nmdsFig2014,
  ncol = 1,
  rel_heights = c(10, 5)
)

# ggsave("Fig1_pretrt.png", width = 10, height = 15, dpi = 300)




# Treatment (2019, 2024) Analysis -----------------------------------------------------------

# Composition - PERMANOVA (qualitatively same results if running 2019 and 2024 in separate models)

abundanceTrt <- abundance %>% 
  mutate(taxa=paste(order, family, sep='_')) %>% 
  select(year, watershed, block, plot, burn_trt, litter, nutrient, taxa, count) %>% 
  filter(year!=2014) %>% 
  pivot_wider(names_from='taxa', values_from = 'count', values_fill = 0)  

permanovaTrt <- vegan::adonis2(formula = abundanceTrt[,8:156] ~ burn_trt + nutrient + litter, data=abundanceTrt[,1:7], permutations=999, method='bray', by='margin')
print(permanovaTrt)


# SIMPER and bar graph

summary(simperTrt <- simper(abundanceTrt[, 8:156], abundanceTrt$burn_trt, permutations = 999))

simperTrtTable <- summary(simperTrt)$Annual_Unburned %>%
  as.data.frame() %>%
  arrange(desc(average)) %>%
  rownames_to_column(var='species')

ggplot(simperTrtTable, aes(x=reorder(species, average), y=average)) +
  geom_col() +
  coord_flip() +
  labs(x='Species', y='Contribution to Dissimilarity')

# sppBC <- metaMDS(abundanceTrt[, 8:156])
# 
# nmds_df <- data.frame(
#   scores(sppBC, display = "sites"),
#   burn_trt = abundanceTrt$burn_trt
# )
# 
# sp_scores <- as.data.frame(scores(sppBC, display='species')) %>%
#   rownames_to_column(var='species')
# 
# nmdsSpecies <- left_join(simperTrtTable, sp_scores) %>%
#   mutate(species=str_to_title(str_replace(species, "_", " "))) %>%
#   separate(species, into=c('order', 'family'), remove=F)
# 
# # function to create ellipse coordinates
# veganCovEllipse <- function(cov, center = c(0,0), scale = 1, npoints = 100) {
#   theta <- seq(0, 2 * pi, length.out = npoints)
#   circle <- cbind(cos(theta), sin(theta))
#   ellipse <- t(center + scale * t(circle %*% chol(cov)))
#   as.data.frame(ellipse)
# }
# 
# ord_ell <- ordiellipse(sppBC, groups = abundanceTrt$burn_trt, display = "sites", kind = "se", conf = 0.95, draw = "none")
# 
# ellipse_df <- bind_rows(lapply(names(ord_ell),
#                                function(g) {
#                                  df <- veganCovEllipse(ord_ell[[g]]$cov, ord_ell[[g]]$center, ord_ell[[g]]$scale)
#                                  colnames(df) <- c("NMDS1", "NMDS2")
#                                  df$burn_trt <- g
#                                  df
#                                }))
# 
# nmdsFigTrt <- ggplot() +
#   geom_point(data=nmds_df, aes(x = NMDS1, y = NMDS2, color = burn_trt), size = 6) + # plot loadings
#   geom_path(data = ellipse_df, aes(x = NMDS1, y = NMDS2, color = burn_trt), linewidth = 1.5) +
#   scale_color_manual(values = c("#de1a24", "#056517")) +
#   labs(color = "Burn Treatment") +
#   # geom_point(data=nmdsSpecies, aes(x=NMDS1, y=NMDS2, size=average), color='darkgrey') +
#   geom_text(data=subset(nmdsSpecies, cumsum<0.8), aes(x=NMDS1, y=NMDS2, label=family), color='black', size=5) + # species loadings
#   annotate("text", x=-Inf, y=Inf, label='(e)', hjust=-0.2, vjust=1.2, size=6) +
#   theme(axis.text = element_text(size = 24, color = "black"),
#         legend.text = element_text(size = 18),
#         legend.position='none')



# Richness
richnessTrt <- lmer(richness ~ burn_trt*nutrient*litter*as.factor(year) + (1 | watershed), data = subset(communityStructure, year!=2014))
summary(richnessTrt)
anova(richnessTrt)
emmeans(richnessTrt, ~ nutrient*litter*as.factor(year))

richnessFigTrt <- ggplot(subset(communityStructure, year!=2014), aes(x = nutrient, y = richness, color=litter)) +
  geom_boxplot() +
  scale_color_manual(values = c("darkgreen", "tan")) +
  xlab("") +
  ylab("Family Richness") +
  annotate("text", x=-Inf, y=Inf, label='(c)', hjust=-0.2, vjust=1.2, size=6) +
  facet_wrap(~year) +
  theme(legend.position=c(0.2,0.9),
        plot.margin = margin(5, 20, 5, 5))


# Evenness
evennessTrt <- lmer(Evar ~ burn_trt*nutrient*litter*as.factor(year) + (1 | watershed), data = subset(communityStructure, year!=2014))
summary(evennessTrt)
anova(evennessTrt)
emmeans(evennessTrt, ~ nutrient*litter*as.factor(year))

evennessFigTrt <- ggplot(subset(communityStructure, year!=2014), aes(x = nutrient, y = Evar, color=litter)) +
  geom_boxplot() +
  scale_color_manual(values = c("darkgreen", "tan")) +
  xlab("") +
  ylab("Family Evenness") +
  annotate("text", x=-Inf, y=Inf, label='(d)', hjust=-0.2, vjust=1.2, size=6) +
  facet_wrap(~year) +
  theme(legend.position=c(0.2,0.9),
        plot.margin = margin(5, 20, 5, 5))


# Abundance
countTrt <- lmer(total_count ~ as.factor(year)*burn_trt*nutrient*litter + (1 | watershed), data = subset(communityStructure, year!=2014))
summary(countTrt)
anova(countTrt)
emmeans(countTrt, ~ burn_trt*as.factor(year))

countFigTrt <- ggplot(subset(communityStructure, year!=2014), aes(x = burn_trt, y = total_count)) +
  geom_boxplot() +
  xlab("") +
  ylab("Total Abundance") +
  annotate("text", x=-Inf, y=Inf, label='(a)', hjust=-0.2, vjust=1.2, size=6) +
  facet_wrap(~year) +
  theme(legend.position=c(0.2,0.9),
        plot.margin = margin(5, 20, 5, 5))


# Biomass
biomassTrt <- lmer(invertebrate_biomass ~ burn_trt*nutrient*litter + (1 | watershed), data = subset(biomass, year==2024))
summary(biomassTrt)
anova(biomassTrt)
emmeans(biomassTrt, ~nutrient)

biomassFigTrt <- ggplot(subset(biomass, year==2024), aes(x = nutrient, y = invertebrate_biomass)) +
  geom_boxplot() +
  xlab("") +
  ylab("Total Biomass (mg)") +
  annotate("text", x=-Inf, y=Inf, label='(b)', hjust=-0.2, vjust=1.2, size=6) +
  theme(legend.position='none',
        plot.margin = margin(5, 5, 5, 20))


# Combined figure
plot_grid(
  countFigTrt, biomassFigTrt,
  richnessFigTrt, evennessFigTrt,
  ncol = 2,
  rel_spacing = 0.1
)

# ggsave("Fig2_trt.png", width = 25, height = 25, dpi = 300)



# # Functional Groups
# functionalGroupTrt <- lmer(funct_count ~ burn_trt*nutrient*litter*as.factor(year) + (1|watershed),
#                            data = subset(functionalStructure, year!=2014 & eco_functional_group=='herbivore' & funct_count<300))
# summary(functionalGroupTrt)
# anova(functionalGroupTrt)
# 
# functionalGroupTrt <- lmer(funct_count ~ burn_trt*nutrient*litter*as.factor(year) + (1|watershed),
#                            data = subset(functionalStructure, year!=2014 & eco_functional_group=='predator' & funct_count<40))
# summary(functionalGroupTrt)
# anova(functionalGroupTrt)


# Path Analysis -----------------------------------------------

allData <- communityStructure %>% 
  left_join(plant) %>% 
  left_join(CN) %>% 
  filter(year!=2014) %>% 
  mutate(burn_trt=as.factor(burn_trt),
         litter=as.factor(litter),
         nutrient=as.factor(nutrient))

SEMdata <- allData %>% 
  mutate(across(c(richness, Evar, total_count, litter_biomass, live_biomass, avg_CN, plant_richness), ~ as.numeric(scale(.x))),
         burn=ifelse(burn_trt=='Annual', 1, 0),
         litter_trt=ifelse(litter=='Present', 1, 0),
         nutrient_trt = c(Carbon = -1, Control = 0, Nitrogen = 1)[as.character(nutrient)])

emmeans::emm_options(lmer.df = "satterthwaite",
                     disable.pbkrtest = TRUE)


# invertebrate abundance
count_sem <- psem(
  lmer(total_count ~ burn + litter_biomass + live_biomass + avg_CN + plant_richness + (1|watershed), data = SEMdata),
  lmer(litter_biomass ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  lmer(live_biomass ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  lmer(avg_CN ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  lmer(plant_richness ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  avg_CN %~~% live_biomass,
  avg_CN %~~% litter_biomass,
  avg_CN %~~% plant_richness,
  live_biomass %~~% litter_biomass,
  live_biomass %~~% plant_richness,
  litter_biomass %~~% plant_richness,
  data = SEMdata
)
summary(count_sem)


# invertebrate richness
div_sem <- psem(
  lmer(richness ~ burn + litter_biomass + live_biomass + avg_CN + plant_richness + (1|watershed), data = SEMdata),
  lmer(litter_biomass ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  lmer(live_biomass ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  lmer(avg_CN ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  lmer(plant_richness ~ burn + litter_trt + nutrient_trt + (1|watershed), data = SEMdata),
  avg_CN %~~% live_biomass,
  avg_CN %~~% litter_biomass,
  avg_CN %~~% plant_richness,
  live_biomass %~~% litter_biomass,
  live_biomass %~~% plant_richness,
  litter_biomass %~~% plant_richness,
  data = SEMdata
)
summary(div_sem)


## treatment plant regressions
summary(lm(total_count ~ plant_richness, data = allData))


rich_trt <- ggplot(abun_trt %>% filter(year != 2014),
                   aes(x = plant_richness, y = total_abun)) +
  geom_point(aes(shape = plot_trt, color = litter_trt), size = 3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Plant Richness") +
  ylab("Arthropod Abundance")

summary(lm(total_count ~ live_biomass, data = allData))


live_trt <- ggplot(allData, aes(x = live_biomass, y = total_count)) +
  geom_point(aes(shape = nutrient, color = litter), size = 3) +  
  geom_smooth(method = "lm", se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Live Biomass") +
  ylab("Arthropod Abundance")


summary(lm(total_count ~ litter_biomass, data = allData))


litter_trt <- ggplot(allData, aes(x = litter_biomass, y = total_count)) +
  geom_point(aes(shape = nutrient, color = litter), size = 3) +  
  geom_smooth(method = "lm",  se = F, color = "black") +
  scale_shape_manual(values = c(15, 19, 17)) +
  scale_color_manual(values = c("#337539", "#dccd7d")) +
  xlab("Litter Biomass") +
  ylab("Arthropod Abundance")
