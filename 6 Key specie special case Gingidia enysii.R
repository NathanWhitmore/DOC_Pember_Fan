library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)
library(arm)
library(sf)
library(performance)

# load modified data frames
pember <- readRDS("Pember.rds")
cover.sf <- readRDS("Cover sf.rds")
cover <- cover.sf

# quick check
# pember %>% filter(NVSSpeciesName == "Raoulia monroi") # absent 
# notice Gingidia enysii has a zero entry
temp <- pember %>% filter(NVSSpeciesName == "Gingidia enysii") # exists
temp

# standardise pembr year
pember$Year <- as.numeric(pember$Year) - 2018

# simplify
cover <- cover[, c("Year", "Plot", "Subplot", "BG", "fence.dist", "Transect")]
pember <- pember[, c("Year", "Plot", "Subplot","NVSSpeciesName",
                     "TaxonBioStatus", "TaxonGrowthForm", "Indigenous", "Proportion")]

# moniker
cover$moniker <- paste(cover$Year, cover$Plot, cover$Subplot)
pember$moniker <- paste(pember$Year, pember$Plot, pember$Subplot)

# remove duplicate column names from cover
cover$Year <- NULL
cover$Plot <- NULL
cover$Subplot <- NULL

# make into dataframe
cover <- as.data.frame(cover)

# species names
sort(unique(pember$NVSSpeciesName))

# species
# my.species <- "Brachyscome pinnata"
# 
# my.species <- "Sonchus novae-zelandiae"

my.species <- "Gingidia enysii"

pember$key.species <- ifelse(pember$NVSSpeciesName == my.species,
                             pember$Proportion, 
                             0 )

# pember$not.key.species <- 1 - pember$key.species 

# simplify
critical <- pember %>%
  group_by(Year, Plot, Subplot, moniker) %>%
  summarise(Species.prop = sum(key.species))

# spatial join
all <- left_join(critical, cover, by = "moniker")
head(all)

# only keep presence data
presences <- all %>% filter(Species.prop != 0)
head(presences)

# correct for bare ground
presences$correct.prop <- (1 - presences$BG) * presences$Species.prop

# make spatial
presences.sf <- st_as_sf(presences)

# start with presence of bare ground
all$correct.prop <- (1 - all$BG) * all$Species.prop
all$hurdle <- ifelse(all$correct.prop == 0, 0, 1)
all$hurdle <- ifelse(is.na(all$correct.prop), 0, all$hurdle)

# add X-Y coordinates
all<- all %>% st_as_sf()
all.coord <- st_coordinates(all)
all <- cbind(all, all.coord )

# check total should be 12000
nrow(all)


# restrict to southern section for Sonchus and Brachyscome
# all <- all %>% filter(Y < 5225800)

# restrict to RR section for Gingidia enysii
all <- all %>% filter(Plot %in% c("RR3", "RR10"))

# map
ggplot()+
  theme_bw()+
 # geom_sf(data = all %>%
 #           filter(Y< 5225800), 
 #         shape = 3, size = 9, colour = "red")+
  geom_sf(data = cover.sf, 
          shape = 1)+
  geom_sf(data = presences.sf  , 
          aes(size = correct.prop), 
          alpha = 0.7,
          colour = "forestgreen")+
  facet_wrap(~Year)+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle(my.species)+
  labs(size = "Corrected proportion")+
  theme(plot.title = element_text(face = "italic"))

ggsave("Gingidia enysii.png", scale =1.1, height = 6, width =8)


# restrict analysis to RR 3
all <- all %>% filter(Plot == "RR3")

# Note Subplot only as we are only looking at a single plot
Cand.models.hurd <- list()

Cand.models.hurd[[1]] <- glmmTMB(hurdle ~ 1 + (1|Subplot), family = "binomial", 
                                 data = all)
Cand.models.hurd[[2]] <- glmmTMB(hurdle ~ Year  + (1|Subplot), family = "binomial",
                                 data = all)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models.hurd), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models.hurd, formula))))

# AIC table to 4 digits
hurdle <- aictab(cand.set = Cand.models.hurd, modnames = Modnames, sort = TRUE)
hurdle

# summary
summary(Cand.models.hurd[[2]])
ranef(Cand.models.hurd[[1]])

# diagnostics (some issues but we are really stretching the modelling here)
res <- simulateResiduals(Cand.models.hurd[[1]])
plot(res)


# PART 2 beta regression (can we determine the % of bare ground when present)
plant.prop <- all %>% filter(Species.prop != 0)



# model selection

Cand.models.prop <- list()

Cand.models.prop[[1]] <- glmmTMB(Species.prop ~ 1 + (1|Subplot), family = beta_family(link = "logit"), 
                                 data = plant.prop)
Cand.models.prop[[2]] <- glmmTMB(Species.prop ~ Year  + (1|Subplot), family = beta_family(link = "logit"),
                                 data = plant.prop)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models.prop), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models.prop, formula))))

# AIC table to 4 digits
plant.prop.aic <- aictab(cand.set = Cand.models.prop, modnames = Modnames, sort = TRUE)
plant.prop.aic

# diagnostics - has issues
res <- simulateResiduals(Cand.models.prop[[2]]) # reasonably good
plot(res)

# summary
summary(Cand.models.prop[[2]])
ranef(Cand.models.prop[[2]])

# model performance
model_performance(Cand.models.prop[[2]])

# Model 9 singular 
performance::check_singularity(Cand.models.prop[[9]])
performance::check_singularity(Cand.models.prop[[2]])

VarCorr(Cand.models.prop[[2]])

# average loss
1-exp(-0.13608 )

# make predictions
plant.prop$predicted <- fitted(Cand.models.prop[[2]])
plant.prop$group <- paste(plant.prop$Plot,  plant.prop$Subplot)

# predictions
ggplot()+
  theme_bw()+
  geom_line(data = plant.prop, aes(x = Year, y = predicted, 
                                   group = group))

ggsave("Expected loss")

# predictive agreement
ggplot()+
  theme_bw()+
  geom_point(data = plant.prop, aes(y = Species.prop, x = predicted))+
  geom_abline(slope = 1, intercept = 0, colour = "red")+
  xlab("\nPredicted") +
  ylab("Species proportion\n") +
  theme(panel.grid = element_blank())







