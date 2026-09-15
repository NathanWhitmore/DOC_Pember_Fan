library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)
library(arm)
library(sf)

# load modified data frames
pember <- readRDS("Pember.rds")
cover.sf <- readRDS("Cover sf.rds")
cover <- cover.sf

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
my.species <- "Sonchus novae-zelandiae"

# my.species <- "Pilosella officinarum"

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

ggplot()+
  geom_point(data = presences, aes(y = correct.prop, x= BG))+
  facet_wrap(~Year)

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

# restrict to southern section
all <- all %>% filter(Y < 5225800)

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

# ggsave("Brachyscome pinnata.png", scale =1.1, height = 6, width =8)
ggsave("Sonchus novae-zelandiae.png", scale =1.1, height = 6, width =8)

# temporal autocorrection

all$PlotSubplot <- interaction(
  all$Plot,
  all$Subplot,
  drop = TRUE
)

# model selection only on southern portion
# note scale(fence.dist) can nearly perfectly predict (1|Plot/Subplot)
Cand.models.hurd <- list()

Cand.models.hurd[[1]] <- glmmTMB(hurdle ~ 1 + (1|Plot), family = "binomial", 
                                 data = all)
Cand.models.hurd[[2]] <- glmmTMB(hurdle ~ Year  + (1|Plot), family = "binomial",
                                 data = all)
Cand.models.hurd[[3]] <- glmmTMB(hurdle ~ as.factor(Year) + (1|Plot), family = "binomial",
                                 data = all)
Cand.models.hurd[[4]] <- glmmTMB(hurdle ~ scale(fence.dist) + (1|Plot), family = "binomial",
                                 data = all)
Cand.models.hurd[[5]] <- glmmTMB(hurdle ~ Year + scale(fence.dist) + (1|Plot), family = "binomial",
                                 data = all)
Cand.models.hurd[[6]] <- glmmTMB(hurdle ~ as.factor(Year) + scale(fence.dist) + (1|Plot), family = "binomial",
                                 data = all)

# temporal auto correlation
# this produces an incredibly small value for ranef
 Cand.models.hurd[[7]] <- glmmTMB(hurdle ~ as.factor(Year) + ar1(as.factor(Year) + 0 | Plot) + (1|Plot), family = "binomial", data = all)

# ranef(Cand.models.hurd[[7]]
# diagnose(Cand.models.hurd[[7]])


# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models.hurd), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models.hurd, formula))))

# AIC table to 4 digits
hurdle <- aictab(cand.set = Cand.models.hurd, modnames = Modnames, sort = TRUE)
hurdle

# summary
summary(Cand.models.hurd[[4]])




# diagnostics - has issues
res <- simulateResiduals(Cand.models.hurd[[5]])
plot(res)

ranef(Cand.models.hurd[[1]])


# basic diagnostics
testOutliers(res, type = "bootstrap") # okay
testUniformity(res)
testDispersion(res) # large under dispersion
plotResiduals(res, all$Year)


# SPATIAL AUTOCORRELATION

# aggregate DHARMa residuals to Plot
res.plot <- recalculateResiduals(
  res,
  group = all$Plot
)

# One coordinate pair per Plot (not subplot like in all)
coords.plot <- all %>%
  group_by(Plot) %>%
  summarise(
    X = first(X),
    Y = first(Y),
  )


# check
nrow(coords.plot)
length(residuals(res.plot))

# test spatial autocorrelation
# no evidence at the plot scale
testSpatialAutocorrelation(
  res.plot,
  x = coords.plot$X,
  y = coords.plot$Y
)






