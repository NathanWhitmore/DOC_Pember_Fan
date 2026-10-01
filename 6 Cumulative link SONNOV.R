library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)
library(arm)
library(sf)
library(performance)
library(ordinal)

# load modified data frames
pember <- readRDS("Pember.rds")
cover.sf <- readRDS("Cover sf.rds")
cover <- cover.sf

# quick check
# pember %>% filter(NVSSpeciesName == "Raoulia monroi") # absent 

# standardise pembr year
pember$Year <- as.numeric(pember$Year) - 2018

# simplify
cover <- cover[, c("Year", "Plot", "Subplot", "BG", "fence.dist", "Transect")]
pember <- pember[, c("Year", "Plot", "Subplot","NVSSpeciesName",
                     "TaxonBioStatus", "TaxonGrowthForm", "Indigenous",  "Proportion", 
                     "Cover")]

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
my.species <- "Sonchus novae-zelandiae"

# my.species <- "Sonchus novae-zelandiae"
pember$key.species <- ifelse(pember$NVSSpeciesName == my.species,
                             pember$Proportion, 
                             0 )

# fix cover class
pember$Cover <- ifelse(is.na(pember$Cover), 0, pember$Cover)
pember$Cover <- factor(
  pember$Cover,
  levels = c("0", "P", "1", "2", "3", "4", "5", "6"),
  ordered = TRUE
)
 
# pember$not.key.species <- 1 - pember$key.species 

# simplify
just.key.species <- pember %>% filter(pember$NVSSpeciesName == my.species)
just.key.species <- just.key.species[, c("moniker", "Cover")]

unique(just.key.species$Cover)

# this fixes any possible accidental absences
critical <- pember %>%
  group_by(Year, Plot, Subplot, moniker) %>%
  summarise(Species.prop = sum(key.species))

critical <- left_join(critical, just.key.species, by = "moniker")

# fix up cover class
critical$Cover[is.na(critical$Cover)] <- "0"

# spatial join
all <- left_join(critical, cover, by = "moniker")

# only keep presence data
presences <- all %>% filter(Species.prop != 0)
head(presences)

# correct for bare ground
presences$correct.prop <- (1 - presences$BG) * presences$Species.prop

# make spatial
presences.sf <- st_as_sf(presences)

# start with presence of bare ground
all$correct.prop <- (1 - all$BG) * all$Species.prop
all$presence <- ifelse(all$correct.prop == 0, 0, 1)
all$presence <- ifelse(is.na(all$correct.prop), 0, all$presence)

# add X-Y coordinates
all<- all %>% st_as_sf()
all.coord <- st_coordinates(all)
all <- cbind(all, all.coord )

# check total should be 12000
nrow(all)

# restrict to southern section for Sonchus and Brachyscome
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
          colour = "purple")+
  facet_wrap(~Year)+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle(my.species)+
  labs(size = "Corrected proportion")+
  theme(plot.title = element_text(face = "italic"))

ggsave("Sonchus novae-zelandiae.png", scale =1.1, height = 6, width =8)



#  make subplot id because clmm can't handle it otherwise
all$SubplotID <- interaction(
  all$Plot,
  all$Subplot,
  drop = TRUE
)

clmm.mod <- list()

clmm.mod[[1]] <- clmm(Cover ~ 1 + (1|Plot) + (1 | Subplot), data = all)
clmm.mod[[2]] <- clmm(Cover ~ Year  + (1|Plot) + (1 | SubplotID), data = all)
clmm.mod[[3]] <- clmm(Cover ~ scale(fence.dist) + (1|Plot) + (1 | SubplotID), data = all)
clmm.mod[[4]] <- clmm(Cover ~ Year + scale(fence.dist) + (1|Plot) + (1 | SubplotID), data = all)
clmm.mod[[5]] <- clmm(Cover ~ as.factor(Year) + scale(fence.dist) + (1|Plot) + (1 | SubplotID),  data = all)
clmm.mod[[6]] <- clmm(Cover ~ as.factor(Year) + (1|Plot) + (1 | SubplotID), data = all)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(clmm.mod), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(clmm.mod, formula))))

# AIC table to 4 digits
clmm.aic <- aictab(cand.set = clmm.mod, modnames = Modnames, sort = TRUE) %>% 
  as.data.frame()

clmm.aic <- clmm.aic %>% dplyr::select(-ModelLik, -Cum.Wt)
clmm.aic[, 3:6] <- round(clmm.aic[, 3:6],3)
clmm.aic

# manual predictions
# summary
summary(clmm.mod[[2]])

# checks - Plot random effect not really supported
ranef(clmm.mod[[2]])

# simplify scale distance
# don't need to use in prediction as mean = 0
fence.mean <- mean(all$fence.dist, na.rm = TRUE)
fence.sd <- sd(all$fence.dist, na.rm = TRUE)

fence.dist.scaled <- (all$fence.dist - fence.mean)/fence.sd



# makes two columns for each random effect
my.ranef <- ranef(clmm.mod[[2]]) %>% as.data.frame()
my.ranef$total.random <- rowSums(my.ranef)

# choose subplots where the species is likely to be found
plot.subplot <- quantile(my.ranef$total.random, 0.97)


# prediction loop

summary(clmm.mod[[2]])

store <- NULL

my.year <- 0:7

for(i in 1:8){
  
  # get betas
  my.year <- 0:7
  
  C0 <- plogis(13.622 - ( -0.08346 * my.year[i] + plot.subplot))
  C1 <- plogis(14.282  - ( -0.08346 * my.year[i] + plot.subplot))
  C2 <- plogis(16.114  - ( -0.08346 * my.year[i] + plot.subplot))
  C3 <- plogis(18.215  - ( -0.08346 * my.year[i] + plot.subplot))

  
  
  store[[i]] <- data.frame(cover.0 = round(C0 - 0, 3),
                           cover.1 = round(C1 - C0, 3),
                           cover.2 = round(C2 - C1, 3),
                           cover.3 = round(C3 - C2, 3),
                           cover.4 = round(1 - C3, 3),
                           cover.5 = 0,
                           cover.6 = 0)
  
}

my.predictions  <- bind_rows(store)
my.predictions$Year <- 2018:2025

my.predictions


my.long <- my.predictions  |>
  pivot_longer(cols = cover.0:cover.6,
               names_to =  "Cover",
               values_to = "Proportion")

# filter out cover class 0
my.long <-  my.long  |> filter(Cover != "cover.0")

my.long <- my.long %>% 
  mutate(Cover = str_replace_all(Cover, "cover.", ""))

my.long$Cover <- as.numeric(as.character(my.long$Cover))

str(my.long)


ggplot()+
  theme_bw()+
  geom_col(data = my.long, aes( x= Cover, y= Proportion),
           fill = "purple", colour = "white") +
  facet_grid(~Year) +
  theme(axis.title = element_text(
    size = 14,
    colour = "grey40"
  )) + theme(axis.title = element_text(
    size = 14,
    colour = "grey40")) +
  theme(strip.text = element_text(
    size = 12,
    colour = "grey40"
  ))+
  theme(panel.grid = element_blank())+
  ylab("Count\n")+
  xlab("\nCover class")

ggsave("Predicted changes in SONNOV.png", scale = 1.1, height = 6, width =8)

