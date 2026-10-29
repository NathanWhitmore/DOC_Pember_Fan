library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)
library(arm)
library(sf)
library(performance)
library(kableExtra)

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
                     "TaxonBioStatus", "TaxonGrowthForm", "Indigenous", "Cover", "Proportion")]

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
my.species <- "Gingidia enysii"


# do this as a precrusor to ensure all plots are represented
pember$key.species <- ifelse(pember$NVSSpeciesName == my.species,
                             pember$Proportion, 
                             0 )

# pember$not.key.species <- 1 - pember$key.species 

# correct pember levels

pember$Cover <- factor(pember$Cover, 
                    levels = c("0", "1", "2", "3", "4", "5", "6"),
                    ordered = TRUE)

# simplify
critical <- pember %>%
  group_by(Year, Plot, Subplot, moniker) %>%
  summarise(Species.prop = sum(key.species))

# filter for species
cover <- pember %>% filter(NVSSpeciesName == my.species)
cover$Plot <- NULL
cover$Subplot <- NULL
cover$Year <- NULL

# spatial join
all <- left_join(critical, cover, by = "moniker")

# check  presence data
# all are in RR3 except for 1 instance in RR#
presences <- all %>% filter(Species.prop != 0)
presences

# restrict to just RR3
names(all)

# legacy data frame for later graphing
legacy <- all

# restrict to RR section for Gingidia enysii
all <- all %>% filter(Plot %in% c("RR3"))

all$Cover <- ifelse(is.na(all$Cover), "0", all$Cover)


# make subplot id because clmm can't handle it otherwise
all$SubplotID <- interaction(
  all$Plot,
  all$Subplot,
  drop = TRUE
)

# modelling  
all <- all |> ungroup()

# annoyingly ordered factors are wiped
all$Cover <- factor(all$Cover, 
                       levels = c("0",  "1", "2", "3", "4", "5", "6"),
                       ordered = TRUE)

clmm.mod <- list()

clmm.mod[[1]] <- clmm(Cover ~ 1 +  (1 | SubplotID), data = all)
clmm.mod[[2]] <- clmm(Cover ~ Year + (1 | SubplotID), data = all)
clmm.mod[[3]] <- clmm(Cover ~ as.factor(Year) +  (1 | SubplotID), data = all)

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

summary(clmm.mod[[1]])

# latex table
kable(clmm.aic, "latex")

# reattached spatial data
spatial <- readRDS("Spatial.rds")

trans.position <- left_join(legacy, spatial, by = "Plot" ) |>
  st_as_sf()


# map
ggplot()+
  theme_bw()+
 # geom_sf(data = all %>%
 #           filter(Y< 5225800), 
 #         shape = 3, size = 9, colour = "red")+
  geom_sf(data = cover.sf, 
          shape = 1)+
  geom_sf(data = trans.position , 
          aes(size = as.numeric(as.character(Cover))), 
          alpha = 0.4,
          colour = "forestgreen")+
  facet_wrap(~Year)+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle(my.species)+
  labs(size = "Cover class")+
  theme(plot.title = element_text(face = "italic"))+
  scale_size_continuous(
    breaks = 1:6,
    limits = c(1, 6)
  )


ggsave("Gingidia enysii.png", scale =1.1, height = 6, width =8)

