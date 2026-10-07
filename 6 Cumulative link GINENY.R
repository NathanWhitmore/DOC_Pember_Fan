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
                    levels = c("0", "P", "1", "2", "3", "4", "5", "6"),
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
                       levels = c("0", "P", "1", "2", "3", "4", "5", "6"),
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

# ggsave("Gingidia enysii.png", scale =1.1, height = 6, width =8)

# makes two columns for each random effect
my.ranef <- ranef(clmm.mod[[5]]) %>% as.data.frame()
my.ranef$total.random <- rowSums(my.ranef)


plot.subplot <- quantile(my.ranef$total.random, 0.95)
store <- NULL

my.year <- 0:7

  for(i in 1:8){
       
         # get betas
         my.year <- c(0, clmm.mod[[5]]$beta)
       
           
         C0 <- plogis(1.4988 - ( my.year[i] + plot.subplot))
         C1 <- plogis(2.4018 - ( my.year[i] + plot.subplot))
         C2 <- plogis(3.3364 - ( my.year[i] + plot.subplot))
         C3 <- plogis(5.1509 - ( my.year[i] + plot.subplot))
         C4 <- plogis(7.2198 - ( my.year[i] + plot.subplot))
         C5 <- plogis(9.9134 - ( my.year[i] + plot.subplot))
         
           
           store[[i]] <- data.frame(cover.0 = round(C0 - 0, 3),
                                        cover.1 = round(C1 - C0, 3),
                                        cover.2 = round(C2 - C1, 3),
                                        cover.3 = round(C3 - C2, 3),
                                        cover.4 = round(C4 - C3, 3),
                                        cover.5 = round(C5 - C4, 3),
                                        cover.6 = round(1 - C5, 3))
           
  }
 
my.predictions  <- bind_rows(store)
my.predictions$Year <- 2018:2025

my.long <- my.predictions  |>
  pivot_longer(cols = cover.0:cover.6,
                   +               names_to =  "Cover",
                                  values_to = "Proportion")

# filter out cover class 0
my.long <-  my.long  |> filter(Cover != "cover.0")

my.long <- my.long %>% 
  mutate(Cover = str_replace_all(Cover, "cover.", ""))

my.long$Cover <- as.numeric(as.character(my.long$Cover))


ggplot()+
  theme_bw()+
  geom_col(data = my.long, aes( x= Cover, y= Proportion),
                        fill = "orange", colour = "white") +
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
 
ggsave("Predicted changes in bare ground cover class.png", scale = 1.1, height = 6, width =8)


