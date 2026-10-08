library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)
library(arm)
library(sf)
library(performance)
library(ordinal)
library(kableExtra)

# load modified data frames
pember <- readRDS("Pember cover.rds")

pember$Year <- as.numeric(pember$Year) -2018
pember$moniker <- paste(pember$Plot, pember$Subplot, pember$Year)
pember$Plot <- NULL
pember$Subplot <- NULL
pember$Year <- NULL

# add in index data
index <- readRDS("Plot index.rds")
index$moniker <- paste(index$Plot, index$Subplot, index$Year)

# add in spatial data
spatial <- readRDS("Spatial.rds")
spatial$Subplot <- NULL

# simplify
spatial <- as.data.frame(spatial)

pember <- pember[, c("NVSSpeciesName",
                     "TaxonBioStatus", "TaxonGrowthForm", "Indigenous",  
                     "Cover", "moniker")]

# species
my.species <- "Pilosella officinarum"
pember <- pember %>% filter(NVSSpeciesName == my.species)

# fix cover class
pember$Cover <- ifelse(is.na(pember$Cover), 0, pember$Cover)
pember$Cover <- factor(
  pember$Cover,
  levels = c("0",  "1", "2", "3", "4", "5", "6"),
  ordered = TRUE
)
 
# join with index
all <- left_join(index, pember, by = "moniker")

# find duplicated moniker value
issues <- all$moniker[duplicated(all$moniker)]
dup <- all %>% filter(moniker %in% issues)

# insure duplicates are removed
all <- all %>%
  distinct(moniker, .keep_all = TRUE)

# fix up cover class
all$Cover[is.na(all$Cover)] <- "0"

# index join
all <- left_join(all, spatial, by = "Plot")


# add X-Y coordinates
all<- all %>% st_as_sf()
all.coord <- st_coordinates(all)
all <- cbind(all, all.coord )

# check total should be 12000
nrow(all)

# original spatial
og.sf <- readRDS("Spatial.rds")

all$Transect.type <- substring(all$Plot, 1, 2)

# map
ggplot()+
  theme_bw()+
 # geom_sf(data = og.sf, 
  #        shape = 1)+
  # geom_sf(data = all %>% filter(as.numeric(Cover) >= 1)  , 
  geom_sf(data = all %>% filter(as.numeric(Cover) == 1) ,
          aes(size = Cover), 
          alpha = 0.05,
          colour = "skyblue")+
  facet_wrap(~Year)+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  ggtitle(my.species)+
  labs(size = "Corrected proportion")+
  theme(plot.title = element_text(face = "italic"))

ggsave("Pilosella officinarum.png", scale =1.1, height = 6, width =8)



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
clmm.mod[[6]] <- clmm(Cover ~ as.factor(Year) + scale(fence.dist) + Transect.type + (1|Plot) + (1 | SubplotID),  data = all)

# clmm.mod[[6]] <- clmm(Cover ~ as.factor(Year) + (1|Plot) + (1 | SubplotID), data = all)

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

# kable
kable(clmm.aic, "latex")


# manual predictions
# summary
summary(clmm.mod[[5]])

# checks - Plot random effect not really supported
ranef(clmm.mod[[6]])

# simplify scale distance
# don't need to use in prediction as mean = 0
fence.mean <- mean(all$fence.dist, na.rm = TRUE)
fence.sd <- sd(all$fence.dist, na.rm = TRUE)

fence.dist.scaled <- (all$fence.dist - fence.mean)/fence.sd

# makes two columns for each random effect
my.ranef <- ranef(clmm.mod[[5]]) %>% as.data.frame()
my.ranef$total.random <- rowSums(my.ranef)

# choose subplots where the species is likely to be found
plot.subplot <- quantile(my.ranef$total.random, 0.5)


# prediction loop

summary(clmm.mod[[6]])

store <- NULL

my.year <- 0:7

for(i in 1:8){
  
  # get betas
  my.year <- c(0, clmm.mod[[6]]$beta[1:7])
  Transect.type <- 1.72792 
  
  C0 <- plogis(-2.2873886 - ( my.year[i] + plot.subplot + Transect.type + -0.97028 * 2))
  C1 <- plogis(-1.7844467 - ( my.year[i] + plot.subplot + Transect.type+ -0.97028 * 2))
  C2 <- plogis(-1.2499843  - ( my.year[i] + plot.subplot + Transect.type+ -0.97028 * 2))
  C3 <- plogis(-0.0005413  - ( my.year[i] + plot.subplot + Transect.type+ -0.97028 * 2))
  C4 <- plogis(1.4563170 - ( my.year[i] + plot.subplot + Transect.type+ -0.97028 * 2))
  C5 <- plogis(4.4427758 - ( my.year[i] + plot.subplot + Transect.type+ -0.97028 * 2))
  
  
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


ggplot()+
  theme_bw()+
  geom_col(data = my.long, aes( x= Cover, y= Proportion),
           fill = "skyblue", colour = "white") +
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
  ylab("Percentage in class \n(excluding absences)\n")+
  xlab("\nCover class")


ggsave("Predicted changes in PILOFF.png", scale = 1.1, height = 6, width =8)

# check original

re <- ranef(clmm.mod[[6]]) %>% as.data.frame()
head(re)
re$total <- rowSums(re)
re$transect <- substring(rownames(re),1,2)

ggplot()+
  geom_histogram(data = re, aes(x = total, fill = transect))

ggplot()+
  theme_bw()+
  geom_histogram(data = all %>% filter(as.numeric(Cover) >= 1),
                 aes(x = as.numeric(Cover)), binwidth = 1, 
                 colour = "white", 
                 fill = "skyblue")+
  facet_grid(~Year)

