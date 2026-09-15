library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)
library(arm)
library(sf)

# read in data
ground <- read.csv("Ground cover2.csv")

# make into data.frame
ground <- as.data.frame(ground)

# make subplot factor
ground$Subplot <- as.factor(ground$Subplot)

# make Cover factor
ground$Cover <- as.factor(ground$Cover)

# renaming categories
ground <- ground %>%
  ungroup() %>%
  mutate(
    GroundCover = fct_recode(
      GroundCover,
      "BG" = "BG2",
      "BG" = "BG3",
      "Li" = "Li3"
    )
  )


# translate values to means
L1 <- mean(0:1)  / 100 # 1
L2 <- mean(1:5) / 100 # 2
L3 <- mean(6:25)  / 100 # 3
L4 <- mean(26:50) / 100 # 4
L5 <- mean(51:75)  / 100 # 5
L6 <- mean(76:100) / 100 # 6

# make low value for "P"
ground <-  ground %>%
  mutate(Perc = as.numeric(recode(Cover, 
                                 # "P" = P, # Not found in cover
                                  "1" = L1,
                                  "2" = L2,
                                  "3" = L3,
                                  "4" = L4,
                                  "5" = L5,
                                  "6" = L6))
  )

# make into dataframe
ground <- as.data.frame(ground)

# make subplot factor
ground$Subplot <- as.factor(ground$Subplot)

# summarise in terms of raw percentage
ground <- ground %>% 
  group_by(Year, Plot, Subplot) %>%
  mutate(Weight = sum(Perc, na.rm = TRUE))

# corrected perc
ground$Proportion <- ground$Perc/ ground$Weight

# make into data.frame
ground <- as.data.frame(ground)

# check unique years
unique(ground$Year)

set.seed(18)

# graph
ggplot()+
  theme_bw()+
  geom_col(data = ground, aes(x = Year, y = Proportion, fill = GroundCover), 
           position = "fill")+
  theme(axis.title = element_text(
    face = 2,
    size = 14,
    colour = "grey40"
  )) +
  theme(strip.text = element_text(
    size = 12,
    colour = "grey40"
  )) +
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0.5))

# modelling

non.vege <-  ground %>% 
  mutate(
  GroundCover = fct_recode(
    GroundCover,
    "BG" = "BG",
    "Not BG" = "BR",  # getting rid of bare rock
    "Not BG" = "Li",
    "Not BG" = "L",
    "Not BG" = "M",
  )
)

# calculate subplot proportions
bare.ground <- non.vege %>%
  ungroup() %>%
  group_by(Year, Plot, Subplot, GroundCover) %>%
  summarise(Prop = sum(Proportion, na.rm = TRUE))


bare.ground.wide <- bare.ground %>% 
  pivot_wider(names_from = GroundCover,
              values_from = Prop,
              values_fill = 0)

# make sure there are no zeros or 1s
# bare.ground.wide$BG <- ifelse(bare.ground.wide$BG == 0, 0.001, bare.ground.wide$BG)
# bare.ground.wide$BG <- ifelse(bare.ground.wide$BG == 1, 0.999, bare.ground.wide$BG)

# read in spatial data
my.coord.sf <- readRDS("my_coord_sf.rds")

# join with spatial data
bare.wide.sf <- left_join(bare.ground.wide, my.coord.sf, by = "Plot") %>% 
  st_as_sf()

# check in case missed results - okay
# anti_join(bare.ground, my.coord.sf, by = "Plot")

# combine spatial data
bare.wide.sf$Year  <- bare.wide.sf$Year - 2018

# add in fence distance
fence <- st_read("fence line/Fence_line.json")

# incorporate fence line 
fence <- fence  %>% st_transform(crs = 2193)

# finalised data set
final.bare.wide <- bare.wide.sf

final.bare.wide$fence.dist <- as.numeric(
  st_distance(bare.wide.sf, fence 
  )
)

# possible that transect type matters
final.bare.wide$Transect <- substr(final.bare.wide$Plot,1,2)

# map version
cover.map.version <- final.bare.wide

saveRDS(cover.map.version, "Cover sf.rds")

# remove geometry
st_geometry(final.bare.wide) <- NULL

##### HURDLE MODEL

# start with presence of bare ground
final.bare.wide$hurdle <- ifelse(final.bare.wide$BG == 0, 0, 1)



Cand.models.hurd <- list()

Cand.models.hurd[[1]] <- glmmTMB(hurdle ~ 1 + (1|Plot/Subplot), family = "binomial", 
                                data = final.bare.wide)
Cand.models.hurd[[2]] <- glmmTMB(hurdle ~ Year  + (1|Plot/Subplot), family = "binomial",
                                data = final.bare.wide)
Cand.models.hurd[[3]] <- glmmTMB(hurdle ~ scale(fence.dist) + (1|Plot/Subplot), family = "binomial",
                                data = final.bare.wide)
Cand.models.hurd[[4]] <- glmmTMB(hurdle ~ Year + scale(fence.dist) + (1|Plot/Subplot), family = "binomial",
                                data = final.bare.wide)
Cand.models.hurd[[5]] <- glmmTMB(hurdle ~ as.factor(Year) + (1|Plot/Subplot), family = "binomial",
                                data = final.bare.wide)
Cand.models.hurd[[6]] <- glmmTMB(hurdle ~ Transect + (1|Plot/Subplot), family = "binomial", 
                                 data = final.bare.wide)
Cand.models.hurd[[7]] <- glmmTMB(hurdle ~ Year  + Transect  +(1|Plot/Subplot), family = "binomial",
                                 data = final.bare.wide)
Cand.models.hurd[[8]] <- glmmTMB(hurdle ~ scale(fence.dist) + Transect + (1|Plot/Subplot), family = "binomial",
                                 data = final.bare.wide)
Cand.models.hurd[[9]] <- glmmTMB(hurdle ~ Year + scale(fence.dist) + Transect + (1|Plot/Subplot), family = "binomial",
                                 data = final.bare.wide)
Cand.models.hurd[[10]] <- glmmTMB(hurdle ~ as.factor(Year) + Transect + (1|Plot/Subplot), family = "binomial",
                                 data = final.bare.wide)


# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models.hurd), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models.hurd, formula))))

# AIC table to 4 digits
hurdle <- aictab(cand.set = Cand.models.hurd, modnames = Modnames, sort = TRUE)
hurdle

# summary
summary(Cand.models.hurd[[5]])

# diagnostics - all good no issues
res <- simulateResiduals(Cand.models.hurd[[5]])
plot(res)


# PART 2 beta regression (can we determine the % of bare ground when present)
some.bg <- final.bare.wide %>% filter(BG != 0)

# model selection

Cand.models.prop <- list()

Cand.models.prop[[1]] <- glmmTMB(BG ~ 1 + (1|Plot/Subplot), family = beta_family(link = "logit"), 
                                data = some.bg)
Cand.models.prop[[2]] <- glmmTMB(BG ~ Year  + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = some.bg)
Cand.models.prop[[3]] <- glmmTMB(BG ~ scale(fence.dist) + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = some.bg)
Cand.models.prop[[4]] <- glmmTMB(BG ~ Year + scale(fence.dist) + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = some.bg)
Cand.models.prop[[5]] <- glmmTMB(BG ~ as.factor(Year) + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = some.bg)

Cand.models.prop[[6]] <- glmmTMB(BG ~ Transect + (1|Plot/Subplot), family = beta_family(link = "logit"), 
                                     data = some.bg)
Cand.models.prop[[7]] <- glmmTMB(BG ~ Year  + Transect + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                     data = some.bg)
Cand.models.prop[[8]] <- glmmTMB(BG ~ scale(fence.dist) + Transect + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                     data = some.bg)
Cand.models.prop[[9]] <- glmmTMB(BG ~ Year + scale(fence.dist) + Transect + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                     data = some.bg)
Cand.models.prop[[10]] <- glmmTMB(BG ~ as.factor(Year) + Transect + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                     data = some.bg)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models.prop), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models.prop, formula))))

# AIC table to 4 digits
BG.presence <- aictab(cand.set = Cand.models.prop, modnames = Modnames, sort = TRUE)
BG.presence

# summary
summary(Cand.models.presence[[7]])

# diagnostics performing very poorly
res <- simulateResiduals(Cand.models.presence[[7]])
plot(res)

# give away multimodal residuals
testOutliers(res, type = "bootstrap") # good

# multimodal distribution of proportions
ggplot()+
  theme_bw()+
  geom_histogram(data = some.bg, aes(x = BG), binwidth =0.02, 
                 fill = "forestgreen")+
  facet_grid(Year~.)+
  theme(aspect.ratio = 0.1)+
  theme(panel.grid = element_blank())+
  labs(y = "Count\n", x = "\nProportion")

# all data
ggplot()+
  theme_bw()+
  geom_sf(data = cover.map.version %>% filter(BG >0.02),
          aes(size = BG, colour = Transect))+
  facet_wrap(~Year)

# plot level
plot.bg <- bare.wide.sf %>% 
  group_by(Year, Plot) %>%
  summarise(mean = mean(BG))

# make transect name
plot.bg$Transect <- substr(plot.bg$Plot,1,2)
  

ggplot()+
  theme_bw()+
  geom_sf(data = plot.bg %>% filter(mean >0.02),
          aes(size = mean, colour = Transect), alpha = 0.7)+
  facet_wrap(~Year)+
  scale_colour_manual(values = c("purple", "forestgreen"))+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())
  













