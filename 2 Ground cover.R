library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)

# read in data
ground <- read.csv("Ground cover2.csv")

# make into data.frame
ground <- as.data.frame(ground)

# make subplot factor
ground$Subplot <- as.factor(ground$Subplot)

# make low value for "P"
ground <-  ground %>%
  mutate(Perc = as.numeric(recode(Cover, 
                                  "P" = P,
                                  "1" = L1,
                                  "2" = L2,
                                  "3" = L3,
                                  "4" = L4,
                                  "5" = L5,
                                  "6" = L6))
  )

# make into dat.frame
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

# exploratory graphing 
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
    "BG" = "BR",
    "V" = "Li",
    "V" = "L",
    "V" = "M",
  )
)

# claulate sub plot proportions
bare.ground <- non.vege %>%
  ungroup() %>%
  group_by(Year, Plot, Subplot, GroundCover) %>%
  summarise(Prop = sum(Proportion, na.rm = TRUE))


bare.ground.wide <- bare.ground %>% 
  pivot_wider(names_from = GroundCover,
              values_from = Prop,
              values_fill = 0.001)

# make sure there are no zeros or 1s
bare.ground.wide$BG <- ifelse(bare.ground.wide$BG == 0, 0.001, bare.ground.wide$BG)
bare.ground.wide$BG <- ifelse(bare.ground.wide$BG == 1, 0.999, bare.ground.wide$BG)

# read in spatial data
my.coord.sf <- readRDS("my_coord_sf.rds")

# join 
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

st_geometry(final.bare.wide) <- NULL

# model selection

Cand.models.sub <- list()

Cand.models.sub[[1]] <- glmmTMB(BG ~ 1 + (1|Plot/Subplot), family = beta_family(link = "logit"), 
                                data = final.bare.wide)
Cand.models.sub[[2]] <- glmmTMB(BG ~ Year  + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = final.bare.wide)
Cand.models.sub[[3]] <- glmmTMB(BG ~ as.factor(Year) + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = final.bare.wide)
Cand.models.sub[[4]] <- glmmTMB(BG ~ fence.dist + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = final.bare.wide)
Cand.models.sub[[5]] <- glmmTMB(BG ~ Year + fence.dist + (1|Plot/Subplot), family = beta_family(link = "logit"),
                                data = final.bare.wide)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models.sub), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models.sub, formula))))

# AIC table to 4 digits
w.subplot <- aictab(cand.set = Cand.models.sub, modnames = Modnames, sort = TRUE)
w.subplot

# summary
summary(Cand.models.sub[[3]])




