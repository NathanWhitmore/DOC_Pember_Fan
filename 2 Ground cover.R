library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)
library(arm)
library(sf)
library(ordinal)

# read in data
ground <- read.csv("Ground cover.csv")

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

# translate values to means for spatial graphing
L1 <- mean(0:1)  / 100 # 1
L2 <- mean(1:5) / 100 # 2
L3 <- mean(6:25)  / 100 # 3
L4 <- mean(26:50) / 100 # 4
L5 <- mean(51:75)  / 100 # 5
L6 <- mean(76:100) / 100 # 6

# make low value for "P"
ground <-  ground %>%
  mutate(Perc = as.numeric(recode(Cover, 
                                  "1" = L1,
                                  "2" = L2,
                                  "3" = L3,
                                  "4" = L4,
                                  "5" = L5,
                                  "6" = L6))
  )

# make into dataframe
ground <- as.data.frame(ground)

# check levels are appropriate
levels(ground$Cover)

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


# make bg for cover analysis
bg <- ground %>% filter(GroundCover == "BG")
bg$Cover <- as.numeric(as.character(bg$Cover))

# check histogram of cover for bare ground
ggplot()+
  theme_bw()+
  geom_histogram(data = bg, aes(x = Cover), binwidth = 1, 
                 colour = "white", 
                 fill = "orange")+
  facet_grid(~Year)+
  theme(axis.title = element_text(
    size = 14,
    colour = "grey40"
  )) +
  theme(strip.text = element_text(
    size = 12,
    colour = "grey40"
  ))+
  theme(panel.grid = element_blank())+
  ylab("Count\n")+
  xlab("\nCover class")

ggsave("Changes in bare ground cover class.png", scale = 1.1, height = 6, width =8)

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

ggsave("Changes in ground cover.png", scale = 1.1, height = 6, width =8)

# modelling

non.vege <- ground %>% 
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

# insert missing values
bare.ground.wide <- bare.ground %>% 
  pivot_wider(names_from = GroundCover,
              values_from = Prop,
              values_fill = 0)


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


# all data
ggplot()+
  theme_bw()+
  geom_sf(data = cover.map.version %>% filter(BG >0.005),
          aes(size = BG, colour = Transect), alpha = 0.3)+
  facet_wrap(~Year)+
  theme(panel.grid = element_blank())+
  theme(axis.ticks = element_blank())+
  theme(axis.title = element_blank())+
  theme(axis.text = element_blank())+
  scale_colour_manual(values =c("purple", "forestgreen"))+
  labs(size = "Corrected proportion",
       colour = "Transect type")

ggsave("Cover proportion.png", scale =1.1, height =6, width =8)

# plot level
plot.bg <- bare.wide.sf %>% 
  group_by(Year, Plot) %>%
  summarise(mean = mean(BG))

# make transect name
plot.bg$Transect <- substr(plot.bg$Plot,1,2)
  

# tidy data and join (no duplicates present)
bg$Year <- bg$Year - 2018
bg$moniker <- paste(bg$Plot, bg$Subplot, bg$Year)

bg$Year <- NULL
bg$Plot <- NULL
bg$Subplot <- NULL

key <- final.bare.wide[, c("Year", "Plot", "Subplot", "fence.dist")]
key$moniker <- paste(key$Plot, key$Subplot, key$Year)
key <- as.data.frame(key)

# length(unique(key$moniker))
# table(duplicated(key$moniker))

# bare ground cover
bgc <- left_join(key, bg, by = "moniker")

# find duplicated moniker value
bgc$moniker[duplicated(bgc$moniker)]
dup <- bgc %>% filter(moniker == "CA9 240 1")

# insure duplicates are removed
bgc <- bgc %>%
  distinct(moniker, .keep_all = TRUE)

# a few plots seem to missing I'll fill them in
index <- readRDS("Plot index.rds")

bgc <- left_join(index, bgc, by = "moniker")
bgc$Plot.y <- NULL
bgc$Subblot.y <- NULL
bgc$Year.y <- NULL
bgc <- bgc %>% rename(
                      Plot = Plot.x,
                      Subplot = Subplot.x,
                      Year = Year.x)

# ensure 0 is inserted as a cover class
bgc$Cover <- ifelse(is.na(bgc$Cover), 0, bgc$Cover)
bgc$Cover <- factor(bgc$Cover, 
                    levels = c("0", "1", "2", "3", "4", "5", "6"),
                    ordered = TRUE)



# add in all transects subject to
# Among transects where bare ground occurs, 
# how does the amount of bare ground change through time?




# make suplot id because glmm can't handle it otherwise
bgc$SubplotID <- interaction(
  bgc$Plot,
  bgc$Subplot,
  drop = TRUE
)


# cumulative link models

clmm.mod <- list()

clmm.mod[[1]] <- clmm(Cover ~ 1 + (1|Plot) + (1 | SubplotID), data = bgc)
clmm.mod[[2]] <- clmm(Cover ~ Year  + (1|Plot) + (1 | SubplotID), data = bgc)
clmm.mod[[3]] <- clmm(Cover ~ scale(fence.dist) + (1|Plot) + (1 | SubplotID), data = bgc)
clmm.mod[[4]] <- clmm(Cover ~ Year + scale(fence.dist) + (1|Plot) + (1 | SubplotID), data = bgc)
clmm.mod[[5]] <- clmm(Cover ~ as.factor(Year) + (1|Plot) + (1 | SubplotID), data = bgc)
clmm.mod[[6]] <- clmm(Cover ~ as.factor(Year) + scale(fence.dist) + (1|Plot) + (1 | SubplotID),  data = bgc)

# temporal auto correlation - cant really cope with a 
#  Cand.models.hurd[[7]] <- glmmTMB(hurdle ~ as.factor(Year) + ar1(as.factor(Year) + 0 | Plot/Subplot) + (1|Plot/Subplot), 
#                                  family = "binomial", 
#                                   data = all)


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

library(kableExtra)

kable(clmm.aic, "latex")

# summary
summary(clmm.mod[[5]])

# checks - Plot random effect are supported
ranef(clmm.mod[[5]])

# makes two columns for each random effect
my.ranef <- ranef(clmm.mod[[5]]) %>% as.data.frame()
my.ranef$total.random <- rowSums(my.ranef)




plot.subplot <- quantile(my.ranef$total.random, 0.95)

store <- NULL

my.year <- 0:7

for(i in 1:8){
  
  # get betas
  my.year <- c(0, clmm.mod[[5]]$beta)

  
C0 <- plogis(1.4989- ( my.year[i] + plot.subplot))
C1 <- plogis(2.4018  - ( my.year[i] + plot.subplot))
C2 <- plogis(3.3363 - ( my.year[i] + plot.subplot))
C3 <- plogis(5.1506 - ( my.year[i] + plot.subplot))
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
               names_to =  "Cover",
               values_to = "Proportion")

# filter out cover class 0
my.long <-  my.long  |> filter(Cover != "cover.0")

my.long <- my.long %>% 
  mutate(Cover = str_replace_all(Cover, "cover.", ""))
my.long$Cover <- as.numeric(as.character(my.long$Cover))

str(my.long$Cover)

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
  ylab("Percentage in class \n(excluding absences)\n")+
  xlab("\nCover class")

ggsave("Predicted changes in bare ground cover class.png", scale = 1.1, height = 6, width =8)


