library(tidyverse)
library(readxl)
library(randomcoloR)
library(performance)
library(gridExtra)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)


# don't use these just yet
df.2019 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "2019 data")

df.2020.1 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                        sheet = "JAN 2020")

# use these
df.2020.12 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                         sheet = "Dec2020")
df.2021 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "Dec2021")
df.2022 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "Dec2022")
df.2023 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "Dec2023")
df.2024 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "Dec2024")
df.2025 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "Dec2025")
df.2025 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "Dec2025")

# keep only standardised columns
# had to add in ""conv ground cover" column to df.2019 and df.2020jan tabs

df.2019 <- df.2019[,1:10]
df.2020.1<- df.2020.1[,1:10]
df.2020.12 <- df.2020.12[,1:10]
df.2021 <- df.2021[,1:10]
df.2022 <- df.2022[,1:10]
df.2023 <- df.2023[,1:10]
df.2024 <- df.2024[,1:10]
df.2025 <- df.2025[,1:10]

# add in year
df.2019$Year <- "2019"
df.2020.1$Year <- "2020.1"
df.2020.12$Year <- "2020.12"
df.2021$Year <- "2021"
df.2022$Year <- "2022"
df.2023$Year <- "2023"
df.2024$Year <- "2024"
df.2025$Year <- "2025"

# combine
df <- rbind(df.2019, df.2020.1, df.2020.12, df.2021, df.2022, df.2023, df.2024, df.2025)

# rename cover
df <- df %>% rename(Cover = `Rooted inside Ring / Cover class`)

# ground cover@ strip out ground cover
ground.cover <- df %>% 
  filter(!(is.na(GroundCover)))

df <- df %>% 
  filter(is.na(GroundCover))

df$GroundCover <- NULL

# change cover values to perc
# : 1 = <1% cover, 2= 1-5%, 3=6-25%, 4=26-50%, 5=51-75%, 6=76-100%).

# 
P <- 0.05 / 100
L1 <- mean(0:1)  / 100 #p
L2 <- mean(1:5) / 100 # 1
L3 <- mean(6:25)  / 100 # 2
L4 <- mean(26:50) / 100 # 3
L5 <- mean(51:75)  / 100 # 4
L6 <- mean(76:100) / 100 # 5

# expedited removal of errors
df <- df %>% filter(!(Cover %in% c("??", "4?")))

# make low value for "P"
df <-  df %>%
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
df <- as.data.frame(df)

# make subplot factor
df$Subplot <- as.factor(df$Subplot)

# summarise in terms of raw percentage
df <- df %>% 
  group_by(Year, Plot, Subplot) %>%
  mutate(Weight = sum(Perc, na.rm = TRUE))

# corrected perc
df$Proportion <- df$Perc/ df$Weight

## overall
set.seed(18)

# remove unknown species
df.no.unknown <- df %>% filter(NVSSpeciesName != "(Unknown)")

# no of species
df.no.unknown$Bio.simple <- ifelse(df.no.unknown$TaxonBioStatus != "Exotic", 
                                   "Indigenous", "Exotic")

no.species <- df.no.unknown %>%
  group_by(Year, Plot, Subplot, Bio.simple ) %>%
  summarise(Richness = n_distinct(NVSSpeciesName))

no.species$Year <- as.factor(no.species$Year)
no.species$Year <- as.numeric(no.species$Year) -1

# simple correlation  graph
no.species.wide <- no.species %>%
  pivot_wider(names_from = Bio.simple,
              values_from = Richness)


# indig vs exotic richness
indig <- no.species %>% 
  filter(Bio.simple != "Exotic") 

exotic <- no.species %>% 
  filter(Bio.simple == "Exotic") 

no.species.wide <- no.species.wide %>%
  drop_na(Indigenous, Exotic, Year, Plot)

# modelling
# genpois

Models <- list()

Models[[1]] <- glmmTMB(Indigenous ~ 1 + (1|Plot/Subplot), family = genpois(), data = no.species.wide)
Models[[2]] <- glmmTMB(Indigenous ~ Year + (1|Plot/Subplot), family = genpois(), data = no.species.wide)
Models[[3]] <- glmmTMB(Indigenous ~ Exotic + (1|Plot/Subplot), family = genpois(), data = no.species.wide)
Models[[4]] <- glmmTMB(Indigenous ~ Year + Exotic + (1|Plot/Subplot), family = genpois(), data = no.species.wide)  
Models[[5]] <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot/Subplot), family = genpois(), data = no.species.wide)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Models), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Models, formula))))

# AIC table to 4 digits
aictab(cand.set = Models, modnames = Modnames, sort = TRUE)

# summary
summary(Models[[2]])

# annual reduction due to year
1- exp(-0.016758)

# check convergence
Models[[2]]$sdr$pdHess

# random effects meaningful
ranef(Models[[2]])

# diagnostics
res <- simulateResiduals(Models[[2]], n = 1000)
plot(res)

# do we need subplot as a random effects quick and dirty
test <- list()

test[[1]] <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot), family = genpois(), data = no.species.wide, REML = TRUE)
test[[2]] <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot/Subplot), family = genpois(), data = no.species.wide, REML = TRUE)

Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(test, formula))))

aictab(cand.set = test, modnames = Modnames, sort = TRUE)

# THere is no practical need to treat subplot as a random effect
ranef(test[[1]])
ranef(test[[2]])



# make graph
no.species.wide$fitted <- predict(Models[[2]], type = "response")

ggplot()+
  theme_bw()+
  geom_point(data = no.species.wide, aes(y = Indigenous, x = Year, colour = Plot), size =2, alpha = 0.2,
             position = position_jitter(width =0.05, height =0))+
  # geom_point(data = no.species.wide, aes(y = fitted, x = Year, colour = Plot), size = 1)+
  facet_wrap(~Plot)+
  theme(legend.position = 'none')

ggplot()+
  theme_bw()+
  geom_point(data = no.species.wide, 
             aes(x = fitted, y  = Indigenous), size = 2, 
             alpha = 0.2,
             position = position_jitter(width =0, height =0.2))+
  ggtitle("Predictive agreement (all data)")+
  xlab("\nFitted\n\n\n")+
  geom_abline(intercept =0, slope =1, colour = "red")+
  theme(aspect.ratio = 1)



xlab("\nExpected Indigenous")+
  theme(aspect_ratio = 1)+
  theme(panel.grid = element_blank())+
  theme(legend.position = "bottom")




