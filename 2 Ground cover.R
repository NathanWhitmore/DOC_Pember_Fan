library(tidyverse)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)
library(glmmTMB)

# read in data
ground <- read.csv("Ground cover.csv")

# make into dat.frame
ground <- as.data.frame(ground)

# make subplot factor
ground$Subplot <- as.factor(ground$Subplot)

# summarise by cover
ground <- ground %>% 
  group_by(Year, Plot, Subplot) %>%
  mutate(Weight = sum(Perc, na.rm = TRUE))

# corrected perc
ground$Proportion <- ground$Perc/ ground$Weight

# just one transect as proof of concept
one <- ground %>% 
  filter(Plot == "RR2") %>%
  as.data.frame()

ggplot()+
  geom_col(data = ground, aes(x = Subplot, y = Proportion,
                           fill = GroundCover))+
  facet_grid(Year~.)

ground <- as.data.frame(ground)


ground <- ground %>%
  mutate(
    GroundCover = fct_recode(
      GroundCover,
      "BG" = "BG2",
      "BG" = "BG3",
      "Li" = "Li3"
    )
  )

# modelling

ground$Year  <- 2019 - ground$Year 

ground <- ground %>% drop_na(Proportion)

ground$Proportion <- ifelse(ground$Proportion == 1, 0.999, ground$Proportion)

Cand.models <- list()

Cand.models[[1]] <- glmmTMB(Proportion ~ 1 + GroundCover + (1|Plot), family = beta_family(link = "logit"), data = ground)
Cand.models[[2]] <- glmmTMB(Proportion ~ Year * GroundCover + (1|Plot), family = beta_family(link = "logit"), data = ground)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models, formula))))

# AIC table to 4 digits
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)


summary(Cand.models[[2]])  
