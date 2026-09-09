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
  group_by(Year ,Plot,  Bio.simple ) %>%
  summarise(Richness = n_distinct(NVSSpeciesName))

no.species$Year <- as.factor(no.species$Year)
no.species$Year <- as.numeric(no.species$Year) -1

# simple correlation  graph
no.species.wide <- no.species %>%
  pivot_wider(names_from = Bio.simple,
              values_from = Richness)

ggplot()+
  theme_bw()+
  #geom_point(data = no.species.wide, aes(y = Exotic, x = Indigenous), size = 5, 
  #           position = position_jitter(width =0.1),
  #          alpha = 0.5)+
  geom_smooth(data = no.species.wide, aes(y = Exotic, x = Indigenous, colour = as.factor(Year)), 
              method = lm, fill = NA)

# ratio 
no.species.wide$indig.ratio <- no.species.wide$Indigenous / (no.species.wide$Indigenous + no.species.wide$Exotic)

ggplot()+
  geom_line(data = no.species.wide, aes(y = indig.ratio , 
                                        x = as.factor(Year)),
            alpha = 0.5)


ggplot()+
  geom_boxplot(data = no.species.wide, aes(y = indig.ratio , 
                                           x = as.factor(Year)),
               alpha = 0.5)

# indig vs exotic richness
indig <- no.species %>% 
  filter(Bio.simple != "Exotic") 


exotic <- no.species %>% 
  filter(Bio.simple == "Exotic") 

no.species.wide <- no.species.wide %>%
  drop_na(Indigenous, Exotic, Year, Plot)



m.comp <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot), family = compois(), data = no.species.wide)
m.pois <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot), family = "poisson", data = no.species.wide)
m.gen.pois <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot), family = genpois(), data = no.species.wide)

summary(m.gen.pois)
summary(m.comp)

# check
AIC(m.pois, m.comp, m.gen.pois)

Cand.models <- list()

Cand.models[[1]] <- glmmTMB(Indigenous ~ 1 + (1|Plot), family = compois(), data = no.species.wide)
Cand.models[[2]] <- glmmTMB(Indigenous ~ Year + (1|Plot), family = compois(), data = no.species.wide)
Cand.models[[3]] <- glmmTMB(Indigenous ~ Exotic + (1|Plot), family = compois(), data = no.species.wide)
Cand.models[[4]] <- glmmTMB(Indigenous ~ Year + Exotic + (1|Plot), family = compois(), data = no.species.wide)  
Cand.models[[5]] <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot), family = compois(), data = no.species.wide)
# Cand.models[[6]] <- glmmTMB(Indigenous ~ Year * poly(Exotic,2) + (1|Plot), family = compois(), data = no.species.wide)

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models, formula))))

# AIC table to 4 digits
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

# summary
summary(Cand.models[[5]])

# check diagnostics
res <- simulateResiduals(Cand.models[[5]], n = 1000)
plot(res)

testDispersion(res)
testUniformity(res)
testZeroInflation(res)

# this looks fine
plotResiduals(res, model.frame(m.gen.pois)$Exotic)
plotResiduals(res, model.frame(m.gen.pois)$Year)
plotResiduals(res, model.frame(m.gen.pois)$Plot)

# model performance
model_performance(Cand.models[[5]])

# predictions

no.species.wide$fitted <- predict(Cand.models[[5]], type = "response")




A <- ggplot()+
  theme_bw()+
  #geom_point(data = no.species.wide, aes(y = Exotic, x = Indigenous), size = 5, 
  #           position = position_jitter(width =0.1),
  #          alpha = 0.5)+
  geom_smooth(data = no.species.wide, aes(y = Exotic, x = Indigenous, colour = as.factor(Year)), 
              method = lm, fill = NA) +
  ggtitle("Real trajectory")+
  scale_x_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 18))+
  scale_y_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 9))


B <- ggplot()+
  theme_bw()+
  #geom_point(data = no.species.wide, aes(y = Exotic, x = Indigenous), size = 5, 
  #           position = position_jitter(width =0.1),
  #          alpha = 0.5)+
  geom_smooth(data = no.species.wide, aes(y = Exotic, x = fitted, colour = as.factor(Year)), 
              method = lm, fill = NA) +
  ggtitle("Modelled trajectory")+
  scale_x_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 18))+
  scale_y_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 9))+
  xlab("Expected Indigenous")



grid.arrange(A, B, ncol = 2)

# make example

summary(Cand.models[[6]])

example <- no.species.wide %>% filter(Plot == "RR2" | Plot == "CA1"  | Plot == "RR7"| Plot == "CA15")

example$fitted <- predict(Cand.models[[5]], newdata = example, type = "response")
manual <- exp(2.383900 + -0.006249 * example$Year + 0.015851 * example$Exotic -0.005358 * (example$Year * example$Exotic))


C <- ggplot()+
  theme_bw()+
  geom_point(data = example, aes(y = Exotic, x = manual, colour = Plot), size = 7)+
  geom_path(data = example, aes(y = Exotic, x = manual, group = Year))+
  geom_text(data = example, aes(y = Exotic, x = manual, label = Year))+
  ggtitle("Modelled trajectory")+
  #scale_x_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 18))+
  #scale_y_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 9))+
  xlab("\nExpected Indigenous")+
  theme(aspect_ratio = 1)+
  theme(panel.grid = element_blank())+
  theme(legend.position = "bottom")


D <- ggplot()+
  theme_bw()+
  geom_point(data = example, aes(y = Exotic, x = Indigenous, colour = Plot), size = 9)+
  geom_text(data = example, aes(y = Exotic, x = Indigenous, label = Year), size =5 )+
  ggtitle("Raw data")+
  #scale_x_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 18))+
  #scale_y_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 9))+
  xlab("\nExpected Indigenous")+
  theme(aspect_ratio = 1)+
  theme(panel.grid = element_blank())+
  theme(legend.position = "bottom")

E <- ggplot()+
  theme_bw()+
  geom_point(data = no.species.wide, 
             aes(x = fitted, y  = Indigenous), size = 5, alpha = 0.2)+
  ggtitle("Predictive agreement (all data)")+
  xlab("\nFitted\n\n\n")+
  scale_x_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 18))+
  scale_y_continuous(breaks = 4:20, labels = 4:20,  limits = c(3, 18))+
  geom_abline(intercept =0, slope =1, colour = "red")

# quick model comparison
plot.a <- grid.arrange(D, C,  ncol = 2)

ggsave(
  "Biodiversity plot raw and modelled.png",
  plot = plot.a,
  scale = 1.1,
  height = 6,
  width = 9,
  dpi = 300
)

plot.b <- grid.arrange(C, E,  ncol = 2)

ggsave(
  "Biodiversity plot modelled and correlation.png",
  plot = plot.b,
  scale = 1.1,
  height = 6,
  width = 9,
  dpi = 300
)

# correlation (predictive agreement)
mod1 <- lm(Indigenous ~ fitted, data = no.species.wide)
summary(mod1)