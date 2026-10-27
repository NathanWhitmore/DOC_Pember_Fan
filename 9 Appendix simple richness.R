library(tidyverse)
library(readxl)
library(randomcoloR)
library(performance)
library(gridExtra)
library(AICcmodavg)
library(glmmTMB)
library(DHARMa)


# load modified data frames
df <- readRDS("Pember cover.rds")

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

# ratio 
no.species.wide$indig.ratio <- no.species.wide$Indigenous / (no.species.wide$Indigenous + no.species.wide$Exotic)

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

Cand.models[[1]] <- glmmTMB(Indigenous ~ 1 + (1|Plot), family = genpois(), data = no.species.wide)
Cand.models[[2]] <- glmmTMB(Indigenous ~ Year + (1|Plot), family = genpois(), data = no.species.wide)
Cand.models[[3]] <- glmmTMB(Indigenous ~ Exotic + (1|Plot), family = genpois(), data = no.species.wide)
Cand.models[[4]] <- glmmTMB(Indigenous ~ Year + Exotic + (1|Plot), family = genpois(), data = no.species.wide)  
Cand.models[[5]] <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot), family = genpois(), data = no.species.wide)
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

# save bioidiversity data
saveRDS(no.species.wide, "no.species.wide.rds")


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

summary(Cand.models[[5]])

example <- no.species.wide %>% filter(Plot == "RR2" | Plot == "CA1"  | Plot == "RR7"| Plot == "CA15")

example$fitted <- predict(Cand.models[[5]], newdata = example, type = "response")
# manual <- exp(2.383900 + -0.006249 * example$Year + 0.015851 * example$Exotic -0.005358 * (example$Year * example$Exotic))


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