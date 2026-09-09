library(tidyverse)
library(readxl)
library(randomcoloR)

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

# summarise by cover
df <- df %>% 
  group_by(Year, Plot, Subplot) %>%
  mutate(Weight = sum(Perc, na.rm = TRUE))

# corrected perc
df$Proportion <- df$Perc/ df$Weight

# filter for types
# Just forbs (as an example)
unique(df$TaxonGrowthForm)
type <- df %>% filter(TaxonGrowthForm == "Forb") 

# palette
set.seed(18)
unique.color <- length(unique(type$NVSSpeciesName))
my.colour <- distinctColorPalette(k = unique.color)

# graph
ggplot()+
  theme_bw()+
  geom_col(data = type, aes(x = as.factor(Year), y = Proportion, fill = NVSSpeciesName), 
           position = "fill")+
  facet_grid(TaxonBioStatus~TaxonGrowthForm)+
  scale_fill_manual(values = my.colour)+
  # scale_x_continuous(breaks = 2020:2025, labels =  2020:2025) +
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
  guides(fill = guide_legend(ncol = 2))+
  labs(fill = "")+
  theme(aspect.ratio = 0.5)+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0.5))


## overall
set.seed(18)

df.no.unknown <- df %>% filter(NVSSpeciesName != "(Unknown)"&
                                 TaxonBioStatus != "Unknown")

unique.color <- length(unique(df.no.unknown$NVSSpeciesName))
my.colour <- distinctColorPalette(k = unique.color)

ggplot()+
  theme_bw()+
  geom_col(data = df.no.unknown, aes(x = Year, y = Proportion, fill = NVSSpeciesName)
           # , position = "fill"
  )+
  scale_fill_manual(values = my.colour)+
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
  guides(fill = guide_legend(ncol = 3))+
  labs(fill = "")+
  facet_grid(.~TaxonBioStatus)

ggsave("Overall composition.png", scale = 1.2, heigh =6, width =8)


# simple biostatus

names(df.no.unknown)

ggplot()+
  theme_bw()+
  geom_col(data = df.no.unknown, aes(x = Year, y = Proportion, fill = TaxonBioStatus)
           #, position = "fill"
  )+
  scale_fill_manual(values = my.colour)+
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
  labs(fill = "Biostatus")


# top 5 in each category

names(df)

top20 <- df.no.unknown  %>% 
  group_by(TaxonBioStatus,TaxonGrowthForm,  NVSSpeciesName) %>%
  summarise(Overall.prop = sum(Proportion, na.rm = TRUE)) %>% 
  slice_max(order_by = Overall.prop, n = 5) %>%
  ungroup() %>%
  arrange(Overall.prop) %>%
  slice_max(order_by = Overall.prop, n = 20)

# get trajectory of major species
major.species <- unique(top20$NVSSpeciesName)

major <- df.no.unknown %>% 
  filter(NVSSpeciesName %in% major.species)

unique.color <- length(unique(major$NVSSpeciesName))
my.colour <- distinctColorPalette(k = unique.color)


ggplot()+
  theme_bw()+
  geom_col(data = major, aes(x = Year, y = Proportion, fill = NVSSpeciesName), 
           position = "fill")+
  scale_fill_manual(values = my.colour)+
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
  guides(fill = guide_legend(ncol = 2))+
  labs(fill = "")+
  facet_grid(.~TaxonBioStatus)

# Focus on "Brachyscome pinnata" 
brac <- df.no.unknown  %>% 
  filter(NVSSpeciesName == "Brachyscome pinnata")

# is this species associated with higher biodiversity?

hist(brac$Weight)

# look for Sonchus novae zealandeae.


ggplot() +
  theme_bw()+
  geom_point(data = brac, aes(x = Year, y = Weight), 
             position = position_jitter(width = 0.1))+
  geom_boxplot(data = brac, aes(x = Year, y = Weight), fill = NA, colour = "purple")


# no of species
names(df.no.unknown)

df.no.unknown$Bio.simple <- ifelse(df.no.unknown$TaxonBioStatus != "Exotic", 
                                   "Indigenous", "Exotic")

unique(df.no.unknown$TaxonBioStatus)

no.species <- df.no.unknown %>%
  group_by(Year,Plot, Subplot, Bio.simple ) %>%
  summarise(Richness = n_distinct(NVSSpeciesName))

no.species$Subplot <- as.factor(no.species$Subplot)
no.species$Year <- as.factor(no.species$Year)
no.species$Year <- as.numeric(no.species$Year)

indig <- no.species %>% 
  filter(Bio.simple != "Exotic") %>% 
  filter(Plot == "RR2") %>%
  ungroup()

str(indig)
ggplot()+
  geom_point(data = indig, aes(x = Subplot, y = Richness),
             position = position_jitter(width =0.01, height =0.05), 
             alpha =0.5)+
  geom_smooth(data = indig, aes(x = as.numeric(Subplot), y = Richness))

library(AICcmodavg)
library(glmmTMB)
library(DHARMa)

Cand.models <- list()

Cand.models[[1]] <- glmmTMB(Richness ~ 1 + (1|Plot/Subplot), family = compois(), data = no.species)
Cand.models[[2]] <- glmmTMB(Richness ~ Year + (1|Plot/Subplot), family = compois(), data = no.species)
Cand.models[[3]] <- glmmTMB(Richness ~ Bio.simple + (1|Plot/ Subplot), family = compois(), data = no.species)
Cand.models[[4]] <- glmmTMB(Richness ~ Year + Bio.simple + (1|Plot/ Subplot), family = compois(), data = no.species)  
Cand.models[[5]] <- glmmTMB(Richness ~ Year * Bio.simple + (1|Plot/Subplot), family = compois(), data = no.species)  

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Cand.models, formula))))

# AIC table to 4 digits
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

# check residuals
res <- simulateResiduals(Cand.models[[5]])
plot(res)


testOverdispersion(res)
testDispersion(res)

# predictions

new.data <- na.omit(no.species)
new.data$pred <- predict(Cand.models[[5]], newdata = new.data, type = "response")

ggplot()+
  geom_point(data = new.data, aes(x= Richness, y= pred))+
  geom_abline(slope = 1, intercept = 0, colour = "red", lwd =1)+
  theme(aspect.ratio = 1)


summary(Cand.models[[5]])
