library(tidyverse)
library(readxl)
library(randomcoloR)
library(plotly)

# load modified data frames
df <- readRDS("Pember.rds")

# stipulate my.theme
my.theme <- theme_bw()+
  theme(axis.title = element_text(
  face = 2,
  size = 14,
  colour = "grey40"
)) +
  theme(strip.text = element_text(
    size = 12,
    colour = "grey40"
  )) +
  theme(aspect.ratio = 0.5)+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust = 0.5))

# Taxon Biostatus

set.seed(18)
unique.color <- length(unique(df$NVSSpeciesName))
my.colour <- distinctColorPalette(k = unique.color)

# graph all vege
ggplot()+
  theme_bw()+
  geom_col(data = df, aes(x = as.factor(Year), y = Proportion, fill = TaxonBioStatus), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")

ggsave("TaxonBiostatus.png", scale = 1.1, height =6, width =8)

# TaxonGrowthForm
ggplot()+
  theme_bw()+
  geom_col(data = df, aes(x = as.factor(Year), y = Proportion, fill = TaxonGrowthForm), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")

ggsave("TaxonGrowthForm.png", scale = 1.1, height =6, width =8)

# make dual carrier
df$GrowthBio <- paste(df$TaxonGrowthForm, df$Indigenous)

# add in transect
df$Transect <- substring(df$Plot,1,2)

ggplot()+
  theme_bw()+
  geom_col(data = df, aes(x = as.factor(Year), y = Proportion, fill = GrowthBio), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  facet_wrap(~Transect)

ggsave("TaxonGrowthForm by Transect type.png", scale = 1.1, height =6, width =8)





ggplot()+
  theme_bw()+
  geom_col(data = df, aes(x = as.factor(Year), y = Proportion, fill = GrowthBio), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  facet_grid(Transect~Indigenous)

ggsave("TaxonGrowthForm by Transect x Biostatus.png", scale = 1.1, height =6, width =8)



############plotly

# set.seed(18)
# unique.color <- length(unique(sum.all$NVSSpeciesName))
# my.colour <- distinctColorPalette(k = unique.color)

library(plotly)

df$Year <- as.factor(df$Year)

sum.all <- df %>%
  group_by(Year, Transect, Indigenous,TaxonGrowthForm) %>%
  summarise(Total = sum(Proportion)) %>%
  ungroup()

sum.all <- sum.all %>%
  group_by(Year, Transect, Indigenous) %>%
  mutate(Denominator = sum(Total))

sum.all <- sum.all %>%
  mutate(Proportion = round(Total/Denominator,3))

# sum.all <- sum.all %>% arrange(TaxonGrowthForm)

my.plot <- ggplot()+
  theme_bw()+
  geom_col(data = sum.all, aes(x = Year, y = Proportion, 
                               fill = TaxonGrowthForm), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  facet_grid(Transect~Indigenous)+
  theme(legend.position = "none")

# scale_fill_manual(values = my.colour)

ggplotly(my.plot ) #, tooltip = "NVSSpeciesName")


# grasses
grass <- df %>% filter(TaxonGrowthForm  == "Graminoid")
 
grass.all <- grass %>%
  group_by(Year, Transect, NVSSpeciesName) %>%
  summarise(Total = sum(Proportion)) %>%
  ungroup()

grass.all  <- grass.all  %>%
  group_by(Year, Transect) %>%
  mutate(Denominator = sum(Total))

grass.all <- grass.all %>%
  mutate(Proportion = round(Total/Denominator,4))


my.grass <- ggplot()+
  theme_bw()+
  geom_col(data = grass.all, aes(x = Year, y = Proportion, 
                               fill = NVSSpeciesName), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  facet_grid(.~Transect)+
  theme(legend.position = "none")+
  ggtitle("Grasses only")

ggplotly(my.grass ) 

# restricted for ggplot

top8 <- grass.all %>%
  group_by(NVSSpeciesName) %>%
  summarise(total = sum(Proportion, na.rm = TRUE)) %>%
  slice_max(total, n = 8) %>%
  pull(NVSSpeciesName)

grass.all.top8 <- grass.all
grass.all.top8$Species <- ifelse(grass.all.top8$NVSSpeciesName %in% top8,
                                  grass.all.top8$NVSSpeciesName, "Other")


ggplot()+
  theme_bw()+
  geom_col(data = grass.all.top8, aes(x = Year, y = Proportion, 
                                 fill = Species), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  facet_grid(.~Transect)


ggsave("Top 8 grass species.png", scale = 1.1, height =6, width =8)

sort(unique(grass.all$NVSSpeciesName))

# forbs

forbs <- df %>% filter(TaxonGrowthForm  == "Forb")

forbs.all <- forbs %>%
  group_by(Year, Transect, NVSSpeciesName) %>%
  summarise(Total = sum(Proportion)) %>%
  ungroup()

forbs.all  <- forbs.all  %>%
  group_by(Year, Transect) %>%
  mutate(Denominator = sum(Total))

forbs.all <- forbs.all %>%
  mutate(Proportion = round(Total/Denominator,4))

my.forbs <- ggplot()+
  theme_bw()+
  geom_col(data = forbs.all, aes(x = Year, y = Proportion, 
                                 fill = NVSSpeciesName), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  facet_grid(.~Transect)+
  theme(legend.position = "none")+
  ggtitle("Forbs only")

ggplotly(my.forbs ) 

# top 8

sort(unique(forbs$NVSSpeciesName))

top8.forbs <- forbs %>%
  group_by(NVSSpeciesName) %>%
  summarise(total = sum(Proportion, na.rm = TRUE)) %>%
  slice_max(total, n = 8) %>%
  pull(NVSSpeciesName)  


forbs.all.top8 <- forbs.all
forbs.all.top8$Species <- ifelse(forbs.all.top8$NVSSpeciesName %in% top8.forbs,
                                 forbs.all.top8$NVSSpeciesName, "Other")

ggplot()+
  theme_bw()+
  geom_col(data = forbs.all.top8, aes(x = Year, y = Proportion, 
                                      fill = Species), 
           position = "fill")+
  ylab("Proportion\n")+
  xlab("\nMonitoring year")+
  facet_grid(.~Transect)


ggsave("Top 8 forb species.png", scale = 1.1, height =6, width =8)

