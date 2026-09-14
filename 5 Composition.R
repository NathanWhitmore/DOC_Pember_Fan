library(tidyverse)
library(readxl)
library(randomcoloR)

# don't use these just yet
df.2018 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                      sheet = "2019 data")

df.2019 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
                        sheet = "JAN 2020")

# use these

df.2020 <- read_excel("Copy of NC_LEESVALLEY_MONITORING_DEC2025.xlsx",
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


df.2018 <- df.2018[,1:10]
df.2019 <- df.2019[,1:10]
df.2020 <- df.2020[,1:10]
df.2021 <- df.2021[,1:10]
df.2022 <- df.2022[,1:10]
df.2023 <- df.2023[,1:10]
df.2024 <- df.2024[,1:10]
df.2025 <- df.2025[,1:10]

# add in year
df.2018$Year <- "2018"
df.2019$Year <- "2019"
df.2020$Year <- "2020"
df.2021$Year <- "2021"
df.2022$Year <- "2022"
df.2023$Year <- "2023"
df.2024$Year <- "2024"
df.2025$Year <- "2025"

# combine
df <- rbind(df.2018, df.2019, df.2020, df.2021, df.2022, df.2023, df.2024, df.2025)

# rename cover
df <- df %>% rename(Cover = `Rooted inside Ring / Cover class`)

# ground cover@ strip out ground cover
df <- df %>% 
  filter(is.na(GroundCover))


# change cover values to perc
# : 1 = <1% cover, 2= 1-5%, 3=6-25%, 4=26-50%, 5=51-75%, 6=76-100%).

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

df <- as.data.frame(df)
df$GroundCover <- NULL
df$`Overhanging Ring` <- NULL

# corrected perc
df$Proportion <- df$Perc/ df$Weight
df$Proportion  <- ifelse(is.na(df$Proportion), 0, df$Proportion)
df$Perc <- NULL
df$Weight <- NULL
df$Cover <- NULL

# Indigenous comparison
df$Indigenous <- ifelse(df$TaxonBioStatus == "Exotic", "Exotic", "Indigenous")

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
  facet_grid(Indigenous~TaxonGrowthForm)+
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
