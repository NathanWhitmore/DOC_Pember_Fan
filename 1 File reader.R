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

# Indigenous comparison
df$Indigenous <- ifelse(df$TaxonBioStatus == "Exotic", "Exotic", "Indigenous")


### Ground cover extraction
# ground cover@ strip out ground cover
ground.cover <- df %>% 
  filter(!(is.na(GroundCover)))

write.csv(ground.cover, "Ground cover.csv", row.names = FALSE)


### Extractions for presence and cover

# ground cover@ strip out ground cover
df <- df %>% 
  filter(is.na(GroundCover))


# change cover values to approximate % solely for graphing purposes
L1 <- mean(0:1)  / 100 #p
L2 <- mean(1:5) / 100 # 1
L3 <- mean(6:25)  / 100 # 2
L4 <- mean(26:50) / 100 # 3
L5 <- mean(51:75)  / 100 # 4
L6 <- mean(76:100) / 100 # 5

# expedited removal of errors
df <- df %>% filter(!(Cover %in% c("??")))
df <- df %>%   mutate(
  Cover = str_replace_all(Cover, fixed("4?"), "4")
)


### Cover extraction
# create presence
df$Presence <- ifelse(df$Cover %in% c("1", "2", "3", "4", "5", "6", "P"),
                      1, 0)

# make distinct cover dataset
df.cover <-  df %>%
  mutate(Perc = as.numeric(recode(Cover,
                                  "P" = 0,
                                  "1" = L1,
                                  "2" = L2,
                                  "3" = L3,
                                  "4" = L4,
                                  "5" = L5,
                                  "6" = L6))
  )

# make subplot factor
df.cover$Subplot <- as.factor(df.cover$Subplot)

# save into memory
saveRDS(df.cover, "Pember cover.rds")


### Quasi-percentage extraction
# summarise by cover
df <- df.cover %>% 
  group_by(Year, Plot, Subplot) %>%
  mutate(Weight = sum(Perc, na.rm = TRUE))

df <- as.data.frame(df)

df$GroundCover <- NULL

# df$`Overhanging Ring` <- NULL

# corrected perc
df$Proportion <- df$Perc/ df$Weight
df$Proportion  <- ifelse(is.na(df$Proportion), 0, df$Proportion)
df$Perc <- NULL
df$Weight <- NULL

# save into memory
saveRDS(df, "Pember presence.rds")

#### Plot index for ensuring zeros aren't missed
Plot <- unique(df$Plot)
Subplot <- unique(df$Subplot)
Year <- as.numeric(unique(df$Year)) - 2018

# expand grid
index <- expand_grid(Plot, Subplot, Year)
index$moniker <- paste(index$Plot, index$Subplot, index$Year)

# save into memory
saveRDS(index, "Plot index.rds")

# spatial data
convenant <- st_read("Lees_Valley_Convenant_Transects_20191220.gpx")
reserve <-st_read("Lees_Valley_Rec_Reserve_transects_20191220.~dbf.gpx")

my.spatial <- rbind(convenant,reserve)
my.spatial <- my.spatial  |> st_transform(crs = 2193)

my.spatial <- my.spatial |> rename(Plot = name)
my.spatial <- my.spatial %>% select(Plot, geometry)

# read river and fence distances 
ashley <- st_read("Ashley river boundary.kml")
ashley <- st_union(ashley) %>% st_transform(crs = 2193)

fence <- st_read("fence line/Fence_line.json")
fence <- fence  %>% st_transform(crs = 2193)

# make sf base
# my.data.sf <- st_as_sf(my.spatial, coords =c("X", "Y"), crs = 2193)

# incorporate fence line 

my.spatial$river.dist <- as.numeric(
  st_distance(my.spatial, ashley )
)

my.spatial$fence.dist <- as.numeric(
  st_distance(my.spatial, fence )
)

# save into memory
saveRDS(my.spatial, "Spatial.rds")


