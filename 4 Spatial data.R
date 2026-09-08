library(tidyverse)
library(sf)
library(ggrepel)
library(ggspatial)

convenant <- st_read("Lees_Valley_Convenant_Transects_20191220.gpx")
reserve <-st_read("Lees_Valley_Rec_Reserve_transects_20191220.~dbf.gpx")


ggplot()+
  theme_void()+
  geom_sf(data = convenant, shape = 21, fill = "forestgreen", colour = "white", size = 3 )+
  geom_sf(data = reserve, shape = 21, fill = "hotpink", colour = "white", size  = 3)+
  geom_sf_text(data = convenant, aes(label = name), vjust = 2)+
  geom_text_repel(
    data = reserve,
    aes(label = name, geometry = geometry),
    stat = "sf_coordinates", 
    point.padding = 2,
  )+
  annotation_scale(location = "bl")+
  annotation_north_arrow(location = "bl",
                         pad_x = unit(0.1, "cm"),
                         pad_y = unit(1, "cm"))+
  ggtitle("Pember Fan: plot layout")+
  theme(plot.title = element_text(size = 16, face =2, colour = "grey"))

ggsave("Pember fan layout.png", scale = 1.3, height =7, width =7,
       bg = "white")
