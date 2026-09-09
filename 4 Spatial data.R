library(tidyverse)
library(sf)
library(ggrepel)
library(ggspatial)
library(sf)
library(glmmTMB)

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

# combine all data

my.spatial <- rbind(convenant,reserve)
my.spatial <- my.spatial  %>% st_transform(crs = 2193)

# spatial covariance
my.coord <- st_coordinates(my.spatial) %>% as.data.frame()
my.coord$Plot <- my.spatial$name

my.coord  <- my.coord  |> 
  mutate(pos = numFactor(X, Y), spatial_group = factor(1))

my.data <- left_join(no.species.wide, my.coord, by = "Plot") %>% as.data.frame()
str(my.data )

m.test1 <- glmmTMB( Indigenous ~ Year * Exotic + (1 | Plot/Subplot) + exp(pos + 0 | spatial_group),
                    family = genpois(),
                    data = my.data
)

m.test2 <- glmmTMB( Indigenous ~ Year + Exotic + (1 | Plot) + exp(pos + 0 | spatial_group),
                    family = genpois(),
                    data = my.data
)


diagnose(m.test2)

m.plot <- glmmTMB(
  Indigenous ~ Year + Exotic +
    (1 | Plot),
  family = genpois(),
  data = my.data
)

m.spatial <- glmmTMB(
  Indigenous ~ Year + Exotic +
    exp(pos + 0 | spatial_group),
  family = genpois(),
  data = my.data
)

m.spatial <- glmmTMB(
  Indigenous ~ Year + Exotic +
    exp(pos + 0 | spatial_group),
  family = genpois(),
  data = my.data
)

m.spatial2 <- glmmTMB(
  Indigenous ~ Year * Exotic +
    exp(pos + 0 | spatial_group)+
    (1 | Plot:Subplot),
  family = genpois(),
  data = my.data
)


VarCorr(m.spatial2)

diagnose(m.plot)
diagnose(m.spatial)

VarCorr(m.plot)
VarCorr(m.spatial)


VarCorr(m.test2)

AIC(m.plot, m.spatial)
logLik(m.plot)
logLik(m.spatial)


# Is position effectively equivalent to Plot?
my.data |>
  distinct(Plot, pos) |>
  count(Plot)


m.spatial2$fit$par

summary(m.spatial2)


m <- glmmTMB(
  Indigenous ~ Year + (1 | Plot/Subplot),
  family = genpois(),
  data = my.data
)


plot.re <- data.frame(
  Plot = rownames(re),
  effect = re[, "(Intercept)"]
)

coords <- my.data |>
  distinct(Plot, X, Y)

# check there are 30 rows
nrow(coords)

# rejoin
plot.re <- left_join(plot.re, coords, by = "Plot")

library(spdep)

# make distnace matrix
xy <- cbind(plot.re$X, plot.re$Y)
d <- as.matrix(dist(xy))

# check non-zero range
range(d[d > 0])

# find nearest neighbour within 500m
nb <- dnearneigh(
  xy,
  d1 = 0,
  d2 = 500
)

# summary of nearest neighbour
summary(nb)

# spatial wieghts for Moran test
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# significant
moran.test(
  plot.re$effect,
  lw,
  zero.policy = TRUE
)

# is this significant across different number of neighbours?

for(k in 2:6) {
  
  nb <- knn2nb(knearneigh(xy, k = k))
  lw <- nb2listw(nb, style = "W")
  
  test <- moran.mc(
    plot.re$effect,
    lw,
    nsim = 9999
  )
  
  cat(
    "k =", k,
    " Moran I =", test$statistic,
    " p =", test$p.value, "\n"
  )
}


# random effects
plot.re

# map random effects

my.coord.sf <- my.coord %>% st_as_sf(coords = c("X", "Y"), crs = 2193)

viz <- left_join(my.coord.sf, plot.re, by = "Plot")
viz$Group <- substr(viz$Plot, 1, 2)

# read in rivers
river.lines <- st_read("river/nz-river-name-lines-pilot.shp")
river.polygons <- st_read("river/nz-river-polygons-topo-150k.shp")

adj.viz <- st_bbox(  c(xmin=  1533588 ,
          ymin = 5224676 ,
          xmax = 1536088 + 600, 
          ymax  = (5227264 + 200)))

river.lines.cropped <- st_crop(river.lines, st_bbox(adj.viz))
river.polygons.cropped <- st_crop(river.polygons, st_bbox(adj.viz))

ggplot()+
  geom_sf(data = river)

ggplot()+
  theme_void()+
  geom_sf(data = river.lines.cropped)+
  geom_sf(data = river.polygons.cropped)+
  geom_sf(data = viz, aes(colour =effect, fill = effect, shape = Group), colour = "grey", size =5)+
  scale_fill_gradient(low = "yellow", high = "purple")+
  scale_shape_manual(values = c(21, 23))+
  annotation_scale(location = "bl")+
  annotation_north_arrow(location = "bl",
                         pad_x = unit(0.1, "cm"),
                         pad_y = unit(1, "cm"))+
  ggtitle("Pember fan: random effecst for biodiversity appear spatially autocorrelated")+
  theme(plot.title = element_text(size = 16, face =2, colour = "grey"))+
  labs(fill = "Random effect size", 
       shape = "Plot type",
       subtitle = "model: indig.biodiv ~ Year + (1|Plot/Subplot), family = genpois()")

ggsave("Pember fan random effect spatial signal.png",  
       scale = 1.2, width =9, height = 6,
       bg = "white")

