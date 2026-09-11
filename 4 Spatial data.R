library(tidyverse)
library(sf)
library(ggrepel)
library(ggspatial)
library(sf)
library(glmmTMB)
library(spdep)

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

# m.test1 <- glmmTMB( Indigenous ~ Year * Exotic + (1 | Plot/Subplot) + exp(pos + 0 | spatial_group),
#                    family = genpois(),
#                    data = my.data
#)


# checking spatial distribution of random effects

m <- glmmTMB(
  Indigenous ~ Year + (1 | Plot/Subplot),
  family = genpois(),
  data = my.data
)

re <- ranef(m)$cond$Plot %>% as.data.frame()

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

# make distance matrix
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

adj.viz <- st_bbox(  c(xmin =  1533588 ,
          ymin = 5224676 ,
          xmax = 1536088 + 600, 
          ymax  = (5227264 + 200)))

river.lines.cropped <- st_crop(river.lines, st_bbox(adj.viz))
river.polygons.cropped <- st_crop(river.polygons, st_bbox(adj.viz))

# read kml

ashley <- st_read("Ashley river boundary.kml")
fence <- st_read("fence line/Fence_line.json")

ggplot()+
  theme_void()+
  geom_sf(data = river.lines.cropped)+
  geom_sf(data = river.polygons.cropped)+
 # geom_sf(data = ashley, colour = "blue", lwd =2)+
 # geom_sf(data = fence, colour = "red")+
  geom_sf(data = viz, aes(colour =effect, fill = effect, shape = Group), colour = "grey", size =5)+
  scale_fill_gradient(low = "yellow", high = "purple")+
  scale_shape_manual(values = c(21, 23))+
  annotation_scale(location = "bl")+
  annotation_north_arrow(location = "bl",
                         pad_x = unit(0.1, "cm"),
                         pad_y = unit(1, "cm"))+
  ggtitle("Pember fan: random effects for biodiversity appear spatially autocorrelated")+
  theme(plot.title = element_text(size = 16, face =2, colour = "grey"))+
  labs(fill = "Random effect size", 
       shape = "Plot type",
       subtitle = "model: Indigenous ~ Year + (1|Plot/Subplot), family = genpois()")

ggsave("Pember fan random effect spatial signal.png",  
       scale = 1.2, width =9, height = 6,
       bg = "white")


# incorporate ashley river distance 
my.data.sf <- st_as_sf(my.data, coords =c("X", "Y"), crs = 2193)
ashley <- st_union(ashley) %>% st_transform(crs = 2193)

# incorporate fence line 
fence <- fence  %>% st_transform(crs = 2193)

my.data$river.dist <- as.numeric(
  st_distance(my.data.sf, ashley )
)

my.data$fence.dist <- as.numeric(
  st_distance(my.data.sf, fence )
)

# model selection

Models <- list()

Models[[1]] <- glmmTMB(Indigenous ~ 1 + (1|Plot/Subplot), family = genpois(), data = my.data)
Models[[2]] <- glmmTMB(Indigenous ~ Year + (1|Plot/Subplot), family = genpois(), data = my.data)
Models[[3]] <- glmmTMB(Indigenous ~ Exotic + (1|Plot/Subplot), family = genpois(), data = my.data)
Models[[4]] <- glmmTMB(Indigenous ~ Year + Exotic + (1|Plot/Subplot), family = genpois(), data = my.data)  
Models[[5]] <- glmmTMB(Indigenous ~ Year * Exotic + (1|Plot/Subplot), family = genpois(), data = my.data)
Models[[6]] <- glmmTMB(Indigenous ~ Year + scale(river.dist) + (1|Plot/Subplot), family = genpois(), data = my.data) 
Models[[7]] <- glmmTMB(Indigenous ~ Year + scale(fence.dist) + (1|Plot/Subplot), family = genpois(), data = my.data)
Models[[8]] <- glmmTMB(Indigenous ~ scale(river.dist) + (1|Plot/Subplot), family = genpois(), data = my.data)
Models[[9]] <- glmmTMB(Indigenous ~ scale(fence.dist) + (1|Plot/Subplot), family = genpois(), data = my.data)
Models[[10]] <- glmmTMB(Indigenous ~ Year * scale(fence.dist) + (1|Plot/Subplot), family = genpois(), data = my.data)
# Models[[7]] <- glmmTMB(Indigenous ~ Year + Exotic + scale(river_dist) + (1|Plot/Subplot), family = genpois(), data = my.data)
#Models[[8]] <- glmmTMB(Indigenous ~ Year * scale(river_dist) + (1|Plot/Subplot), family = genpois(), data = my.data) 

# create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Models), sep = " ")
Modnames <- paste(sub(".*formula =*(.*?) *, .*", "\\1", 
                      unlist(lapply(Models, formula))))

# AIC table to 4 digits
aictab(cand.set = Models, modnames = Modnames, sort = TRUE)

# summary
summary(Models[[10]])

# fitted
my.data$fitted <- fitted(Models[[10]])

ggplot()+ 
  geom_violin(data = my.data, aes(y = Indigenous, x = as.factor(round(fitted ))))+
  geom_abline(intercept = 0, slope = 1, colour = "red", lwd = 2)


# simplified gen vs com
gen <- glmmTMB(Indigenous ~ Year + scale(fence.dist) + (1|Plot), family = genpois(), data = my.data)
com <- glmmTMB(Indigenous ~ Year + scale(fence.dist) + (1|Plot), family = compois(), data = my.data)

# diagnostics
res <- simulateResiduals(Models[[10]])
plot(res)
res.gen <- simulateResiduals(gen)
res.com <- simulateResiduals(com)

# test dispersion compois performs much better
testDispersion(res)
testDispersion(res.gen)
testDispersion(res.com)

# test K-S test (ideally want to see 0) compois performs much better
testUniformity(res)
testUniformity(res.gen)
testUniformity(res.com)

# are subplot random effects warranted
rand.without <- glmmTMB(Indigenous ~ Year * scale(fence.dist) + (1|Plot), family = genpois(), data = my.data, REML = TRUE)
rand.with <- glmmTMB(Indigenous ~ Year * scale(fence.dist) + (1|Plot/Subplot), family = genpois(), data = my.data, REML = TRUE)

# subplots are totally justified
AIC(rand.without, rand.with)


# new random effects
m <- glmmTMB(
  Indigenous ~ Year * scale(fence.dist) + (1 | Plot/Subplot),
  family = genpois(),
  data = my.data
)

re <- ranef(m)$cond$Plot %>% as.data.frame()

plot.re <- data.frame(
  Plot = rownames(re),
  effect = re[, "(Intercept)"]
)

coords <- my.data %>%
  distinct(Plot, X, Y)

# check there are 30 rows
nrow(coords)

# rejoin
plot.re <- left_join(plot.re, coords, by = "Plot")

my.coord.sf <- my.coord %>% st_as_sf(coords = c("X", "Y"), crs = 2193)

viz <- left_join(my.coord.sf, plot.re, by = "Plot")
viz$Group <- substr(viz$Plot, 1, 2)

# COM-Poisson provided somewhat better residual behaviour in simpler models, 
# but was computationally prohibitive for the full hierarchical models; 
# generalized-Poisson models were therefore used for the main analysis, 
# with COM-Poisson fits to reduced models used to assess sensitivity of the principal conclusions.

ggplot()+
  theme_void()+
  geom_sf(data = river.lines.cropped)+
  geom_sf(data = river.polygons.cropped)+
  geom_sf(data = ashley, colour = "blue", lwd =2)+
  geom_sf(data = fence, colour = "red")+
  geom_sf(data = viz , aes(colour =effect, fill = effect, shape = Group),colour = "grey", size =5)+
  scale_fill_gradient(low = "yellow", high = "purple")+
  scale_shape_manual(values = c(21, 23))+
  annotation_scale(location = "bl")+
  annotation_north_arrow(location = "bl",
                         pad_x = unit(0.1, "cm"),
                         pad_y = unit(1, "cm"))+
  ggtitle("Pember fan: random effects for biodiversity appear spatially autocorrelated")+
  theme(plot.title = element_text(size = 16, face =2, colour = "grey"))+
  labs(fill = "Random effect size", 
       shape = "Plot type",
       subtitle = "model: Indigenous ~ Year + scale(fence.dist) + (1|Plot/Subplot), family = genpois()")

ggsave("Pember fan random effect fence spatial.png",  
       scale = 1.2, width =9, height = 6,
       bg = "white")


# trajectory of 

summary(Models[[10]])
ranef(Models[[10]])

model_performance(Models[[10]])

# fence distances 
range(my.data$fence.dist)
fence.mean <- mean(my.data$fence.dist)
fence.sd <- sd(my.data$fence.dist)
fence.dist <- seq(from = 100, to= 1000, by =100)


newdata <- expand_grid(Year = 0:8,
                      fence.dist  = fence.dist)
newdata$fence.scale <-(newdata$fence.dist-fence.mean)/fence.sd 

newdata$fitted <-   exp(0.483491 + (-0.016209 * newdata$Year) +  
  (-0.114952 * newdata$fence.scale) +
  (0.005134 * newdata$fence.scale * newdata$Year))

newdata$fence.dist  <- as.factor(newdata$fence.dist )

str(newdata)

ggplot()+
  theme_bw() +
  geom_line(data = newdata, aes(x = Year, y= fitted, colour = fence.dist), lwd =1)+
  theme(panel.grid = element_blank())+
  labs(colour = "Distance from fence line (m)",
       x = "\nYear",
       y= "Expected indigenous biodiversity \n(per subplot)\n")+
  theme(axis.title = element_text(size = 16)) +
  theme(axis.text = element_text(size = 14))
  
       
ggsave("Pember fan distance from fence.png",  
       scale = 1.2, width =9, height = 6,
       bg = "white")
                       

