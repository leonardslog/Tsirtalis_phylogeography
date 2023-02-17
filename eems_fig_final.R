rm(list = ls())
currentdir <- getwd()
library(easypackages)
setwd(currentdir)
libraries("ggplot2", "rnaturalearth", "raster", "sp", "rgeos", "rgdal")

# map params
n_america.extent <- extent(-130,-59.5,24.5,54) 

# load prepped data (same projection)
load("lakes_robin.robj")
load("rivers_robin.robj")
load("US_rivers_robin.robj")
load("NorthAmer_robin.robj")
NorthAmer_noborders <- gUnaryUnion(NorthAmer) # remove internal borders

# Flint River doesn't show up on plot for some reason, hack fix:
NAm.rs <- readOGR("USA_Rivers_and_Streams-shp/")
n_america.extent <- extent(-130,-59.5,24.5,54) 
NAm.rs <- crop(NAm.rs,n_america.extent)
crs(NAm.rs) <- crs(NorthAmer)
subset_rivers <- c("Coosa River","Flint River", "Apalachicola River",
                   "Chattahoochee River")
# subset rivers
target_rivers <- c("Rio Grande","Mississippi", 
                   "Ohio", "Arkansas", "Alabama",
                   "Mobile", "Coosa", "Flint River",
                   "Apalachicola", "Chattahoochee","Saskatchewan","Columbia","Missouri")

# make graticule and bounding box (longlat)
grat <- readOGR("ne_10m_graticules_all", layer="ne_10m_graticules_15") 
grat_cropped <- crop(grat, n_america.extent)
grat_df <- fortify(grat_cropped)

bbox <- readOGR("ne_10m_graticules_all", layer="ne_10m_wgs84_bounding_box") 
bbox_cropped <- crop(bbox, n_america.extent)
bbox_df<- fortify(bbox_cropped)

crs(lakes)
# reproject graticule (Robin)
grat_robin <- spTransform(grat_cropped, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "))  # reproject graticule
grat_df_robin <- fortify(grat_robin)
bbox_robin <- spTransform(bbox_cropped, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "))  # reproject bounding box
bbox_robin_df <- fortify(bbox_robin)

# load and extract relevant EEMS data
load("run2.Robj")
eemsdata <- run2$mrates01$data

ggplot(data=bbox_df, aes(long,lat, group=group), colour = "black") +
  geom_polygon(data=grat_df, aes(long,lat, group=group), fill="white") +
  geom_polygon(data=NorthAmer_noborders, aes(x=long, lat, group = group), colour = "dark grey", fill = "white", size = 0.25) + 
  geom_raster(data = eemsdata, aes(x = x, y = y, fill = z), inherit.aes = FALSE) +
  scale_fill_gradient2(low="darkorange", mid="white", high="cyan", midpoint = 0,limits=c(-2.5,2.5),name = "log(m)") +
  geom_polygon(data=NorthAmer_noborders, aes(x=long, lat, group = group), colour = "dark grey", fill = "transparent", size = 0.25) +
  geom_polygon(data=lakes, aes(x=long, lat, group = group), colour = "dark grey", fill = "transparent", size = 0.25) +
  geom_path(data = subset(rivers, name %in% target_rivers), aes(x=long, lat, group=group), colour = "dark grey", size = 0.25) +
  geom_path(data = subset(NAm.rs, Name %in% subset_rivers), aes(x=long, lat, group=group), colour = "dark grey", size = 0.25) +
  coord_equal(xlim = c(-130,-59.5), ylim = c(24.5,54), expand=F) +
  theme(panel.grid.major = element_line(colour = "white", size = 0.30),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.85, 0.25),
        plot.margin = unit(c(0,0,0,0)+0.5, "cm")) +
  # panel.border = element_rect(colour = "white")) +
  xlab("Longitude") +
  ylab("Latitude")
