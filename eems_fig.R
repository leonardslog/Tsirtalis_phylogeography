rm(list = ls())
currentdir <- getwd()
setwd(currentdir)
library(easypackages)
libraries("ggplot2","rnaturalearth","raster","sp","rgeos","rgdal","terra","tidyterra","sf","geodata")

# set map projections and extent
crs.universal <- "+proj=robin +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
n_america.extent <- extent(-130,-59.5,24.5,54)

# load prepped data (same projection)
load("lakes_robin.robj")
load("rivers_robin.robj")
load("US_rivers_robin.robj")

# conversion from spatvector to spdf takes too long
# usa <- gadm("usa", resolution=1,level=1,".")
# usa.spdf <-  	as(usa, "Spatial") 
# can <- gadm("canada", resolution=1,level=1,".")
# can.spdf <- as.data.frame(can, "Spatial")
# mex <- gadm("mexico", resolution=1,level=1,".")
# mex.spdf <- as.data.frame(mex, "Spatial")

usa <-getData('GADM', country='USA', level=1)
can <-getData('GADM', country='CAN', level=1)
mex <-getData('GADM', country='MEX', level=1)
northamer <- rbind(can,usa,mex)
northamer.crop <- crop(northamer,n_america.extent)
crs(northamer.crop) <- crs.universal

# read in water bodies and set CRS
riverslakes.centerlines <- rnaturalearth::ne_download(scale=10, category = 'physical', type = "rivers_lake_centerlines")
lakes.centerlines <- rnaturalearth::ne_download(scale=10, category = 'physical', type = "lakes")
rivers <- crop(riverslakes.centerlines,n_america.extent)
lakes <- crop(lakes.centerlines, n_america.extent)
lakes.noGL <-  subset(lakes, name != "Saginaw Bay") # remove border btwn Saginaw Bay and Lake Huron
plot(lakes.noGL)
US.rivers <- readOGR("USA_Rivers_and_Streams-shp/") # from https://hub.arcgis.com/datasets/esri::usa-rivers-and-streams/explore?location=39.980858%2C-119.086063%2C4.31

# subset rivers and lakes
target_rivers1 <- c("Rio Grande","Mississippi River", 
                   "Ohio River", "Arkansas River", "Alabama River",
                   "Mobile River", "Coosa River","Flint River",
                   "Apalachicola River", "Chattahoochee River")
target_rivers2 <- c("Saskatchewan","Columbia","Missouri")


# Flint River doesn't show up on plot for some reason, hack fix:
NAm.rs <- readOGR("USA_Rivers_and_Streams-shp/")
NAm.rs <- crop(NAm.rs,n_america.extent)
crs(NAm.rs) <- crs(northamer.crop)
subset_rivers <- c("Coosa River","Flint River", "Apalachicola River",
                   "Chattahoochee River")
# subset rivers
target_rivers3 <- c("Rio Grande","Mississippi", 
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
crs(lakes.noGL)
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
  geom_polygon(data=northamer.crop, aes(x=long, lat, group = group), colour = "dark grey", fill = "white", size = 0.25) +
  geom_raster(data = eemsdata, aes(x = x, y = y, fill = z), inherit.aes = FALSE) +
  scale_fill_gradient2(low="darkorange", mid="white", high="cyan", midpoint = 0,limits=c(-2.5,2.5),name = "log(m)") +
  geom_polygon(data=northamer.crop, aes(x=long, lat, group = group), colour = "dark grey", fill = "transparent", size = 0.25) +
  # geom_spatvector(data=northamer_crop, aes(x=long, lat), colour = "dark grey", fill = "transparent", size = 0.25) +
  geom_path(data = subset(riverslakes.centerlines, name %in% target_rivers1), aes(x=long, lat, group=group), colour = "dark blue", size = 0.25) +
  geom_path(data = subset(riverslakes.centerlines, name %in% target_rivers3), aes(x=long, lat, group=group), colour = "dark blue", size = 0.25) +
  geom_path(data = subset(riverslakes.centerlines, name %in% target_rivers2), aes(x=long, lat, group=group), colour = "dark blue", size = 0.25) +
  geom_path(data = subset(NAm.rs, Name %in% subset_rivers), aes(x=long, lat, group=group), colour = "dark blue", size = 0.25) +
  geom_polygon(data=lakes.noGL, aes(x=long, lat, group = group), colour = "dark grey", fill = "white", size = 0.25) +
  coord_equal(xlim = c(-130,-59.5), ylim = c(24.5,54), expand=F) +
  theme(panel.grid.major = element_line(colour = "white", size = 0.30),
        panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = c(0.85, 0.25),
        plot.margin = unit(c(0,0,0,0)+0.5, "cm")) +
  # panel.border = element_rect(colour = "white")) +
  xlab("Longitude") +
  ylab("Latitude")
