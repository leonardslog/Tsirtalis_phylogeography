library(easypackages)
libraries("maps","rgdal","raster","RColorBrewer","scales",
          "maptools","gridExtra","rgeos","geodata","terra",
          "rnaturalearthdata","rnaturalearth")

# set map projections and extent
crs.universal <- "+proj=robin +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
n_america.extent <- extent(-130,-59.5,24.5,54)

# read in water bodies and set CRS
riverslakes.centerlines <- rnaturalearth::ne_download(scale=10, category = 'physical', type = "rivers_lake_centerlines")
lakes.centerlines <- rnaturalearth::ne_download(scale=10, category = 'physical', type = "lakes")
rivers <- crop(riverslakes.centerlines,n_america.extent)
lakes <- crop(lakes.centerlines, n_america.extent)
US.rivers <- readOGR("USA_Rivers_and_Streams-shp/") # from https://hub.arcgis.com/datasets/esri::usa-rivers-and-streams/explore?location=39.980858%2C-119.086063%2C4.31

# pull elevation data
srtm.USA <- elevation_30s(country="USA",path = ".")
srtm.CAN <- elevation_30s(country="CAN",path = ".")
srtm.MEX <- elevation_30s(country="MEX",path = ".")
srtm.NA <- mosaic(srtm.CAN,srtm.USA,srtm.MEX)

usa <- gadm("usa", resolution=1,level=1,".")
can <- gadm("canada", resolution=1,level=1,".")
mex <- gadm("mexico", resolution=1,level=1,".")
northamer <- rbind(can,usa,mex)
northamer_crop <- crop(northamer,n_america.extent)
crs(northamer_crop) <- crs.universal

# plot map
terra::plot(srtm.NA,
            col=alpha(over.col(100),0.4),
            ylim=c(24.5, 54),xlim=c(-130,-59.5),
            mar=c(11,8,10,8)-3,legend=F,xlab="longtitude",ylab="Latitude")
plot(northamer_crop,lwd=1, col=alpha("dark grey",0.3), border="dark grey",add=T)
plot(l.cl,xlim=x, ylim=y,col="white",border=alpha("dark grey",0.9),lwd=1,add=T)
plot(SBxSL_noborders,xlim=x, ylim=y,col="white",border=alpha("dark grey",0.9),lwd=1,add=T)
plot(subset(US.rivers, Name %in% target_rivers), col=alpha("blue",0.9),lwd=1, add=T)
plot(subset(rl.cl, name %in% target_rivers2), col=alpha("blue",0.9),lwd=1, add=T)
