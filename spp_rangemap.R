library(easypackages)
libraries("smoothr", "dplyr", "gdistance", "rgdal", "rgeos", "maps", "maptools","sp", "sf",
          "adehabitatHR","scales", "RColorBrewer", "terra", "geodata","base",
          "spdep", "rangemap","alphahull","raster")

###################
#  country shape  #
###################

# base map settings
crs.universal <- "+proj=robin +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
mycrs <- "+proj=robin +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
n_america.extent <- extent(-130,-59.5,24.5,62)

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
crs(srtm.NA) <- crs.universal

# pull in admin border data
usa <- gadm("usa", resolution=1,level=1,".")
can <- gadm("canada", resolution=1,level=1,".")
mex <- gadm("mexico", resolution=1,level=1,".")
northamer <- rbind(can,usa,mex)
northamer_crop <- crop(northamer,n_america.extent)
crs(northamer_crop) <- crs.universal

# subset rivers and lakes
target_rivers <- c("Rio Grande","Mississippi River", 
                   "Ohio River", "Arkansas River", "Alabama River",
                   "Mobile River", "Coosa River","Flint River",
                   "Apalachicola River", "Chattahoochee River")
target_rivers2 <- c("Saskatchewan","Columbia","Missouri")

target_lakes <- c("Lake Michigan", "Lake Huron", 
                  "Lake Erie", "Lake Ontario", "Lake Superior",
                  "Great Salt Lake","Lake Nipigon",
                  "Lake Manitoba","Lake Winnipeg","Lake Winnipegosis")

############################
# sirtalis range shapefile #
############################

sirtalis_range <- readOGR("iucn_rangedata-sirtalis/")
crs(sirtalis_range) <- crs.universal
plot(sirtalis_range,col=alpha("red",0.4),add=T)

##########################
# spp range construction #
##########################

############
# sirtalis #
############
sirt <- read.csv("sirtalis_clean.csv", header=TRUE)
sirt.pts <- na.omit(sirt[,(c("decimalLongitude","decimalLatitude"))])
sirt.pts.nodupes <- sirt.pts[!duplicated(paste(sirt.pts$decimalLongitude,sirt.pts$decimalLatitude)),]

sirt.spts <- SpatialPoints(sirt.pts.nodupes, proj4string=CRS(crs.universal), bbox = NULL)
sirt.spol <- hull_polygon(sirt.spts, hull_type = "concave", concave_distance_lim=5,
                          verbose = TRUE)
sirt.buff <- gBuffer(sirt.spol)
ONT.poly <- as(northamer_crop[northamer_crop$NAME_1=="Ontario"],"Spatial") # convert
crs(ONT.poly) <- crs.universal
sirtalis_ONT <- gIntersection(sirtalis_range,ONT.poly)
sirt.range <- rbind(sirt.buff,sirtalis_ONT,makeUniqueIDs = TRUE)
sirt.range <- gUnaryUnion(sirt.range)
sirt.range <- gIntersection(sirtalis_range,sirt.range)
plot(sirt.range,col=alpha("#b2df8a", 0.5), add=T)

###########
# similis #
###########
simi <- read.csv("similis_clean.csv", header=TRUE)
simi.pts <- na.omit(simi[,(c("decimalLongitude","decimalLatitude"))])
simi.pts.nodupes <- simi.pts[!duplicated(paste(simi.pts$decimalLongitude,simi.pts$decimalLatitude)),]
simi.spts <- SpatialPoints(simi.pts.nodupes, proj4string=CRS(crs.universal), bbox = NULL)
simi.poly <- hull_polygon(simi.spts, hull_type = "concave", concave_distance_lim=2,
                          verbose = TRUE)
simi.buff <- gBuffer(simi.poly,width = 0.25)
simi.range <- smooth(simi.buff, method = "ksmooth", smoothness=3, n=10)
simi.range <- gIntersection(sirtalis_range, simi.range)
plot(sirtalis_range)
plot(simi.range,col=alpha("#FFA7FE", 0.5),add=T)

##############
# pallidulus #
##############

palli <- read.csv("pallidulus_clean.csv", header=TRUE)
palli.pts <- na.omit(palli[,(c("decimalLongitude","decimalLatitude"))])
palli.pts.nodupes <- palli.pts[!duplicated(paste(palli.pts$decimalLongitude,palli.pts$decimalLatitude)),]
palli.spts <- SpatialPoints(palli.pts.nodupes, proj4string=CRS(crs.universal), bbox = NULL)
palli.poly <- hull_polygon(palli.spts, hull_type = "concave", concave_distance_lim=1.75,
                           verbose = TRUE)
palli.poly.smooth <- smooth(palli.poly, method = "ksmooth", smoothness=5, n=10)
# plot(palli.poly.smooth)
palli.buff <- gBuffer(palli.poly.smooth)
palli.buff <- gIntersection(palli.buff, as(northamer_crop,"Spatial"))
# QBC.poly <- as(NAm.prj[NAm.prj$name=="Québec"],"Spatial") # convert
QBC.poly <- as(northamer_crop[northamer_crop$NAME_1=="Québec"],"Spatial") # convert

crs(QBC.poly) <- crs.universal
sirtalis_QBC <- gIntersection(sirtalis_range,QBC.poly)
palli.range <- rbind(palli.buff,sirtalis_QBC,makeUniqueIDs = TRUE)
palli.range <- gBuffer(palli.range,width = 0)
plot(sirtalis_range)
plot(palli.range,col=alpha("#866117", 0.5),add=T)

#############
# annectens #
#############

annect <- read.csv("annectens_clean.csv", header=TRUE)
annect.pts <- na.omit(annect[,(c("decimalLongitude","decimalLatitude"))])
annect.pts.nodupes <- annect.pts[!duplicated(paste(annect.pts$decimalLongitude,annect.pts$decimalLatitude)),]
annect.spts <- SpatialPoints(annect.pts.nodupes, proj4string=CRS(mycrs), bbox = NULL)
annect.poly <- hull_polygon(annect.spts, hull_type = "concave", concave_distance_lim=1.75,
                            verbose = TRUE)
# plot(annect.poly,add=T,col="#33a02c")
annect.buff <- gBuffer(annect.poly,width = 0.5)
annect.buff <- gIntersection(annect.buff,sirtalis_range)
# plot(annect.buff,add=T,col=alpha("#33a02c", 0.5))
annect.range <- annect.buff
# annect.poly.smooth <- smooth(annect.poly, method = "ksmooth", smoothness=5, n=10)
# TX.poly <- as(northamer_crop[northamer_crop$NAME_1=="Texas"],"Spatial") # convert
# crs(TX.poly) <- crs.universal
# sirtalis_TX <- gIntersection(sirtalis_range,TX.poly)
# annect_TX <- gDifference(sirtalis_TX,sirt.range)
# annect.range <- rbind(annect_TX, annect.buff,makeUniqueIDs = TRUE)
# annect.range <- gUnaryUnion(annect.range)
# annect.range <- gDifference(annect.range,dors.range)

#################
# semifasciatus #
#################
semi <- read.csv("semifasciatus_clean.csv", header=TRUE)
semi.pts <- na.omit(semi[,(c("decimalLongitude","decimalLatitude"))])
semi.pts.nodupes <- semi.pts[!duplicated(paste(semi.pts$decimalLongitude,semi.pts$decimalLatitude)),]
semi.spts <- SpatialPoints(semi.pts.nodupes[-c(353,354),], proj4string=CRS(mycrs), bbox = NULL)
semi.poly <- hull_polygon(semi.spts, hull_type = "concave", concave_distance_lim=1.75,
                          verbose = TRUE)
# plot(semi.poly, col=alpha("#E31A1C",0.5),add=T)
semi.buff <- gBuffer(semi.poly,width = 0.5)
semi.range <- gIntersection(semi.buff,sirtalis_range)
# plot(semi.range,add=T,col=alpha("#E31A1C", 0.5))

##############
# parietalis #
##############

parie <- read.csv("parietalis_clean.csv", header=TRUE)
parie.pts <- na.omit(parie[,(c("decimalLongitude","decimalLatitude"))])
parie.pts.nodupes <- parie.pts[!duplicated(paste(parie.pts$decimalLongitude,parie.pts$decimalLatitude)),]
parie.spts <- SpatialPoints(parie.pts.nodupes, proj4string=CRS(mycrs), bbox = NULL)
parie.poly <- hull_polygon(parie.spts, hull_type = "concave", concave_distance_lim=1.75,
                           verbose = TRUE)
# plot(sirtalis_range)
# plot(parie.poly, col=alpha("#ffff99",0.5),add=T)
parie.buff <- gBuffer(parie.poly,width = 1.5)
# plot(parie.buff, col=alpha("#ffff99",0.5))
parie.range <- gIntersection(parie.buff,sirtalis_range)
# parie.range <- smooth(parie.range, method = "ksmooth", smoothness=10, n=10)
# plot(sirtalis_range)
# plot(northamer_crop,add=T)
# plot(parie.range,add=T,col=alpha("#ffff99", 0.5))

############
# dorsalis #
############

dors <- read.csv("dorsalis_clean.csv", header=TRUE)
dors.pts <- na.omit(dors[,(c("decimalLongitude","decimalLatitude"))])
dors.pts.nodupes <- dors.pts[!duplicated(paste(dors.pts$decimalLongitude,dors.pts$decimalLatitude)),]
dors.spts <- SpatialPoints(dors.pts.nodupes, proj4string=CRS(mycrs), bbox = NULL)
dors.spol <- hull_polygon(dors.spts, hull_type = "concave", concave_distance_lim=25,
                          verbose = TRUE)
NM.poly <- as(northamer_crop[northamer_crop$NAME_1=="New Mexico"],"Spatial") # convert
Chi.poly <- as(northamer_crop[northamer_crop$NAME_1=="Chihuahua"],"Spatial") # convert
dorsrange.poly <- rbind(NM.poly,Chi.poly)
dorsrange.buff <- gBuffer(dorsrange.poly)
dorsrange.poly <- smooth(dorsrange.poly, method = "ksmooth", smoothness=10, n=10)
# plot(dorsrange.poly)
# dorsrange.poly <- gIntersection(dorsrange.poly)
dors.range <- gIntersection(sirtalis_range,dorsrange.buff)
# plot(sirtalis_range)
# plot(dors.range, col=alpha("#8F6996",0.5),add=T)

##########
# fitchi #
##########

fitchi <- read.csv("fitchi_clean.csv", header=TRUE)
fitchi.pts <- na.omit(fitchi[,(c("decimalLongitude","decimalLatitude"))])
fitchi.pts.nodupes <- fitchi.pts[!duplicated(paste(fitchi.pts$decimalLongitude,fitchi.pts$decimalLatitude)),]
fitchi.spts <- SpatialPoints(fitchi.pts.nodupes, proj4string=CRS(mycrs), bbox = NULL)
fitchi.poly <- hull_polygon(fitchi.spts, hull_type = "concave", concave_distance_lim=1.75,
                            verbose = TRUE)
# plot(sirtalis_range)
# plot(fitchi.poly, col=alpha("#ff7f00",0.5),add=T)
fitchi.buff <- gBuffer(fitchi.poly,width = .25)
fitchi.range <- gIntersection(fitchi.buff,sirtalis_range)
fitchi.range <- smooth(fitchi.range, method = "ksmooth", smoothness=10, n=10)
# plot(sirtalis_range)
# plot(fitchi.range,col=alpha("#ff7f00",0.5),add=T)

##############
# infernalis #
##############

infer <- read.csv("infernalis_clean.csv", header=TRUE)
infer.pts <- na.omit(infer[,(c("decimalLongitude","decimalLatitude"))])
infer.pts.nodupes <- infer.pts[!duplicated(paste(infer.pts$decimalLongitude,infer.pts$decimalLatitude)),]
infer.pts.nodupes <- data.frame(infer.pts[!duplicated(paste(infer.pts$decimalLongitude,infer.pts$decimalLatitude)),])
infer.spts <- SpatialPoints(infer.pts.nodupes[-c(601,728,730,739,740),], proj4string=CRS(mycrs), bbox = NULL)
infer.poly <- hull_polygon(infer.spts, hull_type = "concave", concave_distance_lim=1.75,
                           verbose = TRUE)
# plot(sirtalis_range)
# plot(infer.poly, col=alpha("#D700FF",0.5),add=T)
infer.buff <- gBuffer(infer.poly,width = 0.25)
infer.range <- gIntersection(infer.buff,sirtalis_range)
infer.range <- smooth(infer.range, method = "ksmooth", smoothness=10, n=10)
# plot(sirtalis_range)
# plot(infer.range,col=alpha("#D700FF",0.5),add=T)


#############
# concinnus #
#############

conci <- read.csv("concinnus_clean.csv", header=TRUE)
conci.pts <- na.omit(conci[,(c("decimalLongitude","decimalLatitude"))])
conci.pts.nodupes <- conci.pts[!duplicated(paste(conci.pts$decimalLongitude,conci.pts$decimalLatitude)),]
conci.spts <- SpatialPoints(conci.pts.nodupes[-c(878:879,884:886,889,890:894,901:908),], proj4string=CRS(crs.universal))
conci.poly <- hull_polygon(conci.spts, hull_type = "concave", concave_distance_lim=1,
                           verbose = TRUE)
# plot(sirtalis_range)
# plot(conci.poly, col=alpha("#fb9a99",0.5),add=T)
conci.buff <- gBuffer(conci.poly,width = .1)
conci.range <- gIntersection(conci.buff,sirtalis_range)
conci.range <- smooth(conci.range, method = "ksmooth", smoothness=10, n=10)
# plot(sirtalis_range)
# plot(conci.range,add=T,col=alpha("#fb9a99",0.5))

###############
# tetrataenia #
###############

tetra <- read.csv("tetrataenia_clean.csv", header=TRUE)
tetra.pts <- na.omit(tetra[,(c("decimalLongitude","decimalLatitude"))])
tetra.pts.nodupes <- tetra.pts[!duplicated(paste(tetra.pts$decimalLongitude,tetra.pts$decimalLatitude)),]
tetra.spts <- SpatialPoints(tetra.pts.nodupes, proj4string=CRS(crs.universal), bbox = NULL)
tetra.poly <- hull_polygon(tetra.spts, hull_type = "concave", concave_distance_lim=1.75,
                           verbose = TRUE)
# plot(sirtalis_range)
# plot(tetra.poly, col=alpha("#04E4FB",0.5),add=T)
tetra.buff <- gBuffer(tetra.poly,width = 0.5)
tetra.range <- gIntersection(tetra.buff,sirtalis_range)
tetra.range <- smooth(tetra.range, method = "ksmooth", smoothness=10, n=10)
# plot(sirtalis_range)
# plot(tetra.range,add=T,col=alpha("#04E4FB",0.5))

###############
# pickeringii #
###############

pick <- read.csv("pickeringii_clean.csv", header=TRUE)
pick.pts  <- na.omit(pick[,(c("decimalLongitude","decimalLatitude"))])
pick.pts.nodupes <- pick.pts[!duplicated(paste(pick.pts$decimalLongitude,pick.pts$decimalLatitude)),]
pick.spts <- SpatialPoints(pick.pts.nodupes, proj4string=CRS(crs.universal), bbox = NULL)
pick.poly <- hull_polygon(pick.spts, hull_type = "concave", concave_distance_lim=1.75,
                          verbose = TRUE)
# plot(sirtalis_range)
# plot(sirtalis_range,add=T)
# plot(pick.poly, col=alpha("#3E00FF",0.5),add=T)
pick.buff <- gBuffer(pick.poly,width = 0.25)
pick.smooth <- smooth(pick.buff, method = "ksmooth", smoothness=10, n=10)
pick.range <- gIntersection(pick.smooth,as(northamer_crop,"Spatial"))
# plot(sirtalis_range)
# plot(pick.range,col=alpha("#3E00FF",0.5),add=T)


sirtalis_range_new <- rbind(sirtalis_range,sirt.range,simi.range,palli.range,semi.range,parie.range,
                            annect.range,dors.range,infer.range,fitchi.range,conci.range,tetra.range,pick.range)

sirt.range.df <- data.frame(sirt.range)
sirt.range.spdf <- SpatialPolygonsDataFrame(sirt.range,data = as.data.frame(sirt.range))

#########################
#  now plot everything  #
#########################

over.col <- colorRampPalette(c("white", "black"))
terra::plot(northamer_crop,
            ylim=c(24.5, 62),xlim=c(-130,-59.5),
            pax=list(xaxt="n", yaxt="n"),
            mar=c(11,8,9.5,7.5)-3,legend=F,
            lwd=1, border="dark grey")
sbar(500, xy="bottomright", type="bar", lonlat=TRUE,below="km", divs=4, cex=.8)
# plot(lakes.centerlines,xlim=x, ylim=y,col="white",border=alpha("dark grey",0.9),lwd=1,add=T)
# plot(subset(US.rivers, Name %in% target_rivers), col=alpha("blue",0.9),lwd=1, add=T)
# plot(subset(riverslakes.centerlines, name %in% target_rivers2), col=alpha("blue",0.9),lwd=1, add=T)

plot(sirt.range,col=alpha("#b2df8a", 0.7), border=NA,add=T)
plot(simi.range,col=alpha("#FFA7FE", 1),border=NA,add=T)
plot(palli.range,col=alpha("#866117", 0.7),border=NA,add=T)
plot(semi.range,col=alpha("#E31A1C", 0.8),border=NA,add=T)
plot(parie.range,col=alpha("#ffff99", 0.8),border=NA,add=T)
plot(annect.range,col=alpha("#33a02c", 0.7),border=NA,add=T)
plot(dors.range,col=alpha("#8F6996",0.7),border=NA,add=T)
plot(infer.range,col=alpha("#D700FF",1),border=NA,add=T)
plot(fitchi.range,col=alpha("#ff7f00",0.7),border=NA,add=T)
plot(conci.range,col=alpha("#3fff00",1),border=NA,add=T) # #3fff00 #fb9a99
plot(tetra.range,col=alpha("#04E4FB",1),border=NA,add=T)
plot(pick.range,col=alpha("#3E00FF",1),border=NA,add=T)

plot(subset(lakes.centerlines, name %in% target_lakes),xlim=x, ylim=y,col="white",border=alpha("dark grey",0.9),lwd=1,add=T)

# add study samples
samples <- read.csv("sirtalis_samples.csv",stringsAsFactors = F)
samples.sp <- SpatialPoints(coords = samples[c("longitude","latitude")],proj4string=CRS(crs.universal))
points(samples.sp,pch=20)
