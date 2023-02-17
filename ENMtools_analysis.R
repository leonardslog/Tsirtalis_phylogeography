rm(list=ls())
currentdir <- getwd()
setwd(currentdir)
options(java.parameters = "-Xmx10g")
library(easypackages)
libraries("dismo", "ENMTools","dplyr", "jsonlite","ade4","gdistance", "adehabitatHR", "adehabitatLT", "adehabitatMA", 
          "rgdal", "rgeos", "lattice", "maps", "maptools", "raster", "sp", "car", "rJava", "rpaleoclim",
          "sm", "rgl", "rpanel", "rms", "virtualspecies", "base", "ecospat", "spdep", "viridis")

#### helper functions ####

# windows_map <- function(){
#   windows(width=8, height=6.5)
#   opt <- par(mar=rep(0, 4), oma=c(0,0,0,0), omi=c(0.5,0.5,0.5,0.5), mai=rep(0, 4))
# }

#### standardize projection and extent ####

bioclim.proj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
EPSG.26907 <- "+proj=utm +zone=7 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" 
NorthAmerica <- extent(-130,-59.5,24.5,54)

#### load and extract occurrence data for T. sirtalis samples ####

samples <- read.csv("sample_summary.csv", header = TRUE)
head(samples)
points.samples <- unique.data.frame(samples[,cbind(7,6)])
head(points.samples)

#### load and extract unique GBIF research grade occurrence data for T. sirtalis ####

gbif.data <- read.csv("gbif.csv")

head(gbif.data)
total.gbif <- unique.data.frame(gbif.data[,cbind(3,2)])
head(total.gbif)
gbif.spdf <- SpatialPoints(coords=total.gbif, proj4string=CRS(bioclim.proj))
# plot(gbif.spdf,pch=19)

#### combine all points

total.data <- merge(data.frame(points.samples,row.names=NULL), data.frame(total.gbif, row.names=NULL),all.x=TRUE,all.y=TRUE)
total.points <- total.data[,c(2,1)] 
head(total.points)

#### subset population points for population polygons and make total

# total polygon
pbuffer.total <- rbind(samples[,c(6,7)]+0.0001, samples[,c(6,7)]-0.0001)
points.total <- rbind(pbuffer.total[,2:1]+0.00001,pbuffer.total[,2:1]-0.00001)
# plot(points.total,pch=19)    
points.total.sp <-SpatialPoints(points.total,proj4string=CRS(bioclim.proj))
poly.total <-mcp(points.total.sp,percent=100)
range.total <- raster(poly.total)
poly.total@data$id <- 1
poly.total@proj4string <- CRS(bioclim.proj)

# Southeast polygon
pbuffer.SE <- rbind((subset(samples[,c(6,7)]+0.0001, samples[,2]=="Southeast")),
                 (subset(samples[,c(6,7)]-0.0001, samples[,2]=="Southeast")))
points.SE <- rbind(pbuffer.SE[,2:1]+0.00001,pbuffer.SE[,2:1]-0.00001)
# plot(points.SE, col="#d7191c",pch=19)    
points.SE.sp <-SpatialPoints(points.SE,proj4string=CRS(bioclim.proj))
poly.SE <-mcp(points.SE.sp,percent=100)
range.SE <- raster(poly.SE)
poly.SE@data$id <- 1
poly.SE@proj4string <- CRS(bioclim.proj)
# plot(poly.SE)

# East polygon
pbuffer.E <- rbind((subset(samples[,c(6,7)]+0.0001, samples[,2]=="East")),
                (subset(samples[,c(6,7)]-0.0001, samples[,2]=="East")))
points.E <- rbind(pbuffer.E[,2:1]+0.00001,pbuffer.E[,2:1]-0.00001)
# plot(points.E, col="#2c7bb6",pch=19)    
points.E.sp <-SpatialPoints(points.E,proj4string=CRS(bioclim.proj))
poly.E <-mcp(points.E.sp,percent=100)
range.E <- raster(poly.E)
poly.E@data$id <- 2
poly.E@proj4string <- CRS(bioclim.proj)
# plot(poly.E)

# Central polygon
pbuffer.C <- rbind((subset(samples[,c(6,7)]+0.0001, samples[,2]=="Central")),
                 (subset(samples[,c(6,7)]-0.0001, samples[,2]=="Central")))
points.C <- rbind(pbuffer.C[,2:1]+0.00001,pbuffer.C[,2:1]-0.00001)
# plot(points.C, col="#fdae61",pch=19)    
points.C.sp <-SpatialPoints(points.C,proj4string=CRS(bioclim.proj))
poly.C <-mcp(points.C.sp,percent=100)
range.C <- raster(poly.C)
poly.C@data$id <- 3
poly.C@proj4string <- CRS(bioclim.proj)
# plot(poly.C)

# West polygon
pbuffer.W <- rbind((subset(samples[,c(6,7)]+0.0001, samples[,2]=="West")),
                (subset(samples[,c(6,7)]-0.0001, samples[,2]=="West")))
points.W <- rbind(pbuffer.W[,2:1]+0.00001,pbuffer.W[,2:1]-0.00001)
# plot(points.W, col="#abd9e9",pch=19)    
points.W.sp <-SpatialPoints(points.W, proj4string=CRS(bioclim.proj))
poly.W <-mcp(points.W.sp,percent=100)
range.W <- raster(poly.W)
poly.W@data$id <- 4
poly.W@proj4string <- CRS(bioclim.proj)
# plot(poly.W)

# samples polygon
pbuffer.samples <- rbind(samples[,c(6,7)]+0.0001, samples[,c(6,7)]-0.0001)
points.samples <- rbind(pbuffer.samples[,2:1]+0.00001, pbuffer.samples[,2:1]-0.00001)
# plot(points.samples, pch=19)    
points.samples.sp <-SpatialPoints(points.samples, proj4string=CRS(bioclim.proj))
poly.samples <-mcp(points.samples.sp,percent=100)
range.samples <- raster(poly.samples)
poly.samples@data$id <- 5
poly.samples@proj4string <- CRS(bioclim.proj)
# plot(poly.samples)

# plot lineage polygons and g
load("NorthAmer.cropped.robj")
NorthAmer.r <- raster(NorthAmer)
# plot(NorthAmer, col="grey")
# points(points.SE, col="#d7191c",pch=19,cex=1)
# points(points.W, col="#abd9e9",pch=19,cex=1)    
# points(points.C, col="#fdae61",pch=19,cex=1)    
# points(points.E, col="#2c7bb6",pch=19,cex=1)   
# plot(poly.SE, add=T)
# plot(poly.E, add=T)
# plot(poly.C, add=T)
# plot(poly.W, add=T)

# assign GBIF+sample points to containing polygon
total.spdf <- SpatialPointsDataFrame(total.data[,c("longitude","latitude")],
                                     total.data[,c("longitude","latitude")],
                                     proj4string = CRS(bioclim.proj))
                                     
poly.df <- rbind(poly.SE, poly.E, poly.C, poly.W, makeUniqueIDs = TRUE)


# LOOP USING gDistance, DISTANCES STORED IN LIST OBJECT

#distance of each point to each polygon
Fdist <- list()
for(i in 1:dim(total.spdf)[1]) {
  pDist <- vector()
  for(j in 1:dim(poly.df)[1]) { 
    pDist <- append(pDist, gDistance(total.spdf[i,],poly.df[j,])) 
  }
  Fdist[[i]] <- pDist
} 
# RETURN POLYGON (NUMBER) WITH THE SMALLEST DISTANCE FOR EACH POINT  
min.dist <- unlist(lapply(Fdist, FUN=function(x) which(x == min(x))[1])) 

# RETURN DISTANCE TO NEAREST POLYGON
polyDist <- unlist(lapply(Fdist, FUN=function(x) min(x)[1])) 

# CREATE POLYGON-ID AND MINIMUM DISTANCE COLUMNS IN POINT FEATURE CLASS
total.spdf@data <- data.frame(total.spdf@data, polyID=min.dist, PDist=polyDist)

# plot results
# plot(NorthAmer, axes=T, col='grey')
# plot(total.spdf,col=total.spdf@data$polyID, pch=20) #polygon assignment

# subset by polygon
total.pts.SE <- total.spdf[total.spdf@data$polyID.1==1,]
total.pts.E <- total.spdf[total.spdf@data$polyID.1==2,]
total.pts.C <- total.spdf[total.spdf@data$polyID.1==3,]
total.pts.W <- total.spdf[total.spdf@data$polyID.1==4,]

# subset and rarify coords by polygon
coords.SE <- unique(as.data.frame(coordinates(total.pts.SE)))
coords.E <- unique(as.data.frame(coordinates(total.pts.E)))
coords.C <- unique(as.data.frame(coordinates(total.pts.C)))
coords.W <- unique(as.data.frame(coordinates(total.pts.W)))

# download bioclim data, remove highly correlated variables
bioclim <- getData('worldclim', var='bio', res=2.5)
# bioclim

## crop bioclim data
bioclim.NA <- crop(bioclim,NorthAmerica) # subset to North America
# writeRaster(x=bioclim.NA$bio1,filename="bio1.asc",format="ascii") # raster template for waterdistance in QGis
# examine correlations btwn climate variables
# raster.cor.plot(bioclim.NA)

#### load water body data, format to bioclim raster dims,
#### and export raster to make distance raster 

# H20 <- readOGR(dsn="data/NA_Lakes_and_Rivers/hydrography_l_rivers_v2.shp")
# H20.spdf <-spTransform(H20,CRS(bioclim.proj))
bio1.r <- bioclim.NA$bio1
# plot(bio1.r)
# lines(H20.spdf)
# bioclim.H20.NA <- bio1.r  # this will be the template
# bioclim.H20.NA[] <- NA  # assigns all values as NA
# summary(bioclim.H20.NA) # shows you what you have: all NA's
# H20.raster <- rasterize(H20.spdf, bioclim.H20.NA, field=1) # this takes a while, so write raster after
# writeRaster(H20.raster, filename="H20.raster.tiff")
# summary(H20.raster)          # pixels crossed by water have "1"
# plot(H20.raster, add=TRUE)
# lines(H20.spdf)

# load water distance raster constructed in QGIS
waterdistance <- raster("waterdistance.tif") # distance raster created in QGIS

# read elevation data to manipulate in QGIS
# elevation.USA <- getData('alt', country='USA')
# elevation.CAN <- getData('alt', country='CAN')
# elevation.MEX <- getData('alt', country='MEX')

# load QGIS-constructed landcover tiff
elevation <- raster("NA_merged_elevation.tif")
elevation.spdf <- rasterToPoints(elevation,spatial=TRUE)
elevation.tmp.r <- bio1.r  # this will be the template
elevation.tmp.r[] <- NA  # assigns all values as NA
elevation.r <- rasterize(elevation.spdf, elevation.tmp.r, field="NA_merged_elevation")

bioclim.NA.reduced <- removeCollinearity(bioclim.NA, multicollinearity.cutoff = .75, 
                                     select.variables = FALSE, sample.points = FALSE,
                                     plot=FALSE)
env.vars <- addLayer(bioclim.NA, waterdistance, elevation.r)
env.vars.reduced <- removeCollinearity(env.vars, multicollinearity.cutoff = .75, 
                                         select.variables = FALSE, sample.points = FALSE,
                                         plot=FALSE)
env <- env.vars[[c("bio1", "bio2", "bio4", "bio5", "bio8", "bio12", "bio14", "bio15","bio18","waterdistance","layer")]]
env <- setMinMax(env)

env.sb <- env.vars[[c("bio1", "bio4", "bio12", "bio14", "bio15", "bio18")]]

#### SDMs for each lineage
sirtalis.samples <- enmtools.species(species.name="samples",
                                  range = range.samples,
                                  presence.points = points.samples)
sirtalis.samples.E <- enmtools.species(species.name = "east",
                                  range = range.samples,
                                  presence.points = points.E)
sirtalis.samples.SE <- enmtools.species(species.name = "southeast",
                                  range = range.samples,
                                  presence.points <- points.SE)
sirtalis.samples.C <- enmtools.species(species.name = "central",
                                  range = range.samples,
                                  presence.points <- points.C)
sirtalis.samples.W <- enmtools.species(species.name = "west",
                                  range = range.samples,
                                  presence.points = points.W)

#### create background points for each group

# making random background points with dismo package
sirtalis.samples$background.points <-  randomPoints(mask = env, n = 200, excludep = TRUE, p = sirtalis.samples$presence.points)
sirtalis.samples.SE$background.points <- randomPoints(mask = env, n = 200, excludep = TRUE, p = sirtalis.samples.SE$presence.points)
sirtalis.samples.E$background.points <- randomPoints(mask = env, n = 200, excludep = TRUE, p = sirtalis.samples.E$presence.points)
sirtalis.samples.C$background.points <- randomPoints(mask = env, n = 200, excludep = TRUE, p = sirtalis.samples.C$presence.points)
sirtalis.samples.W$background.points <- randomPoints(mask = env, n = 200, excludep = TRUE, p = sirtalis.samples.W$presence.points)

#### background raster buffer
# sirtalis.snps.SE$range <- background.raster.buffer(sirtalis.snps.SE$presence.points, 1000, mask=env)
# sirtalis.snps.E$range <- background.raster.buffer(sirtalis.snps.E$presence.points, 1000, mask=env)
# sirtalis.snps.MW$range <- background.raster.buffer(sirtalis.snps.MW$presence.points, 1000, mask=env)
# sirtalis.snps.W$range <- background.raster.buffer(sirtalis.snps.W$presence.points, 1000, mask=env)

#### maxent models

sirtalis.samples.max <- enmtools.maxent(species = sirtalis.samples,
                                           bg.source = "points",
                                           env = env, test.prop = .25)

sirtalis.samples.SE.max <- enmtools.maxent(species = sirtalis.samples.SE,
                                        bg.source = "points",
                                        env = env, test.prop = .25)

sirtalis.samples.E.max <- enmtools.maxent(species = sirtalis.samples.E, 
                                       bg.source = "points",
                                       env = env, test.prop = .25)

sirtalis.samples.C.max <- enmtools.maxent(species = sirtalis.samples.C,
                                        bg.source = "points",
                                        env = env, test.prop = .25)

sirtalis.samples.W.max <- enmtools.maxent(species = sirtalis.samples.W,
                                       bg.source = "points",
                                       env = env, test.prop = .25)

### plot without presence/absence points in the way

plot.enmtools.maxent <- function(x, ...){
  suit.points <- data.frame(rasterToPoints(x$suitability))
  colnames(suit.points) <- c("Longitude", "Latitude", "Suitability")
  suit.plot <- ggplot(data = suit.points, aes_string(y = "Latitude", x = "Longitude")) +
    geom_raster(aes_string(fill = "Suitability")) +
    scale_fill_viridis(option = "B", guide = guide_colourbar(title = "Suitability")) +
    coord_fixed() + theme_classic()
}

total.plot <- plot.enmtools.maxent(sirtalis.samples.max)
SE.plot <- plot.enmtools.maxent(sirtalis.samples.SE.max)
E.plot <- plot.enmtools.maxent(sirtalis.samples.E.max)
C.plot <- plot.enmtools.maxent(sirtalis.samples.C.max)
W.plot <- plot.enmtools.maxent(sirtalis.samples.W.max)

#### niche equivalency test
# id.E.SE <- identity.test(sirtalis.samples.E,sirtalis.samples.SE, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="identity_test/E_SE/")
# id.E.C <- identity.test(sirtalis.samples.E,sirtalis.samples.C, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="identity_test/E_C/")
# id.C.W <- identity.test(sirtalis.samples.C,sirtalis.samples.W, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="identity_test/C_W/")
# id.E.W <- identity.test(sirtalis.samples.E,sirtalis.samples.W, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="identity_test/E_W/")
# id.SE.W <- identity.test(sirtalis.samples.SE,sirtalis.samples.W, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="identity_test/SE_W/")
# id.SE.C <- identity.test(sirtalis.samples.SE,sirtalis.samples.C, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="identity_test/SE_C/")

#### background/similarity test
# bg.E.SE <- background.test(sirtalis.samples.E,sirtalis.samples.SE, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="background_test/E_SE/")
# write.csv(bg.E.SE$reps.overlap,"bg_east_southeast.csv")
# 
# bg.E.C <- background.test(sirtalis.samples.E,sirtalis.samples.C, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="background_test/E_C/")
# write.csv(bg.E.C$reps.overlap,"bg_east_central.csv")
# 
# bg.C.W <- background.test(sirtalis.samples.C,sirtalis.samples.W, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="background_test/C_W/")
# write.csv(bg.C.W$reps.overlap,"bg_central_west.csv")
# 
# bg.E.W <- background.test(sirtalis.samples.E,sirtalis.samples.W, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="background_test/E_W/")
# write.csv(bg.E.W$reps.overlap,"bg_east_west.csv")
# 
# bg.SE.W <- background.test(sirtalis.samples.SE,sirtalis.samples.W, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="background_test/SE_W/")
# write.csv(bg.SE.W$reps.overlap,"bg_southeast_west.csv")
# 
# bg.SE.C <- background.test(sirtalis.samples.SE,sirtalis.samples.C, env=env, type="mx", nreps=100, low.memory = TRUE, rep.dir="background_test/SE_C/")
# write.csv(bg.SE.C$reps.overlap,"bg_southeast_central.csv")

#####################################################################
#####################################################################
#####################################################################

# Pliocene: mid-Pliocene warm period (3.205 Ma), v1.0*
mpwp_2.5m <- load_paleoclim("data/mPWP_v1_r2_5m.zip")
mpwp_2.5m.NA <- crop(mpwp_2.5m,NorthAmerica)
mpwp.vars <- removeCollinearity.terra(mpwp_2.5m.NA, multicollinearity.cutoff = .75, 
                                         select.variables = FALSE, sample.points = FALSE,
                                         plot=FALSE)
mpwp <- mpwp_2.5m.NA[[c("bio_1", "bio_4", "bio_12", "bio_14", "bio_15", "bio_18")]]
mpwp <- setMinMax(mpwp)
names(mpwp) <- c("bio1", "bio4", "bio12", "bio14", "bio15", "bio18")

# Pleistocene: MIS19 (ca. 787 ka), v1.0*
mis19_2.5m <- load_paleoclim("data/MIS19_v1_r2_5m.zip")
mis19_2.5m.NA <- crop(mis19_2.5m,NorthAmerica)
mis19.vars <- removeCollinearity(mis19_2.5m.NA, multicollinearity.cutoff = .75, 
                                select.variables = FALSE, sample.points = FALSE,
                                plot=FALSE)
mis19 <- mis19_2.5m.NA[[c("bio_1", "bio_4", "bio_12", "bio_14", "bio_15", "bio_18")]]
mis19 <- setMinMax(mis19)
names(mis19) <- c("bio1", "bio4", "bio12", "bio14", "bio15", "bio18")

# Pleistocene: Last Interglacial (ca. 130 ka), v1.0
lig_2.5m <- load_paleoclim("data/LIG_v1_2_5m.zip")
lig_2.5m.NA <- crop(lig_2.5m,NorthAmerica)
lig.vars <- removeCollinearity(lig_2.5m.NA, multicollinearity.cutoff = .75, 
                                 select.variables = FALSE, sample.points = FALSE,
                                 plot=FALSE)
lig <- lig_2.5m.NA[[c("bio_1", "bio_4", "bio_12", "bio_14", "bio_15", "bio_18")]]
lig <- setMinMax(lig)
names(lig) <- c("bio1", "bio4", "bio12", "bio14", "bio15", "bio18")

# Pliocene: M2 (ca. 3.3 Ma), v1.0*
m2_2.5m <- load_paleoclim("data/M2_v1_r2_5m.zip")
m2_2.5m.NA <- crop(m2_2.5m,NorthAmerica)
m2.vars <- removeCollinearity(m2_2.5m.NA, multicollinearity.cutoff = .75, 
                               select.variables = FALSE, sample.points = FALSE,
                               plot=FALSE)
m2 <- m2_2.5m.NA[[c("bio_1", "bio_4", "bio_12", "bio_14", "bio_15", "bio_18")]]
m2 <- setMinMax(m2)
names(m2) <- c("bio1", "bio4", "bio12", "bio14", "bio15", "bio18")

# Pleistocene: Last Glacial Maximum (ca. 21 ka), v1.2b**, NCAR CCSM4
lgm_2.5m <- load_paleoclim("data/chelsa_LGM_v1_2B_r2_5m.zip")
lgm_2.5m.NA <- crop(lgm_2.5m,NorthAmerica)
lgm.vars <- removeCollinearity(lgm_2.5m.NA, multicollinearity.cutoff = .75, 
                              select.variables = FALSE, sample.points = FALSE,
                              plot=FALSE)
lgm <- lgm_2.5m.NA[[c("bio_1", "bio_4", "bio_12", "bio_14", "bio_15", "bio_18")]]
lgm <- setMinMax(lgm)
names(lgm) <- c("bio1", "bio4", "bio12", "bio14", "bio15", "bio18")

# SDMs w/ reduced paleoclim variables
sirtalis.samples.SE.max.sb <- enmtools.maxent(species = sirtalis.samples.SE, bg.source = "points", env = env.sb, test.prop = .25)
SE.sb.plot <- plot.enmtools.maxent(sirtalis.samples.SE.max.sb)

sirtalis.samples.E.max.sb <- enmtools.maxent(species = sirtalis.samples.E, bg.source = "points", env = env.sb, test.prop = .25)
E.sb.plot <- plot.enmtools.maxent(sirtalis.samples.E.max.sb)

sirtalis.samples.C.max.sb <- enmtools.maxent(species = sirtalis.samples.C, bg.source = "points", env = env.sb, test.prop = .25)
C.sb.plot <- plot.enmtools.maxent(sirtalis.samples.C.max.sb)

sirtalis.samples.W.max.sb <- enmtools.maxent(species = sirtalis.samples.W, bg.source = "points", env = env.sb, test.prop = .25)
W.sb.plot <- plot.enmtools.maxent(sirtalis.samples.W.max.sb)

sirtalis.samples.max.sb <- enmtools.maxent(species = sirtalis.samples, bg.source = "points", env = env.sb, test.prop = .25)
samples.sb.plot <- plot.enmtools.maxent(sirtalis.samples.max.sb)


# model projections onto historic environments
SE.lgm <- dismo::predict(sirtalis.samples.SE.max.sb$model, lgm)
SE.lgm.df <- as.data.frame(SE.lgm,xy=T)
SE.m2 <- dismo::predict(sirtalis.samples.SE.max.sb$model, m2)
SE.m2.df <- as.data.frame(SE.m2,xy=T)
SE.lig <- dismo::predict(sirtalis.samples.SE.max.sb$model, lig)
SE.lig.df <- as.data.frame(SE.lig,xy=T)
SE.mis19 <- dismo::predict(sirtalis.samples.SE.max.sb$model, mis19)
SE.mis19.df <- as.data.frame(SE.mis19,xy=T)
SE.mpwp <- dismo::predict(sirtalis.samples.SE.max.sb$model, mpwp)
SE.mpwp.df <- as.data.frame(SE.mpwp,xy=T)

E.lgm <- dismo::predict(sirtalis.samples.E.max.sb$model, lgm)
E.lgm.df <- as.data.frame(E.lgm,xy=T)
E.m2 <- dismo::predict(sirtalis.samples.E.max.sb$model, m2)
E.m2.df <- as.data.frame(E.m2,xy=T)
E.lig <- dismo::predict(sirtalis.samples.E.max.sb$model, lig)
E.lig.df <- as.data.frame(E.lig,xy=T)
E.mis19 <- dismo::predict(sirtalis.samples.E.max.sb$model, mis19)
E.mis19.df <- as.data.frame(E.mis19,xy=T)
E.mpwp <- dismo::predict(sirtalis.samples.E.max.sb$model, mpwp)
E.mpwp.df <- as.data.frame(E.mpwp,xy=T)

C.lgm <- dismo::predict(sirtalis.samples.C.max.sb$model, lgm)
C.lgm.df <- as.data.frame(C.lgm,xy=T)
C.m2 <- dismo::predict(sirtalis.samples.C.max.sb$model, m2)
C.m2.df <- as.data.frame(C.m2,xy=T)
C.lig <- dismo::predict(sirtalis.samples.C.max.sb$model, lig)
C.lig.df <- as.data.frame(C.lig,xy=T)
C.mis19 <- dismo::predict(sirtalis.samples.C.max.sb$model, mis19)
C.mis19.df <- as.data.frame(C.mis19,xy=T)
C.mpwp <- dismo::predict(sirtalis.samples.C.max.sb$model, mpwp)
C.mpwp.df <- as.data.frame(C.mpwp,xy=T)

W.lgm <- dismo::predict(sirtalis.samples.W.max.sb$model, lgm)
W.lgm.df <- as.data.frame(W.lgm,xy=T)
W.m2 <- dismo::predict(sirtalis.samples.W.max.sb$model, m2)
W.m2.df <- as.data.frame(W.m2,xy=T)
W.lig <- dismo::predict(sirtalis.samples.W.max.sb$model, lig)
W.lig.df <- as.data.frame(W.lig,xy=T)
W.mis19 <- dismo::predict(sirtalis.samples.W.max.sb$model, mis19)
W.mis19.df <- as.data.frame(W.mis19,xy=T)
W.mpwp <- dismo::predict(sirtalis.samples.W.max.sb$model, mpwp)
W.mpwp.df <- as.data.frame(W.mpwp,xy=T)

total.lgm <- dismo::predict(sirtalis.samples.max.sb$model, lgm)
total.lgm.df <- as.data.frame(total.lgm,xy=T)
total.m2 <- dismo::predict(sirtalis.samples.max.sb$model, m2)
total.m2.df <- as.data.frame(total.m2,xy=T)
total.lig <- dismo::predict(sirtalis.samples.max.sb$model, lig)
total.lig.df <- as.data.frame(total.lig,xy=T)
total.mis19 <- dismo::predict(sirtalis.samples.max.sb$model, mis19)
total.mis19.df <- as.data.frame(total.mis19,xy=T)
total.mpwp <- dismo::predict(sirtalis.samples.max.sb$model, mpwp)
total.mpwp.df <- as.data.frame(total.mpwp,xy=T)


plot.prediction <- function(input, ...){
  ggplot() +
    geom_polygon(data = input, aes(x=x, y), colour = "white", fill = "white", size = 0.25) +
    geom_raster(input,mapping = aes(x = x, y = y, fill = maxent)) +
    scale_fill_viridis(option = "B", guide = guide_colourbar(title = "Suitability"),na.value = "white") +
    coord_fixed() + 
    xlab("Longitude") +
    ylab("Latitude") +
    theme(
      axis.line = element_line(),
      axis.text.y = element_text(colour = "black", size = 10),
      axis.text.x = element_text(colour = "black", size = 10),
      legend.justification = c(0, 1), legend.position = c(.90, .35),
      panel.background = element_blank()
    )
}

plot.prediction(SE.lgm.df)
plot.prediction(SE.m2.df)
plot.prediction(SE.lig.df)
plot.prediction(SE.mis19.df)
plot.prediction(SE.mpwp.df)

plot.prediction(W.lgm.df)
plot.prediction(W.m2.df)
plot.prediction(W.lig.df)
plot.prediction(W.mis19.df)
plot.prediction(W.mpwp.df)

plot.prediction(C.lgm.df)
plot.prediction(C.m2.df)
plot.prediction(C.lig.df)
plot.prediction(C.mis19.df)
plot.prediction(C.mpwp.df)

plot.prediction(E.lgm.df)
plot.prediction(E.m2.df)
plot.prediction(E.lig.df)
plot.prediction(E.mis19.df)
plot.prediction(E.mpwp.df)

plot.prediction(total.lgm.df)
plot.prediction(total.m2.df)
plot.prediction(total.lig.df)
plot.prediction(total.mis19.df)
plot.prediction(total.mpwp.df)
