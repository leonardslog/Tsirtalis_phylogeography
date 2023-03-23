library(easypackages)
libraries("maps","rgdal","raster","RColorBrewer","scales",
          "maptools","gridExtra","rgeos","geodata","terra",
          "rnaturalearthdata","rnaturalearth","vcfR","plotrix")

# import data
# Coordinates file must be in long,lat order without individual names

vcf <- read.vcfR(file = "sirtalis.vcf")
genlight <- vcfR2genlight(vcf,n.cores = 8)
locality.info <- read.csv("sirtalis_samples.csv",stringsAsFactors = F)
head(locality.info)
names <- as.matrix(genlight@ind.names)
tess.q.matrix.k4 <- read.table("tess.q.matrix.k4.q", sep=" ")
head(tess.q.matrix.k4)

# build dataframe w/ genlight names and coordinates
genlight.coords <- setNames(data.frame(matrix(ncol = 3, nrow = length(genlight@ind.names))), c("Sample", "x", "y"))

for (i in 1:length(genlight$ind.names)) {
  genlight.coords[i,1:3] <- c(
    genlight@ind.names[i],
    subset(locality.info$longitude[], locality.info$Sample==genlight@ind.names[i]),
    subset(locality.info$latitude[], locality.info$Sample==genlight@ind.names[i]))
}
head(genlight.coords)

m <- tess.q.matrix.k4
rownames(m) <- names
ca <- round(m, 4)

# coerce coordinates to double vector to sidestep error
coordinates <- as.matrix(cbind(
  as.double(genlight.coords$x), as.double(genlight.coords$y)))
head(coordinates)

# prep qmatrix and coordinate data
sp.qmatrix.k4 <- cbind(ca, coordinates)
head(sp.qmatrix.k4)
nclust <- ncol(tess.q.matrix.k4)

# k4 pop order    C           S         E           W
div.colors <- c("#ABD9E9","#2C7BB6", "#D7191C", "#FDAE61") 
div.colors2 <- c("#FDAE61","#D7191C", "#2C7BB6", "#ABD9E9") # adjusted color order

#########################
# North America basemap #
#########################

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
crs(srtm.NA) <- crs.universal


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

# remove border btwn Saginaw Bay and Lake Huron
SaginawBay <- subset(lakes, name == "Saginaw Bay")
LakeHuron <- subset(lakes, name == "Lake Huron")
SBxSL <- rbind(SaginawBay,LakeHuron)
SBxSL_noborders <- gUnaryUnion(SBxSL)

# plot map
over.col <- colorRampPalette(c("white", "black"))
terra::plot(srtm.NA,
            col=alpha(over.col(100),0.4),
            ylim=c(24.5, 54),xlim=c(-130,-59.5),
            mar=c(11,8,10,8)-3,legend=F,xlab="Longtitude",ylab="Latitude")
plot(northamer_crop,lwd=1, col=alpha("dark grey",0.3), border="dark grey",add=T)
# plot(lakes.centerlines,xlim=x, ylim=y,col="white",border=alpha("dark grey",0.9),lwd=1,add=T)
plot(subset(US.rivers, Name %in% target_rivers), col=alpha("blue",0.9),lwd=1, add=T)
plot(subset(riverslakes.centerlines, name %in% target_rivers2), col=alpha("blue",0.9),lwd=1, add=T)
plot(lakes.centerlines,xlim=x, ylim=y,col="white",border=alpha("dark grey",0.9),lwd=1,add=T)
plot(SBxSL_noborders,xlim=x, ylim=y,col="white",border=alpha("dark grey",0.9),lwd=1,add=T)

# add pie chart: tess.k4
apply(sp.qmatrix.k4,1, function(z) {
  zz <- data.matrix(z[1:nclust])
  index <- which(zz !=0)
  plotrix::floating.pie(xpos = z[nclust+1], ypos = z[nclust+2], x = zz[index],
               radius = .7, col = alpha(div.colors2[index],0.8))})
