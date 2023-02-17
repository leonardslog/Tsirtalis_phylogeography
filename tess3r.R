rm(list=ls())
currentdir <- getwd()
setwd(currentidr)

library(easypackages)
libraries("BiocManager","tess3r","dismo", "adegenet", "devtools","maps","mapproj","maptools",
          "raster","rgdal","sf","scales","dplyr","ggplot2","plotrix","poppr",
          "mapdata","adegenet","RColorBrewer","LEA","vcfR","rnaturalearth",
          "rgeos")

# note: verify shared order of coordinates file and .vcf!
# ipyRAD vcf's are alphanumeric

# Convert from .vcf format to .geno format
vcf2geno("sirtalis_tesspops_p1r80.vcf","sirtalis_usnps_p1r80.geno")
geno2lfmm("sirtalis_usnps_p1r80.geno","sirtalis_usnps_p1r80.lfmm")

#Now you should have a .lfmm file named with your base_file_name string. Missing data = '9'
#from the conversion, but you need to change from '9' to 'NA'
total_p1r80_usnps_data <- read.table("sirtalis_usnps_p1r80.lfmm")
total_p1r80_usnps_data <- as.matrix(total_p1r80_usnps_data)
mydim <- dim(total_p1r80_usnps_data)
total_p1r80_usnps_data <- as.vector(total_p1r80_usnps_data)
total_p1r80_usnps_data[which(total_p1r80_usnps_data==9)] <- NA
total_p1r80_usnps_data <- matrix(data=total_p1r80_usnps_data,nrow=mydim[1],byrow=F)
total_p1r80_usnps_data <- as.data.frame(total_p1r80_usnps_data)

# import data
# Coordinates file must be in long,lat order without individual names

vcf <- read.vcfR(file = "sirtalis_tesspops_p1r80.vcf")
genlight <- vcfR2genlight(vcf,n.cores = 8)
locality.info <- read.csv("sirtalis_samples.csv",stringsAsFactors = F)
head(locality.info)
names <- as.matrix(genlight@ind.names)

# build dataframe w/ genlight names and coordinates

genlight.coords <- setNames(data.frame(matrix(ncol = 3, nrow = length(genlight@ind.names))), c("Sample", "x", "y"))

for (i in 1:length(genlight$ind.names)) {
  genlight.coords[i,1:3] <- c(
    genlight@ind.names[i],
    subset(locality.info$longitude[], locality.info$Sample==genlight@ind.names[i]),
    subset(locality.info$latitude[], locality.info$Sample==genlight@ind.names[i]))
}
head(genlight.coords)

# coerce coordinates to double vector to sidestep error
coordinates <- as.matrix(cbind(
  as.double(genlight.coords$x), as.double(genlight.coords$y)))

# run tess3 
# tess3.obj <- tess3(total_p1r80_usnps_data,XProba = NULL,coord = coordinates, K = 1:10, rep = 20,
                   # keep = "all", openMP.core.num = 8, method = "projected.ls", ploidy = 2)
# save object to file
# saveRDS(tess3.obj, file = "tess3.obj.rds")
# save(tess3.obj,file="tess3.robj")

# if revisiting this script, just reload the tess object
load("tess3.robj")

# cross validation plot
plot(tess3.obj, pch = 16, col = "blue",cex=2,
     xlab = "Number of populations (K)",
     ylab = "Cross-validation score")

# Choose your k value to plot
tess.q.matrix.k4 <- qmatrix(tess3.obj, K = 4)
# write.table(x = tess.q.matrix.k4,"tess.q.matrix.k4.q", sep=" ", row.names=FALSE, col.names=FALSE)

# save cluster assignment %s to file for current K
m <- tess.q.matrix.k4
rownames(m) <- names
ca <- round(m, 4)
# cluster.assignments <- write.csv(ca, file = "k4_cluster_assignments.csv")


############################
###     admixture map    ###
############################

# base map settings
# crs.universal <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs.universal <- "+proj=robin +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
n_america.extent <- extent(-130,-59.5,24.5,54)

# # read in water bodies and set CRS
# riverslakes.centerlines <- readOGR("ne_10m_rivers_lake_centerlines/ne_10m_rivers_lake_centerlines.shp")
# US.rivers <- readOGR("USA_Rivers_and_Streams-shp/")
# save(US.rivers,file="US.rivers.robj")
# lakes.centerlines <- readOGR("ne_10m_lakes/ne_10m_lakes.shp")
# save(lakes.centerlines,file="lakes.centerlines.robj")
# rivers <- crop(riverslakes.centerlines,n_america.extent)
# save(rivers,file="rivers.cropped.robj")
# lakes <- crop(lakes.centerlines, n_america.extent)
# save(lakes,file="lakes.cropped.robj")

# # one way to get North America boundaries and combine + crop
# gadm.USA <- getData('GADM', country='USA',level=0)
# gadm.MEX <- getData('GADM', country='MEX',level=0)
# gadm.CAN <- getData('GADM', country='CAN', level=0)
# gadm.NorthAmer <- rbind(gadm.CAN,gadm.USA,gadm.MEX)
# NorthAmer <- crop(gadm.NorthAmer, n_america.extent)
# save(NorthAmer,file="NorthAmer.cropped.robj")

# # load and crop SRTM elevation data (.tif prepped as virtual raster in QGIS)
# NorthAmer.srtm <- crop(raster("NorthAmer_srtm.tif"), n_america.extent)
# save(NorthAmer.srtm,file="NorthAmer.srtm.cropped.robj")

# load prepped data
load("NorthAmer.srtm.cropped.robj")
load("NorthAmer.cropped.robj")
load("US.rivers.robj")
load("rivers.cropped.robj")
load("lakes.cropped.robj")

crs(NorthAmer.srtm) <- crs.universal
crs(NorthAmer) <- crs.universal
crs(lakes) <- crs.universal
crs(rivers) <- crs.universal
crs(US.rivers) <- crs.universal


# subset rivers and lakes
target_rivers <- c("Rio Grande","Mississippi River", 
                   "Ohio River", "Arkansas River", "Alabama River",
                   "Mobile River", "Coosa River","Flint River",
                   "Apalachicola River", "Chattahoochee River")
target_rivers2 <- c("Saskatchewan","Columbia","Missouri")

# plot base map
over.col <- colorRampPalette(c("white", "black"))
par(mar=c(14,5,5,6)-0.5)
plot(NorthAmer.srtm,legend=F, axes=T, box=T, col=over.col(100),xlab = "Longitude", ylab = "Latitude")
plot(NorthAmer,col=alpha("black",0.5),border=NA,add=T)
plot(lakes,xlim=x, ylim=y,col="white",lty="blank",add=T)
plot(subset(US.rivers, Name %in% target_rivers), col=alpha("dark grey",0.9),lwd=1, add=T)
plot(subset(rivers, name %in% target_rivers2), col=alpha("dark grey",0.9),lwd=1, add=T)
plot(lakes,xlim=x, ylim=y,col="white",lty="blank",add=T) # Saskatchewan River overlaps Lake Winnipeg

# prep qmatrix and coordinate data
sp.qmatrix.k4 <- cbind(ca, coordinates)
nclust <- ncol(tess.q.matrix.k4)

# k4 pop order    C           S         E           W
div.colors <- c("#ABD9E9","#2C7BB6", "#D7191C", "#FDAE61") 
div.colors2 <- c("#FDAE61","#D7191C", "#2C7BB6", "#ABD9E9") # adjusted color order

# add pie chart: tess.k4
apply(sp.qmatrix.k4,1, function(z) {
  zz <- data.matrix(z[1:nclust])
  index <- which(zz !=0)
  floating.pie(xpos = z[nclust+1], ypos = z[nclust+2], x = zz[index],
               radius = .7, col = alpha(div.colors2[index],0.8))})
