rm(list=ls())
currentdir <- getwd()
setwd(currentdir)

library(easypackages)
libraries("tess3r","vcfR","LEA")

# note: verify shared order of coordinates file and .vcf!
# ipyRAD vcfs are alphanumeric

# Convert from .vcf format to .geno format
vcf2geno("sirtalis.vcf","sirtalis.geno")
geno2lfmm("sirtalis.geno","sirtalis.lfmm")

# Now you should have a .lfmm file named with your base_file_name string. Missing data = '9'
# from the conversion, but you need to change from '9' to 'NA'
sirtalis_data <- read.table("sirtalis.lfmm")
sirtalis_data <- as.matrix(sirtalis_data)
mydim <- dim(sirtalis_data)
sirtalis_data <- as.vector(sirtalis_data)
sirtalis_data[which(sirtalis_data==9)] <- NA
sirtalis_data <- matrix(data=sirtalis_data,nrow=mydim[1],byrow=F)
sirtalis_data <- as.data.frame(sirtalis_data)

# import data
# Coordinates file must be in long,lat order without individual names

vcf <- read.vcfR(file = "sirtalis.vcf")
genlight <- vcfR2genlight(vcf,n.cores = 8)
locality.info <- read.csv("samples.csv",stringsAsFactors = F)
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
tess3.obj <- tess3(sirtalis_data,XProba = NULL,coord = coordinates, K = 1:10, rep = 20,
  keep = "all", openMP.core.num = 8, method = "projected.ls", ploidy = 2)
# save object to file
save(tess3.obj,file="tess3_test.robj")

# cross validation plot
plot(tess3.obj, pch = 16, col = "blue",cex=2,
     xlab = "Number of populations (K)",
     ylab = "Cross-validation score")

# Choose your k value to plot
tess.q.matrix.k4 <- qmatrix(tess3.obj, K = 4)
write.table(x = tess.q.matrix.k4,"tess.q.matrix.k4.q", sep=" ", row.names=FALSE, col.names=FALSE)

# save cluster assignment %s to file for current K
m <- tess.q.matrix.k4
rownames(m) <- names
ca <- round(m, 4)
cluster.assignments <- write.csv(ca, file = "k4_cluster_assignments.csv")
