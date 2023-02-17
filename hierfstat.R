rm(list = ls(all.names = TRUE))
setwd("~/Documents/sirtalis/hierfstat/")
library(easypackages)
libraries("rgdal", "raster","poppr","RColorBrewer","maps",
          "mapdata","plotrix","adegenet","calibrate",
          "viridis","mapproj")
library(hierfstat)

gen <- read.genepop("populations.snps.gen")
vcf <- read.VCF("sirtalis_allsites.vcf.gz",BiAllelic = T)
isPoly(gen)
isPoly(vcf)
vcfdat <- with(vcf@snps,table(A1,A2))
genet.dist(vcfdat,diploid=TRUE,method="Fst")
vcf
