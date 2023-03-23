rm(list = ls(all.names = TRUE))
currentdir <- getwd()
setwd(currentdir)

library(easypackages)
libraries("hierfstat","adegenet")

auto <- read.genepop("auto.gen")
reduced <- read.genepop("reduced.gen")
z <- read.genepop("z.gen")
total <- read.genepop("total.gen")
east <- read.genepop("east.gen")
southeast <- read.genepop("southeast.gen")
central <- read.genepop("central.gen")
west <- read.genepop("west.gen")

iqtree.z <- read.genepop("iqtree_z.gen")
iqtree.auto <- read.genepop("iqtree_autosomes.gen")


total.gd <- genet.dist(total,diploid=TRUE,method="WC84")
total.bs <- basic.stats(total)

reduced.gd <- genet.dist(reduced,diploid=TRUE,method="WC84")
reduced.bs <- basic.stats(reduced)

auto.gd <- genet.dist(auto,diploid=TRUE,method="WC84")
auto.bs <- basic.stats(auto)

z.gd <- genet.dist(z,diploid=FALSE,method="WC84")
z.bs <- basic.stats(z)

east.bs <- basic.stats(east)
southeast.bs <- basic.stats(southeast)
central.bs <- basic.stats(central)
west.bs <- basic.stats(west)


iqtree.z.gd <- genet.dist(iqtree.z,diploid=TRUE,method="WC84")
iqtree.z.bs <- basic.stats(iqtree.z)

iqtree.auto.gd <- genet.dist(iqtree.auto,diploid=TRUE,method="WC84")
iqtree.auto.bs <- basic.stats(iqtree.auto)

