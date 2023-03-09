rm(list=ls())
setwd("~/Documents/sirtalis/pca/")
library(easypackages)
libraries("adegenet","RColorBrewer","parallel","RColorBrewer", "ggplot2")


###########################
# import and convert data #
###########################

locality.info <- read.csv("sample_summary.csv",stringsAsFactors = T)
head(locality.info)

total.vcf <- vcfR::read.vcfR("~/Documents/sirtalis/stacks/ele_map/sirtalis/p1r80/populations.snps.vcf")
# check that all samples in vcf and sample info are represented
all(colnames(total.vcf@gt)[-1]== locality.info$sample) 
total.gl <- vcfR::vcfR2genlight(total.vcf)
ploidy(total.gl) <- 2
pop(total.gl) <- locality.info$population

z.vcf <- vcfR::read.vcfR("~/Documents/sirtalis/stacks/ele_map/sirtalis/z_p1r80/populations.snps.vcf")
# check that all samples in vcf and sample info are represented
all(colnames(z.vcf@gt)[-1]== locality.info$sample) 
z.gl <- vcfR::vcfR2genlight(z.vcf)
ploidy(z.gl) <- 2
pop(z.gl) <- locality.info$population

autosomes.vcf <- vcfR::read.vcfR("~/Documents/sirtalis/stacks/ele_map/sirtalis/autosomes_p1r80/populations.snps.vcf")
# check that all samples in vcf and sample info are represented
all(colnames(autosomes.vcf@gt)[-1]== locality.info$sample) 
autosomes.gl <- vcfR::vcfR2genlight(autosomes.vcf)
ploidy(autosomes.gl) <- 2
pop(autosomes.gl) <- locality.info$population


# div colors  "#2C7BB6", "#FDAE61", "#D7191C", "#ABD9E9"
# K = 4 pop order      E          C           S         W
# div.colors <- c("#2C7BB6", "#FDAE61", "#D7191C", "#ABD9E9")
# colors incorrectly assigned in figure
div.colors <- c("#FDAE61", "#2C7BB6", "#D7191C", "#ABD9E9")

#############
# TOTAL PCA #
#############

total.pca <- glPca(total.gl, n.cores = 8, nf = 3, parallel = TRUE)
total.pca.scores <- as.data.frame(total.pca$scores)
total.pca.scores$Population <- pop(total.gl)

set.seed(9)
p <- ggplot(total.pca.scores, aes(x=PC1, y=PC2, colour=Population)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = div.colors) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p <- p + xlab("PC1 (12.11%)") + ylab("PC2 (9.48%)")
p

####################
# Z CHROMOSOME PCA #
####################

z.pca <- glPca(z.gl, n.cores = 8, nf = 3, parallel = TRUE)
z.pca.scores <- as.data.frame(z.pca$scores)
z.pca.scores$Population <- pop(z.gl)

p <- ggplot(z.pca.scores, aes(x=PC1, y=PC2, colour=Population)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = div.colors) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p <- p + xlab("PC1 (1.86%)") + ylab("PC2 (0.97%)")
p

################
# AUTOSOME PCA #
################

autosomes.pca <- glPca(autosomes.gl, n.cores = 8, nf = 3, parallel = TRUE)
autosomes.pca.scores <- as.data.frame(autosomes.pca$scores)
autosomes.pca.scores$Population <- pop(autosomes.gl)

p <- ggplot(autosomes.pca.scores, aes(x=PC1, y=PC2, colour=Population)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = div.colors) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p <- p + xlab("PC1 (9.34%)") + ylab("PC2 (8.11%)")
p
