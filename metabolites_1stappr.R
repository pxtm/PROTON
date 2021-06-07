#############################################################################################################
################################# PROTON metabolites first approach #########################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 26.05.2021
## last edit: 27.05.2021

## libraries
library(tidyverse)
library(factoextra)
library(mdatools)
library(rsample)
library(WGCNA)

## functions
source('lm_function.R')

## data
# lipids <- read.table('Metabolomics/Lipidomics/edited_lipid_trans_scaled_imp.txt', header=T, row.names=1, sep='\t')
# gc <- read.table('Metabolomics/GCGC/edited_PROTON_GCGC_trans_scale_imp.txt', header=T, row.names=1, sep='\t')
# mtdt <- read.table('mtdt_selected.txt', row.names=1, header=T, sep='\t')
# 
# # data edit - sample 36204 needs to be removed because it has not cell count data
# common.samples <- Reduce(intersect, list(row.names(mtdt), row.names(lipids), row.names(gc)))
# mtdt <- mtdt[row.names(mtdt)%in%common.samples,]
# lipids <- lipids[row.names(lipids)%in%common.samples,]
# gc <- gc[row.names(gc)%in%common.samples,]
# 
# identical(row.names(lipids), row.names(mtdt))
# identical(row.names(gc), row.names(mtdt))
# lipids <- lipids[match(row.names(mtdt), row.names(lipids)),]

# saveRDS(gc, 'GCGC_mets.rds')
# saveRDS(lipids, 'Lipids_mets.rds')
# saveRDS(mtdt, 'mtdt_common.rds')

mtdt <- readRDS('mtdt_common.rds')
gc <- readRDS('GCGC_mets.rds')
lipids <- readRDS('Lipids_mets.rds')
identical(row.names(mtdt), row.names(gc))
identical(row.names(mtdt), row.names(lipids))

## metabolomics data correction
gc <- cbind(gc, 'egfr' = mtdt$eGFR)
residuals.df <- list()
for (met in seq(1:398)){
  model <- lm(gc[,met] ~ egfr, data = gc)
  residuals.met <- residuals(model)
  residuals.df[[met]]<- residuals.met
}

gc.res <- do.call(cbind.data.frame, residuals.df)
names(gc.res) <- names(gc)[1:398]

lipids <- cbind(lipids, 'egfr' = mtdt$eGFR)
residuals.df <- list()
for (met in seq(1:7470)){
  model <- lm(lipids[,met] ~ egfr, data = lipids)
  residuals.met <- residuals(model)
  residuals.df[[met]]<- residuals.met
}

lipids.res <- do.call(cbind.data.frame, residuals.df)
names(lipids.res) <- names(lipids)[1:7470]

saveRDS(lipids.res, 'residuals_lipids_egfr.rds')
saveRDS(gc.res, 'residuals_GCGC_egfr.rds')

## lipidomics
pca.lipid <- prcomp(lipids.res)
fviz_screeplot(pca.lipid)
fviz_pca_ind(pca.lipid, geom.ind = 'point', col.ind = as.factor(mtdt$Groups), addEllipses = T, ellipse.level=.95)

# GC (I guess it means targeted?)
pca.gc <- prcomp(gc.res)
fviz_screeplot(pca.gc)
fviz_pca_ind(pca.gc, geom.ind = 'point', col.ind = mtdt$age)
fviz_pca_ind(pca.gc, geom.ind = 'point', col.ind = as.factor(mtdt$Groups), addEllipses = T, ellipse.level=.95)

grouping.annot <- data.frame('Group'= mtdt$Groups, 'DM' = mtdt$dm); row.names(grouping.annot) <- row.names(mtdt)
pheatmap::pheatmap(gc.res, annotation_row = grouping.annot)

gc.avg <- gc.res %>% 
  add_column('group' = as.factor(mtdt$Groups)) %>% 
  group_by(group) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  t() %>% 
  as.data.frame()

names(gc.avg) <- gc.avg[1,]
gc.avg <- gc.avg[-1,]
gc.avg[] <- as.data.frame(sapply(gc.avg, as.numeric))
gc.avg$met <- row.names(gc.avg)
gc.avg.m <- reshape2::melt(gc.avg)
gc.avg.m$variable <- factor(gc.avg.m$variable, levels=c('Controls', 'Normo', 'Micro', 'Macro'))
ggplot(gc.avg.m, aes(variable, met, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low='coral2', high='darkred')

## PLS-DA (supervised)
mtdt$ctl <- plyr::revalue(mtdt$Groups, c('Normo' = 'T1D', 'Micro' = 'T1D', 'Macro' = 'T1D'))
gc.split <- initial_split(gc.res, prop = .7)
gc.train <- training(gc.split)
gc.test <- testing(gc.split)

mtdt.gc.train <- mtdt[row.names(mtdt)%in%row.names(gc.train),]
mtdt.gc.test <- mtdt[row.names(mtdt)%in%row.names(gc.test),]

gc.plsda.t1d <- plsda(gc.train, mtdt.gc.train$ctl, 3, cv = 1)
summary(gc.plsda.t1d)
plotPredictions(gc.plsda.t1d)
res <- predict(gc.plsda.t1d, gc.test, mtdt.gc.test$t1d)
summary(res)
plotPredictions(res)

mtdt.t1d <- mtdt[mtdt$Groups!='Controls',]
mtdt.t1d$normo <- plyr::revalue(mtdt.t1d$Groups, c('Micro' = 'noNormo', 'Macro' = 'noNormo'))
mtdt.t1d$micro <- plyr::revalue(mtdt.t1d$Groups, c('Normo' = 'noMicro', 'Macro' = 'noMicro'))
mtdt.t1d$macro <- plyr::revalue(mtdt.t1d$Groups, c('Normo' = 'noMacro', 'Micro' = 'noMacro'))

gc.res.t1d <- gc.res[row.names(gc.res)%in%row.names(mtdt.t1d),]

# normo
gc.t1d.split <- initial_split(gc.res.t1d, prop = .7)
gc.t1d.train <- training(gc.t1d.split)
gc.t1d.test <- testing(gc.t1d.split)

mtdt.t1d.gc.train <- mtdt.t1d[row.names(mtdt.t1d)%in%row.names(gc.t1d.train),]
mtdt.t1d.gc.test <- mtdt.t1d[row.names(mtdt.t1d)%in%row.names(gc.t1d.test),]

gc.plsda.t1d.normo <- plsda(gc.t1d.train, mtdt.t1d.gc.train$normo, 3, cv = 1)
summary(gc.plsda.t1d.normo)
plotPredictions(gc.plsda.t1d.normo)
getConfusionMatrix(gc.plsda.t1d.normo$calres)
res <- predict(gc.plsda.t1d.normo, gc.t1d.test, mtdt.t1d.gc.test$normo)
summary(res)
plotPredictions(res)
out.normo <- as.data.frame(lm.associations(mtdt.t1d, gc.res.t1d, 'normo', control_by = c('age', 'sex', 'bmi', 'diet')))
out.normo[out.normo$normo.fdr<.1,]


## lipidomics ----
# reduce dimensionality --- cluster them!
allowWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(lipids, powerVector = powers, verbose = 5)
#Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 4
adjacency = adjacency(lipids, power = softPower)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04)
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(lipids, colors = dynamicColors)
MEs = MEList$eigengenes
saveRDS(MEs, 'lipids_clusters.rds')
annot.row <- data.frame('group' = mtdt$Groups, 'sex'=mtdt$sex); row.names(annot.row) <- row.names(mtdt)
pheatmap::pheatmap(MEs, annotation_row = annot.row)
