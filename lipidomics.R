#############################################################################################################
#################################### PROTON lipids Controls vs T1D ##########################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 05.07.2021
## last edit: 

## libraries
library(tidyverse)
library(factoextra)
library(mdatools)
library(rsample)
library(WGCNA)
library(ggpubr)
library(ropls)
library(ggbeeswarm)
library(grid)
library(ggrepel)
library(tayloRswift)
library(OptimalCutpoints)
library(effsize)
library(randomForest)

## functions
source('lm_function.R')

## data
mtdt <- readRDS('mtdt_common.rds')
lipids.res <- readRDS('residuals_lipids_egfr.rds')

row.names(mtdt) == row.names(lipids.res)

##### Controls vs T1D ------------------------------------------------------------------------------------
## PCA for controls vs T1D 
mtdt$T1D <- as.factor(plyr::revalue(mtdt$Groups, c('Normo' = 'T1D', 'Micro' = 'T1D', 'Macro' = 'T1D')))
pca.lipids.res <- prcomp(lipids.res, scale. = T, center = T)
fviz_screeplot(pca.lipids.res)
fviz_pca_ind(pca.lipids.res, geom = 'point', col.ind = mtdt$T1D, addEllipses = T)
scores.lipids.res <- as.data.frame(pca.lipids.res$x)
row.names(scores.lipids.res) == row.names(mtdt)
scores.lipids.res$T1D <- mtdt$T1D
p1 <- ggplot(scores.lipids.res, aes(PC1, PC2, color = T1D)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_taylor() +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 15.1%)", 
       y = "Principal Component 2 (var.expl. 8.3%)")

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## PLS-DA for controls vs T1D 
lip.split <- initial_split(lipids.res, prop = .7)
lip.train <- training(lip.split)
mylist <-list('total' = row.names(lipids.res), 'train' = row.names(lip.train))
common_values <- intersect(row.names(lipids.res), row.names(lip.train))
train.index <- lapply(mylist, function(x) which(x %in% common_values))$total

# initial fit
t1d.pls <- opls(lipids.res, mtdt$T1D, subset = train.index)
summary(t1d.pls)

# check predictions on the training subset
trainVi <- getSubsetVi(t1d.pls)
table(mtdt$T1D[trainVi], fitted(t1d.pls))
table(mtdt$T1D[-trainVi], predict(t1d.pls, lipids.res[-trainVi,]))

## variables
vip.scores <- data.frame('metabolite' = names(sort(getVipVn(t1d.pls))), 'VIP' = sort(getVipVn(t1d.pls)))
vip.scores$metabolite <- factor(vip.scores$metabolite, levels = vip.scores$metabolite)
ggplot(vip.scores, aes(metabolite, VIP)) + geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust =.5, hjust=1))

sig.contributors <- vip.scores[vip.scores$VIP>=1.5,]

## retain only highly contributors 
lipids.res.reduced <- lipids.res[,names(lipids.res)%in%sig.contributors$metabolite]
pca.lipids.res.red <- prcomp(lipids.res.reduced, scale. = T, center = T)
fviz_screeplot(pca.lipids.res.red)
fviz_pca_ind(pca.lipids.res.red, geom = 'point', col.ind = mtdt$T1D, addEllipses = T)
fviz_pca_ind(pca.lipids.res.red, geom = 'point', col.ind = mtdt$Groups, addEllipses = T)
fviz_pca_var(pca.lipids.res.red)
scores.lipids.res.red <- as.data.frame(pca.lipids.res.red$x)
scores.lipids.res.red$T1D <- mtdt$T1D
p1 <- ggplot(scores.lipids.res.red, aes(PC1, PC2, color = T1D)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_taylor() +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 29.5%)", 
       y = "Principal Component 2 (var.expl. 9.8%)", 
       color = 'Clinical\ngroup')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

p1.allgroups <- ggplot(scores.lipids.res.red, aes(PC1, PC2, color = mtdt$Groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_taylor() +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 29.5%)", 
       y = "Principal Component 2 (var.expl. 9.8%)")

ggExtra::ggMarginal(
  p = p1.allgroups,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## univariate approach 
identical(row.names(mtdt), row.names(lipids.res))
lipids.res$T1D <- mtdt$T1D
lipids.res.ctrl <- lipids.res[lipids.res$T1D == 'Controls',]
lipids.res.t1d <- lipids.res[lipids.res$T1D == 'T1D',]

lipids.res$T1D
eff.size <- sapply(lipids.res[,1:7470], function(x){
  cliff.delta(x, lipids.res$T1D, data = lipids.res)$estimate
})

p.vals <- sapply(lipids.res[,1:7470], function(x){
  t.test(x~T1D, data = lipids.res)$p.val
})

volcano.ctl.t1d <- data.frame('Effect.size' = eff.size, 
                              'p.values' = p.vals, 'fdr' = p.adjust(p.vals, method = 'fdr'))

ggplot(lipids.res, aes(T1D, Cer.d41.1.)) +
  geom_boxplot() + geom_quasirandom()

ggplot(volcano.ctl.t1d, aes(Effect.size, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.05, 'red', 'black'))) +
  geom_hline(yintercept = -log10(.05), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = -log10(.001), linetype = 'dashed', color = 'darkolivegreen2') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('T1D', x = .1, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Controls', x = .9, y = .19, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  scale_color_manual(values = c('black', 'red')) + guides(color = F) +
  # geom_text_repel(aes(label = ifelse(fdr<0.05, as.character(row.names(volcano.ctl.t1d)), ''))) +
  theme_bw() +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

dim(volcano.ctl.t1d[volcano.ctl.t1d$fdr<.05,])

row.names(volcano.ctl.t1d[volcano.ctl.t1d$fdr<.1,])
names(lipids.res.reduced)
venn.plot <- VennDiagram::venn.diagram(x = list('Univariate' = row.names(volcano.ctl.t1d[volcano.ctl.t1d$fdr<.1,]),
                                                'PLS-DA' = names(lipids.res.reduced)), filename = NULL, 
                                       fill = c('red', 'blue'))
grid::grid.newpage()
grid::grid.draw(venn.plot)

lipids.res$T1D.num <- plyr::revalue(lipids.res$T1D, replace = c('Controls' = 0, 'T1D' = 1))
met.aucs <- c()
low.aucs <- c()
high.aucs <- c()
for (i in names(lipids.res)[1:7470]){
  opt.test <- optimal.cutpoints(X = i, status= 'T1D.num', data = lipids.res, methods = c('ValueSe', 'ValueSp', 'Youden', 'MaxDOR', 'ValueNPV', 'ValuePPV'), tag.healthy = 0)
  met.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[1]
  low.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[2]
  high.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[3]
}

names(met.aucs) <- names(lipids.res)[1:7470]
names(met.aucs) == names(p.vals)

met.aucs <- data.frame('Metabolite' = names(met.aucs), 'AUC' = met.aucs, 'Low.AUC' = low.aucs, 'High.AUC' = high.aucs, 'p-value' = p.vals, 'FDR' = p.adjust(p.vals, method = 'fdr'))
met.aucs <- met.aucs[order(met.aucs$AUC),]
met.aucs$Metabolite <- factor(met.aucs$Metabolite, levels = met.aucs$Metabolite)
pdf('lipids_AUC_T1DvsCTL')
ggplot(met.aucs, aes(AUC, Metabolite, color = ifelse(FDR<.1, 'red', 'black'))) + 
  geom_point() +
  scale_color_manual(values = c('black', 'red')) +
  # geom_errorbar(aes(xmin = Low.AUC, xmax = High.AUC)) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 0.66, linetype = 'dashed') +
  geom_vline(xintercept = 0.75, linetype = 'dotted') +
  geom_vline(xintercept = 0.33, linetype = 'dashed') +
  geom_vline(xintercept = 0.25, linetype = 'dotted') +
  xlim(0,1) + guides(color=F) + 
  theme_bw() + coord_flip()

volcano.ctl.t1d <- volcano.ctl.t1d[match(row.names(met.aucs), row.names(volcano.ctl.t1d)),]
row.names(met.aucs) == row.names(volcano.ctl.t1d)
met.aucs %>% 
  add_column('cliff.delta' = volcano.ctl.t1d$Effect.size) %>% 
  write.table('CTLvsT1D_univariatelipidomics.txt', sep = '\t', row.names = F)



##### Normo vs macro ------------------------------------------------------------------------------------
## PCA for normo vs macro
ind.normo.macro <- c(row.names(mtdt[mtdt$Groups=='Normo',]), row.names(mtdt[mtdt$Groups=='Macro',]))
lip.res.nor.mac <- lipids.res[row.names(lipids.res)%in%ind.normo.macro,]
mtdt.nor.mac <- mtdt[row.names(mtdt)%in%ind.normo.macro,]

pca.lipids.res.nor.mac <- prcomp(lip.res.nor.mac[,1:7470], scale. = T, center = T)
fviz_screeplot(pca.lipids.res.nor.mac)
fviz_pca_ind(pca.lipids.res.nor.mac, geom = 'point', col.ind = mtdt.nor.mac$Groups, addEllipses = T)
scores.lipids.res <- as.data.frame(pca.lipids.res.nor.mac$x)
row.names(scores.lipids.res) == row.names(mtdt.nor.mac)
scores.lipids.res$groups <- mtdt.nor.mac$Groups
p1 <- ggplot(scores.lipids.res, aes(PC1, PC2, color = groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 15.4%)", 
       y = "Principal Component 2 (var.expl. 9.4%)",
       color = 'Clinical\ngroups')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## PLS-DA for controls vs T1D 
lip.res.nor.mac$T1D <- NULL
lip.split <- initial_split(lip.res.nor.mac, prop = .7)
lip.train <- training(lip.split)
mylist <-list('total' = row.names(lip.res.nor.mac), 'train' = row.names(lip.train))
common_values <- intersect(row.names(lip.res.nor.mac), row.names(lip.train))
train.index <- lapply(mylist, function(x) which(x %in% common_values))$total

# initial fit
nor.mac.pls <- opls(lip.res.nor.mac, mtdt.nor.mac$Groups, subset = train.index, predI = 5)
summary(nor.mac.pls)

# check predictions on the training subset
trainVi <- getSubsetVi(nor.mac.pls)
table(mtdt.nor.mac$Groups[trainVi], fitted(nor.mac.pls))
table(mtdt.nor.mac$Groups[-trainVi], predict(nor.mac.pls, lip.res.nor.mac[-trainVi,]))

## variables
vip.scores <- data.frame('metabolite' = names(sort(getVipVn(nor.mac.pls))), 'VIP' = sort(getVipVn(nor.mac.pls)))
vip.scores$metabolite <- factor(vip.scores$metabolite, levels = vip.scores$metabolite)
ggplot(vip.scores, aes(metabolite, VIP)) + geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust =.5, hjust=1))

sig.contributors <- vip.scores[vip.scores$VIP>=1.5,]

## retain only highly contributors 
lipids.res.nor.mac.red <- lip.res.nor.mac[,names(lip.res.nor.mac)%in%sig.contributors$metabolite]
pca.lipids.res.nor.mac.red <- prcomp(lipids.res.nor.mac.red, scale. = T, center = T)
fviz_screeplot(pca.lipids.res.nor.mac.red)
fviz_pca_ind(pca.lipids.res.nor.mac.red, geom = 'point', col.ind = mtdt.nor.mac$Groups, addEllipses = T)
fviz_pca_var(pca.lipids.res.nor.mac.red)
scores.lipids.res.nor.mac.red <- as.data.frame(pca.lipids.res.nor.mac.red$x)
scores.lipids.res.nor.mac.red$groups <- mtdt.nor.mac$Groups
p1 <- ggplot(scores.lipids.res.nor.mac.red, aes(PC1, PC2, color = groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 18.3%)", 
       y = "Principal Component 2 (var.expl. 12.7%)", 
       color = 'Clinical\ngroup')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## univariate approach 
identical(row.names(mtdt.nor.mac), row.names(lip.res.nor.mac))
lip.res.nor.mac$groups <- factor(mtdt.nor.mac$Groups)

lip.res.nor.mac$groups
eff.size <- sapply(lip.res.nor.mac[,1:7470], function(x){
  cliff.delta(x, lip.res.nor.mac$groups, data = lip.res.nor.mac)$estimate
})

p.vals <- sapply(lip.res.nor.mac[,1:7470], function(x){
  t.test(x~groups, data = lip.res.nor.mac)$p.val
})

volcano.nor.mac <- data.frame('Effect.size' = eff.size, 
                              'p.values' = p.vals, 'fdr' = p.adjust(p.vals, method = 'fdr'))

ggplot(lipids.res, aes(T1D, Cer.d41.1.)) +
  geom_boxplot() + geom_quasirandom()

ggplot(volcano.nor.mac, aes(Effect.size, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.05, 'red', 'black'))) +
  geom_hline(yintercept = -log10(.05), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = -log10(.001), linetype = 'dashed', color = 'darkolivegreen2') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('T1D', x = .1, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Controls', x = .9, y = .19, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  scale_color_manual(values = c('black', 'red')) + guides(color = F) +
  geom_text_repel(aes(label = ifelse(fdr<0.05, as.character(row.names(volcano.nor.mac)), ''))) +
  theme_bw() +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

####

##### Normo vs micro ------------------------------------------------------------------------------------
## PCA for normo vs micro
ind.normo.micro <- c(row.names(mtdt[mtdt$Groups=='Normo',]), row.names(mtdt[mtdt$Groups=='Micro',]))
lip.res.nor.mic <- lipids.res[row.names(lipids.res)%in%ind.normo.micro,]
mtdt.nor.mic <- mtdt[row.names(mtdt)%in%ind.normo.micro,]

pca.lipids.res.nor.mic <- prcomp(lip.res.nor.mic[,1:7470], scale. = T, center = T)
fviz_screeplot(pca.lipids.res.nor.mic)
fviz_pca_ind(pca.lipids.res.nor.mic, geom = 'point', col.ind = mtdt.nor.mic$Groups, addEllipses = T)
scores.lipids.res <- as.data.frame(pca.lipids.res.nor.mic$x)
row.names(scores.lipids.res) == row.names(mtdt.nor.mic)
scores.lipids.res$groups <- mtdt.nor.mic$Groups
p1 <- ggplot(scores.lipids.res, aes(PC1, PC2, color = groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 15.4%)", 
       y = "Principal Component 2 (var.expl. 9.4%)",
       color = 'Clinical\ngroups')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## PLS-DA for controls vs T1D 
lip.res.nor.mic$T1D <- NULL
lip.split <- initial_split(lip.res.nor.mic, prop = .7)
lip.train <- training(lip.split)
mylist <-list('total' = row.names(lip.res.nor.mic), 'train' = row.names(lip.train))
common_values <- intersect(row.names(lip.res.nor.mic), row.names(lip.train))
train.index <- lapply(mylist, function(x) which(x %in% common_values))$total

# initial fit
nor.mic.pls <- opls(lip.res.nor.mic, mtdt.nor.mic$Groups, subset = train.index, predI = 5)
summary(nor.mic.pls)

# check predictions on the training subset
trainVi <- getSubsetVi(nor.mic.pls)
table(mtdt.nor.mic$Groups[trainVi], fitted(nor.mic.pls))
table(mtdt.nor.mic$Groups[-trainVi], predict(nor.mic.pls, lip.res.nor.mic[-trainVi,]))

## variables
vip.scores <- data.frame('metabolite' = names(sort(getVipVn(nor.mic.pls))), 'VIP' = sort(getVipVn(nor.mic.pls)))
vip.scores$metabolite <- factor(vip.scores$metabolite, levels = vip.scores$metabolite)
ggplot(vip.scores, aes(metabolite, VIP)) + geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust =.5, hjust=1))

sig.contributors <- vip.scores[vip.scores$VIP>=1.5,]

## retain only highly contributors 
lipids.res.nor.mic.red <- lip.res.nor.mic[,names(lip.res.nor.mic)%in%sig.contributors$metabolite]
pca.lipids.res.nor.mic.red <- prcomp(lipids.res.nor.mic.red, scale. = T, center = T)
fviz_screeplot(pca.lipids.res.nor.mic.red)
fviz_pca_ind(pca.lipids.res.nor.mic.red, geom = 'point', col.ind = mtdt.nor.mic$Groups, addEllipses = T)
fviz_pca_var(pca.lipids.res.nor.mic.red)
scores.lipids.res.nor.mic.red <- as.data.frame(pca.lipids.res.nor.mic.red$x)
scores.lipids.res.nor.mic.red$groups <- mtdt.nor.mic$Groups
p1 <- ggplot(scores.lipids.res.nor.mic.red, aes(PC1, PC2, color = groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 15.9%)", 
       y = "Principal Component 2 (var.expl. 9.6%)", 
       color = 'Clinical\ngroup')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## univariate approach 
identical(row.names(mtdt.nor.mic), row.names(lip.res.nor.mic))
lip.res.nor.mic$groups <- factor(mtdt.nor.mic$Groups)

lip.res.nor.mic$groups
eff.size <- sapply(lip.res.nor.mic[,1:7470], function(x){
  cliff.delta(x, lip.res.nor.mic$groups, data = lip.res.nor.mic)$estimate
})

p.vals <- sapply(lip.res.nor.mic[,1:7470], function(x){
  t.test(x~groups, data = lip.res.nor.mic)$p.val
})

volcano.nor.mic <- data.frame('Effect.size' = eff.size, 
                              'p.values' = p.vals, 'fdr' = p.adjust(p.vals, method = 'fdr'))

ggplot(lipids.res, aes(T1D, Cer.d41.1.)) +
  geom_boxplot() + geom_quasirandom()

ggplot(volcano.nor.mic, aes(Effect.size, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.05, 'red', 'black'))) +
  geom_hline(yintercept = -log10(.05), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = -log10(.001), linetype = 'dashed', color = 'darkolivegreen2') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('T1D', x = .1, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Controls', x = .9, y = .19, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  scale_color_manual(values = c('black', 'red')) + guides(color = F) +
  geom_text_repel(aes(label = ifelse(fdr<0.05, as.character(row.names(volcano.nor.mic)), ''))) +
  theme_bw() +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))
hist(volcano.nor.mic$fdr); min(volcano.nor.mic$fdr)
####

##### micro vs macro ------------------------------------------------------------------------------------
## PCA for micro vs macro
ind.micro.macro <- c(row.names(mtdt[mtdt$Groups=='Micro',]), row.names(mtdt[mtdt$Groups=='Macro',]))
lip.res.mic.mac <- lipids.res[row.names(lipids.res)%in%ind.micro.macro,]
mtdt.mic.mac <- mtdt[row.names(mtdt)%in%ind.micro.macro,]

pca.lipids.res.mic.mac <- prcomp(lip.res.mic.mac[,1:7470], scale. = T, center = T)
fviz_screeplot(pca.lipids.res.mic.mac)
fviz_pca_ind(pca.lipids.res.mic.mac, geom = 'point', col.ind = mtdt.mic.mac$Groups, addEllipses = T)
scores.lipids.res <- as.data.frame(pca.lipids.res.mic.mac$x)
row.names(scores.lipids.res) == row.names(mtdt.mic.mac)
scores.lipids.res$groups <- mtdt.mic.mac$Groups
p1 <- ggplot(scores.lipids.res, aes(PC1, PC2, color = groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 15.5%)", 
       y = "Principal Component 2 (var.expl. 9.2%)",
       color = 'Clinical\ngroups')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## PLS-DA for controls vs T1D 
lip.res.mic.mac$T1D <- NULL
lip.split <- initial_split(lip.res.mic.mac, prop = .7)
lip.train <- training(lip.split)
mylist <-list('total' = row.names(lip.res.mic.mac), 'train' = row.names(lip.train))
common_values <- intersect(row.names(lip.res.mic.mac), row.names(lip.train))
train.index <- lapply(mylist, function(x) which(x %in% common_values))$total

# initial fit
mic.mac.pls <- opls(lip.res.mic.mac, mtdt.mic.mac$Groups, subset = train.index, predI = 5)
summary(mic.mac.pls)

# check predictions on the training subset
trainVi <- getSubsetVi(mic.mac.pls)
table(mtdt.mic.mac$Groups[trainVi], fitted(mic.mac.pls))
table(mtdt.mic.mac$Groups[-trainVi], predict(mic.mac.pls, lip.res.mic.mac[-trainVi,]))

## variables
vip.scores <- data.frame('metabolite' = names(sort(getVipVn(mic.mac.pls))), 'VIP' = sort(getVipVn(mic.mac.pls)))
vip.scores$metabolite <- factor(vip.scores$metabolite, levels = vip.scores$metabolite)
ggplot(vip.scores, aes(metabolite, VIP)) + geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust =.5, hjust=1))

sig.contributors <- vip.scores[vip.scores$VIP>=1.5,]

## retain only highly contributors 
lipids.res.mic.mac.red <- lip.res.mic.mac[,names(lip.res.mic.mac)%in%sig.contributors$metabolite]
pca.lipids.res.mic.mac.red <- prcomp(lipids.res.mic.mac.red, scale. = T, center = T)
fviz_screeplot(pca.lipids.res.mic.mac.red)
fviz_pca_ind(pca.lipids.res.mic.mac.red, geom = 'point', col.ind = mtdt.mic.mac$Groups, addEllipses = T)
fviz_pca_var(pca.lipids.res.mic.mac.red)
scores.lipids.res.mic.mac.red <- as.data.frame(pca.lipids.res.mic.mac.red$x)
scores.lipids.res.mic.mac.red$groups <- mtdt.mic.mac$Groups
p1 <- ggplot(scores.lipids.res.mic.mac.red, aes(PC1, PC2, color = groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 36.4%)", 
       y = "Principal Component 2 (var.expl. 14.7%)", 
       color = 'Clinical\ngroup')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## univariate approach 
identical(row.names(mtdt.mic.mac), row.names(lip.res.mic.mac))
lip.res.mic.mac$groups <- factor(mtdt.mic.mac$Groups)

lip.res.mic.mac$groups
eff.size <- sapply(lip.res.mic.mac[,1:7470], function(x){
  cliff.delta(x, lip.res.mic.mac$groups, data = lip.res.mic.mac)$estimate
})

p.vals <- sapply(lip.res.mic.mac[,1:7470], function(x){
  t.test(x~groups, data = lip.res.mic.mac)$p.val
})

volcano.mic.mac <- data.frame('Effect.size' = eff.size, 
                              'p.values' = p.vals, 'fdr' = p.adjust(p.vals, method = 'fdr'))

ggplot(lip.res.mic.mac, aes(groups, Cer.d41.1.)) +
  geom_boxplot() + geom_quasirandom()

ggplot(volcano.mic.mac, aes(Effect.size, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.05, 'red', 'black'))) +
  geom_hline(yintercept = -log10(.05), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = -log10(.001), linetype = 'dashed', color = 'darkolivegreen2') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('T1D', x = .1, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Controls', x = .9, y = .19, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  scale_color_manual(values = c('black', 'red')) + guides(color = F) +
  geom_text_repel(aes(label = ifelse(fdr<0.05, as.character(row.names(volcano.nor.mac)), ''))) +
  theme_bw() +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))
hist(volcano.mic.mac$fdr); min(volcano.mic.mac$fdr)
####
