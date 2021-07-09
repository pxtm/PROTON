#############################################################################################################
################################# PROTON metabolites Controls vs T1D ########################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 08.06.2021
## last edit: 24.06.2021

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

## functions
source('lm_function.R')

## data
mtdt <- readRDS('mtdt_common.rds')
gc.res <- readRDS('residuals_GCGC_egfr.rds')

row.names(mtdt) == row.names(gc.res)

## GCGC - controls vs T1D ------------------------------------------------------------------------------------------------------
# normality by Shapiro-Wilk
# hist(sapply(gc.res, function(a)(shapiro.test(a)$p.val)))

qqplots <- list()
for (i in names(gc.res)){
  qqplots[[i]] <- ggqqplot(gc.res[,i], title = i)
}

multi.page <- ggarrange(plotlist = qqplots, nrow=4, ncol=4)
ggexport(multi.page, filename = 'Metabolomics/qqplots_GCGC.pdf')

## PCA for controls vs T1D ----
mtdt$T1D <- as.factor(plyr::revalue(mtdt$Groups, c('Normo' = 'T1D', 'Micro' = 'T1D', 'Macro' = 'T1D')))
pca.gc.res <- prcomp(gc.res, scale. = T, center = T)
fviz_screeplot(pca.gc.res)
fviz_pca_ind(pca.gc.res, geom = 'point', col.ind = mtdt$T1D, addEllipses = T)
scores.gc.res <- as.data.frame(pca.gc.res$x)
row.names(scores.gc.res) == row.names(mtdt)
scores.gc.res$T1D <- mtdt$T1D
p1 <- ggplot(scores.gc.res, aes(PC1, PC2, color = T1D)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_taylor() +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 13.1%)", 
       y = "Principal Component 2 (var.expl. 4.6%)")

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## PLS-DA for controls vs T1D ----
gc.split <- initial_split(gc.res, prop = .7)
gc.train <- training(gc.split)
mylist <-list('total' = row.names(gc.res), 'train' = row.names(gc.train))
common_values <- intersect(row.names(gc.res), row.names(gc.train))
train.index <- lapply(mylist, function(x) which(x %in% common_values))$total

# initial fit
t1d.pls <- opls(gc.res, mtdt$T1D, predI = 5, subset = train.index)
summary(t1d.pls)

# check predictions on the training subset
trainVi <- getSubsetVi(t1d.pls)
table(mtdt$T1D[trainVi], fitted(t1d.pls))
table(mtdt$T1D[-trainVi], predict(t1d.pls, gc.res[-trainVi,]))

## variables
vip.scores <- data.frame('metabolite' = names(sort(getVipVn(t1d.pls))), 'VIP' = sort(getVipVn(t1d.pls)))
vip.scores$metabolite <- factor(vip.scores$metabolite, levels = vip.scores$metabolite)
ggplot(vip.scores, aes(metabolite, VIP)) + geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust =.5, hjust=1))

sig.contributors <- vip.scores[vip.scores$VIP>=1,]

## retain only highly contributors 
gc.res.reduced <- gc.res[,names(gc.res)%in%sig.contributors$metabolite]
pca.gc.res.red <- prcomp(gc.res.reduced, scale. = T, center = T)
fviz_screeplot(pca.gc.res.red)
fviz_pca_ind(pca.gc.res.red, geom = 'point', col.ind = mtdt$T1D, addEllipses = T)
fviz_pca_ind(pca.gc.res.red, geom = 'point', col.ind = mtdt$Groups, addEllipses = T)
fviz_pca_var(pca.gc.res.red)
scores.gc.res.red <- as.data.frame(pca.gc.res.red$x)
scores.gc.res.red$T1D <- mtdt$T1D
p1 <- ggplot(scores.gc.res.red, aes(PC1, PC2, color = T1D)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_taylor() +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 11.2%)", 
       y = "Principal Component 2 (var.expl. 8.7%)")

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

p1.allgroups <- ggplot(scores.gc.res.red, aes(PC1, PC2, color = mtdt$Groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_taylor() +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 11.2%)", 
       y = "Principal Component 2 (var.expl. 8.7%)")

ggExtra::ggMarginal(
  p = p1.allgroups,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## univariate approach ----
identical(row.names(mtdt), row.names(gc.res))
gc.res$T1D <- mtdt$T1D
gc.res.ctrl <- gc.res[gc.res$T1D == 'Controls',]
gc.res.t1d <- gc.res[gc.res$T1D == 'T1D',]

gc.res$T1D
eff.size <- sapply(gc.res[,1:398], function(x){
  cliff.delta(x, gc.res$T1D, data = gc.res)$estimate
})

p.vals <- sapply(gc.res[,1:398], function(x){
  t.test(x~T1D, data = gc.res)$p.val
})

volcano.ctl.t1d <- data.frame('Effect.size' = eff.size, 
                              'p.values' = p.vals, 'fdr' = p.adjust(p.vals, method = 'fdr'))

ggplot(gc.res, aes(T1D, X1.5.Anhydrosorbitol)) +
    geom_boxplot() + geom_quasirandom()

ggplot(volcano.ctl.t1d, aes(Effect.size, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.05, 'red', 'black'))) +
  geom_hline(yintercept = -log10(.05), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = -log10(.001), linetype = 'dashed', color = 'darkolivegreen2') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('T1D', x = .1, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Controls', x = .9, y = .18, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  scale_color_manual(values = c('black', 'red')) + guides(color = F) +
  geom_text_repel(aes(label = ifelse(fdr<0.05, as.character(row.names(volcano.ctl.t1d)), ''))) +
  theme_bw() +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

dim(volcano.ctl.t1d[volcano.ctl.t1d$fdr<.1,])

row.names(volcano.ctl.t1d[volcano.ctl.t1d$fdr<.1,])
names(gc.res.reduced)
venn.plot <- VennDiagram::venn.diagram(x = list('Univariate' = row.names(volcano.ctl.t1d[volcano.ctl.t1d$fdr<.1,]),
                                   'PLS-DA' = names(gc.res.reduced)), filename = NULL, 
                                   fill = c('red', 'blue'))
grid::grid.newpage()
grid::grid.draw(venn.plot)

gc.res$T1D.num <- plyr::revalue(gc.res$T1D, replace = c('Controls' = 0, 'T1D' = 1))
met.aucs <- c()
low.aucs <- c()
high.aucs <- c()
for (i in names(gc.res)[1:398]){
  opt.test <- optimal.cutpoints(X = i, status= 'T1D.num', data = gc.res, methods = c('ValueSe', 'ValueSp', 'Youden', 'MaxDOR', 'ValueNPV', 'ValuePPV'), tag.healthy = 0)
  met.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[1]
  low.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[2]
  high.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[3]
}

names(met.aucs) <- names(gc.res)[1:398]
names(met.aucs) == names(p.vals)

met.aucs <- data.frame('Metabolite' = names(met.aucs), 'AUC' = met.aucs, 'Low.AUC' = low.aucs, 'High.AUC' = high.aucs, 'p-value' = p.vals, 'FDR' = p.adjust(p.vals, method = 'fdr'))
met.aucs <- met.aucs[order(met.aucs$AUC),]
met.aucs$Metabolite <- factor(met.aucs$Metabolite, levels = met.aucs$Metabolite)
ggplot(met.aucs, aes(AUC, Metabolite, color = ifelse(FDR<.1, 'red', 'black'))) + 
  geom_point() +
  scale_color_manual(values = c('black', 'red')) +
  # geom_errorbar(aes(xmin = Low.AUC, xmax = High.AUC)) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 0.66, linetype = 'dashed') +
  geom_vline(xintercept = 0.75, linetype = 'dotted') +
  geom_vline(xintercept = 0.33, linetype = 'dashed') +
  geom_vline(xintercept = 0.25, linetype = 'dotted') +
  xlim(0,1) +
  theme_bw() 

volcano.ctl.t1d <- volcano.ctl.t1d[match(row.names(met.aucs), row.names(volcano.ctl.t1d)),]
row.names(met.aucs) == row.names(volcano.ctl.t1d)
met.aucs %>% 
  add_column('cliff.delta' = volcano.ctl.t1d$Effect.size) %>% 
  write.table('CTLvsT1D_univariatePolarMetabolites.txt', sep = '\t', row.names = F)

## GCGC - Normo vs Macro ------------------------------------------------------------------------------------------------------
gc.res$T1D <- NULL
ind.normo.macro <- c(row.names(mtdt[mtdt$Groups=='Normo',]), row.names(mtdt[mtdt$Groups=='Macro',]))
gc.res.nor.mac <- gc.res[row.names(gc.res)%in%ind.normo.macro,]
mtdt.nor.mac <- mtdt[row.names(mtdt)%in%ind.normo.macro,]

pca.gc.res.nor.mac <- prcomp(gc.res.nor.mac[,1:398], scale. = T, center = T)
fviz_screeplot(pca.gc.res.nor.mac)
fviz_pca_ind(pca.gc.res.nor.mac, geom = 'point', col.ind = mtdt.nor.mac$Groups, addEllipses = T)

## OPLS-DA for controls vs T1D ----
gc.split <- initial_split(gc.res.nor.mac, prop = .7)
gc.train <- training(gc.split)
mylist <-list('total' = row.names(gc.res.nor.mac), 'train' = row.names(gc.train))
common_values <- intersect(row.names(gc.res.nor.mac), row.names(gc.train))
train.index <- lapply(mylist, function(x) which(x %in% common_values))$total

# initial fit
nor.mac.pls <- opls(gc.res.nor.mac[,1:398], mtdt.nor.mac$Groups, predI = 5, subset = train.index)
summary(nor.mac.pls)

# check predictions on the training subset
trainVi <- getSubsetVi(nor.mac.pls)
table(mtdt.nor.mac$Groups[trainVi], fitted(nor.mac.pls))
table(mtdt.nor.mac$Groups[-trainVi], predict(nor.mac.pls, gc.res.nor.mac[,1:398][-trainVi,]))

## variables
vip.scores <- data.frame('metabolite' = names(sort(getVipVn(nor.mac.pls))), 'VIP' = sort(getVipVn(nor.mac.pls)))
vip.scores$metabolite <- factor(vip.scores$metabolite, levels = vip.scores$metabolite)
ggplot(vip.scores, aes(metabolite, VIP)) + geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust =.5, hjust=1))
sig.contributors <- vip.scores[vip.scores$VIP>=1,]

## retain only highly contributors 
gc.res.nor.mac.reduced <- gc.res.nor.mac[,names(gc.res.nor.mac)%in%sig.contributors$metabolite]
pca.gc.res.nor.mac.red <- prcomp(gc.res.nor.mac.reduced, scale. = T, center = T)
fviz_screeplot(pca.gc.res.nor.mac.red)
fviz_pca_ind(pca.gc.res.nor.mac.red, geom = 'point', col.ind = mtdt.nor.mac$Groups, addEllipses = T)
fviz_pca_var(pca.gc.res.nor.mac.red)

## univariate approach ----
identical(row.names(mtdt.nor.mac), row.names(gc.res.nor.mac))
gc.res.nor.mac$alb <- mtdt.nor.mac$Groups
gc.res.nor <- gc.res.nor.mac[gc.res.nor.mac$alb == 'Normo',]
gc.res.mac <- gc.res.nor.mac[gc.res.nor.mac$alb == 'Macro',]

gc.res.nor.mac$alb
eff.size <- sapply(gc.res.nor.mac[,1:398], function(x){
  cliff.delta(x, gc.res.nor.mac$alb, data = gc.res.nor.mac)$estimate
})

p.vals <- sapply(gc.res.nor.mac[,1:398], function(x){
  t.test(x~alb, data = gc.res.nor.mac)$p.val
})

volcano.nor.mac <- data.frame('Effect.size' = eff.size, 
                              'p.values' = p.vals, 'fdr' = p.adjust(p.vals, method = 'fdr'))

ggplot(gc.res, aes(T1D, X1.5.Anhydrosorbitol)) +
  geom_boxplot() + geom_quasirandom()

ggplot(volcano.nor.mac, aes(Effect.size, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<=.1, 'red', 'black'))) +
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_hline(yintercept = -log10(.001), linetype = 'dashed', color = 'darkolivegreen2') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('Normo', x = .1, y = .5, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Macro', x = .9, y = .5, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  scale_color_manual(values = c('black', 'red')) + guides(color = F) +
  geom_text_repel(aes(label = ifelse(fdr<0.1, as.character(row.names(volcano.ctl.t1d)), ''))) +
  theme_bw() +
  labs(x = '-log10(FDR)', y = "Cliff's Delta effect size")

dim(volcano.ctl.t1d[volcano.ctl.t1d$fdr<.05,])
volcano.ctl.t1d %>% 
  filter(fdr <=.1) %>% 
  arrange(Effect.size)

## GCC: macro vs micro-albuminuria ------------------------------------------------------------------------------------------------------------------------------------------
individuals <- row.names(mtdt[mtdt$Groups=='Micro'|mtdt$Groups=='Macro',])
gc.macmic <- gc.res[row.names(gc.res)%in%individuals,]
mtdt.macmic <- mtdt[row.names(mtdt)%in%individuals,]
identical(row.names(gc.macmic), row.names(mtdt.macmic))

## PCA
macmic.pca <- prcomp(gc.macmic[, 1:398], scale. = T, center = T)
fviz_screeplot(macmic.pca)
fviz_pca_biplot(macmic.pca)
macmic.scores <- as.data.frame(macmic.pca$x)
p2 <- ggplot(macmic.scores, aes(PC1, PC2, color = mtdt.macmic$Groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 13.9%)", 
       y = "Principal Component 2 (var.expl. 4.9%)")

ggExtra::ggMarginal(
  p = p2,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## PLS-DA
gc.split <- initial_split(gc.macmic, prop = .7)
gc.train <- training(gc.split)
mylist <-list('total' = row.names(gc.macmic), 'train' = row.names(gc.train))
common_values <- intersect(row.names(gc.macmic), row.names(gc.train))
train.index <- lapply(mylist, function(x) which(x %in% common_values))$total

# initial fit
macmic.pls <- opls(gc.macmic[,1:398], mtdt.macmic$Groups, subset = train.index, predI = 3) ## forcing the machine....
summary(macmic.pls)

# check predictions on the training subset
trainVi <- getSubsetVi(macmic.pls)
table(mtdt.macmic$Groups[trainVi], fitted(macmic.pls))
table(mtdt.macmic$Groups[-trainVi], predict(macmic.pls, gc.macmic[-trainVi,1:398]))

## variables
vip.scores <- data.frame('metabolite' = names(sort(getVipVn(macmic.pls))), 'VIP' = sort(getVipVn(macmic.pls)))
vip.scores$metabolite <- factor(vip.scores$metabolite, levels = vip.scores$metabolite)
ggplot(vip.scores, aes(metabolite, VIP)) + geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust =.5, hjust=1))

sig.contributors <- vip.scores[vip.scores$VIP>=1,]

## retain only highly contributors 
gc.macmic.reduced <- gc.macmic[,names(gc.macmic)%in%sig.contributors$metabolite]
pca.gc.macmic.red <- prcomp(gc.macmic.reduced, scale. = T, center = T)
fviz_screeplot(pca.gc.macmic.red)
fviz_pca_ind(pca.gc.macmic.red, geom = 'point', col.ind = mtdt.macmic$Groups, addEllipses = T)
fviz_pca_var(pca.gc.macmic.red)
scores.gc.macmic.red <- as.data.frame(pca.gc.macmic.red$x)
scores.gc.macmic.red$Groups <- mtdt.macmic$Groups
p3 <- ggplot(scores.gc.macmic.red, aes(PC1, PC2, color = Groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_manual(values = c('#43475b', '#a88f92')) +
  theme_bw() +
  labs(x = "Principal Component 1 (var. expl. 16.3%)", 
       y = "Principal Component 2 (var.expl. 6%)")

ggExtra::ggMarginal(
  p = p3,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

## univariate approach
identical(row.names(mtdt.macmic), row.names(gc.macmic))
gc.macmic$Groups <- mtdt.macmic$Groups

eff.size <- sapply(gc.macmic[,1:398], function(x){
  cliff.delta(x, gc.macmic$Groups, data = gc.macmic)$estimate
})

p.vals <- sapply(gc.macmic[,1:398], function(x){
  t.test(x~Groups, data = gc.macmic)$p.val
})

volcano.macmic <- data.frame('Effect.size' = eff.size, 
                              'p.values' = p.vals, 'fdr' = p.adjust(p.vals, method = 'fdr'))

ggplot(gc.macmic, aes(Groups, Ribitol_2)) +
  geom_boxplot() + geom_quasirandom()

ggplot(volcano.macmic, aes(Effect.size, -log10(p.vals))) +
  geom_point(aes(color = ifelse(p.vals<=.1, 'red', 'black'))) +
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  # geom_hline(yintercept = -log10(.001), linetype = 'dashed', color = 'darkolivegreen2') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('Micro', x = .1, y = .14, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Macro', x = .9, y = .15, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  scale_color_manual(values = c('black', 'red')) + guides(color = F) +
  geom_text_repel(aes(label = ifelse(p.vals<0.1, as.character(row.names(volcano.macmic)), ''))) +
  theme_bw() +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

volcano.macmic %>% 
  filter(p.vals <= .05) %>% 
  arrange(Effect.size)


gc.macmic$Groups.num <- plyr::revalue(gc.macmic$Groups, replace = c('Macro' = 0, 'Micro' = 1))
met.aucs <- c()
low.aucs <- c()
high.aucs <- c()
for (i in names(gc.macmic)[1:398]){
  opt.test <- optimal.cutpoints(X = i, status= 'Groups.num', data = gc.macmic, methods = c('ValueSe', 'ValueSp', 'Youden', 'MaxDOR', 'ValueNPV', 'ValuePPV'), tag.healthy = 0)
  met.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[1]
  low.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[2]
  high.aucs[i] <- opt.test$ValueSe$Global$measures.acc$AUC[3]
}

names(met.aucs) <- names(gc.res)[1:398]
names(met.aucs) == names(p.vals)

met.aucs <- data.frame('Metabolite' = names(met.aucs), 'AUC' = met.aucs, 'Low.AUC' = low.aucs, 'High.AUC' = high.aucs, 'p-value' = p.vals, 'FDR' = p.adjust(p.vals, method = 'fdr'))
met.aucs <- met.aucs[order(met.aucs$AUC),]
met.aucs$Metabolite <- factor(met.aucs$Metabolite, levels = met.aucs$Metabolite)
ggplot(met.aucs, aes(AUC, Metabolite, color = ifelse(p.value<.1, 'red', 'black'))) + 
  geom_point() +
  scale_color_manual(values = c('black', 'red')) +
  # geom_errorbar(aes(xmin = Low.AUC, xmax = High.AUC)) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 0.66, linetype = 'dashed') +
  geom_vline(xintercept = 0.75, linetype = 'dotted') +
  geom_vline(xintercept = 0.33, linetype = 'dashed') +
  geom_vline(xintercept = 0.25, linetype = 'dotted') +
  xlim(0,1) +
  theme_bw() 

volcano.macmic <- volcano.macmic[match(row.names(met.aucs), row.names(volcano.macmic)),]
row.names(met.aucs) == row.names(volcano.macmic)
met.aucs %>% 
  add_column('cliff.delta' = volcano.macmic$Effect.size) %>% 
  write.table('MicrovsMacro_univariatePolarMetabolites.txt', sep = '\t', row.names = F)

## can we get a predictive model?
# linear model
library(glmnet); library(caret); library(MASS); library(caTools); library(ROCR)
gc.macmic.model <- gc.macmic[,1:398]
gc.macmic.model$Group <- gc.macmic$Groups.num

set.seed(101) 
sample <- sample.split(row.names(gc.macmic.model), SplitRatio = .7)
train <- subset(gc.macmic.model, sample == TRUE)
test  <- subset(gc.macmic.model, sample == FALSE)

fullmodel <- lm(Group ~ ., data = train)
formula.full <- paste0('Group ~', paste(names(fullmodel$coefficients[!is.na(fullmodel$coefficients)])[2:76], collapse = '+'))
fullmodel <- lm(formula.full, data = train)
emptymodel <- lm(Group ~ 1, data = train)

model.fw <- step(emptymodel, direction = 'forward', scope = formula(fullmodel))
model.both <- step(emptymodel, direction = 'both', scope = formula(fullmodel))

p <- predict(model.both, newdata = test, type='response')
pr <- prediction(p, test$Group)
prf <- performance(pr, measure='tpr', x.measure = 'fpr')
auc.perf <- performance(pr, measure = 'auc')
auc <- unlist(slot(auc.perf,"y.values"))
aucprint <- paste(c("AUC ="),round(auc,digits=3),sep="")
aucprint
plot(prf,col="darkgreen")
abline(a=0,b=1,col="red")

# random forest
library(randomForest)
train$Group <- as.factor(train$Group)
rf.class <- randomForest(Group~., data = train, importance = T)
