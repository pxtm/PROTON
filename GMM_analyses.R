#############################################################################################################
#################################### PROTON metaG GMMs analyses #############################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 07.06.2021
## last edit: 

## libraries
library(tidyverse)
library(vegan)
library(ape)
library(reshape2)
library(ggrepel)
library(tayloRswift)
library(grid)
library(corrplot)

## functions
source('lm_function.R')

## load data
mtdt.red <- read.table('mtdt_selected.txt', header=T, row.names=1, sep='\t')
gmm <- as.data.frame(readRDS('GMMcounts.rds'))
# gmm <- as.data.frame(readRDS('GMMabd.rds'))
gmm.annotation <- read.table('GMMs.v1.07.names', sep='\t')

mtdt.red$sample <- paste0('s', row.names(mtdt.red))
row.names(gmm) == mtdt.red$sample
gmm <- gmm[match(mtdt.red$sample, row.names(gmm)),]

## diversity
gmm <- gmm[,sort(names(gmm))]
gmm.annotation <- gmm.annotation[gmm.annotation$V1%in%names(gmm),]
gmm.annotation$V1 == names(gmm)
names(gmm) <- gmm.annotation$V2

gmm.dist <- vegdist(gmm, method = 'bray')
gmm.pcoa <- pcoa(gmm.dist)
gmm.pcoa$values
gmm.pcoa.scores <- as.data.frame(gmm.pcoa$vectors)

p1 <- ggplot(gmm.pcoa.scores, aes(Axis.1, Axis.2, color = mtdt.red$Groups)) +
  geom_point() + stat_ellipse() +
  scale_color_taylor(palette = 'taylor1989') +
  theme_bw() + 
  labs(x = 'Principal coordinate 1 (expl. var. 34.50%)', y = 'Principal coordinate 2 (expl. var. 21.16%)',
       color = 'Clinical\ngroup')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

row.names(gmm) == mtdt.red$sample
gmm.mean <- gmm %>% 
  add_column('group' = mtdt.red$Groups) %>% 
  group_by(group) %>% 
  summarise(across(where(is.numeric), mean))

gmm.mean <- as.data.frame(t(gmm.mean))  
names(gmm.mean) <- gmm.mean[1,]
gmm.mean <- gmm.mean[-1,]
gmm.mean.m <- melt(t(gmm.mean))
gmm.mean.m$Var1 <- factor(gmm.mean.m$Var1, levels = c('Controls', 'Normo', 'Micro', 'Macro'))
gmm.mean.m$value <- as.numeric(gmm.mean.m$value)

ggplot(gmm.mean.m, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = 'white', high = 'darkred')

## associations to bioclinical -----
## Control vs T1D
mtdt.red$biclass <- plyr::revalue(mtdt.red$Groups, c('Macro' = 'T1D', 'Micro' = 'T1D', 'Normo' = 'T1D'))
gmm.biclass <- lm.associations(mtdt.red, gmm, variables = 'biclass', control_by = c('age', 'sex', 'race', 'bmi'))$biclass
gmm.biclass$delta <- as.numeric(gmm.biclass$delta)
volcano.ctl.t1d <- ggplot(gmm.biclass, aes(delta*-1, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.1, 'red', ''))) +
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('Controls', x = .1, y = .18, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('T1D', x = .9, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  geom_text_repel(aes(label = ifelse(fdr<.1, mgs, ''))) + theme_bw() +
  scale_color_manual(values = c('black', 'firebrick4')) + guides(color = F) +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

loadings <- gmm.biclass %>% 
  filter(fdr<=.1) %>% 
  arrange(delta) %>% 
  mutate('mgs' = factor(mgs, levels = mgs)) %>% 
  mutate('fill' = ifelse(delta>0, 'Increased', 'Decreased')) %>% 
  ggplot(aes(delta*-1, mgs, fill = fill)) + 
    geom_col() +
    geom_vline(xintercept = 0) +
    scale_fill_manual(values = c('firebrick4', 'darkolivegreen2')) +
    theme_bw() + xlim(-0.5, 0.5) +
    labs(x = "Cliff's Delta effect size", y = 'Gut Metabolic Module (GMM)', fill = 'Effect size\ndirection') 

CTLvst1d.plot <- cowplot::plot_grid(volcano.ctl.t1d, loadings, ncol = 2)#, rel_widths = c(1, .8))
CTLvst1d.plot <- cowplot::add_sub(CTLvst1d.plot, "For the volcano plot, Cliff's Delta has been inverted for visualization purposes. Points are labelled if FDR < 10%.\nInversion of the Cliff's Delta values is maintained for interpretation coherence in the individual MGS effect sizes barplots",
                                    fontface = 'italic', size = 10)
cowplot::ggdraw(CTLvst1d.plot)

## normo vs micro-----
row.names(gmm) <- gsub('s', '', row.names(gmm))
ind.normo.micro <- row.names(mtdt.red[mtdt.red$Groups=='Normo'|mtdt.red$Groups=='Micro',])
mtdt.nor.mic <- mtdt.red[row.names(mtdt.red)%in%ind.normo.micro, ]
mtdt.nor.mic$Groups <- factor(mtdt.nor.mic$Groups)
gmm.nor.mic <- gmm[row.names(gmm)%in%ind.normo.micro,]

gmm.nor.mic <- lm.associations(mtdt.nor.mic, gmm.nor.mic, variables = 'Groups', control_by = c('age', 'sex', 'race', 'bmi'))$Groups
gmm.nor.mic$delta <- as.numeric(gmm.nor.mic$delta)
volcano.nor.mic <- ggplot(gmm.nor.mic, aes(delta*-1, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.1, 'red', ''))) +
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('Normo', x = .1, y = .18, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Micro', x = .9, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  geom_text_repel(aes(label = ifelse(fdr<.1, mgs, ''))) + theme_bw() +
  scale_color_manual(values = c('black', 'firebrick4')) + guides(color = F) +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

loadings <- gmm.nor.mic %>% 
  filter(fdr<=.1) %>% 
  arrange(delta) %>% 
  mutate('mgs' = factor(mgs, levels = mgs)) %>% 
  mutate('fill' = ifelse(delta>0, 'Increased', 'Decreased')) %>% 
  ggplot(aes(delta*-1, mgs, fill = fill)) + 
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c('firebrick4', 'darkolivegreen2')) +
  theme_bw() + xlim(-0.5, 0.5) +
  labs(x = "Cliff's Delta effect size", y = 'Gut Metabolic Module (GMM)', fill = 'Effect size\ndirection') 

normovsmicro.plot <- cowplot::plot_grid(volcano.nor.mic, loadings, ncol = 2)#, rel_widths = c(1, .8))
normovsmicro.plot <- cowplot::add_sub(normovsmicro.plot, "For the volcano plot, Cliff's Delta has been inverted for visualization purposes. Points are labelled if FDR < 10%.\nInversion of the Cliff's Delta values is maintained for interpretation coherence in the individual MGS effect sizes barplots",
                                  fontface = 'italic', size = 10)
cowplot::ggdraw(normovsmicro.plot)

## normo vs macro ----
ind.normo.macro <- row.names(mtdt.red[mtdt.red$Groups=='Normo'|mtdt.red$Groups=='Macro',])
mtdt.nor.mac <- mtdt.red[row.names(mtdt.red)%in%ind.normo.macro, ]
mtdt.nor.mac$Groups <- factor(mtdt.nor.mac$Groups)
gmm.nor.mac <- gmm[row.names(gmm)%in%ind.normo.macro,]

gmm.nor.mac <- lm.associations(mtdt.nor.mac, gmm.nor.mac, variables = 'Groups', control_by = c('age', 'sex', 'race', 'bmi'))$Groups
gmm.nor.mac$delta <- as.numeric(gmm.nor.mac$delta)
volcano.nor.mac <- ggplot(gmm.nor.mac, aes(delta*-1, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.1, 'red', ''))) +
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('Normo', x = .1, y = .18, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('macro', x = .9, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  geom_text_repel(aes(label = ifelse(fdr<.1, mgs, ''))) + theme_bw() +
  scale_color_manual(values = c('black', 'firebrick4')) + guides(color = F) +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

loadings <- gmm.nor.mac %>% 
  filter(fdr<=.1) %>% 
  arrange(delta) %>% 
  mutate('mgs' = factor(mgs, levels = mgs)) %>% 
  mutate('fill' = ifelse(delta>0, 'Increased', 'Decreased')) %>% 
  ggplot(aes(delta*-1, mgs, fill = fill)) + 
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c('firebrick4', 'darkolivegreen2')) +
  theme_bw() + xlim(-0.5, 0.5) +
  labs(x = "Cliff's Delta effect size", y = 'Gut Metabolic Module (GMM)', fill = 'Effect size\ndirection') 

normovsmacro.plot <- cowplot::plot_grid(volcano.nor.mac, loadings, ncol = 2)#, rel_widths = c(1, .8))
normovsmacro.plot <- cowplot::add_sub(normovsmacro.plot, "For the volcano plot, Cliff's Delta has been inverted for visualization purposes. Points are labelled if FDR < 10%.\nInversion of the Cliff's Delta values is maintained for interpretation coherence in the individual MGS effect sizes barplots",
                                      fontface = 'italic', size = 10)
cowplot::ggdraw(normovsmacro.plot)

## micro vs macro ----
ind.micro.macro <- row.names(mtdt.red[mtdt.red$Groups=='Normo'|mtdt.red$Groups=='Macro',])
mtdt.mic.mac <- mtdt.red[row.names(mtdt.red)%in%ind.micro.macro, ]
mtdt.mic.mac$Groups <- factor(mtdt.mic.mac$Groups)
gmm.mic.mac <- gmm[row.names(gmm)%in%ind.micro.macro,]

gmm.mic.mac <- lm.associations(mtdt.mic.mac, gmm.mic.mac, variables = 'Groups', control_by = c('age', 'sex', 'race', 'bmi'))$Groups
gmm.mic.mac$delta <- as.numeric(gmm.mic.mac$delta)
volcano.mic.mac <- ggplot(gmm.mic.mac, aes(delta*-1, -log10(fdr))) +
  geom_point(aes(color = ifelse(fdr<.1, 'red', ''))) +
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  annotation_custom(grobTree(textGrob('micro', x = .1, y = .18, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('macro', x = .9, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  geom_text_repel(aes(label = ifelse(fdr<.1, mgs, ''))) + theme_bw() +
  scale_color_manual(values = c('black', 'firebrick4')) + guides(color = F) +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"))

loadings <- gmm.mic.mac %>% 
  filter(fdr<=.1) %>% 
  arrange(delta) %>% 
  mutate('mgs' = factor(mgs, levels = mgs)) %>% 
  mutate('fill' = ifelse(delta>0, 'Increased', 'Decreased')) %>% 
  ggplot(aes(delta*-1, mgs, fill = fill)) + 
  geom_col() +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c('firebrick4', 'darkolivegreen2')) +
  theme_bw() + xlim(-0.5, 0.5) +
  labs(x = "Cliff's Delta effect size", y = 'Gut Metabolic Module (GMM)', fill = 'Effect size\ndirection') 

microvsmacro.plot <- cowplot::plot_grid(volcano.mic.mac, loadings, ncol = 2)#, rel_widths = c(1, .8))
microvsmacro.plot <- cowplot::add_sub(microvsmacro.plot, "For the volcano plot, Cliff's Delta has been inverted for visualization purposes. Points are labelled if FDR < 10%.\nInversion of the Cliff's Delta values is maintained for interpretation coherence in the individual MGS effect sizes barplots",
                                      fontface = 'italic', size = 10)
cowplot::ggdraw(microvsmacro.plot)

## correlations with genera? ----
genus.log <- readRDS('genusLOG.rds')
tax.gen <- read.table('tax.txt', header=T, sep ='\t')
tax.gen <- tax.gen[tax.gen$X%in%names(genus.log),]
tax.gen$X == names(genus.log)
names(genus.log) <- make.unique(paste0('f_', tax.gen$Family, ';g_', tax.gen$Genus))
row.names(genus.log) == row.names(gmm)

genus.gmm <- cbind(genus.log, gmm)
genus.gmm.c <- cor(genus.gmm, method = 'spearman')
genus.gmm.p <- cor.mtest(genus.gmm, method = 'spearman')

corrplot(genus.gmm.c[1:146, 147:247], p.mat = genus.gmm.p$p[1:146, 147:247], 
         insig = 'blank', tl.cex = .3)
