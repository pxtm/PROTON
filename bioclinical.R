#############################################################################################################
################################## PROTON metadata first approach ###########################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 31.05.2021
## last edit: 01.06.2021

## libraries
library(here)
library(tidyverse)
library(ggbeeswarm)
library(tayloRswift)
library(ggpubr)
library(phyloseq)
library(GUniFrac)
library(ape)
library(factoextra)
library(ggrepel)

## data
mtdt <- read.table('mtdt_selected.txt', header=T, row.names=1, sep='\t')
mtdt[,2:46] <- sapply(mtdt[,2:46], as.numeric)
mtdt$Groups <- factor(mtdt$Groups, levels=c('Controls', 'Normo', 'Micro', 'Macro'))
sapply(mtdt, function(x)(sum(is.na(x))))

## work with complete cases for the PCA
mtdt2 <- mtdt[complete.cases(mtdt),]
pca <- prcomp(mtdt2[,2:46])
fviz_screeplot(pca)
fviz_pca_biplot(pca)

fviz_pca_ind(pca, geom='point', col.ind = mtdt2$Groups, addEllipses = T)
scores.pca <- as.data.frame(pca$x)
scores.plot <- ggplot(data = scores.pca, aes(PC1, PC2, color = mtdt2$Groups)) +
  geom_point() +
  stat_ellipse() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_color_taylor(palette = 'taylor1989') +
  theme_bw() +
  theme(legend.position = c(.9, .2)) +
  labs(x = 'Prin. Component 1 (expl. var. 53.9%)', y = 'Prin. Component 2 (expl. var. 18.9%)', color = 'Groups')

loadings.pca <- as.data.frame(pca$rotation)
load.plot <- ggplot(loadings.pca, aes(PC1, PC2, label=row.names(loadings.pca))) +
  geom_segment(aes(x=0, y=0, xend=PC1, yend=PC2), arrow=arrow(length=unit(.3, 'cm')), alpha = .7) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_text_repel(aes(label = ifelse(abs(PC1)>0.3|abs(PC2)>0.3, as.character(row.names(loadings.pca)), ''))) +
  theme_bw() +
  labs(x = 'Prin. Component 1 (expl. var. 53.9%)', y = 'Prin. Component 2 (expl. var. 18.9%)')

cowplot::plot_grid(ggExtra::ggMarginal(
  p = scores.plot,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
), load.plot, rel_widths = c(1.6, 1), labels=c('A', 'B'))

## anovas
anova.tukey <- function(dataframe, comparison, variables){
  if(missing(comparison)){
    stop('Comparisor missing!')
  } else {
    output <- list()
    for(i in variables){
      f <- paste0(i, '~', comparison)
      model <- lm(f, data = dataframe)
      model.aov <- aov(model)
      p.model <- summary(model.aov)[[1]]$`Pr(>F)`[1]
      p.pair <- TukeyHSD(model.aov)[[1]][,4]
      output[[i]] <- c('p.model'  = p.model, p.pair)
    }
    return(output)
  }
}

bioclin.comp <- anova.tukey(mtdt, comparison = 'Groups', variables = names(mtdt)[2:46])
bioclin.comp <- as.data.frame(t(do.call(cbind.data.frame, bioclin.comp)))
p.adjust(bioclin.comp$p.model, method ='fdr')
diff.vars <- bioclin.comp %>% 
  add_column('fdr' = p.adjust(bioclin.comp$p.model, method ='fdr'), .after = 'p.model') %>% 
  filter(fdr<=.1) %>% 
  row.names()
diff.vars <- diff.vars[!(diff.vars%in%c('dm'))]

mtdt.diff <- mtdt[, names(mtdt)%in%c(diff.vars, 'Groups')]
mtdt.diff <- reshape2::melt(mtdt.diff)
ggplot(mtdt.diff, aes(Groups, value, color = Groups)) +
  geom_boxplot(color='black', outlier.colour = NA) +
  geom_quasirandom(size=3, alpha=.3) +
  facet_wrap(~variable, scales='free') +
  stat_compare_means(comparisons = list(c('Controls', 'Normo'), c('Controls', 'Micro'), c('Controls', 'Macro'),
                                        c('Normo', 'Micro'), c('Normo', 'Macro'),
                                        c('Micro', 'Macro'))) +
  scale_color_taylor(palette='taylor1989') + 
  theme_bw() + theme(legend.position = 'none')

## heatmap
is.na(mtdt) %>% table()
pheatmap::pheatmap(mtdt2[,2:46], scale = 'column', cluster_rows = F, cluster_cols = F)
