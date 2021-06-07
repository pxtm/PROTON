#############################################################################################################
#################################### PROTON metaG first approach ############################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 25.05.2021
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
library(corrplot)

## load functions
source('lm_function.R'); source('QMP.R'); source('parwise.adonis.r')

## data
# mtdt <- xlsx::read.xlsx('Cohort_overview/list_group.xlsx', header=T, row.names = 1, sheetIndex = 1)
mtdt <- read.table('Meta_data_new_2.csv', header=T, row.names=1, sep=',')
mtdt$id <- paste0('s', mtdt$id)
mtdt$group_name <- factor(mtdt$group_name, levels=c('Controls', 'Normo', 'Micro', 'Macro'))
mtdt$diet <- interaction(factor(mtdt$PC1quant), factor(mtdt$PC2quant), factor(mtdt$PC3quant))
mtdt.red <- read.table('mtdt_selected.txt', header=T, row.names=1, sep='\t')
# mtdt.red <- mtdt.red[match(row.names(mtdt), row.names(mtdt.red)),]
# row.names(mtdt) == row.names(mtdt.red)
# mtdt.red$diet <- mtdt$diet

tax <- read.table('tax.txt', header=T, row.names = 1, sep='\t')
load('Shotgun_QC/Shotgun_QC_processed/data_QC_CM/kuhtod.mgsCounts.RData')
mgs.counts <- as.data.frame(kuhtod.mgscounts); rm(kuhtod.mgscounts)
load('Shotgun_QC/Shotgun_QC_processed/data_QC_CM/kuhtod.kocounts.RData')
ko.counts <- as.data.frame(kuhtod.koCounts); rm(kuhtod.koCounts)

cell.count <- read.table('Cell_count_gutMicrobiome/cellcount_final190910.CSV', header=T, sep=';')

## gene richness
gene.rich <- rowSums(ko.counts!=0)
names(gene.rich) == mtdt$id
gene.rich <- data.frame('sample' = mtdt$id, 'group' = mtdt$group, 'group_name'=mtdt$group_name, 'GR'=unname(gene.rich))

ggplot(gene.rich, aes(group_name, GR, group = interaction(group_name))) +
  geom_violin(aes(color=group_name, fill=group_name), alpha=.3) + 
  geom_quasirandom(aes(color=group_name), size = 3) +
  geom_boxplot(width = .1, alpha = .5, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  stat_compare_means() + 
  stat_compare_means(comparisons = list(c('Controls','Normo'), c('Controls', 'Micro'),
                                        c('Controls', 'Macro'), c('Normo', 'Micro'),
                                        c('Normo', 'Macro'), c('Micro', 'Macro'))) + 
  scale_color_taylor(palette = 'taylor1989') +
  scale_fill_taylor(palette = 'taylor1989') +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, colour = "black"), 
        panel.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        legend.position = 'none') +
  labs(x = "Patient's group", y = "Number of different KOs", color = 'Group', fill ='Group')

hist(gene.rich$GR)
gene.rich$tertile <- as.factor(ntile(gene.rich$GR, 3))
ggplot(gene.rich, aes(group, fill = tertile)) + 
  geom_bar(position = 'dodge')
gene.rich %>% 
  group_by(group, tertile) %>% 
  tally()

# GR associated to bioclinical?
mtdt.red2 <- mtdt.red[complete.cases(mtdt.red),]
gene.rich <- gene.rich[gene.rich$sample%in%paste0('s', row.names(mtdt.red2)),]
gene.rich$sample == paste0('s', row.names(mtdt))
cor.test(gene.rich$GR, as.numeric(mtdt.red2$diet))
cor(gene.rich$GR, sapply(mtdt.red2[,2:47], as.numeric))
lm(gene.rich$GR~mtdt.red2$diet)
mtdt.red2$GR <- gene.rich$GR
ggplot(mtdt.red2, aes(diet, GR)) + 
  geom_boxplot() +
  geom_

## QMPs computation - MGS with samples in rows | cell counts samples as rows -----
# mgs.counts.ds <- as.data.frame(Rarefy(mgs.counts)$otu.tab.rff)
# # saveRDS(mgs.counts.ds, 'rarefied_MGS.rds')
# cell.count$ID <- paste0('s', cell.count$ID)
# cell.count <- cell.count[cell.count$ID%in%row.names(mgs.counts.ds),]
# cell.count <- cell.count[match(row.names(mgs.counts.ds), cell.count$ID),]
# row.names(cell.count) <- cell.count$ID; cell.count$ID <- NULL
# row.names(cell.count) == row.names(mgs.counts.ds)
# qmp <- rarefy_even_sampling_depth(mgs.counts.ds, cell.count)
# qmp <- as.data.frame(qmp)
# qmp <- qmp[,1:1273]
# saveRDS(qmp, 'qmp.rds')

## diversity analyses ----
qmp <- readRDS('qmp.rds')
MGS.not0 <- colSums(qmp!=0)
qmp.ha <- qmp[,names(qmp)%in%names(MGS.not0[MGS.not0>=21])]
qmp.ha <- as.data.frame(t(qmp.ha)) ## put samples in columns
qmp2 <- as.data.frame(sapply(qmp.ha, as.integer))
qmp.m <- as.matrix(qmp2); colnames(qmp.m) <- names(qmp.ha); rownames(qmp.m) <- row.names(qmp.ha)
tax <- tax[row.names(tax)%in%rownames(qmp.m),]
tax.m <- as.matrix(tax); colnames(tax.m) <- names(tax); rownames(tax.m) <- row.names(tax)
row.names(mtdt) <- mtdt$id
physeq <-  phyloseq(otu_table(qmp.m, taxa_are_rows = T), tax_table(tax.m), sample_data(mtdt))

## basic numbers -----
qmp.ha
tax.ha <- tax[row.names(tax)%in%row.names(qmp.ha),]
sapply(tax.ha, function(x)(length(unique(x))))
100-(nrow(tax.ha[tax.ha$Kingdom=='unclassified',])/990)*100
100-(nrow(tax.ha[tax.ha$Phylum=='unclassified',])/990)*100
100-(nrow(tax.ha[tax.ha$Class=='unclassified',])/990)*100
100-(nrow(tax.ha[tax.ha$Order=='unclassified',])/990)*100
100-(nrow(tax.ha[tax.ha$Family=='unclassified',])/990)*100
100-(nrow(tax.ha[tax.ha$Genus=='unclassified',])/990)*100
100-(nrow(tax.ha[tax.ha$Species=='unclassified',])/990)*100

## alpha-div ====
alphas.qmp <- estimate_richness(physeq)
alphas.qmp <- alphas.qmp[names(alphas.qmp)%in%c('Observed', 'Shannon', 'Simpson', 'InvSimpson', 'Fisher')]
alphas.qmp$class <- mtdt$group_name
alphas.qmp <- reshape2::melt(alphas.qmp)
ggplot(alphas.qmp, aes(class, value, group = interaction(class))) + 
  geom_violin(aes(color=class, fill=class), alpha=.3) + 
  geom_quasirandom(aes(color=class), size = 3) +
  geom_boxplot(width = .1, alpha = .5, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  stat_compare_means() + 
  stat_compare_means(comparisons = list(c('Controls','Normo'), c('Controls', 'Micro'),
                                        c('Controls', 'Macro'), c('Normo', 'Micro'),
                                        c('Normo', 'Macro'), c('Micro', 'Macro'))) + 
  facet_wrap(~variable, scales = 'free') + 
  scale_color_taylor(palette = 'taylor1989') +
  scale_fill_taylor(palette = 'taylor1989') +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, colour = "black"), 
        panel.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        legend.position = 'none', strip.background = element_rect(fill=NA)) +
  labs(x = "Patient's group", color = 'Group', fill ='Group')

## beta-div ====
byc.qmp <- distance(physeq, method = 'bray')
byc.pcoa.qmp <- pcoa(byc.qmp)
byc.qmp.scores <- as.data.frame(byc.pcoa.qmp$vectors)
ggplot(byc.qmp.scores, aes(Axis.1, Axis.2, color = mtdt$group_name)) +
  geom_point() + stat_ellipse() +
  scale_color_taylor(palette='taylor1989')

jsd.qmp <- distance(physeq, method = 'jsd')
jsd.pcoa.qmp <- pcoa(jsd.qmp)
jsd.qmp.scores <- as.data.frame(jsd.pcoa.qmp$vectors)
jsd.qmp.vars <- (jsd.pcoa.qmp$values$Relative_eig)*100
p1 <- ggplot(jsd.qmp.scores, aes(Axis.1, Axis.2, color = mtdt$group_name)) +
  geom_point(size=3) + stat_ellipse() +
  scale_color_taylor(palette='taylor1989') +
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  theme_bw() +
  labs(x = paste0('Principal Coordinate 1 (expl. var ', round(jsd.qmp.vars[1], 2), '%)'),
       y = paste0('Principal Coordinate 2 (expl. var ', round(jsd.qmp.vars[2], 2), '%)'),
       color = 'Albuminuria')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

adonis2(t(qmp2)~mtdt$group_name)
pairwise.adonis(t(qmp2), mtdt$group_name)

ggplot(jsd.qmp.scores, aes(Axis.1, Axis.2, color = mtdt$diet)) +
  geom_point(size=3) + 
  scale_color_taylor(palette='taylorRed') +
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  theme_bw() +
  labs(x = paste0('Principal Coordinate 1 (expl. var ', round(jsd.qmp.vars[1], 2), '%)'),
       y = paste0('Principal Coordinate 2 (expl. var ', round(jsd.qmp.vars[2], 2), '%)'),
       color = 'Albuminuria')

## cluster to genus 
ph.genus <- tax_glom(physeq, taxrank = rank_names(physeq)[6])
genus.qmp <- as.data.frame(t(ph.genus@otu_table))
names(genus.qmp) <- make.unique(paste0(as.data.frame(ph.genus@tax_table)$Family, '_', as.data.frame(ph.genus@tax_table)$Genus))
# for enterotype
# saveRDS(genus.qmp, 'genusqmp.rds')
#
sort(colSums(genus.qmp), decreasing = T)
top15 <- c('Bacteroidaceae_Bacteroides', 'Ruminococcaceae_Ruminococcus', 'Ruminococcaceae_Faecalibacterium', 'Rikenellaceae_Alistipes', 'Lachnospiraceae_Blautia',
           'Lachnospiraceae_Roseburia', 'Clostridiaceae_Clostridium', 'Eubacteriaceae_Eubacterium', 'Lachnospiraceae_Coprococcus', 'Bifidobacteriaceae_Bifidobacterium',
           'Prevotellaceae_Prevotella', 'Ruminococcaceae_Subdoligranulum', 'Lachnospiraceae_Dorea', 'Akkermansiaceae_Akkermansia', 'Veillonellaceae_Dialister')
genus.top <- genus.qmp[,names(genus.qmp)%in%top15]
genus.low <- genus.qmp[,!(names(genus.qmp)%in%top15)]

names(genus.top) <- unlist(strsplit(names(genus.top), '_'))[c(F,T)]
genus.plot <- cbind(genus.top, 'Other' = rowSums(genus.low))
genus.plot$individual <- row.names(genus.plot)
genus.plot$group <- mtdt$group_name
genus.plot.m <- reshape2::melt(genus.plot)

genus.plot.m$group <- factor(genus.plot.m$group, levels=c('Controls', 'Normo', 'Micro', 'Macro'))
genus.plot.m$individual <- factor(genus.plot.m$individual, levels=names(sort(rowSums(genus.plot[,1:16]), decreasing = T)))

ggplot(genus.plot.m, aes(individual, value/(10^8), fill=variable)) + 
  geom_col() +
  scale_fill_taylor(palette='lover') +
  facet_wrap(~group, scales='free_x') +
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = .5),
        legend.text = element_text(face = 'italic'),
        panel.border = element_rect(fill=NA), panel.background = element_blank(),
        strip.background = element_rect(fill=NA)) +
  labs(x = 'Individual', y = 'QMP counts (at genus level)/10^8', fill = 'Genus')

## MGS differences between groups ----
# control vs DM
mtdt$biclass <- plyr::revalue(mtdt$group_name, c('Normo' = 'DM1', 'Micro' = 'DM1', 'Macro' = 'DM1'))
qmp.ha2 <- as.data.frame(t(qmp.ha))
biclass <- lm.associations(mtdt, qmp.ha2, variables = 'biclass', control_by = c('age', 'sex', 'race', 'bmi', 'diet'))$biclass
biclass[biclass$fdr<.1,]$mgs
tax[row.names(tax)%in%biclass[biclass$fdr<.1,]$mgs, ]

# differences within DM
biclass[biclass$fdr<.1,]$mgs
biclass.plot <- data.frame(qmp.ha2[,names(qmp.ha2)%in%biclass[biclass$fdr<.1,]$mgs], 'group' = mtdt$group_name)
biclass.plot <- reshape2::melt(biclass.plot)
ggplot(biclass.plot, aes(group, value/(10^8))) +
  geom_violin(aes(color=group, fill=group), alpha=.3) + 
  geom_quasirandom(aes(color=group), size = 3) +
  geom_boxplot(width = .1, alpha = .5, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  stat_compare_means() + 
  stat_compare_means(comparisons = list(c('Controls','Normo'), c('Controls', 'Micro'),
                                        c('Controls', 'Macro'), c('Normo', 'Micro'),
                                        c('Normo', 'Macro'), c('Micro', 'Macro'))) + 
  facet_wrap(~variable, scales = 'free') + 
  scale_color_taylor(palette = 'taylor1989') +
  scale_fill_taylor(palette = 'taylor1989') +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, colour = "black"), 
        panel.background = element_rect(fill = NA), panel.border = element_rect(fill=NA)) +
  labs(x = "Patient's group",color = 'Group', fill ='Group')

## samples.lists
ctl.nor <- mtdt[mtdt$group_name=='Controls'|mtdt$group_name=='Normo',]$id
ctl.micro <- mtdt[mtdt$group_name=='Controls'|mtdt$group_name=='Micro',]$id
ctl.macro <- mtdt[mtdt$group_name=='Controls'|mtdt$group_name=='Macro',]$id
nor.micro <- mtdt[mtdt$group_name=='Normo'|mtdt$group_name=='Micro',]$id
nor.macro <- mtdt[mtdt$group_name=='Normo'|mtdt$group_name=='Macro',]$id
micro.macro <- mtdt[mtdt$group_name=='Micro'|mtdt$group_name=='Macro',]$id

## Control vs Normo ====
ctl.nor.lm <- lm.associations(mtdt[mtdt$id%in%ctl.nor, ],  qmp.ha2[row.names(qmp.ha2)%in%ctl.nor,], variables = 'group_name', control_by = c('age', 'sex', 'race', 'bmi', 'diet'))$group_name
ctl.nor.lm[ctl.nor.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.nor.lm[ctl.nor.lm$fdr<.1,]$mgs, ]

## Control vs Micro ====
ctl.micro.lm <- lm.associations(mtdt[mtdt$id%in%ctl.micro, ],  qmp.ha2[row.names(qmp.ha2)%in%ctl.micro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
ctl.micro.lm[ctl.micro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.micro.lm[ctl.micro.lm$fdr<.1,]$mgs, ]

## Control vs Macro ====
ctl.macro.lm <- lm.associations(mtdt[mtdt$id%in%ctl.macro, ],  qmp.ha2[row.names(qmp.ha2)%in%ctl.macro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
ctl.macro.lm[ctl.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.macro.lm[ctl.macro.lm$fdr<.1,]$mgs, ]

## Normo vs Micro ====
nor.micro.lm <- lm.associations(mtdt[mtdt$id%in%nor.micro, ],  qmp.ha2[row.names(qmp.ha2)%in%nor.micro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
nor.micro.lm[nor.micro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%nor.micro.lm[nor.micro.lm$fdr<.1,]$mgs, ]

## Normo vs Macro ====
nor.macro.lm <- lm.associations(mtdt[mtdt$id%in%nor.macro, ],  qmp.ha2[row.names(qmp.ha2)%in%nor.macro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
nor.macro.lm[nor.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%nor.macro.lm[nor.macro.lm$fdr<.1,]$mgs, ]

## Micro vs Macro ====
micro.macro.lm <- lm.associations(mtdt[mtdt$id%in%micro.macro, ],  qmp.ha2[row.names(qmp.ha2)%in%micro.macro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
micro.macro.lm[micro.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%micro.macro.lm[micro.macro.lm$fdr<.1,]$mgs, ]

