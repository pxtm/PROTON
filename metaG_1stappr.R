#############################################################################################################
#################################### PROTON metaG first approach ############################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 25.05.2021
## last edit: 01.06.2021

## libraries
library(here)
library(microbiomeutilities)
library(tidyverse)
library(ggbeeswarm)
library(tayloRswift)
library(ggpubr)
library(phyloseq)
library(GUniFrac)
library(ape)
library(corrplot)
library(DESeq2)
library(grid)
library(ggrepel)
library(RColorBrewer)
# library(ggtern)

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
cor.test(gene.rich$GR, as.numeric(as.factor(mtdt.red2$diet)))
mtdt.red2$diet <-  as.numeric(as.factor(mtdt.red2$diet))
cor(gene.rich$GR, sapply(mtdt.red2[,2:47], as.numeric))

lm(gene.rich$GR~mtdt.red2$diet)
mtdt.red2$GR <- gene.rich$GR
ggplot(mtdt.red2, aes(diet, GR)) + 
  geom_boxplot() 

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
byc.qmp <- phyloseq::distance(physeq, method = 'bray')
byc.pcoa.qmp <- pcoa(byc.qmp)
byc.qmp.scores <- as.data.frame(byc.pcoa.qmp$vectors)
ggplot(byc.qmp.scores, aes(Axis.1, Axis.2, color = mtdt$group_name)) +
  geom_point() + stat_ellipse() +
  scale_color_taylor(palette='taylor1989')

jsd.qmp <- phyloseq::distance(physeq, method = 'jsd')
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
# qmp.ha2 <- compositions::clr(qmp.ha2)
qmp.log <- log(qmp.ha2+1) # log transform counts for modelling
biclass <- lm.associations(mtdt, qmp.log, variables = 'biclass', control_by = c('age', 'sex', 'race', 'bmi', 'diet'))$biclass
# biclass[biclass$fdr<.1,]$mgs
# tax[row.names(tax)%in%biclass[biclass$fdr<.1,]$mgs, ]
biclass$delta <- as.numeric(biclass$delta)
biclass$delta.low <- as.numeric(biclass$delta.low)
biclass$delta.up <- as.numeric(biclass$delta.up)
# add variables for plotting purposes
row.names(tax) == biclass$mgs
biclass$bacteria <- tax$Name
for (i in seq(1:nrow(biclass))){
  if (biclass[i,]$bacteria == 'unclassified sp.'){
    biclass[i,]$bacteria <- biclass[i,]$mgs
  } else {}
}

for (i in seq(1:nrow(biclass))){
  if (biclass[i,]$bacteria == 'Bacteria sp.'){
    biclass[i,]$bacteria <- paste0(biclass[i,]$mgs, ': Bacteria sp.')
  } else {}
}

biclass$phylum <- factor(tax$Phylum, levels = c('Actinobacteria', 'Bacteroidetes', 'Candidatus Melainabacteria', 'Euryarchaeota',
                                                'Firmicutes', 'Lentisphaerae', 'Proteobacteria', 'Spirochaetes', 'Stramenopiles', 
                                                'Synergistetes', 'Verrucomicrobia', 'unclassified'))
biclass$prevalence <- (colSums(!(qmp.ha2==0))/nrow(qmp.ha2))*100
biclass$qmp.relabd <- colSums(qmp.ha2)/sum(colSums(qmp.ha2))
write.table(biclass, 'diffs_CTLT1D_QMP.txt', sep = '\t')

# add colors for phylum
colourCount = length(levels(biclass$phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#
volcano.all <- ggplot(biclass, aes(delta*-1, -log10(fdr))) +
  geom_point(aes(size = prevalence, color = phylum, alpha=qmp.relabd)) + 
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_alpha(range = c(0.4, 0.8)) +
  annotation_custom(grobTree(textGrob('Controls', x = .1, y = .15, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('T1D', x = .9, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  geom_text_repel(aes(label = ifelse(fdr<.1&abs(delta)>=.3, bacteria, '')), fontface = 'italic') +
  theme_bw() + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = getPalette(colourCount)) +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"), 
       title = 'Differential MGS between healthy individuals (n=50) and T1D-diagnosed patients (n=161)',
       color = 'Phylum', size = 'MGS prevalence\nin the cohort', alpha = 'Mean rel. abundance\nin the total cohort')

ctl.mgs <- biclass %>% 
  filter(fdr<=.1) %>% 
  filter(delta > 0) %>% 
  mutate('Bact.name' = paste0(mgs, ': ', bacteria)) %>% 
  arrange(bacteria) %>% 
  mutate('Bact.name' = factor(Bact.name, levels = rev(Bact.name))) %>% 
  ggplot(aes(delta*-1, Bact.name, fill='darkolivegreen3')) +
    geom_col() + xlim(-0.5, 0) +
    theme_bw() + theme(plot.title = element_text(size = 10), axis.text.y = element_text(face = 'italic')) +
    scale_fill_manual(values  = 'darkolivegreen2') + guides(fill = F) +
    labs(x = "Cliff's Delta effect size", y = 'Metagenomic species (MGS)', title = 'Decreased MGS\nin T1D patients')

t1d.mgs <- biclass %>% 
  filter(fdr<=.1) %>% 
  filter(delta < 0) %>% 
  mutate('Bact.name' = paste0(mgs, ': ', bacteria)) %>% 
  arrange(bacteria) %>% 
  mutate('Bact.name' = factor(Bact.name, levels = rev(Bact.name))) %>% 
  ggplot(aes(delta*-1, Bact.name, fill ='red')) +
    geom_col() +  xlim(0,0.5) +
    theme_bw() + theme(plot.title = element_text(size = 10),  axis.text.y = element_text(face = 'italic')) +
    scale_fill_manual(values  = 'firebrick4') + guides(fill = F) +
    scale_y_discrete(position = "right") +
    labs(x = "Cliff's Delta effect size", y = NULL, title = 'Increased MGS\nin T1D patients')

CTLvsT1D.plot <- cowplot::plot_grid(volcano.all, ctl.mgs, t1d.mgs, ncol = 3, rel_widths = c(1.5, .8, .8))
CTLvsT1D.plot <- cowplot::add_sub(CTLvsT1D.plot, "For the volcano plot, Cliff's Delta has been inverted for visualization purposes. Points are labelled if FDR < 10% and Cliff's Delta absolute value > 0.3. Inversion of the Cliff's Delta values is maintained for interpretation coherence in the individual MGS effect sizes barplots",
                                  fontface = 'italic', size = 10)
cowplot::ggdraw(CTLvsT1D.plot)

## are the QMP counts of the differential MGS associated to years diagnosed?
mtdt.t1d <- mtdt[mtdt$biclass=='DM1',]
qmp.log.t1d <- qmp.log[row.names(qmp.log)%in%row.names(mtdt.t1d),]
qmp.log.t1d <- qmp.log.t1d[,names(qmp.log.t1d)%in%biclass[biclass$fdr<=.1,]$mgs]
cor.p <- sapply(qmp.log.t1d, function(a){
  cor.test(a, mtdt.t1d$dm_duration)$p.val
})
cor.rho <- sapply(qmp.log.t1d, function(a){
  cor.test(a, mtdt.t1d$dm_duration)$estimate
})

cor.qmp.duration.t1d <- data.frame('MGS' = gsub('.cor', '', names(cor.rho)), 'rho' = cor.rho, 'pval' = cor.p, 'fdr' = p.adjust(cor.p, method='fdr'))
mgs.to.plot.duration <- cor.qmp.duration.t1d[cor.qmp.duration.t1d$fdr<=.1,]$MGS
mgs.dur.plot <- data.frame(qmp.log.t1d[,names(qmp.log.t1d)%in%mgs.to.plot.duration], 'duration' = mtdt.t1d$dm_duration)
mgs.dur.plot <- reshape2::melt(mgs.dur.plot, id.vars = 'duration')
ggplot(mgs.dur.plot, aes(value, duration)) + 
  geom_point() +
  geom_smooth(method='lm', se=F) +
  stat_cor(method='spearman') +
  facet_wrap(~variable, scales = 'free')

# differences within DM
biclass[biclass$fdr<.1,]$mgs
biclass.plot <- data.frame(qmp.log[,names(qmp.log)%in%biclass[biclass$fdr<.1,]$mgs], 'group' = mtdt$group_name)
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
mtdt.ctl.normo <- mtdt[mtdt$id%in%ctl.nor,]
qmp.ctl.normo <- qmp.ha2[row.names(qmp.ha2)%in%ctl.nor,]
MGS.not0 <- colSums(qmp.ctl.normo!=0)
qmp.ctl.normo <- qmp.ctl.normo[,names(qmp.ctl.normo)%in%names(MGS.not0[MGS.not0>=10])]
qmp.log.ctl.normo <- log(qmp.ctl.normo+1)
mtdt.ctl.normo$group <- factor(mtdt.ctl.normo$group)
ctl.nor.lm <- lm.associations(mtdt.ctl.normo,  qmp.log.ctl.normo, variables = 'group', control_by = c('age', 'sex', 'race', 'bmi', 'diet'))$group
ctl.nor.lm[ctl.nor.lm$pval<.05,]$mgs
ctl.nor.lm[ctl.nor.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.nor.lm[ctl.nor.lm$fdr<=.1,]$mgs, ]

## Control vs Micro ====
mtdt.ctl.micro <- mtdt[mtdt$id%in%ctl.micro,]
qmp.ctl.micro <- qmp.ha2[row.names(qmp.ha2)%in%ctl.micro,]
MGS.not0 <- colSums(qmp.log.ctl.micro!=0)
qmp.ctl.micro <- qmp.ctl.micro[,names(qmp.ctl.micro)%in%names(MGS.not0[MGS.not0>=10])]
qmp.log.ctl.micro <- log(qmp.ctl.micro+1)
mtdt.ctl.micro$group <- factor(mtdt.ctl.micro$group)
ctl.micro.lm <- lm.associations(mtdt.ctl.micro, qmp.log.ctl.micro, variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
ctl.micro.lm[ctl.micro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.micro.lm[ctl.micro.lm$fdr<.1,]$mgs, ]
ctl.micro.lm$delta <- as.numeric(ctl.micro.lm$delta)
ctl.micro.lm$delta.low <- as.numeric(ctl.micro.lm$delta.low)
ctl.micro.lm$delta.up <- as.numeric(ctl.micro.lm$delta.up)

# add variables for plotting purposes
tax.ctl.micro <- tax[row.names(tax)%in%ctl.micro.lm$mgs,]
row.names(tax.ctl.micro) == ctl.micro.lm$mgs
ctl.micro.lm$bacteria <- tax.ctl.micro$Name
for (i in seq(1:nrow(ctl.micro.lm))){
  if (ctl.micro.lm[i,]$bacteria == 'unclassified sp.'){
    ctl.micro.lm[i,]$bacteria <- ctl.micro.lm[i,]$mgs
  } else {}
}

for (i in seq(1:nrow(ctl.micro.lm))){
  if (ctl.micro.lm[i,]$bacteria == 'Bacteria sp.'){
    ctl.micro.lm[i,]$bacteria <- paste0(ctl.micro.lm[i,]$mgs, ': Bacteria sp.')
  } else {}
}

ctl.micro.lm$phylum <- factor(tax.ctl.micro$Phylum, levels = c('Actinobacteria', 'Bacteroidetes', 'Candidatus Melainabacteria', 'Euryarchaeota',
                                                'Firmicutes', 'Lentisphaerae', 'Proteobacteria', 'Spirochaetes', 'Stramenopiles', 
                                                'Synergistetes', 'Verrucomicrobia', 'unclassified'))
ctl.micro.lm$prevalence <- (colSums(!(qmp.log.ctl.micro==0))/nrow(qmp.log.ctl.micro))*100
ctl.micro.lm$qmp.relabd <- colSums(qmp.log.ctl.micro)/sum(colSums(qmp.log.ctl.micro))
write.table(ctl.micro.lm, 'diffs_CTLMicro_QMP.txt', sep = '\t')

# add colors for phylum
colourCount = length(levels(ctl.micro.lm$phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#
volcano.ctlmicro <- ggplot(ctl.micro.lm, aes(delta*-1, -log10(fdr))) +
  geom_point(aes(size = prevalence, color = phylum, alpha=qmp.relabd)) + 
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_alpha(range = c(0.4, 0.8)) +
  annotation_custom(grobTree(textGrob('Controls', x = .1, y = .15, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Micro', x = .9, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  geom_text_repel(aes(label = ifelse(fdr<.1&abs(delta)>=.3, bacteria, '')), fontface = 'italic') +
  theme_bw() + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = getPalette(colourCount)) +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"), 
       title = 'Differential MGS between healthy individuals (n=50) and T1D-diagnosed patients with micro-albuminuria (n=50)',
       color = 'Phylum', size = 'MGS prevalence\nin the cohort', alpha = 'Mean rel. abundance\nin the total cohort')

ctl.micmgs <- ctl.micro.lm %>% 
  filter(fdr<=.1) %>% 
  filter(delta > 0) %>% 
  mutate('Bact.name' = paste0(mgs, ': ', bacteria)) %>% 
  arrange(bacteria) %>% 
  mutate('Bact.name' = factor(Bact.name, levels = rev(Bact.name))) %>% 
  ggplot(aes(delta*-1, Bact.name, fill='darkolivegreen3')) +
  geom_col() + xlim(-0.5, 0) +
  theme_bw() + theme(plot.title = element_text(size = 10), axis.text.y = element_text(face = 'italic')) +
  scale_fill_manual(values  = 'darkolivegreen2') + guides(fill = F) +
  labs(x = "Cliff's Delta effect size", y = 'Metagenomic species (MGS)', title = 'Decreased MGS\nin micro-albuminuria patients')

micro.mgs <- ctl.micro.lm %>% 
  filter(fdr<=.1) %>% 
  filter(delta < 0) %>% 
  mutate('Bact.name' = paste0(mgs, ': ', bacteria)) %>% 
  arrange(bacteria) %>% 
  mutate('Bact.name' = factor(Bact.name, levels = rev(Bact.name))) %>% 
  ggplot(aes(delta*-1, Bact.name, fill ='red')) +
  geom_col() +  xlim(0,0.5) +
  theme_bw() + theme(plot.title = element_text(size = 10),  axis.text.y = element_text(face = 'italic')) +
  scale_fill_manual(values  = 'firebrick4') + guides(fill = F) +
  scale_y_discrete(position = "right") +
  labs(x = "Cliff's Delta effect size", y = NULL, title = 'Increased MGS\nin micro-albuminuria patients')

CTLvsmicro.plot <- cowplot::plot_grid(volcano.ctlmicro, ctl.micmgs, micro.mgs, ncol = 3, rel_widths = c(1.5, .8, .8))
CTLvsmicro.plot <- cowplot::add_sub(CTLvsmicro.plot, "For the volcano plot, Cliff's Delta has been inverted for visualization purposes. Points are labelled if FDR < 10% and Cliff's Delta absolute value > 0.3. Inversion of the Cliff's Delta values is maintained for interpretation coherence in the individual MGS effect sizes barplots",
                                  fontface = 'italic', size = 10)
cowplot::ggdraw(CTLvsmicro.plot)

## Control vs Macro ====
mtdt.ctl.macro <- mtdt[mtdt$id%in%ctl.macro,]
qmp.ctl.macro <- qmp.ha2[row.names(qmp.ha2)%in%ctl.macro,]
MGS.not0 <- colSums(qmp.ctl.macro!=0)
qmp.ctl.macro <- qmp.ctl.macro[,names(qmp.ctl.macro)%in%names(MGS.not0[MGS.not0>=10])]
qmp.log.ctl.macro <- log(qmp.ctl.macro+1)
mtdt.ctl.macro$group <- factor(mtdt.ctl.macro$group)
ctl.macro.lm <- lm.associations(mtdt.ctl.macro,  qmp.log.ctl.macro, variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
ctl.macro.lm[ctl.macro.lm$fdr<.1,]$mgs

ctl.macro.lm$delta <- as.numeric(ctl.macro.lm$delta)
ctl.macro.lm$delta.low <- as.numeric(ctl.macro.lm$delta.low)
ctl.macro.lm$delta.up <- as.numeric(ctl.macro.lm$delta.up)

# add variables for plotting purposes
tax.ctl.macro <- tax[row.names(tax)%in%ctl.macro.lm$mgs,]
row.names(tax.ctl.macro) == ctl.macro.lm$mgs
ctl.macro.lm$bacteria <- tax.ctl.macro$Name
for (i in seq(1:nrow(ctl.macro.lm))){
  if (ctl.macro.lm[i,]$bacteria == 'unclassified sp.'){
    ctl.macro.lm[i,]$bacteria <- ctl.macro.lm[i,]$mgs
  } else {}
}

for (i in seq(1:nrow(ctl.macro.lm))){
  if (ctl.macro.lm[i,]$bacteria == 'Bacteria sp.'){
    ctl.macro.lm[i,]$bacteria <- paste0(ctl.macro.lm[i,]$mgs, ': Bacteria sp.')
  } else {}
}

ctl.macro.lm$phylum <- factor(tax.ctl.macro$Phylum, levels = c('Actinobacteria', 'Bacteroidetes', 'Candidatus Melainabacteria', 'Euryarchaeota',
                                                               'Firmicutes', 'Lentisphaerae', 'Proteobacteria', 'Spirochaetes', 'Stramenopiles', 
                                                               'Synergistetes', 'Verrucomicrobia', 'unclassified'))
ctl.macro.lm$prevalence <- (colSums(!(qmp.log.ctl.macro==0))/nrow(qmp.log.ctl.macro))*100
ctl.macro.lm$qmp.relabd <- colSums(qmp.log.ctl.macro)/sum(colSums(qmp.log.ctl.macro))
write.table(ctl.macro.lm, 'diffs_CTLMacro_QMP.txt', sep = '\t')

# add colors for phylum
colourCount = length(levels(ctl.macro.lm$phylum))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#
volcano.ctlmacro <- ggplot(ctl.macro.lm, aes(delta*-1, -log10(fdr))) +
  geom_point(aes(size = prevalence, color = phylum, alpha=qmp.relabd)) + 
  geom_hline(yintercept = -log10(.1), linetype = 'dashed', color = 'red') +
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'orange') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  scale_alpha(range = c(0.4, 0.8)) +
  annotation_custom(grobTree(textGrob('Controls', x = .1, y = .15, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Macro', x = .9, y = .1, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  geom_text_repel(aes(label = ifelse(fdr<.1&abs(delta)>=.3, bacteria, '')), fontface = 'italic') +
  theme_bw() + theme(plot.title = element_text(size = 10)) +
  scale_color_manual(values = getPalette(colourCount)) +
  labs(x = "Cliff's Delta effect size", y = bquote("-"~log[10]~"(FDR)"), 
       title = 'Differential MGS between healthy individuals (n=50) and T1D-diagnosed patients with macro-albuminuria (n=50)',
       color = 'Phylum', size = 'MGS prevalence\nin the cohort', alpha = 'Mean rel. abundance\nin the total cohort')

ctl.macmgs <- ctl.macro.lm %>% 
  filter(fdr<=.1) %>% 
  filter(delta > 0) %>% 
  mutate('Bact.name' = paste0(mgs, ': ', bacteria)) %>% 
  arrange(bacteria) %>% 
  mutate('Bact.name' = factor(Bact.name, levels = rev(Bact.name))) %>% 
  ggplot(aes(delta*-1, Bact.name, fill='darkolivegreen3')) +
  geom_col() + xlim(-0.6, 0) +
  theme_bw() + theme(plot.title = element_text(size = 10), axis.text.y = element_text(face = 'italic')) +
  scale_fill_manual(values  = 'darkolivegreen2') + guides(fill = F) +
  labs(x = "Cliff's Delta effect size", y = 'Metagenomic species (MGS)', title = 'Decreased MGS\nin macro-albuminuria patients')

macro.mgs <- ctl.macro.lm %>% 
  filter(fdr<=.1) %>% 
  filter(delta < 0) %>% 
  mutate('Bact.name' = paste0(mgs, ': ', bacteria)) %>% 
  arrange(bacteria) %>% 
  mutate('Bact.name' = factor(Bact.name, levels = rev(Bact.name))) %>% 
  ggplot(aes(delta*-1, Bact.name, fill ='red')) +
  geom_col() +  xlim(0,0.6) +
  theme_bw() + theme(plot.title = element_text(size = 10),  axis.text.y = element_text(face = 'italic')) +
  scale_fill_manual(values  = 'firebrick4') + guides(fill = F) +
  scale_y_discrete(position = "right") +
  labs(x = "Cliff's Delta effect size", y = NULL, title = 'Increased MGS\nin macro-albuminuria patients')

CTLvsmacro.plot <- cowplot::plot_grid(volcano.ctlmacro, ctl.macmgs, macro.mgs, ncol = 3, rel_widths = c(1.5, .8, .8))
CTLvsmacro.plot <- cowplot::add_sub(CTLvsmacro.plot, "For the volcano plot, Cliff's Delta has been inverted for visualization purposes. Points are labelled if FDR < 10% and Cliff's Delta absolute value > 0.3. Inversion of the Cliff's Delta values is maintained for interpretation coherence in the individual MGS effect sizes barplots",
                                    fontface = 'italic', size = 10)
cowplot::ggdraw(CTLvsmacro.plot)

## Normo vs Micro ====
mtdt.nor.micro <- mtdt[mtdt$id%in%nor.micro,]
qmp.nor.micro <- qmp.ha2[row.names(qmp.ha2)%in%nor.micro,]
MGS.not0 <- colSums(qmp.nor.micro!=0)
qmp.nor.micro <- qmp.nor.micro[,names(qmp.nor.micro)%in%names(MGS.not0[MGS.not0>=10])]
qmp.log.nor.micro <- log(qmp.nor.micro+1)
mtdt.nor.micro$group <- factor(mtdt.nor.micro$group)
nor.micro.lm <- lm.associations(mtdt.nor.micro,  qmp.log.nor.micro, variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
nor.micro.lm[nor.micro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%nor.micro.lm[nor.micro.lm$fdr<.1,]$mgs, ]

## Normo vs Macro ====
mtdt.nor.macro <- mtdt[mtdt$id%in%nor.macro,]
qmp.nor.macro <- qmp.ha2[row.names(qmp.ha2)%in%nor.macro,]
MGS.not0 <- colSums(qmp.nor.macro!=0)
qmp.nor.macro <- qmp.nor.macro[,names(qmp.nor.macro)%in%names(MGS.not0[MGS.not0>=10])]
qmp.log.nor.macro <- log(qmp.nor.macro+1)
mtdt.nor.macro$group <- factor(mtdt.nor.macro$group)
nor.macro.lm <- lm.associations(mtdt.nor.macro, qmp.log.nor.macro, variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
nor.macro.lm[nor.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%nor.macro.lm[nor.macro.lm$fdr<.1,]$mgs, ]

## Micro vs Macro ====
mtdt.micro.macro <- mtdt[mtdt$id%in%micro.macro,]
qmp.micro.macro <- qmp.ha2[row.names(qmp.ha2)%in%micro.macro,]
MGS.not0 <- colSums(qmp.micro.macro!=0)
qmp.micro.macro <- qmp.micro.macro[,names(qmp.micro.macro)%in%names(MGS.not0[MGS.not0>=10])]
qmp.log.micro.macro <- log(qmp.micro.macro+1)
mtdt.micro.macro$group <- factor(mtdt.micro.macro$group)
micro.macro.lm <- lm.associations(mtdt.micro.macro,  qmp.log.micro.macro, variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
micro.macro.lm[micro.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%micro.macro.lm[micro.macro.lm$fdr<.1,]$mgs, ]

## common differences between controls and micro-macro
ctl.macro.lm[ctl.macro.lm$fdr<=.1,]$mgs
ctl.micro.lm[ctl.micro.lm$fdr<=.1,]$mgs
venn.plot <- VennDiagram::venn.diagram(x= list('Macro-albuminuria' = ctl.macro.lm[ctl.macro.lm$fdr<=.1,]$mgs,
                                   'Micro-albuminuria' = ctl.micro.lm[ctl.micro.lm$fdr<=.1,]$mgs), 
                          filename = NULL, fill = c('red', 'blue'))

grid::grid.newpage()
grid::grid.draw(venn.plot)

##################################################################################################################################################################################################################################################################################################
## Species level ----
spps.ph <- tax_glom(physeq, taxrank = rank_names(physeq)[7])
spps <- as.data.frame(spps.ph@otu_table)
row.names(spps) <- make.unique(spps.ph@tax_table[,'Species'])
spps <- as.data.frame(t(spps))

## Control vs Normo ====
ctl.nor.lm <- lm.associations(mtdt[mtdt$id%in%ctl.nor, ],  spps[row.names(spps)%in%ctl.nor,], variables = 'group_name', control_by = c('age', 'sex', 'race', 'bmi', 'diet'))$group_name
ctl.nor.lm[ctl.nor.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.nor.lm[ctl.nor.lm$fdr<.1,]$mgs, ]

## Control vs Micro ====
ctl.micro.lm <- lm.associations(mtdt[mtdt$id%in%ctl.micro, ],  spps[row.names(spps)%in%ctl.micro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
ctl.micro.lm[ctl.micro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.micro.lm[ctl.micro.lm$fdr<.1,]$mgs, ]

## Control vs Macro ====
ctl.macro.lm <- lm.associations(mtdt[mtdt$id%in%ctl.macro, ],  spps[row.names(spps)%in%ctl.macro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
ctl.macro.lm[ctl.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%ctl.macro.lm[ctl.macro.lm$fdr<.1,]$mgs, ]

## Normo vs Micro ====
nor.micro.lm <- lm.associations(mtdt[mtdt$id%in%nor.micro, ],  spps[row.names(spps)%in%nor.micro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
nor.micro.lm[nor.micro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%nor.micro.lm[nor.micro.lm$fdr<.1,]$mgs, ]

## Normo vs Macro ====
nor.macro.lm <- lm.associations(mtdt[mtdt$id%in%nor.macro, ],  spps[row.names(spps)%in%nor.macro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
nor.macro.lm[nor.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%nor.macro.lm[nor.macro.lm$fdr<.1,]$mgs, ]

## Micro vs Macro ====
micro.macro.lm <- lm.associations(mtdt[mtdt$id%in%micro.macro, ],  spps[row.names(spps)%in%micro.macro,], variables = 'group', control_by = c('age', 'sex', 'race', 'bmi'))$group
micro.macro.lm[micro.macro.lm$fdr<.1,]$mgs
tax[row.names(tax)%in%micro.macro.lm[micro.macro.lm$fdr<.1,]$mgs, ]


## Ternary plots -- only T1D ------
physeq.t1d <- subset_samples(physeq, group_name!='Controls')
tern.plot <- microbiomeutilities::prep_ternary(physeq.t1d, group = 'group_name', level = 'lowest')

library(plotly)
tern.plot$Phylum <- factor(tern.plot$Phylum, levels = c("Actinobacteria", "Bacteroidetes", "Candidatus Melainabacteria", "Euryarchaeota", "Firmicutes",
                                                           "Lentisphaerae", "Proteobacteria", "Spirochaetes", "Synergistetes", "Verrucomicrobia", "unclassified"))
qmp.t1d <- as.data.frame(t(physeq.t1d@otu_table))
qmp.t1d <- qmp.t1d[,names(qmp.t1d)%in%tern.plot$OTUID]
qmp.relabd <- colSums(qmp.t1d)
qmp.relabd <- qmp.relabd/sum(qmp.relabd)
names(qmp.relabd) == tern.plot$OTUID
tern.plot$total.abd <- qmp.relabd

p <- plot_ly(tern.plot, a = ~Normo, b=~Micro, c=~Macro,mode = "markers", 
             type = "scatterternary", fillcolor = ~Phylum, asrc = 'Normo',
             colors = 'Set3', size = ~total.abd)

layout <- list(
  margin = list(
    b = 40, 
    l = 60, 
    r = 10, 
    t = 25
  ), 
  ternary = list(
    aaxis = list(title = list(text = "Normo-albuminuria")), 
    baxis = list(title = list(text = "Micro-albuminuria")), 
    caxis = list(title = list(text = "Macro-albuminuria"))
  ), 
  hovermode = "closest", 
  showlegend = TRUE
)

p <- layout(p, margin=layout$margin, ternary=layout$ternary, hovermode=layout$hovermode, showlegend=layout$showlegend)
p


legend.sizes <- seq(0, 1, 0.2)
ax = list(zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
mk = list(sizeref=0.1, sizemode="area")
p.legend = plot_ly() %>%
  add_markers(x = 1, y = legend.sizes, size = legend.sizes, showlegend = F, marker = mk) %>%
  layout(xaxis = ax, yaxis = list(showgrid = FALSE))

subplot(p.legend, p, widths = c(0.1, 0.9))

api_create(p, filename = 'TernaryT1D')
