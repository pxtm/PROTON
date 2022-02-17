#############################################################################################################
##################################### PROTON metaG-metabolites ##############################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 27.05.2021
## last edit: 15.07.2021

## libraries
library(MOFA2); library(compositions)

## data
qmp <- readRDS('qmp.rds')
qmp <- as.data.frame(t(qmp)); qmp.names <- row.names(qmp)
qmp <- as.data.frame(sapply(qmp, function(x)(x/sum(x))))
qmp <- as.data.frame(t(qmp))
names(qmp) <- qmp.names
gc <- readRDS('GCGC_mets.rds')
# lipids <- readRDS('Lipids_mets.rds')
lipids <- readRDS('lipids_clusters.rds')
lipids.annot <- data.table::fread('lipidclusters_annotation.txt', header = T, sep ='\t', select = c(1:3))
lipids.annot$name <- make.unique(lipids.annot$name)
lipids.annot$cluster <- paste0('ME', lipids.annot$cluster)
lipids.annot$cluster == names(lipids)
names(lipids) <- lipids.annot$name
rm(lipids.annot)

mtdt <- readRDS('mtdt_common.rds')
gmm <- as.data.frame(readRDS('GMMabd.rds'))
  
row.names(qmp) <- gsub('s', '', row.names(qmp))
qmp <- qmp[row.names(qmp)%in%row.names(mtdt),]
row.names(gmm) <- gsub('s', '', row.names(gmm))
gmm <- gmm[row.names(gmm)%in%row.names(mtdt),]
identical(row.names(qmp), row.names(mtdt))
identical(row.names(qmp), row.names(gc))
identical(row.names(qmp), row.names(lipids))
gmm <- gmm[match(row.names(qmp), row.names(gmm)),]
identical(row.names(qmp), row.names(gmm))

# metaG annotations
gmm.annotation <- read.table('GMMs.v1.07.names', sep='\t')
gmm.annotation <- gmm.annotation[gmm.annotation$V1%in%names(gmm),]
gmm <- gmm[,match(gmm.annotation$V1, names(gmm))]
gmm.annotation$V1 == names(gmm)
names(gmm) <- paste(names(gmm), "|", gmm.annotation$V2)

tax <- read.table('tax.txt', header=T, sep = '\t', row.names=1)
identical(row.names(tax),names(qmp))
names(qmp) <-  paste(names(qmp), "|", tax$Name)

## MOFA analyses  (based on W.Haak&Argelaguet et al, mSystems, 2021) ---- 
gmm <- as.data.frame(clr(gmm))
qmp <- as.data.frame(clr(qmp))
qmp.m <- t(qmp); rownames(qmp.m) <- names(qmp); colnames(qmp.m) <- row.names(qmp)
gc.m <- t(gc); rownames(gc.m) <- names(gc); colnames(gc.m) <- row.names(gc)
lipids <- lipids[,!(names(lipids)%in%names(gc))]
# lipids <- sapply(lipids, log)
lipids.m <- t(lipids); rownames(lipids.m) <- names(lipids); colnames(lipids.m) <- row.names(lipids)
# lipids.names <- data.frame('met' = make.unique(rep('met', 7215)), 'lipid' =rownames(lipids.m))
# rownames(lipids.m) <- lipids.names$met
gmm.m <- t(gmm); rownames(gmm.m) <- names(gmm); colnames(gmm.m) <- row.names(gmm)

biochem <- mtdt; biochem$Groups <- NULL; biochem$dm <- NULL
biochem$diet <- as.numeric(as.factor(biochem$diet))
biochem <- as.data.frame(sapply(biochem, as.numeric))
biochem <- log(biochem+1)
biochem.m <- t(biochem); rownames(biochem.m) <- names(biochem); colnames(biochem.m) <- row.names(mtdt)
rownames(biochem.m)
mofa.data <- list('QMP' = qmp.m, 
                  'GMM' = gmm.m,
                  'GCGC' = gc.m,
                  'Lipidomics' = lipids.m,
                  'Biochem' = biochem.m)
lapply(mofa.data, dim)

MOFAobject <- create_mofa(mofa.data)#, groups = as.vector(mtdt$Groups))
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE
head(data_opts)

model_opts <- get_default_model_options(MOFAobject) # ideallu Gaussian
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts, 
)

outfile = file.path(getwd(),"model_ALL_lipidsclustered.hdf5")
# outfile = file.path(getwd(),"model_ALL_lipids.hdf5")

MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)

## downstream analysis
model <- MOFAobject.trained
plot_data_overview(model)
# add metadata
red.mtdt <- data.frame('sample' = c(unlist(samples_names(model))), 'group' = mtdt$Groups)
red.mtdt$group <- factor(red.mtdt$group, levels = c('Controls', 'Normo', 'Micro', 'Macro'))

# plots
plot_variance_explained(model, plot_total = T)[[2]]
plot_variance_explained(model)

plot_factor(model, 
            factor = 1:15,
            color_by = red.mtdt$group,
            group_by = red.mtdt$group,
            dot_alpha = .8) + 
  scale_fill_taylor() +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.border = element_rect(fill=NA))

plot_factors(model, 
             factors = 1:9,
             color_by = red.mtdt$group) + 
  scale_fill_manual(values=alpha(c("#B1532A", "#B39B8C", "#BDADAB", "#43475B"),.7)) + scale_color_taylor()

plot_factors(model, 
             factors = c(2,3), 
             color_by = red.mtdt$group, 
             dot_size = 4
) + scale_color_taylor() + scale_fill_taylor() +stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 
plot_factors(model, 
             factors = c(3,8), 
             color_by = red.mtdt$group, 
             dot_size = 4
) + scale_color_taylor() + scale_fill_taylor() +stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 


plot_factors(model, 
                  factors = c(1,9), 
                  color_by = red.mtdt$group, 
                  dot_size = 4
) + scale_color_taylor() + scale_fill_taylor() +stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 

plot_factors(model, 
             factors = c(1,7), 
             color_by = red.mtdt$group, 
             dot_size = 4
) + scale_color_taylor() + scale_fill_taylor() +stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 

plot_factors(model, 
             factors = c(1,2), 
             color_by = red.mtdt$group, 
             dot_size = 4
) + scale_color_taylor() + scale_fill_taylor() +stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 



ggpubr::ggarrange(plot_top_weights(model, factor=2, view = 'Biochem', nfeatures =10),
                  plot_top_weights(model, factor=2, view = 'GCGC', nfeatures=10))
ggpubr::ggarrange(plot_top_weights(model, factor=3, view = 'Biochem', nfeatures=10),
                  plot_top_weights(model, factor=3, view = 'GCGC', nfeatures=10))
ggpubr::ggarrange(plot_top_weights(model, factor=8, view = 'QMP', nfeatures =10),
          plot_top_weights(model, factor=8, view = 'GMM', nfeatures=10))

plot_data_heatmap(model,
                  view = "Biochem",         # view of interest
                  factor = 2,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

plot_data_heatmap(model,
                  view = "Lipidomics",         # view of interest
                  factor = 2,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

plot_weights_scatter()
## repeat the MOFA approach removing the bioclinical dataset ----
mofa.data2 <- list('QMP' = qmp.m, 
                  'GMM' = gmm.m,
                  'GCGC' = gc.m,
                  'Lipidomics' = lipids.m)
lapply(mofa.data2, dim)

MOFAobject2 <- create_mofa(mofa.data2)#, groups = as.vector(mtdt$Groups))
plot_data_overview(MOFAobject2)

data_opts <- get_default_data_options(MOFAobject2)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject2) # ideallu Gaussian
head(model_opts)

train_opts <- get_default_training_options(MOFAobject2)
head(train_opts)

MOFAobject2 <- prepare_mofa(
  object = MOFAobject2,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts, 
)

outfile = file.path(getwd(),"model2_noBioclin_lipidsclustered.hdf5")
# outfile = file.path(getwd(),"model_ALL_lipids.hdf5")

MOFAobject.trained2 <- run_mofa(MOFAobject2, outfile, use_basilisk = F)

## downstream analysis
model2 <- MOFAobject.trained2
plot_data_overview(model2)

# plots
plot_variance_explained(model2, plot_total = T)[[2]]
plot_variance_explained(model2)

plot_factor(model2, 
            factor = 1:15,
            color_by = red.mtdt$group,
            group_by = red.mtdt$group,
            dot_alpha = .8) + 
  scale_fill_taylor() +
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.border = element_rect(fill=NA))

plot_factors(model2, 
             factors = 1:9,
             color_by = red.mtdt$group) + 
  scale_fill_manual(values=alpha(c("#B1532A", "#B39B8C", "#BDADAB", "#43475B"),.7)) + scale_color_taylor()

plot_factors(model2, 
             factors = c(7,8), 
             color_by = red.mtdt$group, 
             dot_size = 4
) + scale_color_taylor() + scale_fill_taylor() +stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 

plot_factors(model, 
             factors = c(7,9), 
             color_by = red.mtdt$group, 
             dot_size = 4
) + scale_color_taylor() + scale_fill_taylor() +stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 

plot_weights(model, factor=1, nfeatures=8) #MGS.hg0769
plot_top_weights(model, factor =1, view = 'GCGC')

group <- data.frame('group' = as.numeric(red.mtdt$group)); row.names(group) <- red.mtdt$sample
correlate_factors_with_covariates(model, covariates = group)
