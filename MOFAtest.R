#############################################################################################################
##################################### PROTON metaG-metabolites ##############################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 27.05.2021
## last edit: 

## libraries
library(MOFA2); library(compositions)

## set reticulate
reticulate::use_condaenv(condaenv = "mofa2", conda = "/Users/pmc959/opt/anaconda3/condabin/conda", required = TRUE)

## data
qmp <- readRDS('qmp.rds')
qmp <- as.data.frame(t(qmp)); qmp.names <- row.names(qmp)
qmp <- as.data.frame(sapply(qmp, function(x)(x/sum(x))))
qmp <- as.data.frame(t(qmp))
names(qmp) <- qmp.names
gc <- readRDS('GCGC_mets.rds')
# lipids <- readRDS('Lipids_mets.rds')
lipids <- readRDS('lipids_clusters.rds')
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

## MOFA analyses  (based on W.Haak&Argelaguet et al, mSystems, 2021) ---- 
gmm <- as.data.frame(clr(gmm))
qmp <- as.data.frame(clr(qmp))
qmp.m <- t(qmp); rownames(qmp.m) <- names(qmp); colnames(qmp.m) <- row.names(qmp)
gc.m <- t(gc); rownames(gc.m) <- names(gc); colnames(gc.m) <- row.names(gc)
lipids <- lipids[,!(names(lipids)%in%names(gc))]
# lipids <- sapply(lipids, log)
lipids.m <- t(lipids); rownames(lipids.m) <- names(lipids); colnames(lipids.m) <- row.names(lipids)
# lipids.names <- data.frame('met' = make.unique(rep('met', 7215)), 'lipid' =rownames(lipids.m))
rownames(lipids.m) <- lipids.names$met
gmm.m <- t(gmm); rownames(gmm.m) <- names(gmm); colnames(gmm.m) <- row.names(gmm)

biochem <- mtdt; biochem$Groups <- NULL; biochem$dm <- NULL
biochem$diet <- as.numeric(as.factor(biochem$diet))
biochem <- as.data.frame(sapply(biochem, as.numeric))
biochem.m <- t(biochem); rownames(biochem.m) <- names(biochem); colnames(biochem.m) <- row.names(mtdt)
rownames(biochem.m)
mofa.data <- list('QMP' = qmp.m, 
                  'GMM' = gmm.m,
                  'GCGC' = gc.m,
                  'Lipidomics' = lipids.m,
                  'Biochem' = biochem.m)
lapply(mofa.data, dim)

MOFAobject <- create_mofa(mofa.data, groups = as.vector(mtdt$Groups))
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
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
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)

## downstream analysis
model <- MOFAobject.trained
plot_data_overview(model)
# add metadata
red.mtdt <- data.frame('sample' = c(unlist(samples_names(model))), 'group' = mtdt$Groups); row.names(red.mtdt) <- rep('group1', 209)
samples_metadata(model) <- red.mtdt
head(model@samples_metadata, n=3)

# plots
plot_variance_explained(model, plot_total = T)[[2]]
plot_variance_explained(model, max_r2=15)

plot_factor(model, 
            factor = 1:3,
            color_by = "group")

plot_factors(model, 
             factors = 1:5,
             color_by = "group")

p <- plot_factors(model, 
                  factors = c(3,2), 
                  color_by = "group", 
                  dot_size = 4
) #+ scale_fill_manual(values=category.colors)

p + 
  # geom_density_2d(aes_string(color="color_by")) +
  stat_ellipse(aes(color=color_by), geom = "polygon", alpha=0.25) 
  #scale_color_manual(values=category.colors)

plot_dimred(model,
            method = "UMAP",  # method can be either "TSNE" or "UMAP"
            color_by = "group"
)

plot_weights(model, factor=3, view="Lipidomics", nfeatures=8) #MGS.hg0769
qmp %>% 
  add_column('group' = mtdt$group_name) %>% 
  select(c('MGS.hg0769', 'group')) %>% 
  ggplot(aes(group, MGS.hg0769/(10^8), color = group)) +
    geom_boxplot() + geom_quasirandom()
tax[which(row.names(tax)=='MGS.hg0769'),]
