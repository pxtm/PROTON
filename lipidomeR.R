#############################################################################################################
#################################### PROTON lipidomics lipidomeR ############################################
#############################################################################################################

## libraries
library(here)
library(lipidomeR)

## data
lipids <- read.table('lipids_annotated.txt', header=T, sep = '\t')
mtdt <- read.table('Meta_data_new_2.csv', header=T, row.names=1, sep=',')
mtdt$id <- paste0('s', mtdt$id)
mtdt$group_name <- factor(mtdt$group_name, levels=c('Controls', 'Normo', 'Micro', 'Macro'))
mtdt$diet <- interaction(factor(mtdt$PC1quant), factor(mtdt$PC2quant), factor(mtdt$PC3quant))
mtdt$biclass <- plyr::revalue(mtdt$group_name, c('Normo' = 'DM1', 'Micro' = 'DM1', 'Macro' = 'DM1'))

## mapping lipids
mapping <- map_lipid_names(unique(lipids$Name))

## model lipids
mtdt <- mtdt[mtdt$id%in%names(lipids),]
mtdt <- mtdt[match(names(lipids)[7:217], mtdt$id),]
names(lipids)[7:217] == mtdt$id
lipids2 <- data.frame(t(lipids[,7:217]))
lipids2 <- cbind('sample' = mtdt$id, 'group' = mtdt$group_name, 'diet' = mtdt$diet, 'age' = as.character(mtdt$age), 
                 'sex' = as.factor(mtdt$sex), 'BMI' = as.character(mtdt$bmi), 'biclass' = mtdt$biclass, lipids2)
names(lipids2)[8:479] <- lipids$Name
lipids2$age <- as.numeric(lipids2$age)
lipids2$sex <- as.factor(lipids2$sex)
lipids2$BMI <- as.numeric(lipids2$BMI)
lipids2$diet <- as.numeric(lipids2$diet)
#
## result limma
result.limma <- compute_models_with_limma(x = lipids2, dependent.variables = mapping$Name, independent.variables = c('group'))#, 
figure.output <-
  heatmap_lipidome_from_limma(
    x = result.limma$"model",
    names.mapping = mapping,
    axis.x.carbons = FALSE,
    class.facet = "row",
    plot.all = TRUE,
    plot.individual = FALSE,
    print.figure = TRUE,
    # scales = "free_y"
    space = "free"
  )

figure.output <-
  heatmap_lipidome_from_limma(
    x = result.limma$"model",
    names.mapping = mapping,
    axis.x.carbons = FALSE,
    class.facet = "wrap",
    class.subset = 'TG',
    plot.all = FALSE,
    plot.individual = TRUE,
    print.figure = FALSE,
    space = "free",
    p.adj.method = 'fdr'
  )
cowplot::plot_grid(figure.output[['groupNormo']],
                   figure.output[['groupMicro']],
                   figure.output[['groupMacro']], ncol=3)
