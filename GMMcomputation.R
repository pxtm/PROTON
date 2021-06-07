#############################################################################################################
#################################### PROTON metaG GMM computation ###########################################
#############################################################################################################

## by: Marc Clos-Garcia, PhD
## initial edit: 03.06.2021
## last edit: 

## libraries
library(omixerRpm)
library(GUniFrac)
library(data.table)

# data
load('Shotgun_QC/Shotgun_QC_processed/data_QC_CM/kuhtod.kocounts.RData')
ko.counts <- as.data.frame(kuhtod.koCounts); rm(kuhtod.koCounts)
ko.counts <- as.data.frame(t(ko.counts))

load('Shotgun_QC/IGC.keggOrtholog2geneIndex.RData')
load('Shotgun_QC/ModulDef_boolean_sum.RData')

# rarefy
ko.counts <- as.data.frame(t(ko.counts)) # pass it transposed! KO counts in columns
ko.counts.rar <- Rarefy(ko.counts)$otu.tab.rff
saveRDS(ko.counts.rar, file = 'KOcounts.ds.rds')


##normalize
ko.counts.rar.nor <- t(apply(t(ko.counts), 2, function(i) i/sum(i)))


#####
#KEGG modules ----
#####
koCounts_remove<-ko.counts[, colnames(ko.counts) %in% KEGGset]
#koCounts matrix from CM is missing some KOs (0 counts on all samples) -> find missing KOs
diff_set<-setdiff(KEGGset,colnames(koCounts_remove))

#make 0 value matrix with rownames of all samples and colnames missing KOs found previously + cbin (combine both dataframes)
add<-as.data.frame(matrix(0, ncol=length(diff_set), nrow = length(rownames(koCounts_remove))))
rownames(add)<-rownames(koCounts_remove)
colnames(add)<-diff_set
koCounts_all<-cbind(koCounts_remove,add)

#TRUE and FALSE dataframe of samples and if KO counts are above 0 
sampleHasKo <- as.data.frame(koCounts_all > 0)
#ModulDef_boolean_sum contains the Module KO composition inclding alternative pathways, the overall number of TRUE KOs found in TRUE and FALSE dataframe
ModuleStepCounts <- sapply(ModulDef_boolean_sum, eval, envir = sampleHasKo)
#percentage of actual module count and complete module count
ModuleCompletenessPct <- round(100 * t(t(ModuleStepCounts) / CompleteModuleCount))
#filter with 2/3 of KOs present (can be adjusted)
ModuleCompletenessPct_filter<-as.data.frame(ModuleCompletenessPct > 66)

#KOs abundances
Module_abundance <- sapply(ModulDef_boolean_sum, eval, envir = koCounts_all)
rownames(Module_abundance) <- rownames(ko.counts)

#find abundance of KOs 
Module_abundance[ModuleCompletenessPct_filter == FALSE] <- 0
saveRDS(Module_abundance, file='KEGGcounts.rds')

#rarefy module abundance
Module_abundance.rare<-Rarefy(Module_abundance)$otu.tab.rff

#normalize to module length
Module_abundance.rarelength <- round(t(t(Module_abundance) / CompleteModuleCount))

saveRDS(Module_abundance.rarelength, file='KEGGcounts.ds.rds')


#normalize
Module_abundance.rarelength.nor <- t(apply(t(Module_abundance.rarelength), 2, function(i) i/sum(i)))
saveRDS(Module_abundance.rarelength.nor, file='KEGGabd.rds')



######
#GMM ----
######
tmp<-t(ko.counts)
tmp<-as.data.frame(tmp)

setDT(tmp, keep.rownames = "entry")
#mods2 <- rpm(tmp, minimum.coverage=0.33, annotation = 1)
mods <- rpm(tmp, minimum.coverage=0.66, annotation = 1)
# Load the default mapping database
db <- loadDefaultDB()

# get the name of the first predicted module
getNames(db, mods@annotation[1,])

# get the abundance|coverage as a data.frame with module id and description
annotation<-as.data.frame(mods@annotation)
abundance<-as.data.frame(mods@abundance)

#binding and ready for normalising:
#cbind
abundancerownames<-cbind(annotation,abundance)

abundancerownames2 <- abundancerownames[,-1]
rownames(abundancerownames2) <- abundancerownames[,1]

abundancerownames2 <- t(abundancerownames2)

saveRDS(abundancerownames2, file='GMMcounts.rds')

##normalizing
abundancerownames2nor <- apply(abundancerownames2, 2, function(i) i/sum(i))
saveRDS(abundancerownames2nor, file='GMMabd.rds')
