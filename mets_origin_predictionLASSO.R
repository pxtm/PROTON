###############################################################
################# integration ##########
###############################################################

## libraries
library(plyr)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(repr)
library(glmnet)
library(here)
library(tidyverse)

## data
here()
qmp <- readRDS('qmp.rds')
qmp <- as.data.frame(t(qmp)); qmp.names <- row.names(qmp)
qmp <- as.data.frame(sapply(qmp, function(x)(x/sum(x))))
qmp <- as.data.frame(t(qmp))
names(qmp) <- qmp.names
gc <- readRDS('GCGC_mets.rds')
lipids <- readRDS('Lipids_mets.rds')
lipids.c <- readRDS('lipids_clusters.rds')
mtdt <- readRDS('mtdt_common.rds')
gmm <- as.data.frame(readRDS('GMMabd.rds'))

row.names(qmp) <- gsub('s', '', row.names(qmp))
row.names(gmm) <- gsub('s', '', row.names(gmm))
qmp <- qmp[row.names(qmp)%in%row.names(gc),]
gmm <- gmm[row.names(gmm)%in%row.names(gc),]

identical(row.names(qmp), row.names(gc))
identical(row.names(qmp), row.names(lipids))
identical(row.names(qmp), row.names(mtdt))
gmm <- gmm[match(row.names(qmp), row.names(gmm)),]
identical(row.names(qmp), row.names(gmm))

tax <- read.table('tax.txt', header=T, sep='\t')
## format and impute mtdt -----
sapply(mtdt, class)
mtdt.num <- mtdt[,!(names(mtdt)%in%c('Groups', 'born_by_csection'))]
mtdt.num$diet <- as.numeric(as.factor(mtdt$diet))
sapply(mtdt.num, summary)
preProcValues <- preProcess(mtdt.num %>% 
                              select(age, sex, bmi, height, weight, waist, hip, dm, 
                                     dm_duration, eGFR, sys_24h, dia_24h, office_sys, office_dia, office_hr, 
                                     hba1c_ifcc, hba1c_pct, hscrp, Celiaki, Faeces_frequency, Faeces_regular, bristol_scale, 
                                     bg, hgb, wbc, basophils, eosinophils, lymphocytes, monocytes, neutrophils, bloodplt, 
                                     HCT, p_potassium, p_sodium, p_crea, p_alb, ALAT, p_chol, p_HDL, p_LDL, 
                                     p_VLDL, p_trig, smoker, alcohol_yesno, diet),
                            method = c('knnImpute'),
                            k = 25,
                            knnSummary = mean)
imputed.mtdt.num <- predict(preProcValues, mtdt.num, na.action=na.pass)
procNames <- data.frame(col = names(preProcValues$mean), mean = preProcValues$mean, sd = preProcValues$std)
for(i in procNames$col){
  imputed.mtdt.num[i] <- imputed.mtdt.num[i]*preProcValues$std[i]+preProcValues$mean[i] 
}

mtdt <- cbind(imputed.mtdt.num, 'groups' = as.numeric(as.factor(mtdt$Groups)), 
              'born_by_csection' = as.numeric(as.factor(mtdt$born_by_csection)))
lifestyle.antro <- mtdt[,names(mtdt)%in%c('age', 'sex', 'bmi', 'height', 'weight', 'waist', 'hip', 'office_hr', 'Faeces_frequency', 'Faeces_regular', 'bristol_scale', 'smoker', 'alcohol_yesno', 'born_by_csection', 'diet', 'Celiaki')]
biochem <- mtdt[,names(mtdt)%in%c('eGFR', 'sys_24h', 'dia_24h', 'office_sys', 'office_dia', 'office_hr', 'hba1c_ifcc', 'hba1c_pct', 'hscrp', 'bg', 'hgb', 'wbc', 'basophils', 'eosinophils', 
                                  'lymphocytes', 'monocytes', 'neutrophils', 'bloodplt', 'HCT', 'p_potassium', 'p_sodium', 'p_crea', 'p_alb', 'ALAT', 'p"chol', 'p_HDL', 'p_LDL', 'p_VLDL', 'p_trig')]

## bacterial metabolites: we will define potential bacterial metabolites those metabolite clusters that are associated, at some extent, to microbial community -----
# GC-GC ====
names(gc) <- gsub('\\.', '-', names(gc))
data.model <- cbind(qmp, gc)

set.seed(100) 
index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
train <- data.model[index,] # Create the training data 
test <- data.model[-index,] # Create the test data


lambdas <- 10^seq(2, -3, by = -.1)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lasso.metaG <- lapply(gc, function(x){
  set.seed(100) 
  data.model <- cbind(qmp, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:1273]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:1273]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.metaG.df <- do.call(rbind.data.frame, lasso.metaG)
lasso.metaG.df$fdr <- p.adjust(lasso.metaG.df$pval, method = 'fdr')
lasso.metaG.df.taxa <- lasso.metaG.df[complete.cases(lasso.metaG.df),]
lasso.metaG.df$MGS <- NULL

lasso.metaG.df.taxa2 <- lasso.metaG.df.taxa$MGS
names(lasso.metaG.df.taxa2) <- row.names(lasso.metaG.df.taxa) 
df2<-  list() 
for(i in names(lasso.metaG.df.taxa2)){
  df <- lasso.metaG.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
metaG.mets.relation <- do.call(rbind.data.frame, df2)
# write.table(metaG.mets.relation, 'MGScoeffsLASSOmetaboliteCLUSTERS.txt', sep ='\t', row.names = F)

# GMMs
lasso.GMM <- lapply(gc, function(x){
  set.seed(106) 
  data.model <- cbind(gmm, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:101]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:101]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.GMM.df <- do.call(rbind.data.frame, lasso.GMM)
lasso.GMM.df$fdr <- p.adjust(lasso.GMM.df$pval, method = 'fdr')
lasso.GMM.df.taxa <- lasso.GMM.df[complete.cases(lasso.GMM.df),]
lasso.GMM.df$MGS <- NULL

lasso.GMM.df.taxa2 <- lasso.GMM.df.taxa$MGS
names(lasso.GMM.df.taxa2) <- row.names(lasso.GMM.df.taxa) 
df2<-  list() 
for(i in names(lasso.GMM.df.taxa2)){
  df <- lasso.GMM.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('gmm' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
GMM.mets.relation <- do.call(rbind.data.frame, df2)

# lifestyle.antro
lasso.lifestyle <- lapply(gc, function(x){
  set.seed(102) 
  data.model <- cbind(lifestyle.antro, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:16]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:16]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.lifestyle.df <- do.call(rbind.data.frame, lasso.lifestyle)
lasso.lifestyle.df$fdr <- p.adjust(lasso.lifestyle.df$pval, method = 'fdr')
lasso.lifestyle.df.taxa <- lasso.lifestyle.df[complete.cases(lasso.lifestyle.df),]
lasso.lifestyle.df$MGS <- NULL

lasso.lifestyle.df.taxa2 <- lasso.lifestyle.df.taxa$MGS
names(lasso.lifestyle.df.taxa2) <- row.names(lasso.lifestyle.df.taxa) 
df2<-  list() 
for(i in names(lasso.lifestyle.df.taxa2)){
  df <- lasso.lifestyle.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
lifestyle.mets.relation <- do.call(rbind.data.frame, df2)

# biochem
lasso.biochem <- lapply(gc, function(x){
  set.seed(102) 
  data.model <- cbind(biochem, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:28]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:28]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.biochem.df <- do.call(rbind.data.frame, lasso.biochem)
lasso.biochem.df$fdr <- p.adjust(lasso.biochem.df$pval, method = 'fdr')
lasso.biochem.df.taxa <- lasso.biochem.df[complete.cases(lasso.biochem.df),]
lasso.biochem.df$MGS <- NULL

lasso.biochem.df.taxa2 <- lasso.biochem.df.taxa$MGS
names(lasso.biochem.df.taxa2) <- row.names(lasso.biochem.df.taxa) 
df2<-  list() 
for(i in names(lasso.biochem.df.taxa2)){
  df <- lasso.biochem.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
biochem.mets.relation <- do.call(rbind.data.frame, df2)


lasso.models.GC <- list('QMP' = lasso.metaG.df,
                        'GMM' = lasso.GMM.df,
                        'biochem' = lasso.biochem.df,
                        'lifestyle' = lasso.lifestyle.df
)

lasso.models.GC <- lapply(lasso.models.GC, function(a){
  a <- a[!is.na(a$fdr),]
  a <- a[a$fdr<.1,]
})

lasso.models.GC <- do.call(rbind.data.frame, lasso.models.GC)
lasso.models.GC$data.type <- unlist(strsplit(row.names(lasso.models.GC), '\\.'))[c(T,F)]
lasso.models.GC$metabolite <- unlist(strsplit(row.names(lasso.models.GC), '\\.'))[c(F,T)]
lasso.models.GC$data.type <- factor(lasso.models.GC$data.type, levels = c('QMP', 'GMM', 'biochem', 'lifestyle'))


summary(lasso.models.GC$data.type)
lasso.models.GC$data.type <- factor(lasso.models.GC$data.type, levels=c('biochem', 'QMP', 'GMM', 'lifestyle'))
lasso.models.GC$data.type <- plyr::revalue(lasso.models.GC$data.type, 
                                        c('QMP' = 'Microbiome (QMP)',
                                          'GMM' = 'Gut Metabolic\nModule (GMM)',
                                          'biochem' = 'Bioclinical',
                                          'lifestyle' = 'Lifestyle and\nanthropometrics'))
p1 <- lasso.models.GC %>% 
  # filter(data.type!='lifestyle') %>% 
  ggplot(aes(data.type, var, group = data.type)) +
  geom_violin(alpha = .1, aes(fill=data.type, color=data.type)) +
  ggbeeswarm::geom_quasirandom(alpha = .3, dodge.width = .9, aes(color = data.type)) +
  geom_boxplot(width = .1, alpha = .5, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) + ylim(0,1) +
  labs(x = NULL, y = "Explained variance (/1)", 
       colour = NULL, title = 'Polar metabolites potential origins\nidentified by LASSO modelling')

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'y',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

lasso.models.GC %>% 
  dplyr::mutate('data.type' = factor(data.type)) %>% 
  dplyr::select('data.type') %>% 
  summary()

#stats variance
lasso.models.GC %>% 
  filter(data.type == 'Bioclinical') %>% 
  dplyr::select("var") %>% 
  summary()

lasso.models.GC %>% 
  filter(data.type == 'Microbiome (QMP)') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models.GC %>% 
  filter(data.type == 'Gut Metabolic\nModule (GMM)') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models.GC %>% 
  filter(data.type == 'Lifestyle\nand anthropometrics') %>% 
  dplyr::select('var') %>% 
  summary()

p2 <- lasso.models.GC %>% 
  ggplot(aes(data.type, metabolite, fill = var)) + 
  geom_tile() +
  scale_fill_gradient2(low='red', high='darkred') +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) +
  labs(x=NULL, y='Expl. var. of the metabolite clusters')

cowplot::plot_grid(p1, p2)

## bacterial polar metabolites
# metaG.mets.relation
# met.order <- sort(table(metaG.mets.relation$met))
# metaG.mets.relation$met <- factor(metaG.mets.relation$met, levels(names(met.order)))
# metaG.mets.relation$tax <- sapply(metaG.mets.relation$mgs, function(a){
#   tax[grep(paste0('^', a, '$'), tax$X),]$Name
# })
# tax.met <- tax[tax$X%in%metaG.mets.relation$mgs,]
# metaG.mets.relation$mgs <- factor(metaG.mets.relation$mgs,levels = tax.met$mgs)
# metaG.mets.relation$tax <- factor(metaG.mets.relation$tax, levels = unique(tax.met$Name))
# metaG.mets.relation$tax <-  factor(metaG.mets.relation$tax, levels = rev(levels(metaG.mets.relation$tax)))
# metaG.mets.relation2 <- metaG.mets.relation[metaG.mets.relation$met%in%names(tail(met.order, n=30)),]
# 
# ggplot(metaG.mets.relation2, aes(tax, met, fill ='red')) + 
#   geom_tile() +
#   scale_fill_manual(values = 'darkred') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = .5)) +
#   labs(x = 'Bacteria', y = 'Polar metabolite') + guides(fill=F)

# lipids ALL ====
names(lipids) <- gsub('\\.', '-', names(lipids))
data.model <- cbind(qmp, lipids)

set.seed(100) 
index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
train <- data.model[index,] # Create the training data 
test <- data.model[-index,] # Create the test data


lambdas <- 10^seq(2, -3, by = -.1)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lasso.metaG.lip <- lapply(lipids, function(x){
  set.seed(101) 
  data.model <- cbind(qmp, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:1273]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:1273]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.metaG.lip.df <- do.call(rbind.data.frame, lasso.metaG.lip)
lasso.metaG.lip.df$fdr <- p.adjust(lasso.metaG.lip.df$pval, method = 'fdr')
lasso.metaG.lip.df.taxa <- lasso.metaG.lip.df[complete.cases(lasso.metaG.lip.df),]
lasso.metaG.lip.df$MGS <- NULL

lasso.metaG.lip.df.taxa2 <- lasso.metaG.lip.df.taxa$MGS
names(lasso.metaG.lip.df.taxa2) <- row.names(lasso.metaG.lip.df.taxa) 
df2<-  list() 
for(i in names(lasso.metaG.lip.df.taxa2)){
  df <- lasso.metaG.lip.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('qmp' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
metaG.lipids.relation <- do.call(rbind.data.frame, df2)
# write.table(metaG.lipids.relation, 'MGScoeffsLASSOmetaboliteCLUSTERS.txt', sep ='\t', row.names = F)
# GMMs
lasso.GMM <- lapply(lipids, function(x){
  set.seed(106) 
  data.model <- cbind(gmm, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:101]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:101]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.GMM.df <- do.call(rbind.data.frame, lasso.GMM)
lasso.GMM.df$fdr <- p.adjust(lasso.GMM.df$pval, method = 'fdr')
lasso.GMM.df.taxa <- lasso.GMM.df[complete.cases(lasso.GMM.df),]
lasso.GMM.df$MGS <- NULL

lasso.GMM.df.taxa2 <- lasso.GMM.df.taxa$MGS
names(lasso.GMM.df.taxa2) <- row.names(lasso.GMM.df.taxa) 
df2<-  list() 
for(i in names(lasso.GMM.df.taxa2)){
  df <- lasso.GMM.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('gmm' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
GMM.lipids.relation <- do.call(rbind.data.frame, df2)

# lifestyle.antro
lasso.lifestyle <- lapply(lipids, function(x){
  set.seed(102) 
  data.model <- cbind(lifestyle.antro, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:16]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:16]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.lifestyle.df <- do.call(rbind.data.frame, lasso.lifestyle)
lasso.lifestyle.df$fdr <- p.adjust(lasso.lifestyle.df$pval, method = 'fdr')
lasso.lifestyle.df.taxa <- lasso.lifestyle.df[complete.cases(lasso.lifestyle.df),]
lasso.lifestyle.df$MGS <- NULL

lasso.lifestyle.df.taxa2 <- lasso.lifestyle.df.taxa$MGS
names(lasso.lifestyle.df.taxa2) <- row.names(lasso.lifestyle.df.taxa) 
df2<-  list() 
for(i in names(lasso.lifestyle.df.taxa2)){
  df <- lasso.lifestyle.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
lifestyle.lipids.relation <- do.call(rbind.data.frame, df2)

# biochem
lasso.biochem <- lapply(lipids, function(x){
  set.seed(102) 
  data.model <- cbind(biochem, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:28]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:28]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.biochem.df <- do.call(rbind.data.frame, lasso.biochem)
lasso.biochem.df$fdr <- p.adjust(lasso.biochem.df$pval, method = 'fdr')
lasso.biochem.df.taxa <- lasso.biochem.df[complete.cases(lasso.biochem.df),]
lasso.biochem.df$MGS <- NULL

lasso.biochem.df.taxa2 <- lasso.biochem.df.taxa$MGS
names(lasso.biochem.df.taxa2) <- row.names(lasso.biochem.df.taxa) 
df2<-  list() 
for(i in names(lasso.biochem.df.taxa2)){
  df <- lasso.biochem.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
biochem.lipids.relation <- do.call(rbind.data.frame, df2)


lasso.models.lipids <- list('QMP' = lasso.metaG.lip.df,
                        'GMM' = lasso.GMM.df,
                        'biochem' = lasso.biochem.df,
                        'lifestyle' = lasso.lifestyle.df
)

lasso.models.lipids <- lapply(lasso.models.lipids, function(a){
  a <- a[!is.na(a$fdr),]
  a <- a[a$fdr<.1,]
})

lasso.models.lipids <- do.call(rbind.data.frame, lasso.models.lipids)
lasso.models.lipids$data.type <- unlist(strsplit(row.names(lasso.models.lipids), '\\.'))[c(T,F)]
lasso.models.lipids$metabolite <- unlist(strsplit(row.names(lasso.models.lipids), '\\.'))[c(F,T)]
lasso.models.lipids$data.type <- factor(lasso.models.lipids$data.type, levels = c('QMP', 'GMM', 'biochem', 'lifestyle'))


summary(lasso.models.lipids$data.type)
lasso.models.lipids$data.type <- factor(lasso.models.lipids$data.type, levels=c('biochem', 'QMP', 'GMM', 'lifestyle')) # order by mean variance
lasso.models.lipids$data.type <- plyr::revalue(lasso.models.lipids$data.type, 
                                           c('QMP' = 'Microbiome (QMP)',
                                             'GMM' = 'Gut Metabolic\nModule (GMM)',
                                             'biochem' = 'Bioclinical',
                                             'lifestyle' = 'Lifestyle and\nanthropometrics'))
p3 <- lasso.models.lipids %>% 
  # filter(data.type!='lifestyle') %>% 
  ggplot(aes(data.type, var, group = data.type)) +
    geom_violin(alpha = .1, aes(fill=data.type, color=data.type)) +
    ggbeeswarm::geom_quasirandom(alpha = .3, dodge.width = .9, aes(color = data.type)) +
    geom_boxplot(width = .1, alpha = .5, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) + ylim(0,1) +
  labs(x = NULL, y = "Explained variance (/1)", 
       colour = NULL, title = 'Lipids potential origins\nidentified by LASSO modelling')

ggExtra::ggMarginal(
  p = p3,
  type = 'density',
  margins = 'y',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

lasso.models.lipids %>% 
  dplyr::mutate('data.type' = factor(data.type)) %>% 
  dplyr::select('data.type') %>% 
  summary()

#stats variance
lasso.models.lipids %>% 
  filter(data.type == 'Bioclinical') %>% 
  dplyr::select("var") %>% 
  summary()

lasso.models.lipids %>% 
  filter(data.type == 'Microbiome (QMP)') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models.lipids %>% 
  filter(data.type == 'Gut Metabolic\nModule (GMM)') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models.lipids %>% 
  filter(data.type == 'Lifestyle\nand anthropometrics') %>% 
  dplyr::select('var') %>% 
  summary()

p4 <- lasso.models.lipids %>% 
  ggplot(aes(data.type, metabolite, fill = var)) + 
  geom_tile() +
  scale_fill_gradient2(low='red', high='darkred') +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) +
  labs(x=NULL, y='Expl. var. of the metabolite clusters')

cowplot::plot_grid(p3, p4)

##
# metaG.lipids.relation
# met.order <- sort(table(metaG.lipids.relation$met))
# metaG.lipids.relation$met <- factor(metaG.lipids.relation$met, levels=names(met.order))
# metaG.lipids.relation$tax <- sapply(metaG.lipids.relation$qmp, function(a){
#   tax[grep(paste0('^', a, '$'), tax$X),]$Name
# })
# tax.met <- tax[tax$X%in%metaG.lipids.relation$qmp,]
# metaG.lipids.relation$mgs <- factor(metaG.lipids.relation$qmp,levels = tax.met$X)
# metaG.lipids.relation$tax <- factor(metaG.lipids.relation$tax, levels = unique(tax.met$Name))
# metaG.lipids.relation$tax <-  factor(metaG.lipids.relation$tax, levels = rev(levels(metaG.lipids.relation$tax)))
# metaG.lipids.relation2 <- metaG.lipids.relation[metaG.lipids.relation$met%in%names(tail(met.order, n=30)),]
# 
# ggplot(metaG.lipids.relation2, aes(tax, met, fill ='red')) + 
#   geom_tile() +
#   scale_fill_manual(values = 'darkred') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = .5)) +
#   labs(x = 'Bacteria', y = 'Lipids') + guides(fill=F)
##
# lipids clusters ====
names(lipids.c) <- gsub('\\.', '-', names(lipids.c))
data.model <- cbind(qmp, lipids.c)

set.seed(100) 
index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
train <- data.model[index,] # Create the training data 
test <- data.model[-index,] # Create the test data


lambdas <- 10^seq(2, -3, by = -.1)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lasso.metaG.lip <- lapply(lipids.c, function(x){
  set.seed(101) 
  data.model <- cbind(qmp, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:1273]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:1273]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.metaG.lip.df <- do.call(rbind.data.frame, lasso.metaG.lip)
lasso.metaG.lip.df$fdr <- p.adjust(lasso.metaG.lip.df$pval, method = 'fdr')
lasso.metaG.lip.df.taxa <- lasso.metaG.lip.df[complete.cases(lasso.metaG.lip.df),]
lasso.metaG.lip.df$MGS <- NULL

lasso.metaG.lip.df.taxa2 <- lasso.metaG.lip.df.taxa$MGS
names(lasso.metaG.lip.df.taxa2) <- row.names(lasso.metaG.lip.df.taxa) 
df2<-  list() 
for(i in names(lasso.metaG.lip.df.taxa2)){
  df <- lasso.metaG.lip.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('qmp' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
metaG.lipids.c.relation <- do.call(rbind.data.frame, df2)
# write.table(metaG.lipids.c.relation, 'MGScoeffsLASSOmetaboliteCLUSTERS.txt', sep ='\t', row.names = F)
# GMMs
lasso.GMM <- lapply(lipids.c, function(x){
  set.seed(106) 
  data.model <- cbind(gmm, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:101]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:101]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.GMM.df <- do.call(rbind.data.frame, lasso.GMM)
lasso.GMM.df$fdr <- p.adjust(lasso.GMM.df$pval, method = 'fdr')
lasso.GMM.df.taxa <- lasso.GMM.df[complete.cases(lasso.GMM.df),]
lasso.GMM.df$MGS <- NULL

lasso.GMM.df.taxa2 <- lasso.GMM.df.taxa$MGS
names(lasso.GMM.df.taxa2) <- row.names(lasso.GMM.df.taxa) 
df2<-  list() 
for(i in names(lasso.GMM.df.taxa2)){
  df <- lasso.GMM.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('gmm' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
GMM.lipids.c.relation <- do.call(rbind.data.frame, df2)

# lifestyle.antro
lasso.lifestyle <- lapply(lipids.c, function(x){
  set.seed(102) 
  data.model <- cbind(lifestyle.antro, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:16]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:16]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.lifestyle.df <- do.call(rbind.data.frame, lasso.lifestyle)
lasso.lifestyle.df$fdr <- p.adjust(lasso.lifestyle.df$pval, method = 'fdr')
lasso.lifestyle.df.taxa <- lasso.lifestyle.df[complete.cases(lasso.lifestyle.df),]
lasso.lifestyle.df$MGS <- NULL

lasso.lifestyle.df.taxa2 <- lasso.lifestyle.df.taxa$MGS
names(lasso.lifestyle.df.taxa2) <- row.names(lasso.lifestyle.df.taxa) 
df2<-  list() 
for(i in names(lasso.lifestyle.df.taxa2)){
  df <- lasso.lifestyle.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
lifestyle.lipids.c.relation <- do.call(rbind.data.frame, df2)

# biochem
lasso.biochem <- lapply(lipids.c, function(x){
  set.seed(102) 
  data.model <- cbind(biochem, 'polar.met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:28]), train$polar.met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:28]), train$polar.met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$polar.met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.biochem.df <- do.call(rbind.data.frame, lasso.biochem)
lasso.biochem.df$fdr <- p.adjust(lasso.biochem.df$pval, method = 'fdr')
lasso.biochem.df.taxa <- lasso.biochem.df[complete.cases(lasso.biochem.df),]
lasso.biochem.df$MGS <- NULL

lasso.biochem.df.taxa2 <- lasso.biochem.df.taxa$MGS
names(lasso.biochem.df.taxa2) <- row.names(lasso.biochem.df.taxa) 
df2<-  list() 
for(i in names(lasso.biochem.df.taxa2)){
  df <- lasso.biochem.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
biochem.lipids.c.relation <- do.call(rbind.data.frame, df2)


lasso.models.lipids.c <- list('QMP' = lasso.metaG.lip.df,
                            'GMM' = lasso.GMM.df,
                            'biochem' = lasso.biochem.df,
                            'lifestyle' = lasso.lifestyle.df
)

lasso.models.lipids.c <- lapply(lasso.models.lipids.c, function(a){
  a <- a[!is.na(a$fdr),]
  a <- a[a$fdr<.1,]
})

lasso.models.lipids.c <- do.call(rbind.data.frame, lasso.models.lipids.c)
lasso.models.lipids.c$data.type <- unlist(strsplit(row.names(lasso.models.lipids.c), '\\.'))[c(T,F)]
lasso.models.lipids.c$metabolite <- unlist(strsplit(row.names(lasso.models.lipids.c), '\\.'))[c(F,T)]
lasso.models.lipids.c$data.type <- factor(lasso.models.lipids.c$data.type, levels = c('QMP', 'GMM', 'biochem', 'lifestyle'))


summary(lasso.models.lipids.c$data.type)
lasso.models.lipids.c$data.type <- factor(lasso.models.lipids.c$data.type, levels=c('biochem', 'QMP', 'GMM', 'lifestyle')) # order by mean variance
lasso.models.lipids.c$data.type <- plyr::revalue(lasso.models.lipids.c$data.type, 
                                               c('QMP' = 'Microbiome (QMP)',
                                                 'GMM' = 'Gut Metabolic\nModule (GMM)',
                                                 'biochem' = 'Bioclinical',
                                                 'lifestyle' = 'Lifestyle and\nanthropometrics'))
p5 <- lasso.models.lipids.c %>% 
  # filter(data.type!='lifestyle') %>% 
  ggplot(aes(data.type, var, group = data.type)) +
  geom_violin(alpha = .1, aes(fill=data.type, color=data.type)) +
  ggbeeswarm::geom_quasirandom(alpha = .3, dodge.width = .9, aes(color = data.type)) +
  geom_boxplot(width = .1, alpha = .5, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) +  ylim(0,1) +
  labs(x = NULL, y = "Explained variance (/1)", 
       colour = NULL, title = 'Lipid clusters potential origins\nidentified by LASSO modelling')

ggExtra::ggMarginal(
  p = p5,
  type = 'density',
  margins = 'y',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

lasso.models.lipids.c %>% 
  dplyr::mutate('data.type' = factor(data.type)) %>% 
  dplyr::select('data.type') %>% 
  summary()

#stats variance
lasso.models.lipids.c %>% 
  filter(data.type == 'Bioclinical') %>% 
  dplyr::select("var") %>% 
  summary()

lasso.models.lipids.c %>% 
  filter(data.type == 'Microbiome (QMP)') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models.lipids.c %>% 
  filter(data.type == 'Gut Metabolic\nModule (GMM)') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models.lipids.c %>% 
  filter(data.type == 'Lifestyle\nand anthropometrics') %>% 
  dplyr::select('var') %>% 
  summary()

p6 <- lasso.models.lipids.c %>% 
  ggplot(aes(data.type, metabolite, fill = var)) + 
  geom_tile() +
  scale_fill_gradient2(low='red', high='darkred') +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) +
  labs(x=NULL, y='Expl. var. of the metabolite clusters')

cowplot::plot_grid(p5, p6)

##
# metaG.lipids.c.relation
# met.order <- sort(table(metaG.lipids.c.relation$met))
# metaG.lipids.c.relation$met <- factor(metaG.lipids.c.relation$met, levels=names(met.order))
# metaG.lipids.c.relation$tax <- sapply(metaG.lipids.c.relation$qmp, function(a){
#   tax[grep(paste0('^', a, '$'), tax$X),]$Name
# })
# tax.met <- tax[tax$X%in%metaG.lipids.c.relation$qmp,]
# metaG.lipids.c.relation$mgs <- factor(metaG.lipids.c.relation$qmp,levels = tax.met$X)
# metaG.lipids.c.relation$tax <- factor(metaG.lipids.c.relation$tax, levels = unique(tax.met$Name))
# metaG.lipids.c.relation$tax <-  factor(metaG.lipids.c.relation$tax, levels = rev(levels(metaG.lipids.c.relation$tax)))
# metaG.lipids.c.relation2 <- metaG.lipids.c.relation[metaG.lipids.c.relation$met%in%names(tail(met.order, n=30)),]
# 
# ggplot(metaG.lipids.c.relation, aes(tax, met, fill ='red')) + 
#   geom_tile() +
#   scale_fill_manual(values = 'darkred') +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = .5)) +
#   labs(x = 'Bacteria', y = 'Lipid clusters') + guides(fill=F)
##
# plot all ----
cowplot::plot_grid(p1, p3, p5, ncol = 3)

# save run
save.image('metabolitOrigins.RData')


rev(sort(table(metaG.mets.relation$tax)))
length(sort(table(metaG.mets.relation$tax)))
rev(sort(table(metaG.lipids.relation$tax)))
length(sort(table(metaG.lipids.relation$tax)))
rev(sort(table(metaG.lipids.c.relation$tax)))
length(sort(table(metaG.lipids.c.relation$tax)))

gmm.annot <- read.table('GMMs.v1.07.names', header = F, sep='\t')
GMM.mets.relation$gmm.annot <- sapply(GMM.mets.relation$gmm, function(a){
  gmm.annot[grep(paste0('^', a, '$'), gmm.annot$V1),]$V2
  })
rev(sort(table(GMM.mets.relation$gmm.annot)))
length(sort(table(GMM.mets.relation$gmm.annot)))
GMM.lipids.relation$gmm.annot <- sapply(GMM.lipids.relation$gmm, function(a){
  gmm.annot[grep(paste0('^', a, '$'), gmm.annot$V1),]$V2
})
rev(sort(table(GMM.lipids.relation$gmm.annot)))
length(sort(table(GMM.lipids.relation$gmm.annot)))
GMM.lipids.c.relation$gmm.annot <- sapply(GMM.lipids.c.relation$gmm, function(a){
  gmm.annot[grep(paste0('^', a, '$'), gmm.annot$V1),]$V2
})
rev(sort(table(GMM.lipids.c.relation$gmm.annot)))
length(sort(table(GMM.lipids.c.relation$gmm.annot)))


