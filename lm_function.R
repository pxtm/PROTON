## LM function ----
# df1: metadata
# df2: covariates
# variables: variables to assess
# control_by: variables to control for, coming from metadata

lm.associations <- function(df1, df2, variables, control_by) {
  library(effsize)
  library(caret)
  library(dplyr)
  progression.model <- data.frame(df1, df2)
  if(missing(control_by)){
    print('You are not controlling for any variable')
  } else {
    print(paste0('You are controlling for:', paste(control_by, collapse = '+')))
  }
  lm.outputs <- list()
  for (i in variables){
    print(i)
    result.model <- list()
    for (x in seq(from = ncol(df1)+1, to = ncol(progression.model))) {
      model <-
        data.frame('mgs' = progression.model[, x], progression.model[, 1:ncol(df1)])
      if(missing(control_by)){
        f <- paste('mgs~',i)
      } else {
        f <- paste('mgs~', i, '+', paste(control_by, sep='+'))
      }
      mod.a <-
        lm(
          f,
          data = model,
          na.action = na.omit
        )
      summ.mod.a <- summary(mod.a)
      p.val <-
        summ.mod.a$coefficients[2, 4]
      cl.delta <- cliff.delta(model$mgs, model[,i])
      delta <- cl.delta$estimate
      delta.low <- cl.delta$conf.int[1] #lower
      delta.up <- cl.delta$conf.int[2] #upper
      result <-
        c(names(progression.model)[x],
          p.val,
          delta,
          delta.low,
          delta.up)
      result.model[[x]] <- result
    }
    model.lm <- do.call(rbind.data.frame, result.model)
    names(model.lm) <-
      c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
    model.lm$fdr <- p.adjust(model.lm$pval, method = 'fdr')
    lm.outputs[[i]] <- model.lm
  }
  return(lm.outputs)
}
