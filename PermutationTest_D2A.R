#Permutation Test functions
library(limma)
eBayesArrays <- function(arrayDF, factors, correctMedian = FALSE, nonDataCols = c('Accession','Description','HGNC','Known.GlcNAc',
                                                                                  'PhosphoSitePlus','Wang')){
  if(correctMedian){
    ExpMatrix <- log2(arrayDF[,!(colnames(arrayDF) %in% nonDataCols)])
    for(f in unique(factors[factors!='control'])){
      ExpMatrix[,factors == f] <- ExpMatrix[,factors == f] - 
        median(unlist(ExpMatrix[,factors == f]),na.rm = T) + 
        median(unlist(ExpMatrix[,factors == 'control']), na.rm = T)
    }
  }else{
    ExpMatrix <- log2(arrayDF[,!(colnames(arrayDF) %in% nonDataCols)])
  }
  factors <- factor(factors)
  designArray <- model.matrix(~ 0 + factors)
  colnames(designArray) <- sub('factors(.*$)', '\\1', colnames(designArray))
  fit <- lmFit(ExpMatrix, design = designArray, method = 'robust', maxit = 400)
  contrast.matrix <- makeContrasts(wt-control, D2A-control, wt-D2A, levels=designArray)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  eBayesOut <- eBayes(fit2, robust = TRUE)
  eBayesOut
}
permuteBackgroundwtcont <- function(p.value,tval,inputData,arrayOrder,toCompare, alpha = 0.1, oneway = TRUE,
                                    N = 100){
  #generates a random permutation distribution of eBayes-derived t tests. Only does two way
  #comparison at a time. Need to give single-column of p value and tval. oneway assumes only
  #looking for greater if used. CANNOT be used with array weights.
  if(oneway){
    workingData <- inputData[p.value > (alpha*2) | tval < 0,arrayOrder %in% toCompare]
  }else{
    workingData <- inputData[p.value > alpha,
                             arrayOrder %in% toCompare]
  }
  newOrder <- arrayOrder[arrayOrder %in% toCompare]
  back.tvals <- list()
  for(i in 1:N){
    workingOrder <- sample(newOrder)
    workingOrder <- factor(workingOrder)
    work.design <- model.matrix(~ 0 + workingOrder )
    colnames(work.design) <- sub('workingOrder(.*$)', '\\1',
                                 colnames(work.design))
    currFit <- lmFit(workingData, design = work.design)
    contrast.matrix <-makeContrasts(wt - control, levels=work.design)
    
    currFit2 <- contrasts.fit(currFit, contrast.matrix)
    currFit2 <- eBayes(currFit2, robust = TRUE)
    back.tvals[[i]] <- currFit2$t
  }
  back <- NULL
  back <- unlist(back.tvals)
  pvalout <- NULL
  if(oneway){
    pvalout <- sapply(tval, function(i,back){sum(back>i)/length(back)},back=back)
    # for(i in tval){
    #   pvalout <- c(pvalout,sum(back>i)/length(back))
    # }
  }else{
    pvalout <- sapply(tval, function(i,back){
      if(i>=median(back)){
        2*sum(back>i)/length(back)
      }else{
        2*sum(back<i)/length(back)
      }
    }, back=back)
  }  
  pvalout <- data.frame('p.value' = pvalout)
  rownames(pvalout) <- rownames(inputData)
  output <- list('BackgroundDist' = back, 'pvals' = pvalout)
  output
}

permuteBackgroundD2Acont <- function(p.value,tval,inputData,arrayOrder,toCompare, alpha = 0.1, oneway = TRUE,
                                     N = 100){
  #generates a random permutation distribution of eBayes-derived t tests. Only does two way
  #comparison at a time. Need to give single-column of p value and tval. oneway assumes only
  #looking for greater if used. CANNOT be used with array weights.
  if(oneway){
    workingData <- inputData[p.value > (alpha*2) | tval < 0,arrayOrder %in% toCompare]
  }else{
    workingData <- inputData[p.value > alpha,
                             arrayOrder %in% toCompare]
  }
  newOrder <- arrayOrder[arrayOrder %in% toCompare]
  back.tvals <- list()
  for(i in 1:N){
    workingOrder <- sample(newOrder)
    workingOrder <- factor(workingOrder)
    work.design <- model.matrix(~ 0 + workingOrder )
    colnames(work.design) <- sub('workingOrder(.*$)', '\\1',
                                 colnames(work.design))
    currFit <- lmFit(workingData, design = work.design)
    contrast.matrix <-makeContrasts(D2A - control, levels=work.design)
    
    currFit2 <- contrasts.fit(currFit, contrast.matrix)
    currFit2 <- eBayes(currFit2, robust = TRUE)
    back.tvals[[i]] <- currFit2$t
  }
  back <- NULL
  back <- unlist(back.tvals)
  pvalout <- NULL
  if(oneway){
    pvalout <- sapply(tval, function(i,back){sum(back>i)/length(back)},back=back)
  }else{
    pvalout <- sapply(tval, function(i,back){
      if(i>=median(back)){
        2*sum(back>i)/length(back)
      }else{
        2*sum(back<i)/length(back)
      }
    }, back=back)
  }  
  pvalout <- data.frame('p.value' = pvalout)
  rownames(pvalout) <- rownames(inputData)
  output <- list('BackgroundDist' = back, 'pvals' = pvalout)
  output
}
permuteBackgroundwtvD2A <- function(p.value,tval,inputData,arrayOrder,toCompare, alpha = 0.1, oneway = FALSE,
                                     N = 100){
  #generates a random permutation distribution of eBayes-derived t tests. Only does two way
  #comparison at a time. Need to give single-column of p value and tval. oneway assumes only
  #looking for greater if used. CANNOT be used with array weights.
  if(oneway){
    workingData <- inputData[p.value > (alpha*2) | tval < 0,arrayOrder %in% toCompare]
  }else{
    workingData <- inputData[p.value > alpha,
                             arrayOrder %in% toCompare]
  }
  newOrder <- arrayOrder[arrayOrder %in% toCompare]
  back.tvals <- list()
  for(i in 1:N){
    workingOrder <- sample(newOrder)
    workingOrder <- factor(workingOrder)
    work.design <- model.matrix(~ 0 + workingOrder )
    colnames(work.design) <- sub('workingOrder(.*$)', '\\1',
                                 colnames(work.design))
    currFit <- lmFit(workingData, design = work.design)
    contrast.matrix <-makeContrasts(wt - D2A, levels=work.design)
    
    currFit2 <- contrasts.fit(currFit, contrast.matrix)
    currFit2 <- eBayes(currFit2, robust = TRUE)
    back.tvals[[i]] <- currFit2$t
  }
  back <- NULL
  back <- unlist(back.tvals)
  pvalout <- NULL
  if(oneway){
    pvalout <- sapply(tval, function(i,back){sum(back>i)/length(back)},back=back)
  }else{
    pvalout <- sapply(tval, function(i,back){
      if(i>=median(back)){
        2*sum(back>i)/length(back)
      }else{
        2*sum(back<i)/length(back)
      }
    }, back=back)
  }  
  pvalout <- data.frame('p.value' = pvalout)
  rownames(pvalout) <- rownames(inputData)
  output <- list('BackgroundDist' = back, 'pvals' = pvalout)
  output
}
#Median Correct
medianCorrect <- function(arrayDF, controlCol, listOGTCol, nonDataCols = c('Accession','Description','HGNC','Known.GlcNAc',
                                                                             'PhosphoSitePlus','Wang')){
  #works on data columns as listed and adjusts medians in log2, then converts back.
  log2arrayDF <- arrayDF
  log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)] <- log2(log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)])
  medianControl <- median(unlist(log2arrayDF[,controlCol]), na.rm = TRUE)
  OGTMedians <- list()
  for(i in names(listOGTCol)){
    OGTMedians[[i]] <- median(unlist(log2arrayDF[,listOGTCol[[i]]]), na.rm = TRUE)
    for(j in listOGTCol[[i]]){
      log2arrayDF[,j] <- log2arrayDF[,j] - OGTMedians[[i]] + medianControl
    }
    arrayDF[,listOGTCol[[i]]] <- 2^log2arrayDF[,listOGTCol[[i]]]
  }
  arrayDF
}

medianSignal <- function(arrayDF, listColCategory, nonDataCols = c('Accession','Description','HGNC','Known.GlcNAc',
                                                                   'PhosphoSitePlus','Wang')){
  #gives the by-category log2 median (converted back to non-log2) as columns. Uses names from list
  log2arrayDF <- arrayDF
  log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)] <- log2(log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)])
  CatMedians <- list()
  for(i in names(listColCategory)){
    CatMedians[[paste(i,'Median')]] <- 2^apply(log2arrayDF[,listColCategory[[i]]], 1,median, na.rm = T)
  }
  do.call(cbind.data.frame,CatMedians)
}

meanSignal <- function(arrayDF, listColCategory, nonDataCols = c('Accession','Description','HGNC','Known.GlcNAc',
                                                                   'PhosphoSitePlus','Wang')){
  #gives the by-category log2 median (converted back to non-log2) as columns. Uses names from list
  log2arrayDF <- arrayDF
  log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)] <- log2(log2arrayDF[,!(colnames(log2arrayDF) %in% nonDataCols)])
  CatMeans <- list()
  for(i in names(listColCategory)){
    CatMeans[[paste(i,'Mean')]] <- 2^apply(log2arrayDF[,listColCategory[[i]]], 1,mean, na.rm = T)
  }
  do.call(cbind.data.frame,CatMeans)
}
