#this file contains functions important to the reading in of data and normalizing of arrays
library(MASS)
library(dplyr)
library(reshape2)
#adds concentration for biotinylated V5 control spots 

#normalizes for staining using the biotinylated V5 controls
normStain <- function(dataSet){
  #normalizes data based upon robust linear fit of 
  fitsum <- list()
  for(i in unique(dataSet$genes$Block)){
    fitsum[[i]] <- list()
    for(j in 1:dim(dataSet$R)[2]){
      #now have isolated single block of single array
      x <- dataSet$genes$V5conc[(dataSet$genes$V5conc>0) & (dataSet$genes$Block == i)]
      y <- dataSet$R[(dataSet$genes$V5conc>0) & (dataSet$genes$Block == i),j]
      fit <- rlm(y ~ x, maxit = 200)
      fitsum[[i]][[j]] <- summary(fit)
      #this will allow extraction of slope
    }
  }
  for(i in unique(dataSet$genes$Block)){
    for(j in 1:dim(dataSet$R)[2]){
      tempdata <- dataSet$R[dataSet$genes$Block==i,j]
      if(!is.na(fitsum[[i]][[j]]$coef[1,2]) & (fitsum[[i]][[j]]$coef[2]>0)){
        #If the data is good for a given block, use that fit
        blockintercept <- fitsum[[i]][[j]]$coef[1]
        blockslope <- fitsum[[i]][[j]]$coef[2]
      }else{#if it isn't, use an average of the fits on the array
        blockintercepts <- NULL
        blockslopes <- NULL
        for(k in 1:i){
          blockintercepts <- c(blockintercepts, fitsum[[k]][[j]]$coef[1])
          blockslopes <- c(blockslopes, fitsum[[k]][[j]]$coef[2])
        }
        blockintercept <- mean(blockintercepts)
        blockslope <- mean(blockslopes)
        if(blockslope < 0){
          stop(paste('array: ', as.character(colnames(dataSet$R)[j]), ' block: ', as.character(i), '\n'))
        }
      }
      #ditto slopes
      tempdata <-  tempdata / blockslope
      dataSet$R[dataSet$genes$Block==i,j] <- tempdata
      if(sum(is.na(tempdata))>0){print(c(i,j))}
    }
  }
  dataSet
}
#normalizes for OGT activity based upon relative activity on estrogen receptor alpha
NormERalpha <- function(DataSet,colFactors,relativeActivity,flagInactive = FALSE){
  #set up data as desired into data frame ordered by tidy data principles
  DataSet.df <- data.frame(DataSet$R)
  DataSet.df$Block <- DataSet$genes$Block
  DataSet.df$Description <- DataSet$genes$Description
  DataSet.df$Accession <- DataSet$genes$Name
  DS.tidy <- melt(DataSet.df, id.vars = c('Description','Accession','Block'))
  colnames(DS.tidy)[colnames(DS.tidy)=='variable'] <- 'Array'
  #figure out which arrays are controls, etc.
  contArrays <- make.names(colnames(DataSet$R)[colFactors==0])
  #list OGT-treated arrays by type of array
  OGTArrays <- list()
  for(i in unique(colFactors[colFactors>0])){
    OGTArrays[[i]] <- make.names(colnames(DataSet$R)[colFactors==i])
  }
  #find control median of ERalpha
  contERmedian <- median(DS.tidy$value[(DS.tidy$Accession == 'ERa') &
                                         (DS.tidy$Array %in% contArrays)], na.rm = TRUE)
  #find block ERalpha medians
  blockERmedians <- group_by(DS.tidy[(DS.tidy$Accession == 'ERa') &
                                       !(DS.tidy$Array %in% contArrays),], Array, Block) %>% 
    summarise('ERaMedian' = median(value))
  #find block Activities against ERalpha
  blockERmedians$ERActivity <- blockERmedians$ERaMedian - contERmedian
  #find gene control medians
  geneContMedians <- group_by(DS.tidy[(DS.tidy$Array %in% contArrays),], Accession) %>% 
    summarise('ControlMedian' = median(value))
  #and wt median
  wtERmedian <- median(DS.tidy$value[(DS.tidy$Accession == 'ERa') &
                                       (DS.tidy$Array %in% OGTArrays[[1]])], na.rm = TRUE)
  #now "average" activity for wild-type
  avgActivity <- wtERmedian - contERmedian
  #now go OGT type by OGT type
  for(i in 1:length(OGTArrays)){
    #find by-arrayxBlock accession numbers where there is a sign of activity
    #first join in the control medians, then get gene activities
    geneActivities <- group_by(DS.tidy[(DS.tidy$Array %in% OGTArrays[[i]]),], Array, Block, Accession) %>% 
      summarise('GeneMedian' = median(value)) %>%
      left_join(geneContMedians, by = 'Accession') %>% mutate('Activity' = GeneMedian-ControlMedian) %>%
      left_join(blockERmedians, by = c('Array','Block')) %>%
      mutate('CorrectedActivity' = Activity*avgActivity/ERActivity*relativeActivity[i])
    geneActivities$Activity[geneActivities$Activity<0] <- 0
    geneActivities$CorrectedActivity[geneActivities$CorrectedActivity<0] <- 0
    #now cycle through all the arrays and adjust the single spots, simply subtracting by the Activity
    #and adding the Corrected Activity
    for(j in OGTArrays[[i]]){
      jj <- quo(j)
      toCorrect <- filter(geneActivities, Array == j) %>%
        right_join(DataSet.df[,c('Accession','Block',j)], by = c('Accession','Block')) %>%
        mutate(ValueCorrected = .data[[!!jj]] - .data$Activity + .data$CorrectedActivity)
      DataSet.df[,colnames(DataSet.df)==j] <- toCorrect$ValueCorrected
    }
  }
  DataSet$R <- DataSet.df[,make.names(colnames(DataSet$R))]
  return(DataSet)
}



Collapse2Means <- function(uArray){
  require(dplyr)
  #Takes in uArray (with elements $R, $genes$Description) and collapses to Description, with first modifying the names to have the controls be sensible.
  uArray$Description <- modifyDescription(uArray$genes)
  genesR <- cbind(data.frame(uArray$Description),data.frame(uArray$R))
  colnamesgenes <- c('Description', colnames(uArray$R))
  colnames(genesR) <- colnamesgenes
  
  byGene <- group_by(genesR, Description) %>% summarise_all(funs(mean))
  data.frame(byGene)
}