#shared functions for the analysis of data by either the Calling Score or permutation-test based methods
library(dplyr)
Collapse2Means <- function(uArray, Symbol = FALSE, bySymbol = FALSE, otherCols = NULL){
  #Takes in uArray (with elements $R, $genes$Description) and 
  #collapses to Description (or gene Symbol if the bySymbol argument is TRUE),
  #with first modifying the names to have the controls be sensible.
  #symbol argument is only to be used if HGNC symbols have been added prior to this to the gene
  #data frame of the microarray object
  uArray$Description <- modifyDescription(uArray$genes)
  uArray
  if(Symbol){
    uArray$Symbol <- uArray$genes$Symbol
    if(bySymbol){
      colnamesgenes <- c('Symbol',colnames(uArray$R))
      genesR <- cbind(data.frame(uArray$Symbol),data.frame(uArray$R),data.frame(uArray$genes[,otherCols]))
      colnames(genesR) <- c(colnamesgenes,otherCols)
      byGene <- group_by(genesR, Symbol) %>% summarise_all(funs(mean))
    }else{
      colnamesgenes <- c('Description','Accession','Symbol',colnames(uArray$R))
      genesR <- cbind(data.frame(uArray$Description),data.frame(uArray$genes$Name),data.frame(uArray$Symbol),data.frame(uArray$R),
                      data.frame(uArray$genes[,otherCols]))
      colnames(genesR) <- c(colnamesgenes,otherCols)
      byGene <- group_by(genesR, Accession, Symbol) %>% summarise_all(funs(mean))
      x <- group_by(genesR, Accession, Symbol) %>% summarise('Description' = first(Description))
      byGene$Description <- x$Description
    }
  }else{
    genesR <- cbind(data.frame(uArray$Description),data.frame(uArray$genes$Name),
                    data.frame(uArray$R),data.frame(uArray$genes[,otherCols]))
    colnamesgenes <- c('Description', 'Accession', colnames(uArray$R))
    colnames(genesR) <- c(colnamesgenes,otherCols)
    byGene <- group_by(genesR, Accession) %>% summarise_all(funs(mean)) #this is problematic and should be
    #fixed, but is fine for now...
    #note: gives warnings because you can't take mean of a character vector,
    #but still fastest implementation
    x <- group_by(genesR, Accession) %>% summarise('Description' = first(Description))
    byGene$Description <- x$Description
  }
  data.frame(byGene)
}

addHGNC <- function(Data, masterList, masterListDescCol = "Protein.Header", descriptionCol = "Description",
                    knownOGlcNAc = TRUE, indivDB = FALSE){
  #reads through Data "Description" and uses the dataset listed in masterList to add a column of 
  #HGNC symbols. If "knownOGlcNAc" is TRUE, a second column listing if this is in one of the mined databases
  #(in TRUE/FALSE format) is added. if "indivDB" is TRUE, then 2 columns will be added breaking out by 
  #database from with data is obtained
  Data$HGNC <- ''
  if(knownOGlcNAc){
    Data$Known.GlcNAc <- FALSE
  }
  if(indivDB){
    Data$PhosphoSitePlus <- FALSE
    Data$Wang <- FALSE
  }
  Descs <- unique(Data[,descriptionCol])
  for(i in Descs){
    currData <- masterList[(masterList[,masterListDescCol] == i) &
                             !is.na(masterList[,masterListDescCol]),]

    Data$HGNC[
      Data[,descriptionCol] == i] <- currData$HGNC[1]
    if(knownOGlcNAc){
      Data$Known.GlcNAc[
        Data[,descriptionCol] == i] <- currData$KnownGlyco[1]
    }
    if(indivDB){
      Data$PhosphoSitePlus[
        Data[,descriptionCol] == i] <- currData$PhosphoSitePlus[1]
      Data$Wang[
        Data[,descriptionCol] == i] <- currData$Wang[1]
    }
  }
  Data
}