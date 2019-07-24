#Functions for reading in microarray data
library(limma)
ReadInArrays <- function(targets){
  #reads in the arrays listed in targets by their batch, corrects their backgrounds, and outputs
  #a list by batch of array data in maimage format
  Raw.Data <- list()
  spottypes <- readSpotTypes()
  for(i in unique(targets$Batch)){
    Raw.Data[[i]] <- read.maimages(targets$FileName[targets$Batch==i],
                                   source='genepix.median', 
                                   annotation = c('Block', 'Row', 'Column', 'ID', 'Name',
                                                  'Description'),
                                   wt.fun = wtflags(weight=0,cutoff=-50))
    Raw.Data[[i]]$genes$Status <- controlStatus(spottypes, Raw.Data[[i]])
  }
  Data <- list()
  for(i in unique(targets$Batch)){
    Data[[i]] <- backgroundCorrect(Raw.Data[[i]], method = 'normexp')
  }
  Data
}
modifyAccession <- function(genesFrame){
  Names <- sub('^Hs~Ref:([[:alnum:]]+)-[1-2]~CAT.*$','\\1',genesFrame$Name)
  Names <- sub('^Hs~IVGN:(PM_[0-9]{4})~Ext:([[:alnum:]_\\)\\([:blank:]-]+)~.*$','Invitrogen: \\1, \\2', Names)
  Names <- sub('^Hs~IVGN:([[:alnum:]_-]+)~(CAT_[[:alnum:]-]+)~.*$','Invitrogen: \\1, \\2', Names)
  Names <- sub('^Hs~Abcam:([a-zA-Z-]+)~(CAT_[[:alnum:]]+)~.*$','Abcam: \\1, \\2', Names)
  Names <- sub('^Hs~[[:alpha:]/]+:([[:alnum:]_\\.]*)~[[:alpha:]]+.*$','\\1', Names)
  Names <- sub('^([[:alnum:]-]+)~N/A$','\\1', Names)
  Names <- sub('^([[:alnum:]]+)~RFU:[0-9\\.]+$','\\1', Names)
  Names <- sub('^Invitrogen: (NM_[0-9]+), CAT_.*$', '\\1', Names)
  Names
}
#Changes Controls to contain info on which control in the description line
modifyDescription <- function(genesFrame){
  #generates a modified Description list with 'Control' Swapped for 
  #'Control: (Name Of Control)'
  for(i in which(grepl('Control',genesFrame$Description))){
    genesFrame$Description[i] <- paste('Control: ',
                                       sub('^(.*)~.*$','\\1',
                                           genesFrame$Name[i]),
                                       sep = ''
    )
    
  }
  genesFrame$Description
}
genControlFiles <- function(Data){
  #generates a list of data files for use with addV5conc to add concentrations of V5-biotin control
  #proteins
  ControlFiles <- list()
  for(i in names(Data)){
    ControlFiles[[i]]<-paste('./lotspecific/',as.character(i),'/',as.character(i),'_Control.txt',sep='')
    ControlFiles[[i]]<-read.table(ControlFiles[[i]], header = T, sep='\t', row.names = NULL)
    colnames(ControlFiles[[i]]) <-
      c(colnames(ControlFiles[[i]])[2:length(colnames(ControlFiles[[i]]))], 'X')
    ControlFiles[[i]] <- ControlFiles[[i]][,colnames(ControlFiles[[i]])!='X']#manipulations needed because there is a trailing \t in the data files...
  }
  ControlFiles
}

addV5Conc <- function(geneList, ContFile){
  #uses the appropriate ContFile to add the V5-biotin control spot concentrations to 
  #the geneList. All Concs are in nM
  V5prots <- grep('V5', ContFile$Control.Proteins, value=T)
  geneList[,'V5conc'] <- 0
  for(i in V5prots){
    geneList[grepl(i,geneList$Name),'V5conc'] <-
      ContFile[ContFile$Control.Proteins == i, 'Concentration..nM.']
  }
  geneList
}

#below requires ArrayContent from LitComparison, WRITE THIS UP NICELY AND SAVE DATA FILE
addGeneSymbol <- function(uArray,ArrayContent,locInfo = FALSE){
  uArray$genes[,'Symbol'] <- ''
  if(locInfo){
    uArray$genes[,c('OGTCompartment','NonOGTCompartment')] <- FALSE
    uArray$genes[,c('WeakLocInfo')] <- TRUE
  }
  for(i in 1:length(uArray$genes$Description)){
    ID <- sub('.*~(.*$)','\\1',uArray$genes$ID[i])
    whichGene <- which(grepl(ID, ArrayContent$Spot.ID))
    if(length(whichGene) < 1){#so, if it isn't found in ArrayContent
      if(i > 1){#and its not the first one
        if(uArray$genes$Description[i-1]==uArray$genes$Description[i]){#if it's not the first of its kind
          uArray$genes$Symbol[i] <- uArray$genes$Symbol[i-1] #then make it match
          if(locInfo){
            uArray$genes[i,c('OGTCompartment','NonOGTCompartment','WeakLocInfo')] <- uArray$genes[
              i-1, c('OGTCompartment','NonOGTCompartment','WeakLocInfo')]
          } 
          next()
        }else{#otherwise if it is the first of its kind, who cares
          next()}#so skip it
      }else{# if its the first one, skip it
        next()
      }}
    if(ArrayContent$SymbolOfficial2[whichGene]!=''){
      uArray$genes$Symbol[i] <- paste(ArrayContent$SymbolOfficial[whichGene],
                                      ArrayContent$SymbolOfficial[whichGene],
                                      sep = ' ')
    }else{
      uArray$genes$Symbol[i] <- ArrayContent$SymbolOfficial[whichGene]
    }
    if(locInfo){
      uArray$genes[i,c('OGTCompartment','NonOGTCompartment','WeakLocInfo')] <- ArrayContent[
        whichGene,c('SeesOGT','DoesntSeeOGT','WeakLocInfo')]
    } 
  }
  uArray
}
