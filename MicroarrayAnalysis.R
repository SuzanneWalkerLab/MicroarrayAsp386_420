#Analysis of OGT WT-treated Microarrays and OGT D2A-treated microarrays
#both relative to each other and relative to no-enzyme control


#read in important functions
source('NormalizationFunctions.R')
source('SharedFunctions.R')
source('ReadInFunctions.R')
source('PermutationTest_D2A.R')
#read in essential libraries
library(limma)
library(qvalue)
library(ggplot2)

#read in data
targets <- readTargets('Targets.txt', row.names = 'ArrayID')

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

#correct background on all arrays using normal-exponential model
Data <- list()
for(i in unique(targets$Batch)){
  Data[[i]] <- backgroundCorrect(Raw.Data[[i]], method = 'normexp')
}

#find control files and use them to add concentrations for V5-tagged biotinylated controls
ControlFiles <- list()
for(i in names(Data)){
  ControlFiles[[i]]<-paste('./',as.character(i),'/',as.character(i),'_control.txt',sep='')
  ControlFiles[[i]]<-read.table(ControlFiles[[i]], header = T, sep='\t', row.names = NULL)
  colnames(ControlFiles[[i]]) <-
    c(colnames(ControlFiles[[i]])[2:length(colnames(ControlFiles[[i]]))], 'X')
  ControlFiles[[i]] <- ControlFiles[[i]][,colnames(ControlFiles[[i]])!='X']#manipulations needed because there is a trailing \t in the data files...
}

for(i in names(Data)){
  Data[[i]]$genes <- addV5Conc(Data[[i]]$genes,ControlFiles[[i]])
}

#pull together GlcNAz data into a single list rather than separating by batches
Data.GlcNAz <- Data$HA20358
Data.GlcNAz$R <- cbind(Data.GlcNAz$R,Data$HA20446$R)
Data.GlcNAz$G <- cbind(Data.GlcNAz$G,Data$HA20446$G)
Data.GlcNAz$weights <- cbind(Data.GlcNAz$weights,Data$HA20446$weights)
Data.GlcNAz$targets <- c(Data.GlcNAz$targets,Data$HA20446$targets)
Data.GlcNAz$genes$Description <- modifyDescription(Data.GlcNAz$genes)
Data.GlcNAz$genes$Name <- modifyAccession(Data.GlcNAz$genes)
#normalize for staining
Data.GlcNAz.norm <- normStain(Data.GlcNAz)


colFactors <- c(0,1,2,0,1,2,0,1,2)



Data.GlcNAz.ERa.Norm <- NormERalpha(Data.GlcNAz.norm,colFactors,relativeActivity = c(1,1))
MeansArray <- Collapse2Means(Data.GlcNAz.ERa.Norm)
ArrayContent <- read.csv('ArrayProteinContent.csv', #necessary to appropriately label proteins, file attached; 
                         #derived from Thermo Fisher information and manually annotated, as in Levine et al. JACS 2018
                         stringsAsFactors = FALSE)
MeansArray <- addHGNC(MeansArray, ArrayContent, indivDB = TRUE)
MeansArray.Prot <- MeansArray[!grepl('LINC',MeansArray$HGNC) & !grepl('-AS',MeansArray$HGNC) &
                                !(MeansArray$HGNC=='') & !(is.na(MeansArray$HGNC))&
                                !(grepl('Empty',MeansArray$Description)),]
#dim(MeansArray.Prot)#8968 x 27
#length(unique(MeansArray.Prot$HGNC))#7427

#Now to permutation test
eBayesResult <- eBayesArrays(
  MeansArray.Prot, 
  c('control','wt','D2A','control','wt','D2A','control','wt','D2A')
)
eBayesResultCorr <- eBayesArrays(MeansArray.Prot, 
                                 c('control','wt','D2A','control','wt','D2A','control','wt','D2A'), TRUE)
wtcontPermuted <- permuteBackgroundwtcont(
  eBayesResultCorr$p.value[,1],
  eBayesResultCorr$t[,1],
  MeansArray.Prot[,!(colnames(MeansArray.Prot) %in%
                       c('Accession','Description','HGNC','Known.GlcNAc',
                         'PhosphoSitePlus','Wang'))],
  c('control','wt','D2A','control','wt','D2A','control','wt','D2A'),
  c('wt','control'), alpha = 0.25,oneway = TRUE, N = 100
)
D2AcontPermuted <- permuteBackgroundD2Acont(
  eBayesResultCorr$p.value[,2],
  eBayesResultCorr$t[,2],
  MeansArray.Prot[,!(colnames(MeansArray.Prot) %in%
                       c('Accession','Description','HGNC','Known.GlcNAc',
                         'PhosphoSitePlus','Wang'))],
  c('control','wt','D2A','control','wt','D2A','control','wt','D2A'),
  c('D2A','control'), alpha = 0.25,oneway = TRUE, N = 100
)



library(qvalue)
wtcont.qval <- qvalue(wtcontPermuted$pvals[,1],pi0.method = 'bootstrap')
D2A.qval <- qvalue(D2AcontPermuted$pvals[,1],pi0.method = 'bootstrap')


#make table of final results
FinalResults <- MeansArray.Prot[,c(1,2,12:15, 3,6,9,4,7,10,5,8,11)]

#add Median-subtracted data
FinalResults.medCorrect <- medianCorrect(FinalResults, 7:9, list('wt' = 10:12, 'D2A' = 13:15))
colnames(FinalResults.medCorrect) <- paste(colnames(FinalResults.medCorrect), 'MedianCorrected')

FinalResults <- cbind(FinalResults, FinalResults.medCorrect[,10:15])
FinalResultsMeans <- meanSignal(FinalResults, list('control' = 7:9, 
                                                    'wt'= 10:12,
                                                    'D2A' = 13:15,
                                                    'wtMedCor' = 16:18,
                                                    'D2AMedCor' = 19:21))
FinalResults.sum <- cbind(FinalResults[,1:6],FinalResultsMeans) #note- this excludes individual arrays. Can be obtained from NCBI GEO and

#use of analysis scripts.
FinalResults.sum[,c('wt fold signal','D2A fold signal','wt fold signal (Median Corrected)','D2A fold signal (Median Corrected)',
                    'wt vs. D2A fold signal', 'wt vs. D2A fold signal (Median Corrected)',
                    'wt - D2A signal', 'wt - D2A signal MedCor')] <- cbind(
                      FinalResults.sum[,'wt Mean']/FinalResults.sum[,'control Mean'],
                      FinalResults.sum[,'D2A Mean']/FinalResults.sum[,'control Mean'],
                      FinalResults.sum[,'wtMedCor Mean']/FinalResults.sum[,'control Mean'],
                      FinalResults.sum[,'D2AMedCor Mean']/FinalResults.sum[,'control Mean'],
                      FinalResults.sum[,'wt Mean']/FinalResults.sum[,'D2A Mean'],
                      FinalResults.sum[,'wtMedCor Mean']/FinalResults.sum[,'D2AMedCor Mean'],
                      FinalResults.sum[,'wt Mean'] - FinalResults.sum[,'D2A Mean'],
                      FinalResults.sum[,'wtMedCor Mean'] - FinalResults.sum[,'D2AMedCor Mean']
                    )
FinalResults.sum[,c('wt glycosylation qvalue','D2A glycosylation qvalue')] <- cbind(wtcont.qval$qvalues,
                                                                                    D2A.qval$qvalues)
FinalResults.sum[,c('wt glycosylation pvalue','D2A glycosylation pvalue')] <- cbind(wtcont.qval$pvalues,
                                                                                    D2A.qval$pvalues)
FinalResults.sum[,c('wt negative log pvalue',
                    'D2A negative log pvalue')] <- cbind(-log10(FinalResults.sum$`wt glycosylation pvalue`),
                                                                -log10(FinalResults.sum$`D2A glycosylation pvalue`))
FinalResults.sum[,c('log wt fold signal',
                    'log D2A fold signal',
                    'log wt v. D2A fold signal',
                    'log wt v. D2A fold signal (MedCor)',
                    'log wt fold signal (MedCor)',
                    'log D2A fold signal (MedCor)')] <- cbind(log10(FinalResults.sum$`wt fold signal`),
                                                              log10(FinalResults.sum$`D2A fold signal`),
                                                              log10(FinalResults.sum$`wt vs. D2A fold signal`),
                                                              log10(FinalResults.sum$`wt vs. D2A fold signal (Median Corrected)`),
                                                              log10(FinalResults.sum$`wt fold signal (Median Corrected)`),
                                                              log10(FinalResults.sum$`D2A fold signal (Median Corrected)`))
FinalResults.sum[,c('wt Hit','D2A Hit')] <- cbind(#178 wt hits, 162 D2A hits, 224 total
  FinalResults.sum$`wt fold signal` > 10 & FinalResults.sum$`wt glycosylation qvalue` < 0.05,
  FinalResults.sum$`D2A fold signal` > 10 & FinalResults.sum$`D2A glycosylation qvalue` < 0.05
)

#now adjust to compare directly just those that are a hit in one of the control comparisons
whichData <- FinalResults.sum$`wt Hit` | FinalResults.sum$`D2A Hit`
wtvD2AcontPermuted <- permuteBackgroundwtvD2A(
  eBayesResult$p.value[whichData,3],
  eBayesResult$t[whichData,3],
  MeansArray.Prot[whichData,!(colnames(MeansArray.Prot) %in%
                                c('Accession','Description','HGNC','Known.GlcNAc',
                                  'PhosphoSitePlus','Wang'))],
  c('control','wt','D2A','control','wt','D2A','control','wt','D2A'),
  c('wt', 'D2A'), alpha = 0.25,oneway = FALSE, N = 100
)
wtvD2A.qval <- qvalue(wtvD2AcontPermuted$pvals[,1],pi0.method = 'bootstrap')
#looking at p-value distribution suggests p <0.05 as reasonable, despite 20% FDR- 
#data is noisy for direct comparison, but there is clear evidence of bias if limited confidence 
#in individual hits


FinalResults.sum[,'wt vs. D2A glycosylation qvalue'] <- as.numeric(rep(NA,#note- I'll redo this below for plotting purposes with this as 1
                                                                       length(FinalResults.sum$Accession)))
FinalResults.sum[whichData,'wt vs. D2A glycosylation qvalue'] <- wtvD2A.qval$qvalues
FinalResults.sum[,'wt vs. D2A glycosylation qvaluePlotting'] <- as.numeric(rep(1,#for plotting/selecting hits I don't want missing values
                                                                               length(FinalResults.sum$Accession)))
FinalResults.sum[whichData,'wt vs. D2A glycosylation qvaluePlotting'] <- wtvD2A.qval$qvalues
FinalResults.sum[,'wt vs. D2A glycosylation pvalue'] <- as.numeric(rep(NA,
                                                                       length(FinalResults.sum$Accession)))

FinalResults.sum[whichData,'wt vs. D2A glycosylation pvalue'] <- wtvD2A.qval$pvalues





#cutoffs 1.8 fold, 13 in favor WT, 17 in favor of D2A, 30 total
#max(wtvD2A.qval$qvalues[wtvD2A.qval$pvalues<0.05]) #19% FDR at this point

FinalResults.sum[,'wt v. D2A Hit'] <- ((FinalResults.sum$`wt Hit` | FinalResults.sum$`D2A Hit`) &
                                       (FinalResults.sum$`wt vs. D2A glycosylation pvalue` < 0.05) &
                                         ((FinalResults.sum$`wt vs. D2A fold signal (Median Corrected)`>1.8)|
                                            (FinalResults.sum$`wt vs. D2A fold signal (Median Corrected)`<(1/1.8))))


#Save data
write.csv(FinalResults.sum, row.names=FALSE, file="D2AFinalResultsSummary.csv")
write.csv(FinalResults, row.names=FALSE, file="ForCompositePlatform.csv")

onlyHits <- FinalResults.sum[FinalResults.sum$`wt Hit`|FinalResults.sum$`D2A Hit`,]
onlyHits <- onlyHits[order(onlyHits$`wt fold signal`),]
write.csv(onlyHits, row.names = F, file = 'HitsOnly.csv')

#plot data
library(ggplot2)
#histograms and qqplots 
hist(FinalResults.sum$`log wt fold signal`, breaks = 20, xlim = c(-0.5,2.1))
abline(v=1,lty=2)
qqnorm(FinalResults.sum$`log wt fold signal`, ylim = c(-0.5,2.1))
qqline(FinalResults.sum$`log wt fold signal`, ylim = c(-0.5,2.1))
abline(h=1,lty=2)

hist(FinalResults.sum$`log D2A fold signal`, breaks = 20, xlim = c(-0.5,2.1))
abline(v=1,lty=2)
qqnorm(FinalResults.sum$`log D2A fold signal`, ylim = c(-0.5,2.1))
qqline(FinalResults.sum$`log D2A fold signal`, ylim = c(-0.5,2.1))
abline(h=1,lty=2)

hist(wtvD2A.qval$pvalues,breaks = 20)
#scatterplot
ggplot(FinalResults.sum[(!(FinalResults.sum$`wt Hit`)) |(!(FinalResults.sum$`D2A Hit`)),],
       aes(x = `log wt fold signal`, y = `log D2A fold signal`))+
  geom_point(col = 'grey') + 
  geom_point(data = FinalResults.sum[(FinalResults.sum$`wt Hit` |
                                       FinalResults.sum$`D2A Hit`)&
                                       (!(FinalResults.sum$`wt v. D2A Hit`)),
                                     ], aes(x = `log wt fold signal`, y = `log D2A fold signal`),
             col = 'blue') +
  geom_point(data = FinalResults.sum[((FinalResults.sum$`wt Hit`) |
                                        (FinalResults.sum$`D2A Hit`))&
                                       ((FinalResults.sum$`wt v. D2A Hit`)&(FinalResults.sum$`wt vs. D2A fold signal`>1)),
                                     ], aes(x = `log wt fold signal`, y = `log D2A fold signal`),
             col = 'red') + xlim(-.5, 2.15) + ylim(-.5,2.15) +
  geom_point(data = FinalResults.sum[((FinalResults.sum$`wt Hit`) |
                                        (FinalResults.sum$`D2A Hit`))&
                                       ((FinalResults.sum$`wt v. D2A Hit`)&(FinalResults.sum$`wt vs. D2A fold signal`<1)),
                                     ], aes(x = `log wt fold signal`, y = `log D2A fold signal`),
             col = 'orange') + xlim(-.5, 2.15) + ylim(-.5,2.15) +
  theme_classic()

#correlation analysis Pearson
Hits <- FinalResults.sum$wt.Hit | FinalResults.sum$D2A.Hit
cor(FinalResults.sum$log.wt.fold.signal,FinalResults.sum$log.D2A.fold.signal) #0.778
cor(FinalResults.sum$log.wt.fold.signal[Hits],FinalResults.sum$log.D2A.fold.signal[Hits]) #0.708
cor(FinalResults.sum$log.wt.fold.signal[!Hits],FinalResults.sum$log.D2A.fold.signal[!Hits]) #0.728
