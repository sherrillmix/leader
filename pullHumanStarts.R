library(parallel)
library(abind)
source("~/scripts/R/dna.R")
source("functions.R")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("Homo.sapiens")

if(!exists('knownGenes'))knownGenes<-read.table("ref/knownGene.gtf.gz",stringsAsFactors=FALSE)
dataDir<-'work/align/'
allRnaFiles<-list.files(dataDir,'\\.bam$')
names(allRnaFiles)<-sub('.bam$','',basename(allRnaFiles))
sample<-sub('^[^_]+_','',names(allRnaFiles))
targetFiles<-allRnaFiles[grepl('_HIV_',allRnaFiles)]#&!grepl('^Total',allRnaFiles)])

tmp<-getCdFives(cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,'tx',use.names=TRUE),fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,use.names=TRUE))
cds<-tmp[['cds']]
fives<-tmp[['fives']]
rm(tmp)

windowWidth<-100
startCounts<-cacheOperation('work/humanStartCounts.Rdat',mcmapply,pullCdFiveRegion,fives,cds,SIMPLIFY=FALSE,mc.cores=12,MoreArgs=list(windowWidth=windowWidth))


allCounts<-do.call(abind,c(startCounts,list(along=3)))
treatCounts<-apply(allCounts,1,sum,na.rm=TRUE)
goodCounts<-startCounts[sapply(startCounts,sum,na.rm=TRUE)>50]
meanProp<-getMeanProp(goodCounts)
meanProp236<-apply(allProps[grep('CH0236',rownames(allProps)),,],c(1,2),mean,na.rm=TRUE)
meanProp694<-apply(allProps[grep('CH0694',rownames(allProps)),,],c(1,2),mean,na.rm=TRUE)

pdf('out/meanRiboProps.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  plotTreats(meanProp694)
  title(main='HIV CH0694')
dev.off()

pdf('out/meanRiboByTreatment.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  for(ii in 1:nrow(allProps)){
    message(ii)
    thisMean<-getMeanProp(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)])[ii,,,drop=FALSE]
    plotTreats(thisMean)
    title(main=sprintf('%s (%s reads)',rownames(allProps)[ii],format(treatCounts[ii],big.mark=',')))
  }
dev.off()

