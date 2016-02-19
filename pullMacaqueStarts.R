library(parallel)
library(abind)
library(GenomicFeatures)
source("~/scripts/R/dna.R")
source("functions.R")


dataDir<-'work/align/'
allRnaFiles<-list.files(dataDir,'\\.bam$')
names(allRnaFiles)<-sub('.bam$','',basename(allRnaFiles))
sample<-sub('^[^_]+_','',names(allRnaFiles))
targetFiles<-allRnaFiles[grepl('_SIV_',allRnaFiles)]#&!grepl('^Total',allRnaFiles)])

#getChromInfoFromUCSC('rheMac3')
if(!exists('rheRef'))rheRef <- makeTxDbFromUCSC( genome = "rheMac3" , tablename = "refGene")

tmp<-getCdFives(cdsBy(rheRef,'tx'),fiveUTRsByTranscript(rheRef))
cds<-tmp[['cds']]
fives<-tmp[['fives']]
rm(tmp)


windowWidth<-100
startCounts<-cacheOperation('work/macaqueStartCounts.Rdat',mcmapply,pullCdFiveRegion,fives,cds,SIMPLIFY=FALSE,mc.cores=12,MoreArgs=list(targetFiles=targetFiles,windowWidth=windowWidth))

allCounts<-do.call(abind,c(startCounts,list(along=3)))
treatCounts<-apply(allCounts,1,sum,na.rm=TRUE)
goodCounts<-startCounts[sapply(startCounts,sum,na.rm=TRUE)>50]
meanProp<-getMeanProp(goodCounts)

pdf('out/macaque_meanRiboProps.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  plotTreats(meanProp)
  title(main='SIV 766')
dev.off()

pdf('out/macque_meanRiboByTreatment.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  for(ii in 1:nrow(meanProp)){
    message(ii)
    thisMean<-getMeanProp(lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)],function(x)x[ii,,drop=FALSE]))
    plotTreats(thisMean)
    title(main=sprintf('%s (%s reads)',rownames(meanProp)[ii],format(treatCounts[ii],big.mark=',')))
  }
dev.off()

