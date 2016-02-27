library(parallel)
library(abind)
source("~/scripts/R/dna.R")
source("functions.R")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("Homo.sapiens")

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
targets<-sprintf('%s/%s',dataDir,targetFiles)
names(targets)<-names(targetFiles)
startCounts<-cacheOperation('work/humanStartCounts.Rdat',mcmapply,pullCdFiveRegion,fives,cds,SIMPLIFY=FALSE,mc.cores=8,MoreArgs=list(targetFiles=targets,windowWidth=windowWidth))


allCounts<-do.call(abind,c(startCounts,list(along=3)))
treatCounts<-apply(allCounts,1,sum,na.rm=TRUE)
goodCounts<-startCounts[sapply(startCounts,sum,na.rm=TRUE)>50]
meanProp<-getMeanProp(goodCounts)
meanProp236<-meanProp[grep('CH0236',rownames(meanProp)),]
meanProp694<-meanProp[grep('CH0694',rownames(meanProp)),]

pdf('out/meanRiboProps.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  plotTreats(meanProp694)
  title(main='HIV CH0694')
  plotTreats(meanProp236)
  title(main='HIV CH0236')
dev.off()

pdf('out/meanRiboByTreatment.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  for(ii in 1:nrow(meanProp)){
    message(ii)
    thisMean<-getMeanProp(lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)],function(x)x[ii,,drop=FALSE]))
    plotTreats(thisMean)
    title(main=sprintf('%s (%s reads)',rownames(meanProp)[ii],format(treatCounts[ii],big.mark=',')))
  }
dev.off()

for(ii in 1:nrow(meanProp)){
  png(sprintf('out/pileups/human_%02d.png',ii),height=1000,width=2000,res=200)
    par(mar=c(3.5,3.5,1.1,.1))
    message(ii)
    thisCounts<-lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>100)],function(x)x[ii,,drop=FALSE])
    plot(1,1,type='n',ylim=c(0,1),xlim=c(-99,100),main=rownames(meanProp)[ii],las=1,mgp=c(2.4,.8,0),xlab='Position relative to TIS',ylab='Proportions of reads starting')
    for(thisGene in thisCounts){
      lines(c(-100:99),thisGene/sum(thisGene,na.rm=TRUE),col='#00000001')
    }
    abline(v=-12,lty=2,col='#00000033')
    abline(v=seq(-9,99,3),lty=3,col='#00000033')
  dev.off()
}

