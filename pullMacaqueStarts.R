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

targets<-sprintf('%s/%s',dataDir,targetFiles)
names(targets)<-names(targetFiles)
startCounts<-cacheOperation('work/macaqueStartCounts.Rdat',mcmapply,pullCdFiveRegion,fives,cds,SIMPLIFY=FALSE,mc.cores=12,MoreArgs=list(targetFiles=targets,windowWidth=windowWidth))

allCounts<-do.call(abind,c(startCounts,list(along=3)))
colnames(allCounts)<-c(-100:-1,1:100)
treatCounts<-apply(allCounts,1,sum,na.rm=TRUE)
selector<-sapply(startCounts,sum,na.rm=TRUE)>50
goodCounts<-startCounts[selector]
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
    thisCounts<-lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)],function(x)x[ii,,drop=FALSE])
    thisMean<-getMeanProp(thisCounts)
    plotTreats(thisMean)
    title(main=sprintf('%s (%s reads)',rownames(meanProp)[ii],format(treatCounts[ii],big.mark=',')))
  }
dev.off()


pdf('out/macque_geneExamples.pdf',height=4,width=7)
  par(mar=c(3.5,3.5,1.1,.1))
  for(ii in 1:nrow(meanProp)){
    message(ii)
    thisCounts<-sample(lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>100)],function(x)x[ii,,drop=FALSE]),10)
    for(thisGene in thisCounts){
      plot(c(-100:99),thisGene,main=rownames(meanProp)[ii],type='l',las=1,mgp=c(2.4,.8,0),xlab='Position relative to TIS',ylab='Reads starting count')
      abline(v=-12,lty=2,col='#00000033')
      abline(v=seq(-9,99,3),lty=3,col='#00000033')
    }
  }
dev.off()

for(ii in 1:nrow(meanProp)){
  png(sprintf('out/pileups/macque_%02d.png',ii),height=1000,width=2000,res=200)
    par(mar=c(3.5,3.5,1.1,.1))
    message(ii)
    thisCounts<-lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>100)],function(x)x[ii,,drop=FALSE])
    plot(1,1,type='n',ylim=c(0,1),xlim=c(-99,100),main=rownames(meanProp)[ii],las=1,mgp=c(2.4,.8,0),xlab='Position relative to TIS',ylab='Proportions of reads starting')
    for(thisGene in thisCounts){
      lines(c(-100:99),thisGene/sum(thisGene,na.rm=TRUE),col='#00000005')
    }
    abline(v=-12,lty=2,col='#00000033')
    abline(v=seq(-9,99,3),lty=3,col='#00000033')
  dev.off()
}


#summary(allCounts['DMSO1_SIV_766WT',,2122]==startCounts[[2122]]['DMSO1_SIV_766WT',])
#x<-system('samtools view work/align/DMSO1_SIV_766WT.bam chr9:13812594-13812627',intern=TRUE)
#table(as.numeric(sapply(strsplit(x,'\t'),'[[',4)))
#tail(startCounts[[2122]]['DMSO1_SIV_766WT',],20)
