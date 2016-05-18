library(parallel)
library(abind)
source("~/scripts/R/dna.R")
source("functions.R")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
#library("Homo.sapiens")

dataDir<-'work/align/'
allRnaFiles<-list.files(dataDir,'\\.bam$')
names(allRnaFiles)<-sub('.bam$','',basename(allRnaFiles))
sample<-sub('^[^_]+_','',names(allRnaFiles))
targetFiles<-allRnaFiles[grepl('_HIV_',allRnaFiles)]#&!grepl('^Total',allRnaFiles)])

tmp<-cacheOperation('work/fives.Rdat',getCdFives,cdsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,'tx',use.names=TRUE),fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE))
cds<-tmp[['cds']]
fives<-tmp[['fives']]
rm(tmp)

windowWidth<-100
targets<-sprintf('%s/%s',dataDir,targetFiles)
names(targets)<-names(targetFiles)
startCounts<-cacheOperation('work/humanStartCounts.Rdat',mcmapply,pullCdFiveRegion,fives,cds,SIMPLIFY=FALSE,mc.cores=18,MoreArgs=list(targetFiles=targets,windowWidth=windowWidth))

fiveCounts<-cacheOperation('work/fiveStartCounts.Rdat',cleanMclapply,fives,12,pullFiveRegion,targetFiles=targets,extraCode='library(dnar);library(GenomicRanges);source("functions.R");')
cdsCounts<-cacheOperation('work/cdsStartCounts.Rdat',cleanMclapply,cds,12,pullFiveRegion,targetFiles=targets,extraCode='library(dnar);library(GenomicRanges);source("functions.R");')


allCounts<-do.call(abind,c(startCounts,list(along=3)))
treatCounts<-apply(allCounts,1,sum,na.rm=TRUE)
startSelector<-sapply(startCounts,sum,na.rm=TRUE)>50
goodCounts<-startCounts[startSelector]
meanProp<-getMeanProp(goodCounts)
meanProp236<-meanProp[grep('CH0236',rownames(meanProp)),]
meanProp694<-meanProp[grep('CH0694',rownames(meanProp)),]
geneCounts<-apply(allCounts,c(1,3),sum)
geneCounts[is.na(geneCounts)]<-0
totalCounts<-geneCounts[grep('Total',rownames(geneCounts)),]


treatments<-sub('[0-9]_.*$','',rownames(goodCounts[[1]]))
viruses<-sub('^[^_]+_','',rownames(goodCounts[[1]]))
countBreakout<-lapply(unique(viruses),function(virus){
  out<-lapply(unique(treatments),function(treat){
    thisSelector<-treatments==treat & viruses==virus
    message(virus,' ',treat)
    out<-do.call(rbind,lapply(goodCounts,function(x,selector){
      if(runif(1)<.001)cat('.')
      x[is.na(x)]<-0
      out<-apply(x[selector,],2,sum)
    },thisSelector))
    colnames(out)<-as.character(-100:99) #magic number
    return(out)
  })
  names(out)<-unique(treatments)
  return(out)
})
names(countBreakout)<-unique(viruses)
startCoords<-sapply(fives[startSelector],function(x){
  strand<-strand(x)@values[1]
  coord<-ifelse(strand=='-',end(x)[1],start(x)[1])
  chr<-as.character(seqnames(x)@values[1])
  out<-sprintf('%s:%s%s',chr,coord,strand)
})

drugData<-lapply(countBreakout,function(counts){
  totals<-do.call(cbind,lapply(counts,function(x)apply(x,1,sum)))
  colnames(totals)<-sprintf('total_%s',colnames(totals))
  props<-lapply(counts,function(count)t(apply(count,1,function(x)x/ifelse(sum(x)==0,1,sum(x)))))
  propSubset<-do.call(cbind,mapply(function(prop,name){
    out<-prop[,as.character(-20:30)]
    colnames(out)<-sprintf('%s_%s',name,colnames(out))
    return(out)
  },props,names(props),SIMPLIFY=FALSE))
  sum3<-do.call(cbind,lapply(props,function(count)apply(count[,as.character(seq(-9,30,3))],1,sum)))
  colnames(sum3)<-sprintf('sum3_%s',colnames(sum3))
  out<-cbind('start'=startCoords,propSubset,totals,sum3)
  return(out)
})
for(ii in names(drugData)){
  write.csv(cbind(drugData[[ii]]),sprintf('out/data_%s.csv',ii))
}


fiveTotals<-sapply(fiveCounts,function(x)sum(x[grepl('Total',rownames(x)),]))


