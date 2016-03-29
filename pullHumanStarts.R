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


allCounts<-do.call(abind,c(startCounts,list(along=3)))
treatCounts<-apply(allCounts,1,sum,na.rm=TRUE)
goodCounts<-startCounts[sapply(startCounts,sum,na.rm=TRUE)>50]
meanProp<-getMeanProp(goodCounts)
meanProp236<-meanProp[grep('CH0236',rownames(meanProp)),]
meanProp694<-meanProp[grep('CH0694',rownames(meanProp)),]
geneCounts<-apply(allCounts,c(1,3),sum)
geneCounts[is.na(geneCounts)]<-0
totalCounts<-geneCounts[grep('Total',rownames(geneCounts)),]

meanReorder694<-meanProp694[order(grepl('CHX',rownames(meanProp694))),]
meanReorder236<-meanProp236[order(grepl('CHX',rownames(meanProp236))),]
pdf('out/meanRiboProps.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  plotTreats(meanReorder694,treatCols=ranjitColors)
  title(main='HIV CH0694')
  plotTreats(meanReorder236,treatCols=ranjitColors)
  title(main='HIV CH0236')
dev.off()

pdf('out/meanRiboByTreatment.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  for(ii in 1:nrow(meanProp)){
    message(ii)
    thisMean<-getMeanProp(lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)],function(x)x[ii,,drop=FALSE]))
    thisTreat<-sub('[0-9]_.*$','',rownames(meanProp)[ii])
    plotTreats(thisMean,treatCols=ranjitColors[thisTreat])
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
      lines(c(-100:99),thisGene/sum(thisGene,na.rm=TRUE),col='#00000002')
    }
    abline(v=-12,lty=2,col='#00000033')
    abline(v=seq(-9,99,3),lty=3,col='#00000033')
  dev.off()
}

ltmCols<-sort(rownames(startCounts[[1]])[grep('LTM',rownames(startCounts[[1]]))])
chxCols<-sub('LTM','CHX',ltmCols)
ltmChx<-lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)],function(x)t(apply(x[ltmCols,],1,function(y)y/ifelse(sum(y)==0,1,sum(y)))-apply(x[chxCols,],1,function(y)y/ifelse(sum(y)==0,1,sum(y)))))
ltmProps<-apply(do.call(abind,c(ltmChx,list(along=3))),c(1,2),mean,na.rm=TRUE)
rownames(ltmProps)<-sub('LTM([0-9])','LTM-CHX\\1',rownames(ltmProps))

pdf('out/ltm-chx.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  for(ii in c('CH0694','CH0236')){
    selector<-grep(ii,rownames(ltmProps))
    if(ii=='CH0694')ylim<-c(0,max(meanProp694)*100)
    else if(ii=='CH0236')ylim<-c(0,max(meanProp236)*100)
    else stop(simpleError('Unknown virus'))
    plotTreats(ltmProps[selector,,drop=FALSE],ylab='Difference between LTM and CHX proportions',ylim=ylim,treatCols=c('LTM-CHX'='orange'))
    title(main=sprintf('%s',sub('LTM-CHX[0-9]_','',rownames(ltmProps)[selector][1])))
  }
dev.off()


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
  out<-cbind(propSubset,totals,sum3)
  return(out)
})
for(ii in names(drugData)){
  write.csv(drugData[[ii]],sprintf('out/data_%s.csv',ii))
}


if(FALSE){
  ltm1<-do.call(rbind,lapply(goodCounts,function(x)x['LTM1_HIV_CH0236',]))
  ltm1[is.na(ltm1)]<-0
  ltmProp<-t(apply(ltm1,1,function(x)x/ifelse(sum(x)==0,1,sum(x))))
  colnames(ltmProp)<-c(-100:99)
  prData<-ltmProp[,as.character(-20:30)]
  prData<-cbind(
    prData,
    'sum3'=apply(prData[,as.character(seq(-9,30,3))],1,sum),
    'readCount'=apply(ltm1,1,sum)
  )
  xx<-prcomp(prData)
  pdf('test.pdf');plot(xx);dev.off()
  png('test.png');biplot(xx);dev.off()
}

