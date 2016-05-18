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
if(!exists('rheRef'))rheRef <- cacheOperation('work/rheMac3_refGene.Rdat',makeTxDbFromUCSC, genome = "rheMac3" , tablename = "refGene")
#downloaded from ucsc (should be a cleaner way to do this)
idLookup<-id2name(rheRef,'tx')
geneLookup<-read.table('ref/rheMac3_refSeq.tsv.gz',stringsAsFactors=FALSE)
colnames(geneLookup)<-c('nm','symbol')

tmp<-getCdFives(cdsBy(rheRef,'tx'),fiveUTRsByTranscript(rheRef))
cds<-tmp[['cds']]
fives<-tmp[['fives']]
nameLookup<-id2name(rheRef,'tx')

rm(tmp)


windowWidth<-100

targets<-sprintf('%s/%s',dataDir,targetFiles)
names(targets)<-names(targetFiles)
startCounts<-cacheOperation('work/macaqueStartCounts.Rdat',mcmapply,pullCdFiveRegion,fives,cds,SIMPLIFY=FALSE,mc.cores=8,MoreArgs=list(targetFiles=targets,windowWidth=windowWidth))

allCounts<-do.call(abind,c(startCounts,list(along=3)))
colnames(allCounts)<-c(-100:-1,1:100)
treatCounts<-apply(allCounts,1,sum,na.rm=TRUE)
selector<-sapply(startCounts,sum,na.rm=TRUE)>50
goodCounts<-startCounts[selector]
meanProp<-getMeanProp(goodCounts)
geneCounts<-apply(allCounts,c(1,3),sum)
geneCounts[is.na(geneCounts)]<-0
totalCounts<-geneCounts[grep('Total',rownames(geneCounts)),]
sivApproxCount<-20000  #SIV seems to have 15-25k within +-100 from gag from pullHivStarts.R
sivLikeGenes<-apply(totalCounts,2,function(x)all(abs(x-sivApproxCount)<15000))

meanReorder<-meanProp[order(grepl('CHX',rownames(meanProp))),]

pdf('out/macaque_meanRiboProps.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  plotTreats(meanReorder,treatCols=ranjitColors)
  title(main='SIV 766')
dev.off()

pdf('out/macque_meanRiboByTreatment.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  for(ii in 1:nrow(meanProp)){
    message(ii)
    thisCounts<-lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)],function(x)x[ii,,drop=FALSE])
    thisMean<-getMeanProp(thisCounts)
    thisTreat<-sub('[0-9]_.*$','',rownames(meanProp)[ii])
    plotTreats(thisMean,treatCols=ranjitColors[thisTreat])
    title(main=sprintf('%s (%s reads)',rownames(meanProp)[ii],format(treatCounts[ii],big.mark=',')))
  }
dev.off()


pdf('out/macque_geneExamples.pdf',height=4,width=7)
  par(mar=c(3.5,3.5,1.1,.1))
  thisGenes<-startCounts[sivLikeGenes]
  for(geneName in names(thisGenes)){
    thisGene<-thisGenes[[geneName]]
    for(ii in 1:nrow(meanProp)){
      #thisCounts<-sample(lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>100)],function(x)x[ii,,drop=FALSE]),10)
      thisId<-idLookup[geneName]
      thisGeneName<-paste(unique(geneLookup[geneLookup$nm==thisId,'symbol']),collapse=', ')
      thisCounts<-thisGene[ii,,drop=FALSE]
      plot(c(-100:99),thisCounts/1000,main=sprintf('%s %s',thisGeneName,rownames(meanProp)[ii]),type='l',las=1,mgp=c(2.4,.8,0),xlab='Position relative to TIS',ylab='Reads starting count (1000x)')
      abline(v=-12,lty=2,col='#00000033')
      abline(v=seq(-9,99,3),lty=3,col='#00000033')
    }
  }
dev.off()

treatments<-sub('[0-9]_.*$','',rownames(meanProp))
pdf('out/macque_geneExamples_singlePerDrug.pdf',height=8,width=7)
  par(mar=c(3.5,3.5,1.1,.1),mfrow=c(length(unique(treatments)),1))
  thisGenes<-startCounts[sivLikeGenes]
  for(geneName in names(thisGenes)){
    thisGene<-t(apply(thisGenes[[geneName]],1,function(x)x/sum(x,na.rm=TRUE)))
    for(ii in names(ranjitColors)){
      #thisCounts<-sample(lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>100)],function(x)x[ii,,drop=FALSE]),10)
      thisId<-idLookup[geneName]
      thisGeneName<-paste(unique(geneLookup[geneLookup$nm==thisId,'symbol']),collapse=', ')
      thisCounts<-thisGene[treatments==ii,,drop=FALSE]
      plot(1,1,xlim=c(-100,99),ylim=c(0,max(thisGene)*100),main=sprintf('%s %s',thisGeneName,ii),type='n',las=1,mgp=c(2.4,.8,0),xlab='',ylab='Percent of reads in region')
      title(xlab='Position relative to TIS',mgp=c(1.7,1,0))
      apply(thisCounts,1,function(x)lines(-100:99,x*100,col=ranjitColors[ii]))
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

ltmCols<-sort(rownames(startCounts[[1]])[grep('LTM',rownames(startCounts[[1]]))])
chxCols<-sub('LTM','CHX',ltmCols)
ltmChx<-lapply(startCounts[sapply(startCounts,function(x)sum(x[ii,],na.rm=TRUE)>50)],function(x)t(apply(x[ltmCols,],1,function(y)y/ifelse(sum(y)==0,1,sum(y)))-apply(x[chxCols,],1,function(y)y/ifelse(sum(y)==0,1,sum(y)))))
ltmProps<-apply(do.call(abind,c(ltmChx,list(along=3))),c(1,2),mean,na.rm=TRUE)
rownames(ltmProps)<-sub('LTM([0-9])','LTM-CHX\\1',rownames(ltmProps))

pdf('out/macaque_ltm-chx.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  #for(ii in 1:nrow(ltmProps)){
  plotTreats(ltmProps,ylab='Difference between LTM and CHX proportions',ylim=c(0,max(meanProp)*100))
  title(main=sprintf('%s',sub('LTM[0-9]-CHX[0-9]_','',rownames(ltmProps)[1])))
  #}
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


startCoords<-sapply(fives[selector],function(x){
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
  write.csv(drugData[[ii]],sprintf('out/data_%s.csv',ii))
}
