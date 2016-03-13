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


pdf('out/macaque_meanRiboProps.pdf',height=4,width=7)
  par(mar=c(3,3,1.1,.1))
  plotTreats(meanProp,treatCols=ranjitColors)
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
  for(ii in 1:nrow(ltmProps)){
    plotTreats(ltmProps[ii,,drop=FALSE],ylab='Difference between LTM and CHX proportions')
    title(main=sprintf('%s',sub('LTM[0-9]-CHX[0-9]_','',rownames(ltmProps)[ii])))
  }
dev.off()


#summary(allCounts['DMSO1_SIV_766WT',,2122]==startCounts[[2122]]['DMSO1_SIV_766WT',])
#x<-system('samtools view work/align/DMSO1_SIV_766WT.bam chr9:13812594-13812627',intern=TRUE)
#table(as.numeric(sapply(strsplit(x,'\t'),'[[',4)))
#tail(startCounts[[2122]]['DMSO1_SIV_766WT',],20)
