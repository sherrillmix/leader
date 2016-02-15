library(parallel)
source("~/scripts/R/dna.R")

getHuman<-function(bamFile,regions,strand=c('+','-'),sizeRange=26:30){
  startEnd<-ifelse(all(strand=='-'),'end','start')
  counts<-lapply(regions,function(region){
    reg<-parseRegion(region)
    cmd<-sprintf('~/.local/bin/getstartends %s --region %s --maxGaps 3',bamFile,region)
    output<-system(cmd,intern=TRUE)
    inCon<-textConnection(output)
    out<-read.csv(inCon)
    close(inCon)
    outTable<-unlist(table(out[out$end-out$start+1 %in% sizeRange & out[,startEnd]>=reg$start&out[,startEnd]<=reg$end & out$strand %in% strand,startEnd]))
    outVector<-rep(0,reg$end-reg$start+1)
    names(outVector)<-reg$start:reg$end
    outVector[names(outTable)]<-outTable
    if(all(strand=='-'))outVector<-rev(outVector)
    return(outVector)
  })
  names(counts)<-regions
  return(counts)
}


if(!exists('knownGenes'))knownGenes<-read.table("ref/knownGene.gtf.gz",stringsAsFactors=FALSE)
dataDir<-'work/align/'
allRnaFiles<-list.files(dataDir,'\\.bam$')
names(allRnaFiles)<-sub('.bam$','',basename(allRnaFiles))
sample<-sub('^[^_]+_','',names(allRnaFiles))

if(FALSE){
regs<-c( "ATF1"="uc001rww.5", "CD4"="uc001qqv.3", "ACTB"="uc003sot.5", "GAPDH"="uc001qop.4", "GPI"="uc002nvg.3")
pdf('out/human.pdf',width=16)
for(geneName in names(regs)){ 
  message(geneName)
  ucscId<-regs[geneName]
  target<-knownGenes[knownGenes$V10==ucscId&knownGenes$V3=='start_codon',]
  gene<-knownGenes[knownGenes$V10==ucscId,]
  reg<-sprintf('%s:%d-%d',target$V1,ifelse(target$V7=='+',target$V4,target$V5)-100,ifelse(target$V7=='+',target$V4,target$V5)+100)
  targetFiles<-sub('_virus','',allRnaFiles[sample=='HIV_CH0694'])
  humanStarts<-mclapply(sprintf("%s/%s",dataDir,targetFiles),function(x,...)getHuman(x,reg,strand=target$V7,sizeRange=if(grepl('Total',x))1:100 else 27:29)[[1]],mc.cores=12,mc.preschedule=FALSE)
  names(humanStarts)<-sub('_.*$','',targetFiles)
  treats<-sub('[0-9]$','',names(humanStarts))
  treatCols<-rainbow.lab(length(unique(treats)),lightScale=0,lightMultiple=.8,alpha=.7)
  names(treatCols)<-unique(treats)
  #plotting
    ylim<-range(unlist(humanStarts))
    plot(1,1,type='n',ylim=ylim,xlim=c(-100,100),main=sprintf('%s (%s)',geneName,ucscId),xlab='Position',ylab='Read count')
    abline(v=-12,col='#00000033',lty=2)
    mapply(function(x,col)lines(-100:100,x,col=col),humanStarts,treatCols[treats])
    legend('top',legend=names(treatCols),col=treatCols,lty=1)
    abline(v=ifelse(target$V7=='+',1,-1)*(unlist(gene[gene$V3=='exon',c('V4','V5')])-ifelse(target$V7=='+',target$V4,target$V5)+1),col='#FF000044',lty=3,lwd=3)
  #done plotting
}
dev.off()
}


lee<-read.csv('ref/sd01.csv',comment='#',stringsAsFactors=FALSE)
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("Homo.sapiens")
#x<-transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,'gene')
cds<-cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,'tx',use.names=TRUE)
fives<-fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene,use.names=TRUE)
#some cds dont have 5' UTR (assuming junk)
cds<-cds[names(fives)]
selector<-grepl('^chr[0-9XYM]+$',sapply(fives,function(x)as.character(seqnames(x)@values)[1]))
cds<-cds[selector]
fives<-fives[selector]
coordToReg<-function(coord){
  #may need to adjust start
  sprintf('%s:%d-%d',seqnames(coord),start(coord),end(coord))
}
targetFiles<-allRnaFiles[grepl('_HIV_',allRnaFiles)]#&!grepl('^Total',allRnaFiles)])
dataDir<-'work/align/'

naFillRegs<-function(file,regs,naNs,...){
  outs<-lapply(regs,function(reg,...)do.call(c,getHuman(file,reg,...)),...)
  out<-mapply(function(dat,n){
    if(n==0)return(dat)
    if(n>0)return(c(dat,rep(NA,n)))
    else return(c(rep(NA,-n),dat))
  },outs,naNs,SIMPLIFY=FALSE)
  return(out)
}

windowWidth<-100
startCounts<-cacheOperation('work/humanStartCounts.Rdat',mcmapply,function(five,cd,windowWidth=40,...){
  cat('.')
  strand<-strand(cd)@values[1]
  if(any(strand(cd)!=strand)||any(strand(five)!=strand))stop(simpleError('Strand mismatch'))
  closeFive<-five[c(rev(cumsum(rev(width(five))))[-1],0)<windowWidth]
  closeCd<-cd[c(0,cumsum(width(cd))[-length(cd)])<windowWidth]
  if(strand=='+'){
    newEnd<-start(closeCd)[length(closeCd)]+windowWidth-sum(width(closeCd)[-length(closeCd)])-1
    end(closeCd)[length(closeCd)]<-min(newEnd,end(closeCd)[length(closeCd)])
    newStart<-end(closeFive)[1]-(windowWidth-sum(width(closeFive)[-1]))+1
    start(closeFive)[1]<-max(newStart,start(closeFive)[1])
  }else{
    newEnd<-start(closeFive)[1]+windowWidth-sum(width(closeFive)[-1])-1
    end(closeFive)[1]<-min(newEnd,end(closeFive)[1])
    newStart<-end(closeCd)[length(closeCd)]-(windowWidth-sum(width(closeCd)[-length(closeCd)]))+1
    start(closeCd)[length(closeCd)]<-max(newStart,start(closeCd)[length(closeCd)])
  }
  fiveN<-sum(width(closeFive))
  cdN<-sum(width(closeCd))
  regs<-list(sapply(closeFive,coordToReg),sapply(closeCd,coordToReg))
  #if(strand=='+')revFunc<-function(y)lapply(y,function(z)z) else revFunc<-function(y)lapply(y,function(z)rev(z))
  #humanStarts<-do.call(rbind,lapply(sprintf("%s/%s",dataDir,targetFiles),function(x,revFunc,...)do.call(c,revFunc(getHuman(x,regs,strand=strand,sizeRange=27:29))),revFunc=revFunc))
  naNs<-c(fiveN-windowWidth,windowWidth-cdN)
  humanStarts<-do.call(rbind,lapply(sprintf("%s/%s",dataDir,targetFiles),function(x,...)do.call(c,naFillRegs(x,regs,naNs,strand=strand,sizeRange=27:29))))
  rownames(humanStarts)<-names(targetFiles)
  return(humanStarts)
},fives,cds,SIMPLIFY=FALSE,mc.cores=12,MoreArgs=list(windowWidth=windowWidth))


goodCounts<-startCounts[sapply(startCounts,sum,na.rm=TRUE)>50]
#goodCounts<-startCounts[sapply(startCounts,sum,na.rm=TRUE)>50&sapply(fives,function(x)strand(x)@values[1])=='-']
goodProp<-lapply(goodCounts,function(x)t(apply(x,1,function(y)if(sum(y,na.rm=TRUE)>4)ifelse(is.na(y),0,y)/sum(y,na.rm=TRUE) else rep(NA,length(y)))))
library(abind)
allProps<-do.call(abind,c(goodProp,list(along=3)))
meanProp<-apply(allProps,c(1,2),mean,na.rm=TRUE)

treats<-sub('[0-9]_.*$','',rownames(meanProp))
treatCols<-rainbow.lab(length(unique(treats)),lightScale=0,lightMultiple=.8,alpha=.7)
names(treatCols)<-unique(treats)

pdf('out/meanProps.pdf',height=5,width=7)
  par(mar=c(4,4,.1,.1))
  plot(1,1,type='n',xlim=c(-windowWidth+1,windowWidth),ylim=c(0,max(meanProp))*100,xlab='Offset from TIS',ylab='Mean percent of reads',las=1)
  for(ii in 1:nrow(meanProp)){
    lines((-windowWidth+1):windowWidth,meanProp[ii,]*100,col=treatCols[treats[ii]],lwd=2)
  }
  legend('topright',legend=names(treatCols),col=treatCols,lty=1,inset=.01,lwd=2)
  #abline(v=1,lty=2)
  abline(v=-11,lty=2)
  abline(v=-11+seq(3,36,3),lty=3,col='#00000055')
dev.off()

