
source("~/scripts/R/dna.R")
library(parallel)
library(dnaplotr)
getBaseCounts<-function(bamFile){
  outFile<-sub('bam$','baseCount.gz',bamFile)
  message(bamFile,"->",outFile)
  cmd<-sprintf('~/.local/bin/countbases --strand %s|gzip>%s',bamFile,outFile)
  if(!file.exists(outFile))system(cmd)
  else message("Using cached ",outFile)
  output<-readLines(outFile)
  output[1]<-gsub('\\-','n',gsub('\\+','p',output[1]))
  con<-textConnection(output)
  out<-read.csv(con,stringsAsFactors=FALSE)
  out$pos<-out$pos+1
  return(out)
}

dataDir<-'work/virusAlign'
bulkRnaFiles<-list.files(dataDir,'Total.*bam$')
if(!exists('baseCounts'))print(system.time(baseCounts<-mclapply(sprintf("%s/%s",dataDir,bulkRnaFiles),getBaseCounts,mc.cores=10)))
names(baseCounts)<-sub('_virus.bam','',basename(bulkRnaFiles))

bases<-c('A','C','G','T')
pBases<-sprintf('%sp',bases)
nBases<-sprintf('%sn',bases)

#sapply(sub('Total[0-9]_','',names(baseCounts)),function(x)x)
allData<-list.files('data',full.names=TRUE,recursive=TRUE)
allData<-data.frame('dir'=dirname(allData),'file'=basename(allData),stringsAsFactors=FALSE)
dirs<-allData[sapply(names(baseCounts),grep,allData$file),'dir']
refFiles<-sapply(dirs,function(x)paste(allData[grepl('fasta$',allData$file)&allData$dir==x,c('dir','file')],collapse='/'))
juncFiles<-sapply(dirs,function(x)paste(allData[grepl('sjdbFile.tsv$',allData$file)&allData$dir==x,c('dir','file')],collapse='/'))
names(refFiles)<-names(juncFiles)<-names(baseCounts)

baseCounts<-mapply(function(baseCount,refFile,fileName){
  ref<-read.fa(refFile)$seq
  baseCount$refBase<-substring(ref,baseCount$pos,baseCount$pos)
  baseCount$maxBase<-bases[apply(baseCount[,pBases],1,which.max)]
  baseCount$diversity<-apply(baseCount[,pBases],1,shannon,base=2)
  baseCount$maxProp<-apply(baseCount[,pBases],1,function(x)max(x)/sum(x))
  baseCount$seqContext<-substring(ref,baseCount$pos-10,baseCount$pos+10)
  baseCount$file<-fileName
  return(baseCount)
},baseCounts,refFiles,names(baseCounts),SIMPLIFY=FALSE)


diverseRefPos<-unique(do.call(rbind,lapply(baseCounts,function(x){
  x[x$maxProp<.8,c('ref','pos')]
})))
allBaseCounts<-do.call(rbind,baseCounts)
interestingIds<-unlist(apply(diverseRefPos,1,function(x)which(allBaseCounts$ref==x[1]&allBaseCounts$pos==as.numeric(x[2]))))

diverseBases<-allBaseCounts[interestingIds,!colnames(allBaseCounts) %in% c('n','diversity',nBases)]
diffBases<-do.call(rbind,lapply(baseCounts,function(x){
  x[x$refBase!=x$maxBase&x$diversity<.5,!colnames(x) %in% c('n','diversity',nBases)]
}))
diverseBases<-diverseBases[order(diverseBases$ref,diverseBases$pos,diverseBases$file),]
diverseBases$isDiverse<-diverseBases$maxProp<.8
write.csv(diverseBases,'out/diverseBases.csv',row.names=FALSE)
write.csv(diffBases,'out/diffBases.csv',row.names=FALSE)


treatOrder<-order(sub('(Total[1-2])(.*)','\\2\\1',names(baseCounts)))
pdf('out/bulkRnaCoverage.pdf',width=12)
  for(ii in names(baseCounts)[treatOrder]){
    par(mar=c(2.5,5,1.1,.5))
    pCov<-apply(baseCounts[[ii]][,pBases],1,sum)
    nCov<-apply(baseCounts[[ii]][,nBases],1,sum)
    pos<-baseCounts[[ii]]$pos
    yLim<-range(c(pCov,nCov))
    xLim<-range(pos)
    scale<-1000
    plot(1,1,type='n',xlim=xLim,ylim=yLim/scale,ylab='',xlab='Virus base position',las=1,mgp=c(1.5,.5,0),main=ii)
    title(ylab=sprintf('Coverage (x%d)',scale),mgp=c(3.5,1,0))
    lines(pos,pCov/scale,col='#4daf4a',lwd=2)
    lines(pos,nCov/scale,col='#ff7f00',lwd=2)
    legend('topright',c('Positive','Negative'),col=c("#4daf4a",'#ff7f00'),lwd=2,inset=.02)
    juncs<-read.table(juncFiles[ii])
    abline(v=unique(juncs$V2),col='#0000FF66',lty=2)
    abline(v=unique(juncs$V3),col='#FF000066',lty=2)
    plot(1,1,type='n',xlim=xLim,ylim=log10(yLim+1),ylab='',xlab='Virus base position',las=1,mgp=c(1.5,.5,0),main=ii,yaxt='n')
    title(ylab=sprintf('Coverage',scale),mgp=c(3.5,1,0))
    lines(pos,log10(pCov+1),col='#4daf4a',lwd=2)
    lines(pos,log10(nCov+1),col='#ff7f00',lwd=2)
    legend('topright',c('Positive','Negative'),col=c("#4daf4a",'#ff7f00'),lwd=2,inset=.02)
    abline(v=unique(juncs$V2),col='#0000FF66',lty=2)
    abline(v=unique(juncs$V3),col='#FF000066',lty=2)
    prettyY<-pretty(log10(yLim+1))
    axis(2,c(0,log10(10^prettyY+1)),c(0,sapply(prettyY,function(x)ifelse(x<=2,10^x,as.expression(bquote(10^.(x)))))),las=1)
    #plot(pos,nCov,type='l',xlim=xLim,ylim=range(nCov),ylab='',xlab='Virus base position',las=1,mgp=c(1.5,.5,0),main=ii,lwd=2)
    #title(ylab='Negative strand coverage',mgp=c(3.5,1,0))
  }
dev.off()


pullRegionReads<-function(bamFile,region,onlyUngapped=TRUE){
  output<-strsplit(system(sprintf('samtools view %s %s',bamFile,region),intern=TRUE),'\t')
  names<-sapply(output,'[[',1)
  cigars<-sapply(output,'[[',6)
  starts<-sapply(output,'[[',4)
  seqs<-sapply(output,'[[',10)
  unspliced<-grepl('^[0-9MS]+$',cigars)
  badSpliced<-grepl('[0-9][0-9]+[IDN]',cigars)
  seqs<-seqs[unspliced]
  starts<-as.numeric(starts[unspliced])
  names<-names[unspliced]
  start<-min(starts)
  ends<-starts+nchar(seqs)-1
  end<-max(ends)
  out<-rep(paste(rep('.',end-start+1),collapse=''),length(starts))
  substring(out,starts-start+1,ends-start+1)<-seqs
  return(list('reads'=out,'start'=start,'names'=names))
}
message("Pulling reads for changed bases")
diffReads<-lapply(split(diffBases,1:nrow(diffBases)),function(x){
  message(x$ref,' ',x$pos)
  pullRegionReads(sprintf('work/virusAlign/%s_virus.bam',x$file),sprintf('%s:%d-%d',x$ref,x$pos-25,x$pos+25))
})
diverseReads<-mclapply(split(diverseBases,1:nrow(diverseBases)),function(x){
  message(x$ref,' ',x$pos)
  pullRegionReads(sprintf('work/virusAlign/%s_virus.bam',x$file),sprintf('%s:%d-%d',x$ref,x$pos-25,x$pos+25))
},mc.cores=10)

mapply(function(x,reads){
  message(x$ref,' ',x$pos)
  outFile<-sprintf('out/weirdBases/%s_%d_%s.png',x$ref,x$pos,x$file) 
  png(outFile,width=2000,height=2000)
    plotDNA(reads[['reads']],xStart=reads[['start']],main=sprintf('%s_%d_%s.png',x$ref,x$pos,x$file))
    abline(v=x$pos+c(-.5,.5),lty=2,lwd=2)
  dev.off()
},c(split(diffBases,1:nrow(diffBases)),split(diverseBases,1:nrow(diverseBases))),c(diffReads,diverseReads))

#library(levenR)
#library(dnaplotr)
#ref<-paste(baseCounts[[2]]$refBase,collapse='')
#obs<-paste(baseCounts[[2]]$maxBase,collapse='')
#xx<-levenAlign(ref,obs)
#png('test.png',res=300,width=4000,height=1000);plotDNA(unlist(xx));dev.off()
#png('test.png',res=300,width=4000,height=1000);plotDNA(c(ref,obs),groups=c('ref','obs'));dev.off()
