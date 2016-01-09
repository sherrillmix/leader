
library(parallel)
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

treatOrder<-order(sub('(Total[1-2])(.*)','\\2\\1',names(baseCounts)))
juncFiles<-list.files('data','.*sjdbFile.tsv',recursive=TRUE,full.names=TRUE)
sapply(sub('Total[0-9]_','',names(baseCounts)),function(x)x)

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
  lines(pos,pCov/scale+1,col='blue',lwd=2)
  lines(pos,nCov/scale+1,col='red',lwd=2)
  legend('topright',c('Positive','Negative'),col=c("blue",'red'),lwd=2,inset=.02)
  #plot(pos,nCov,type='l',xlim=xLim,ylim=range(nCov),ylab='',xlab='Virus base position',las=1,mgp=c(1.5,.5,0),main=ii,lwd=2)
  #title(ylab='Negative strand coverage',mgp=c(3.5,1,0))
}
dev.off()
