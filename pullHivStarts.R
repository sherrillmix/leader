source("~/scripts/R/dna.R")
library(parallel)
library(levenR)
library(dnaplotr)

getStarts<-function(bamFile,fillMatrix=TRUE){
  outFile<-sub('bam$','starts.gz',bamFile)
  message(bamFile,"->",outFile)
  cmd<-sprintf('~/.local/bin/getstartends %s --maxGaps 3|gzip>%s',bamFile,outFile)
  if(!file.exists(outFile))system(cmd)
  else message("Using cached ",outFile)
  starts<-read.csv(outFile,stringsAsFactors=FALSE)
  starts<-starts[starts$strand=='+',]
  out<-tapply(starts$start,list(starts$start,starts$end),length)
  if(fillMatrix){
    maxEnd<-max(starts$end)
    message(maxEnd)
    #fill out the matrix so all pos present
    allPos<-1:maxEnd
    missingCols<-!allPos %in% as.numeric(colnames(out))
    fillCols<-do.call(cbind,lapply(1:sum(missingCols),function(x)out[,1]))
    colnames(fillCols)<-as.character(allPos[missingCols])
    fillCols[,]<-NA
    out<-cbind(out,fillCols)
    missingRows<-!allPos %in% as.numeric(rownames(out))
    fillRows<-do.call(rbind,lapply(1:sum(missingRows),function(x)out[1,]))
    rownames(fillRows)<-as.character(allPos[missingRows])
    fillRows[,]<-NA
    out<-rbind(out,fillRows)
    out<-out[as.character(allPos),as.character(allPos)]
  }
  out[is.na(out)]<-0
  return(out)
}

dataDir<-'work/virusAlign'
allRnaFiles<-list.files(dataDir,'\\.bam$')
if(!exists('rnaStarts'))print(system.time(rnaStarts<-mclapply(sprintf("%s/%s",dataDir,allRnaFiles),getStarts,mc.cores=12,mc.preschedule=FALSE)))
names(rnaStarts)<-sub('_virus.bam','',basename(allRnaFiles))


allData<-list.files('data',full.names=TRUE,recursive=TRUE)
allData<-data.frame('dir'=dirname(allData),'file'=basename(allData),stringsAsFactors=FALSE)
dirs<-allData[sapply(names(rnaStarts),grep,allData$file),'dir']
refFiles<-sapply(dirs,function(x)paste(allData[grepl('fasta$',allData$file)&allData$dir==x,c('dir','file')],collapse='/'))
protFiles<-sapply(dirs,function(x)paste(allData[grepl('^prots.csv$',allData$file)&allData$dir==x,c('dir','file')],collapse='/'))
names(refFiles)<-names(protFiles)<-names(rnaStarts)

for(ii in names(rnaStarts)){
  message(ii)
  if(FALSE){
    #not particularly useful since everything is on the diagonal of course
    outFile<-sprintf('out/startEnds/%s.png',ii)
    png(outFile,height=2000,width=2000,res=250)
      cols<-rev(heat.colors(100))
      image(log10(rnaStarts[[ii]]+1),col=cols,xlab='Start',ylab='End')
    dev.off()
  }
  prots<-read.csv(protFiles[ii])
  prots<-prots[order(prots$tss),]
  outFile<-sprintf('out/startEnds/prots_%s.pdf',ii)
  cols<-rev(heat.colors(100))
  pdf(outFile)
    for(jj in 1:nrow(prots)){
      nearbyBases<-prots[jj,'tss']+-25:75
      selectStarts<-rnaStarts[[ii]][as.character(nearbyBases),as.character(nearbyBases)]
      image(nearbyBases,nearbyBases,log10(selectStarts+1),col=cols,xlab='Start',ylab='End',main=prots[jj,'name'],las=1)
      abline(v=prots[jj,'tss']-12,h=prots[jj,'tss']+16,col='#00000033',lty=2)
      abline(v=prots[-jj,'tss']-12,h=prots[-jj,'tss']+16,col='#00000033',lty=3)
    }
    nearbyBases<-1:(prots[1,'tss']+18)
    selectStarts<-rnaStarts[[ii]][as.character(nearbyBases),as.character(nearbyBases)]
    image(nearbyBases,nearbyBases,log10(selectStarts+1),col=cols,xlab='Start',ylab='End',main="Leader",las=1)
    abline(v=prots[1,'tss']-12,h=prots[1,'tss']+16,col='#00000033',lty=2)
  dev.off()
}



