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

plotStarts<-function(rna,prot,nearbyBases=prot[,'tss']+-50:25,nearbyBases2=prot[,'tss']+-25:50,main=prot[,'name']){
  selectStarts<-rna[as.character(nearbyBases),as.character(nearbyBases2)]
  logData<-log10(selectStarts+1)
  breaks<-seq(-1e-9,max(logData)+1e-9,length.out=101)
  prettyLabs<-unique(round(pretty(range(logData))))
  prettyLabs<-prettyLabs[prettyLabs<=max(breaks)]
  image(nearbyBases,nearbyBases2,logData,col=cols,xlab='',ylab='End position',main=main,las=1,breaks=breaks,mgp=c(2.7,1,0))
  title(xlab='Start position',mgp=c(2.2,1,0))
  abline(v=prot[,'tss']-12,h=prot[,'tss']+16,col='#00000033',lty=2)
  abline(v=prot[,'tss']-12-28,h=prot[,'tss']+16-28,col='#00000033',lty=3)
  abline(28,1,col='#00000022')
  #abline(v=prots[-jj,'tss']-12,h=prots[-jj,'tss']+16,col='#00000033',lty=3)
  insetPos<-c(grconvertX(0.01,'nfc','user'),grconvertY(0.025,'nfc','user'),grconvertX(0.4,'nfc','user'),grconvertY(0.045,'nfc','user'))
  breakPos<-(breaks-min(breaks))/max(breaks-min(breaks))*(insetPos[3]-insetPos[1])+insetPos[1]
  rect(breakPos[-1]+.1,insetPos[2],breakPos[-length(breakPos)],insetPos[4],col=cols,xpd=NA,border=NA)
  rect(insetPos[1],insetPos[2],insetPos[3],insetPos[4],xpd=NA)
  prettyPos<-(ifelse(prettyLabs==0,0,log10(10^prettyLabs+1))-min(breaks))/max(breaks-min(breaks))*(insetPos[3]-insetPos[1])+insetPos[1]
  segments(prettyPos,insetPos[2],prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.1,xpd=NA)
  text(prettyPos,insetPos[2]-diff(insetPos[c(2,4)])*.175,ifelse(prettyLabs==0,0,10^prettyLabs),xpd=NA,adj=c(.5,1),cex=.85)
  text(mean(insetPos[c(1,3)]),insetPos[4]+diff(insetPos[c(2,4)])*.1,"Read count",xpd=NA,adj=c(.5,0))
}


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
      plotStarts(rnaStarts[[ii]],prots[jj,],main=sprintf("%s %s",ii,prots[jj,'name']))
    }
    nearbyBases<-1:(prots[1,'tss']+18)
    plotStarts(rnaStarts[[ii]],prots[1,],nearbyBases,nearbyBases,sprintf("%s leader",ii))
  dev.off()
}

bp28<-lapply(rnaStarts[!grepl('Total',names(rnaStarts))],function(rs)sapply(1:(nrow(rs)-28),function(x)rs[x,x+27]))
sampleNames<-sub('^[^_]+_','',names(bp28))
treats<-sub('[0-9]_.*$','',names(bp28))
treatCols<-rainbow.lab(length(unique(treats)),lightScale=0,lightMultiple=.8,alpha=.7)
names(treatCols)<-unique(treats)
pdf('out/28bp.pdf',width=16)
  for(ii in sampleNames){
    theseTreats<-treats[sampleNames==ii]
    theseStarts<-bp28[sampleNames==ii]
    ylim<-range(unlist(theseStarts))
    prots<-read.csv(protFiles[names(theseStarts)[1]])
    plot(1,1,type='n',ylim=ylim,xlim=c(1,length(theseStarts[[1]])),main=ii,xlab='Position',ylab='Read count')
    abline(v=prots$tss-12,col='#00000033',lty=2)
    mapply(function(x,col)lines(1:length(x),x,col=col),theseStarts,treatCols[theseTreats])
    legend('top',legend=names(treatCols),col=treatCols,lty=1)
  }
dev.off()

pdf('out/28bp_prot.pdf',width=8)
  for(ii in sampleNames){
    theseTreats<-treats[sampleNames==ii]
    theseStarts<-bp28[sampleNames==ii]
    prots<-read.csv(protFiles[names(theseStarts)[1]])
    prots<-prots[order(prots$tss),]
    for(jj in 1:nrow(prots)){
      startSubset<-lapply(theseStarts,function(x)x[prots[jj,'tss']+-40:40])
      ylim<-range(unlist(startSubset)+1)
      plot(1,1,type='n',ylim=ylim,xlim=c(-40,40),log='y',main=sprintf('%s %s',ii,prots[jj,'name']),xlab='Position',ylab='Read count',las=1)
      abline(v=-12,lty=2,col='#00000033')
      mapply(function(x,col)lines(-40:40,x+1,col=col),startSubset,treatCols[theseTreats])
      legend('top',legend=names(treatCols),col=treatCols,lty=1,inset=.01)
    }
  }
dev.off()



