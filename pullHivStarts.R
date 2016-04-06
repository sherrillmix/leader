source("~/scripts/R/dna.R")
library(parallel)
library(levenR)
library(dnaplotr)
source("functions.R")

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
spliceFiles<-sapply(dirs,function(x)paste(allData[grepl('^sjdbFile.tsv$',allData$file)&allData$dir==x,c('dir','file')],collapse='/'))
names(refFiles)<-names(protFiles)<-names(spliceFiles)<-names(rnaStarts)

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

bp28<-lapply(rnaStarts,function(rs)sapply(1:(nrow(rs)-28),function(x)sum(rs[x,x+26:28])))
sampleNames<-sub('^[^_]+_','',names(bp28))
treats<-sub('[0-9]_.*$','',names(bp28))
treatCols<-rainbow.lab(length(unique(treats)),lightScale=0,lightMultiple=.8,alpha=.7)
names(treatCols)<-unique(treats)
pdf('out/28bp.pdf',width=16)
  for(ii in unique(sampleNames)){
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

window<-100
pdf('out/prot_drugs.pdf',width=8)
  for(ii in unique(sampleNames)){
    theseTreats<-treats[sampleNames==ii]
    theseStarts<-bp28[sampleNames==ii]
    prots<-read.csv(protFiles[names(theseStarts)[1]])
    prots<-prots[order(prots$tss),]
    thisSplices<-read.table(spliceFiles[names(theseStarts)[1]])
    thisDonors<-unique(thisSplices$V2)
    thisAcceptors<-unique(thisSplices$V3)
    thisRef<-read.fa(refFiles[names(theseStarts)[1]])
    thisAUG<-gregexpr('ATG',thisRef$seq)[[1]]
    thisOneOff<-gregexpr('[CTG]TG|A[ACG]G|AT[CAT]',thisRef$seq)[[1]]
    for(jj in 1:nrow(prots)){
      targetCoords<-prots[jj,'tss']+-window:window
      startSubset<-lapply(theseStarts,function(x)x[targetCoords])
      nearestAcceptor<-max(c(-Inf,thisAcceptors[thisAcceptors<=prots[jj,'tss']]))
      nearestDonor<-min(c(Inf,thisDonors[thisDonors>=prots[jj,'tss']]))
      #if(nearestAcceptor<100)startSubset<-lapply(startSubset,function(x){x[1:(100-nearestAcceptor)]<-NA;return(x)})
      ylim<-range(unlist(startSubset)+1,na.rm=TRUE)
      plot(1,1,type='n',ylim=ylim,xlim=range(targetCoords),main=sprintf('%s %s',ii,prots[jj,'name']),xlab='Position',ylab='Read count',las=1)
      abline(v=prots[jj,'tss']-12,lty=2,col='#00000033')
      abline(v=prots[jj,'tss']+seq(-9,window,3),lty=3,col='#00000033')
      mapply(function(x,col)lines(targetCoords,x,col=col),startSubset,ranjitColors[theseTreats])
      if(nearestAcceptor>par('usr')[1])rect(par('usr')[1],par('usr')[3],nearestAcceptor,par('usr')[4],col='#00000055',border=NA)
      if(nearestDonor<par('usr')[2])rect(par('usr')[2],par('usr')[3],nearestDonor,par('usr')[4],col='#00000022',border=NA)
      for(tcl in c(.3))axis(1,thisOneOff-12,rep('',length(thisOneOff)),col='purple',lwd=0,tcl=tcl,lwd.tick=1)
      for(tcl in c(.5))axis(1,thisAUG-12,rep('',length(thisAUG)),col='orange',lwd=0,tcl=tcl,lwd.tick=1)
      legend('topright',legend=names(ranjitColors[names(ranjitColors)!='Total']),col=ranjitColors[names(ranjitColors)!='Total'],lty=1,inset=.01,bg='white')
    }
  }
dev.off()

pdf('out/leader.pdf',width=8)
  for(ii in unique(sampleNames)){
    theseTreats<-treats[sampleNames==ii]
    theseStarts<-bp28[sampleNames==ii]
    prots<-read.csv(protFiles[names(theseStarts)[1]])
    prots<-prots[order(prots$tss),]
    thisSplices<-read.table(spliceFiles[names(theseStarts)[1]])
    thisDonors<-unique(thisSplices$V2)
    thisAcceptors<-unique(thisSplices$V3)
    thisRef<-read.fa(refFiles[names(theseStarts)[1]])
    thisAUG<-gregexpr('ATG',thisRef$seq)[[1]]
    thisOneOff<-gregexpr('[CTG]TG|A[ACG]G|AT[CAT]',thisRef$seq)[[1]]
    targetCoords<-1:(min(prots$tss)+50)
    startSubset<-lapply(theseStarts,function(x)x[targetCoords])
    ylim<-range(unlist(startSubset)+1,na.rm=TRUE)
    plot(1,1,type='n',ylim=ylim/1000,xlim=range(targetCoords),main=sprintf('%s leader',ii),xlab='Position',ylab='Read count (1000x)',las=1)
    abline(v=prots[,'tss']-12,lty=2,col='#00000033')
    abline(v=thisDonors,lty=3,col='#00000033')
    mapply(function(x,col)lines(targetCoords,x/1000,col=col),startSubset,ranjitColors[theseTreats])
    for(tcl in c(.3))axis(1,thisOneOff-12,rep('',length(thisOneOff)),col='purple',lwd=0,tcl=tcl,lwd.tick=1)
    for(tcl in c(.5))axis(1,thisAUG-12,rep('',length(thisAUG)),col='orange',lwd=0,tcl=tcl,lwd.tick=1)
    legend('topright',legend=names(ranjitColors[names(ranjitColors)!='Total']),col=ranjitColors[names(ranjitColors)!='Total'],lty=1,inset=.01,bg='white')
  }
dev.off()



pdf('out/prot_ltmChx.pdf',width=10)
  for(ii in unique(sampleNames)){
    theseTreats<-treats[sampleNames==ii]
    theseStarts<-bp28[sampleNames==ii]
    theseLtm<-names(theseStarts)[grep('LTM',names(theseStarts))]
    theseChx<-sub('LTM','CHX',theseLtm)
    theseLtmChx<-mapply(function(x,y){
      maxN<-max(length(x),length(y))
      x<-x[1:maxN]
      y<-y[1:maxN]
      x[is.na(x)]<-0
      y[is.na(y)]<-0
      x/sum(x)-y/sum(y)
    },theseStarts[theseLtm],theseStarts[theseChx],SIMPLIFY=FALSE)
    prots<-read.csv(protFiles[names(theseStarts)[1]])
    prots<-prots[order(prots$tss),]
    for(jj in 1:nrow(prots)){
      startSubset<-lapply(theseLtmChx,function(x)x[prots[jj,'tss']+-100:99])
      ylim<-c(0,max(unlist(startSubset))*1000)
      plot(1,1,type='n',ylim=ylim,xlim=c(-100,99),main=sprintf('%s %s',ii,prots[jj,'name']),xlab='Position',ylab='Difference in LTM and CHX proportions (x1000)',las=1)
      abline(v=-12,lty=2,col='#00000033')
      mapply(function(x,col)lines(-100:99,1000*ifelse(x<0,0,x)),startSubset)
      legend('topright',legend='LTM-CHX',col='black',lty=1,inset=.01)
    }
  }
dev.off()


ltmChx<-lapply(unique(sampleNames),function(sample){
  ltms<-names(bp28)[sampleNames==sample&treats=='LTM']
  chxs<-sub('LTM','CHX',ltms)
  ltmChx<-mapply(function(ltm,chx){
    n<-max(length(ltm),length(chx))
    ltm<-ltm[1:n]
    chx<-chx[1:n]
    #ltm/movingStat(ltm,max,100)-chx/movingStat(chx,max,100)
    ltm/sum(ltm)-chx/sum(chx)
  },bp28[ltms],bp28[chxs])
  return(ltmChx)
})
names(ltmChx)<-unique(sampleNames)

pdf('out/virus_ltmChx.pdf',width=10)
  for(ii in names(ltmChx)){
    message(ii)
    thisLtmChx<-ltmChx[[ii]]
    thisLtmChx[is.na(thisLtmChx)]<-0
    thisLtmChx<-ifelse(thisLtmChx<.0001,.0001,thisLtmChx)
    ylim=range(thisLtmChx)
    xlim=c(1,nrow(thisLtmChx))
    plot(1,1,type='l',ylim=ylim,xlim=xlim,xlab='Genome position',ylab='Difference between LTM and CHX proportions',main=ii)
    thisProts<-read.csv(protFiles[colnames(thisLtmChx)[1]])
    thisProts<-thisProts[tolower(thisProts$name)!='pol',]
    thisSplices<-read.table(spliceFiles[colnames(thisLtmChx)[1]])
    thisDonors<-unique(thisSplices$V2)
    thisAcceptors<-unique(thisSplices$V3)
    abline(v=thisProts$tss-12,col='#FF000033')
    axis(1,thisDonors,rep('',length(thisDonors)),col='#0000FF')
    axis(1,thisAcceptors,rep('',length(thisAcceptors)),col='#00FF00')
    box()
    apply(thisLtmChx,2,lines,col='#00000099')
  }
dev.off()

#read.csv(protFiles[ltms][1])
#zz<-apply(ltmChx[[length(ltmChx)]]>0.001,2,which)
#(read.csv(protFiles[ltms][1])$tss-12) %in% intersect(zz[[1]],zz[[2]])

totals<-lapply(rnaStarts[grepl('Total',names(rnaStarts))],function(rs)apply(rs,2,sum))
gagCounts<-sapply(names(totals),function(name){
  prots<-read.csv(protFiles[name])
  gag<-prots[prots$name=='gag','tss']
  sum(totals[[name]][gag-100:100])
})
print(cbind('gag'=gagCounts,'total'=sapply(totals,sum),'propGagx1000'=round(gagCounts/sapply(totals,sum)*1000,3)))
pdf('out/totalVirusCoverage.pdf')
  for(ii in names(totals)){
    plot(totals[[ii]],type='l',main=ii,log='y',ylab='Total RNA start counts',xlab='Genome position')
  }
dev.off()



pdf('out/virusReps.pdf')
par(mfrow=c(2,2))
for(ii in unique(sampleNames)){
  for(jj in c('LTM','CHX','PatA','DMSO')){
      thisReps<-bp28[sampleNames==ii&treats==jj]
      maxN<-max(sapply(thisReps,length))
      thisReps<-lapply(thisReps,function(x)c(x,rep(0,maxN-length(x))))
      plot(thisReps[[1]],thisReps[[2]],main=sprintf('%s %s',ii,jj),xlab='Replicate 1',ylab='Replicate 2')
  }
}
dev.off()


startLike<-unlist(lapply(1:3,function(pos)lapply(c('A','T','C','G'),function(base){out<-'ATG';substring(out,pos,pos)<-base;out})))
startLike<-startLike[startLike!='ATG']

genomeFiles<-list.files('data/','fasta$',recursive=TRUE)
genomes<-sapply(sprintf('data/%s',genomeFiles),function(x)read.fa(x)$seq)
#names(genomes)<-sub('data/','',dirname(genomeFiles))
names(genomes)<- sapply(file.path('data',dirname(genomeFiles)),function(x)sub('_trim.fastq.gz','',sub('^[^_]+_','',list.files(x,'Total1.*fastq.gz'))))

genomeStarts<-lapply(genomes,function(dna){
  pos<-seq(1,nchar(dna)-2)
  nmers<-substring(dna,pos,pos+2)
  return(data.frame('pos'=pos,'nmers'=nmers,'start'=nmers=='ATG','startLike'=nmers %in% startLike,stringsAsFactors=FALSE))
})

viruses<-sub('^[^_]+_','',names(bp28))
treatments<-sub('[0-9]_.*$','',names(bp28))

drugSums<-tapply(bp28,list(viruses,treatments),function(counts){
  nBase<-max(sapply(counts,length))
  counts<-lapply(counts,function(x)c(x,rep(0,nBase-length(x))))
  apply(do.call(rbind,counts),2,sum)
})



startCuts<-lapply(names(genomeStarts),function(virus){
  thisTreats<-treatments[viruses==virus]
  thisData<-drugSums[virus,]
  nBases<-min(sapply(thisData,length))
  thisGenome<-genomeStarts[[virus]]
  #throw out early starts #a bit dangerous if something in first or last few bases is interesting 
  thisGenome<-thisGenome[thisGenome$pos>30&thisGenome$pos<nBases-30,]
  thisStarts<-thisGenome[thisGenome$start|thisGenome$startLike,]
  cuts<-lapply(thisStarts$pos,function(pos){
      plusMinus<--100:99
      selectPos<-pos+plusMinus
      posSelect<-selectPos>0&selectPos<nBases
      selectPos<-selectPos[posSelect]
      out<-do.call(rbind,lapply(thisData,function(x)x[selectPos]))
      colnames(out)<-plusMinus[posSelect]
      return(out)
  })
  names(cuts)<-paste(virus,thisStarts$pos,sep='_pos')
  return(cuts)
})
names(startCuts)<-names(genomeStarts)

drugData<-lapply(startCuts,function(counts){
  totals<-do.call(rbind,lapply(counts,function(x)apply(x,1,sum)))
  colnames(totals)<-sprintf('total_%s',colnames(totals))
  props<-lapply(counts,function(count)t(apply(count,1,function(x)x/ifelse(sum(x)==0,1,sum(x)))))
  propSubset<-do.call(rbind,lapply(props,function(prop){
    do.call(c,lapply(rownames(prop),function(name){
      out<-prop[name,as.character(-20:30)]
      names(out)<-sprintf('%s_%s',name,names(out))
      return(out)
    }))
  }))
  sum3<-do.call(rbind,lapply(props,function(count)apply(count[,as.character(seq(-9,30,3))],1,sum)))
  colnames(sum3)<-sprintf('sum3_%s',colnames(sum3))
  out<-cbind(propSubset,totals,sum3)
  return(out)
})
names(drugData)<-names(startCuts)


for(ii in names(drugData)){
  write.csv(drugData[[ii]],sprintf('out/dataVirus_%s.csv',ii))
}




