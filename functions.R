
getRegionStarts<-function(bamFile,regions,strand=c('+','-'),sizeRange=26:30){
  startEnd<-ifelse(all(strand=='-'),'end','start')
  counts<-lapply(regions,function(region){
    reg<-parseRegion(region)
    cmd<-sprintf('~/.local/bin/getstartends %s --region %s --maxGaps 3',bamFile,region)
    #output<-system(cmd,intern=TRUE)
    #inCon<-textConnection(output)
    out<-read.csv(pipe(cmd),comment.char='',colClasses=c('character','integer','integer','character'),header=TRUE)
    #close(inCon)
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

naFillRegs<-function(file,regs,naNs,...){
  outs<-lapply(regs,function(reg,...)do.call(c,getRegionStarts(file,reg,...)),...)
  out<-mapply(function(dat,n){
    if(n==0)return(dat)
    if(n>0)return(c(dat,rep(NA,n)))
    else return(c(rep(NA,-n),dat))
  },outs,naNs,SIMPLIFY=FALSE)
  return(out)
}

coordToReg<-function(coord){
  #may need to adjust start
  sprintf('%s:%d-%d',seqnames(coord),start(coord),end(coord))
}

getCdFives<-function(cds,fives){
  #some cds dont have 5' UTR (assuming junk)
  cds<-cds[names(fives)]
  selector<-grepl('^chr[0-9XYM]+$',sapply(fives,function(x)as.character(seqnames(x)@values)[1]))
  cds<-cds[selector]
  fives<-fives[selector]
  return(list('cds'=cds,'fives'=fives))
}

pullCdFiveRegion<-function(five,cd,targetFiles,windowWidth=40,...){
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
  #humanStarts<-do.call(rbind,lapply(sprintf("%s/%s",dataDir,targetFiles),function(x,revFunc,...)do.call(c,revFunc(getRegionStarts(x,regs,strand=strand,sizeRange=27:29))),revFunc=revFunc))
  naNs<-c(fiveN-windowWidth,windowWidth-cdN)
  humanStarts<-do.call(rbind,lapply(targetFiles,function(x,...)do.call(c,naFillRegs(x,regs,naNs,strand=strand,sizeRange=27:29))))
  rownames(humanStarts)<-names(targetFiles)
  return(humanStarts)
}

plotTreats<-function(meanProp){
  treats<-sub('[0-9]_.*$','',rownames(meanProp))
  treatCols<-rainbow.lab(length(unique(treats)),lightScale=0,lightMultiple=.8,alpha=.7)
  if(length(treatCols)==1)treatCols<-'black'
  names(treatCols)<-unique(treats)
  plot(1,1,type='n',xlim=c(-windowWidth+1,windowWidth),ylim=c(0,max(meanProp))*100,xlab='Offset from TIS',ylab='Mean percent of reads',las=1,mgp=c(2,.8,0))
  for(ii in 1:nrow(meanProp)){
    lines((-windowWidth+1):windowWidth,meanProp[ii,]*100,col=treatCols[treats[ii]],lwd=2)
  }
  #abline(v=1,lty=2)
  abline(v=-11,lty=2)
  abline(v=-11+seq(3,12+100,3),lty=3,col='#00000033')
  legend('topright',legend=names(treatCols),col=treatCols,lty=1,inset=.01,lwd=2,bg='white')
}

getMeanProp<-function(goodCounts,singleRowLimit=5){
  goodProp<-lapply(goodCounts,function(x)t(apply(x,1,function(y)if(sum(y,na.rm=TRUE)>=singleRowLimit)ifelse(is.na(y),0,y)/sum(y,na.rm=TRUE) else rep(NA,length(y)))))
  allProps<-do.call(abind,c(goodProp,list(along=3)))
  meanProp<-apply(allProps,c(1,2),mean,na.rm=TRUE)
  return(meanProp)
}

nullOmit<-function(x){
  x[!is.null(x)]
}
