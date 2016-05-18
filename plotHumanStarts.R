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

