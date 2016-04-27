if(!exists('allRnaFiles'))source('pullHivStarts.R')

total1s<-allRnaFiles[grep('^Total1',allRnaFiles)]

pdf('test.pdf');
for(ii in names(total1s)){
  message(ii)
  juncs<-read.table(spliceFiles[ii],header=FALSE,stringsAsFactors=FALSE)
  reg<-pasteRegion(juncs[1,'V1'],1,min(juncs$V2)+400)
  ii2<-sub('Total1','Total2',ii)
  x<-pullRegion(reg,file.path(dataDir,allRnaFiles[c(ii,ii2)]),'~/installs/bedCount/bam2depth')
  plot(x$pos,x$counts1,type='l',ylim=range(x[,c("counts1","counts2")]),main=sub('Total1_','',ii),xlab='Position',ylab='Read coverage')
  lines(x$pos,x$counts2)
  abline(v=min(juncs$V2),lty=2)
}
dev.off()

totalReads<-lapply(names(total1s),function(ii){
  message(ii)
  juncs<-read.table(spliceFiles[ii],header=FALSE,stringsAsFactors=FALSE)
  d1<-min(juncs$V2)
  reg<-pasteRegion(juncs[1,'V1'],d1-1,d1+1)
  ii2<-sub('Total1','Total2',ii)
  out<-samView(file.path(dataDir,allRnaFiles[c(ii,ii2)]),sprintf('%s',reg))
  names(out)<-c(ii,ii2)
  return(out)
})
