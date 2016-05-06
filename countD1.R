if(!exists('allRnaFiles'))source('pullHivStarts.R')

total1s<-allRnaFiles[grep('^Total1',allRnaFiles)]

d1Cover<-cacheOperation('work/d1Cover.Rdat',mclapply,names(total1s),function(ii){
  juncs<-read.table(spliceFiles[ii],header=FALSE,stringsAsFactors=FALSE)
  reg<-pasteRegion(juncs[1,'V1'],1,min(juncs$V2)+400)
  ii2<-sub('Total1','Total2',ii)
  x<-pullRegion(reg,file.path(dataDir,allRnaFiles[c(ii,ii2)]),'~/installs/bedCount/bam2depth')
  return(x)
},mc.cores=4)
names(d1Cover)<-names(total1s)

pdf('out/D1.pdf');
for(ii in names(total1s)){
  message(ii)
  juncs<-read.table(spliceFiles[ii],header=FALSE,stringsAsFactors=FALSE)
  x<-d1Cover[[ii]
  plot(x$pos,x$counts1,type='l',ylim=range(x[,c("counts1","counts2")]),main=sub('Total1_','',ii),xlab='Position',ylab='Read coverage')
  lines(x$pos,x$counts2)
  abline(v=min(juncs$V2),lty=2)
}
dev.off()

totalReads<-cacheOperation('work/d1reads.Rdat',mclapply,names(total1s),function(ii){
  message(ii)
  juncs<-read.table(spliceFiles[ii],header=FALSE,stringsAsFactors=FALSE)
  d1<-min(juncs$V2)
  reg<-pasteRegion(juncs[1,'V1'],d1-50,d1+50)
  ii2<-sub('Total1','Total2',ii)
  out<-samView(file.path(dataDir,allRnaFiles[c(ii,ii2)]),sprintf('%s',reg))
  names(out)<-c(ii,ii2)
  return(out)
},mc.cores=4)
names(totalReads)<-names(total1s)
totalReads<-lapply(totalReads,function(x){
  colnames(x)<-c('qName','flag','tName','pos','cigar','seq','file')
  return(x)
})
upDownProp<-sapply(names(totalReads),function(ii){
  x<-totalReads[[ii]]
  juncs<-read.table(spliceFiles[ii],header=FALSE,stringsAsFactors=FALSE)
  d1<-min(juncs$V2)
  selector<-grepl('^[0-9]+M',x$cigar)
  nDistance<-5
  singleSplice<-grepl('^[0-9]+M[0-9]+N[0-9]+',x$cigar)
  nFirst<-as.numeric(sapply(strsplit(x$cigar,'[A-Z]'),'[[',1))
  #unspliced reads ending before d1 - nDistance
  upStream<-sum(selector&x$pos+nchar(x$seq)-1<d1-nDistance)
  #single splice starting before d1 - nDistance and splicing in right location
  upStreamSplice<-sum(abs(x$pos+nFirst-1-d1)<2 & singleSplice & pos<d1-nDistance)
  #reads starting after d1 + nDistance
  downStream<-sum(selector&x$pos>d1+nDistance)
  return(1-downStream/(upStream+upStreamSplice))
})

spliceUnspliceProp<-sapply(names(totalReads),function(ii){
  x<-totalReads[[ii]]
  juncs<-read.table(spliceFiles[ii],header=FALSE,stringsAsFactors=FALSE)
  d1<-min(juncs$V2)
  unspliced<-grepl('^[0-9]+M',x$cigar)
  nOverlap<-10
  readOver<-sum(unspliced&x$pos+nchar(x$seq)-1+nOverlap> d1 +nOverlap&x$pos<d1-nOverlap)
  singleSplice<-grepl('^[0-9]+M[0-9]+N[0-9]+',x$cigar)
  nFirst<-as.numeric(sapply(strsplit(x$cigar,'[A-Z]'),'[[',1))
  nLast<-as.numeric(sapply(strsplit(x$cigar,'[A-Z]'),tail,1))
  spliceOver<-sum(abs(x$pos+nFirst-1-d1)<2 & singleSplice&nFirst>nOverlap & nLast>nOverlap)
  return(1-spliceOver/readOver)
})

print(cbind(spliceUnspliceProp,upDownProp))

