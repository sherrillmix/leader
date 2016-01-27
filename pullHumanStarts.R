
getHuman<-function(bamFile,regions,strand=c('+','-'),sizeRange=26:30){
  counts<-lapply(regions,function(region){
    reg<-parseRegion(region)
    cmd<-sprintf('~/.local/bin/getstartends %s --region %s --maxGaps 3',bamFile,region)
    output<-system(cmd,intern=TRUE)
    out<-read.csv(textConnection(output))
    outTable<-unlist(table(out[out$end-out$start+1 %in% sizeRange & out$start>=reg$start&out$start<=reg$end & out$strand %in% strand,'start']))
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

regs<-c( "ATF1"="uc001rww.5", "CD4"="uc001qqv.3", "ACTB"="uc003sot.5", "GAPDH"="uc001qop.4", "GPI"="uc002nvg.3")
pdf('out/human.pdf',width=16)
for(geneName in names(regs)){ 
  message(geneName)
  ucscId<-regs[geneName]
  target<-knownGenes[knownGenes$V10==ucscId&knownGenes$V3=='start_codon',]
  gene<-knownGenes[knownGenes$V10==ucscId,]
  reg<-sprintf('%s:%d-%d',target$V1,ifelse(target$V7=='+',target$V4,target$V5)-100,ifelse(target$V7=='+',target$V4,target$V5)+100)
  targetFiles<-sub('_virus','',allRnaFiles[sample=='HIV_CH0694'])
  dataDir<-'work/align/'
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



