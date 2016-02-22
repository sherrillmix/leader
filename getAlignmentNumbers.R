source('~/scripts/R/dna.R')
dataDir<-'data'
fastq<-list.files(dataDir,'fastq.gz$',recursive=TRUE)
alignDir<-'work/align'
bam<-list.files(alignDir,'bam$',recursive=TRUE)
readCounts<-cacheOperation('work/readCounts.Rdat',sapply,fastq,function(x)as.numeric(system(sprintf('zgrep -c "^@" %s/%s',dataDir,x),intern=TRUE)))
uniqReads<-cacheOperation('work/alignReadCounts.Rdat',sapply,bam,function(x){
  outFile<-sprintf('%s/%s',alignDir,sub('bam$','uniq',x))
  if(!file.exists(outFile)){
    cmd<-sprintf('samtools view %s/%s|cut -f1|sort|uniq>%s',alignDir,x,outFile)
    message(cmd)
    system(cmd)
  }
  out<-as.numeric(strsplit(system(sprintf('wc -l %s',outFile),intern=TRUE),' ')[[1]][1])
  return(out)
})
