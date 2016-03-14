library(parallel)
source('~/scripts/R/dna.R')
dataDir<-'data'
fastqs<-list.files(dataDir,'fastq.gz$',recursive=TRUE)
names(fastqs)<-sub('\\_trim.fastq.gz$','',basename(fastqs))
alignDir<-'work/align'
bams<-list.files(alignDir,'bam$',recursive=TRUE)
names(bams)<-sub('\\.bam$','',bams)
readCounts<-cacheOperation('work/readCounts.Rdat',sapply,fastqs,function(x)as.numeric(system(sprintf('zcat %s/%s|wc -l',dataDir,x),intern=TRUE))/4)
uniqReads<-cacheOperation('work/alignReadCounts.Rdat',mclapply,bams,function(x){
  outFile<-sprintf('%s/%s',alignDir,sub('bam$','uniq',x))
  if(!file.exists(outFile)){
    cmd<-sprintf('samtools view %s/%s|cut -f1|sort|uniq>%s',alignDir,x,outFile)
    message(cmd)
    system(cmd)
  }
  out<-as.numeric(strsplit(system(sprintf('wc -l %s',outFile),intern=TRUE),' ')[[1]][1])
  return(out)
},mc.cores=3)


if(FALSE){
if(!all(names(bams) %in% names(fastqs))||!all(names(fastqs) %in% names(bams)))stop(simpleError('Mismatch between bam and fastq'))
bamReadIds<-lapply(bams,function(bam){
  headBam<-system(sprintf('samtools view %s/%s|head -100000|cut -f1',alignDir,bam),intern=TRUE)
  uniqBam<-unique(sapply(strsplit(headBam,':'),'[[',1))
  return(uniqBam)
})
fastqReadIds<-cacheOperation('work/fastqIds.Rdat',lapply,fastqs[names(bams)],function(fastq){
  message(fastq)
  headFast<-system(sprintf('zcat %s/%s|grep "^@"',dataDir,fastq),intern=TRUE)
  uniqFast<-sub('^@','',unique(sapply(strsplit(headFast,':'),'[[',1)))
})

if(!all(sapply(bamReadIds,length)==1)||!all(sapply(fastqReadIds,length)==1)||any(unlist(bamReadIds)!=unlist(fastqReadIds)))stop(simpleError('Read ids do not match in fastq and bam'))
}
