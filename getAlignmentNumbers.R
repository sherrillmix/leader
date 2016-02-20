library(parallel)
dataDir<-'data'
fastq<-list.files(dataDir,'fastq.gz',recursive=TRUE)
readCounts<-cacheOperation('work/readCounts.Rdat',mclapply,fastq,function(x)as.numeric(system(sprintf('zgrep -c "^@" %s/%s',dataDir,x),intern=TRUE)))
