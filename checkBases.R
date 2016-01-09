
library(parallel)
getBaseCounts<-function(bamFile){
  message(bamFile)
  outFile<-sub(bamFile,'bam$','baseCount')
  cmd<-sprintf('~/.local/bin/countbases --strand %s|gzip>%s',bamFile,outFile)
  if(!file.exists(outFile))system(cmd)
  else message("Using cached")
  output<-readLines(outFile)
  output[1]<-gsub('-','n',gsub('+','p',output[1]))
  con<-textConnection(output)
  out<-read.csv(con,stringsAsFactors=FALSE)
  return(out)
}

dataDir<-'work/virusAlign'
bulkRnaFiles<-list.files(dataDir,'Total.*bam$')
print(system.time(baseCounts<-mclapply(sprintf("%s/%s",dataDir,bulkRnaFiles),getBaseCounts,mc.cores=10)))
