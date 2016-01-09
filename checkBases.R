
library(parallel)
getBaseCounts<-function(bamFile){
  outFile<-sub('bam$','baseCount.gz',bamFile)
  message(bamFile,"->",outFile)
  cmd<-sprintf('~/.local/bin/countbases --strand %s|gzip>%s',bamFile,outFile)
  if(!file.exists(outFile))system(cmd)
  else message("Using cached ",outFile)
  output<-readLines(outFile)
  output[1]<-gsub('-','n',gsub('+','p',output[1]))
  con<-textConnection(output)
  out<-read.csv(con,stringsAsFactors=FALSE)
  return(out)
}

dataDir<-'work/virusAlign'
bulkRnaFiles<-list.files(dataDir,'Total.*bam$')
if(!exists('baseCounts'))print(system.time(baseCounts<-mclapply(sprintf("%s/%s",dataDir,bulkRnaFiles),getBaseCounts,mc.cores=10)))
