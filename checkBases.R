getBaseCounts<-function(bamFile){
  cmd<-sprintf('countbases --strand %s',bamFile)
  output<-system(cmd,intern=TRUE)
  con<-textConnection(output)
  read.csv(con,stringsAsFactors=FALSE)
}

bulkRnaFiles<-list.files('work/virusAlign/','Total.*bam')
baseCounts<-lapply(bulkRnaFiles,getBaseCounts)
