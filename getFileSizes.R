new<-list.files('work/align/','\\.bam$',full.name=TRUE)
old<-list.files('work/alignBak','\\.bam$',full.name=TRUE)
names(new)<-basename(new)
names(old)<-basename(old)
if(any(names(new)!=names(old)))stop(simpleError('Deal with out of order names'))
check<-data.frame('file'=names(new),'path1'=new,'path2'=old,stringsAsFactors=FALSE)
check$size1<-sapply(check$path1,file.size)
check$size2<-sapply(check$path2,file.size)
round(sort(check$size1-check$size2)/1e6)




