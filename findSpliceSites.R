library('levenR')
source("~/scripts/R/dna.R")

#da: dataframe with columns pos (position of last exonic base before GT or first exonic base after AG from splice donor/acceptor) and type (acceptor or donoor)
findSplicePairings<-function(hiv896,da,target){
  window<-15
  daSurround<-substring(hiv896,da$pos-window,da$pos+window)

  aligns<-mapply(function(x,pos)levenAlign(x,substring(target,pos-1000,pos+1000),trimOuterGaps=TRUE,substring2=TRUE),daSurround,da$pos,SIMPLIFY=FALSE,USE.NAMES=FALSE)
  dists<-sapply(aligns,function(x)mean(strsplit(x[[1]],'')[[1]]==strsplit(x[[2]],'')[[1]]))
  alignStart<-mapply(function(align,pos)gap2NoGap(align[[1]],noGap2Gap(align[[2]],pos)),aligns,window+ifelse(da$type=='acceptor',-1,2))
  da$starts<-sapply(da$pos-1000,function(x)max(1,x))+sapply(aligns,attr,'start')-1+alignStart-1
  checkBases<-substring(target,da$starts,da$starts+1)
  goodDa<-da[(checkBases=="AG"&da$type=='acceptor')|(checkBases=="GT"&da$type=='donor'),]
  goodDa$firstLast<-goodDa$starts+ifelse(goodDa$type=='acceptor',1,0)

  nDa<-table(goodDa$type)
  pairings<-data.frame('donor'=rep(which(goodDa$type=='donor'),each=nDa['acceptor']),'acceptor'=rep(which(goodDa$type=='acceptor'),nDa['donor']))
  pairings<-pairings[goodDa[pairings$donor,'pos']<goodDa[pairings$acceptor,'pos'],]
  pairings$donor<-goodDa[pairings$donor,'firstLast']
  pairings$acceptor<-goodDa[pairings$acceptor,'firstLast']
  return(pairings)
}

hivFasta<-sapply(list.files('data','HIV',full.names=TRUE),list.files,'fasta',full.names=TRUE)
for(fasta in hivFasta){
  message(fasta)
  hiv896<-read.fa('hiv896/inferRef.fa')$seq
  target<-read.fa(fasta)
  da<-read.csv('hiv896/donorAcceptors.csv',stringsAsFactors=FALSE)
  da<-da[da$class %in% c('Good Cryptic','Cannonical'),]
  pairings<-findSplicePairings(hiv896,da,target$seq)
  pairings$chr<-target$name
  pairings$strand<-'+'
  pairings<-pairings[,c('chr','donor','acceptor','strand')]
  outFile<-sprintf('%s/sjdbFile.tsv',dirname(fasta))
  print(pairings)
  write.table(pairings,outFile,row.names=FALSE,quote=FALSE,col.names=FALSE)
}

sivFasta<-sapply(list.files('data','SIV',full.names=TRUE),list.files,'fasta',full.names=TRUE)
for(fasta in sivFasta){
  message(fasta)
  siv<-read.fa('data/SIV_766_WT_Hiseq/SIV_766_WT_R_polyA.fasta')$seq
  target<-read.fa(fasta)
  annots<-read.csv('SIV_766_WT_Splice_Locations.csv',stringsAsFactors=FALSE)
  sivDa<-annots[grep('S[DA][0-9]+',annots$Name),c('Name','Minimum','Maximum')]
  colnames(sivDa)<-c('name','start','end')
  sivDa$type<-ifelse(grepl('SD',sivDa$name),'donor','acceptor')
  sivDa$pos<-ifelse(sivDa$type=='donor',sivDa$end,sivDa$end+1)
  pairings<-findSplicePairings(siv,sivDa,target$seq)
  pairings$chr<-target$name
  pairings$strand<-'+'
  pairings<-pairings[,c('chr','donor','acceptor','strand')]
  outFile<-sprintf('%s/sjdbFile.tsv',dirname(fasta))
  print(pairings)
  write.table(pairings,outFile,row.names=FALSE,quote=FALSE,col.names=FALSE)
}



#--sjdbFileChrStartEnd 
#Chr tab Start tab End tab Strand=+/-   first and last bases of intron 1-based


#checkBases2<-do.call(rbind,mapply(function(align,pos){
  #start<-noGap2Gap(align[[2]],pos)
  #end<-noGap2Gap(align[[2]],pos+1)
  #seqs<-substring(unlist(align),start,end)
  #return(data.frame('ref'=seqs[1],'align'=seqs[2],start,end,stringsAsFactors=FALSE))
#},aligns,window+ifelse(da$type=='acceptor',-1,2),SIMPLIFY=FALSE))

#if(any(checkBases[da$type=='acceptor']!='AG')||any(checkBases[da$type=='donor']!='GT'))stop(simpleError('Problem finding splice sites'))


