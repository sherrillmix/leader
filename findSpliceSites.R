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


findLocations<-function(hiv896,pos,target,window=15){
  surrounds<-substring(hiv896,pos-window,pos+window)
  aligns<-mapply(function(x,pos)levenAlign(x,target,trimOuterGaps=TRUE,substring2=TRUE),surrounds,pos,SIMPLIFY=FALSE,USE.NAMES=FALSE)
  dists<-sapply(aligns,function(x)mean(strsplit(x[[1]],'')[[1]]==strsplit(x[[2]],'')[[1]]))
  alignStart<-mapply(function(align,pos)gap2NoGap(align[[1]],noGap2Gap(align[[2]],pos)),aligns,window+1)
  starts<-sapply(aligns,attr,'start')+alignStart-1
  return(data.frame('start'=starts,'dist'=dists,'ref'=sapply(aligns,'[[',1),'target'=sapply(aligns,'[[',1)))
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
  #work on start sites
  prots<-read.csv('hiv896/prots.csv',stringsAsFactors=FALSE)
  prots<-prots[!prots$isPolyProtein,]
  protStarts<-findLocations(hiv896,prots$start,target$seq,window=30)
  protStarts$codon<-substring(target$seq,protStarts$start,protStarts$start+2)
  prots$tss<-protStarts$start
  if(any(protStarts$codon!='ATG'&prots$name!='Pol'))stop(simpleError('Problem finding protein'))
  outFile<-sprintf('%s/prots.csv',dirname(fasta))
  write.csv(prots[,c('name','desc','tss')],outFile,row.names=FALSE)
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
  annot<-read.csv(list.files(dirname(fasta),'Annotations.csv$',full.names=TRUE),stringsAsFactors=FALSE)
  prots<-annot[grepl('START',annot$Name)|grepl('ORF',annot$Type),]
  prots$tss<-as.numeric(gsub(',','',prots$Minimum))
  prots<-prots[!grepl('exon2',prots$Name)&prots$Name!='Protease',]
  prots$codon<-substring(siv,prots$tss,prots$tss+2)
  if(any(protStarts$codon!='ATG'&prots$name!='Pol'))stop(simpleError('Problem finding protein'))
  outFile<-sprintf('%s/prots.csv',dirname(fasta))
  prots$name<-sub('_.*$','',tolower(prots$Name))
  prots<-prots[!duplicated(prots$name),]
  prots$desc<-prots$Name
  write.csv(prots[,c('name','desc','tss')],outFile,row.names=FALSE)
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


