library(parallel)
library(e1071)
library(dnar)
library(vipor)

if(!exists('features')){
  featureFiles<-list.files('out','data_.*.csv',full.names=TRUE)
  features<-lapply(featureFiles,read.csv,stringsAsFactors=FALSE,row.names=1,check.names=FALSE)
  names(features)<-sub('out/data_([^.]+).csv','\\1',featureFiles)
}

if(!exists('virus')){
  virusFiles<-list.files('out','dataVirus_.*.csv',full.names=TRUE)
  virus<-lapply(virusFiles,read.csv,stringsAsFactors=FALSE,row.names=1,check.names=FALSE)
  names(virus)<-sub('out/dataVirus_([^.]+).csv','\\1',virusFiles)
}
if(any(colnames(virus[[1]][1,])!=colnames(features[[1]][1,])))stop(simpleError('colnames mismatch'))

virus<-lapply(virus,function(x){x[,'ltmChx']<-x[,'LTM_-12']-x[,'CHX_-12'];x})
features<-lapply(features,function(x){x[,'ltmChx']<-x[,'LTM_-12']-x[,'CHX_-12'];x})

selectColumns<-colnames(virus[[1]])
selectColumns<-selectColumns[!grepl('total',selectColumns)]
selectColumns<-selectColumns[grepl('_-12|sum3|ltmChx',selectColumns)]

set.seed(1234)
trainTest<-lapply(features,function(x){
  #x[,'LTM_-12']-x[,'CHX_-12']>.03&
  selector<-x[,'total_Total']>1000&x[,'total_LTM']>100&x[,'total_CHX']>100 &x[,'LTM_-12']>.1 & x[,'CHX_-12']>.1# &x[,'ltmChx']>.03
  print(summary(selector))
  x<-x[selector,]
  trainIds<-sample(nrow(x),round(sum(nrow(x)*.5)))
  list('train'=x[trainIds,selectColumns],'test'=x[-trainIds,selectColumns])
})

svms<-mclapply(names(trainTest),function(datName){
  #message(datName)
  dat<-trainTest[[datName]]
  message('Training ',datName, ' model')
  model<-svm(dat[['train']],as.factor(rep('TIS',nrow(dat[['train']]))),type='one-classification',nu=0.1)
  #print(model)
  #summary(model)
  message('Predicting ',datName, ' test set')
  predicts<-predict(model,dat[['test']])
  message(datName, ' test set:')
  print(summary(predicts))
  message('Predicting ',datName, ' viral data')
  virusPredicts<-predict(model,virus[[datName]][grepl('_start|_oneOff',rownames(virus[[datName]])),selectColumns])
  message(datName, ' viral data:')
  print(summary(virusPredicts))
  return(list(model,predicts,virusPredicts))
},mc.cores=6)
names(svms)<-names(trainTest)
humanQuantileCut<-.25
for(ii in names(svms)){
  thisNames<-rownames(virus[[ii]])[svms[[ii]][[3]]]
  #print(thisNames)
  writeLines(thisNames,sprintf('out/svm_%s.txt',ii))
  thisNames<-names(svms[[ii]][[2]])[svms[[ii]][[2]]]
  writeLines(thisNames,sprintf('out/svmHost_%s.txt',ii))
  pdf(sprintf('out/svm/%s.pdf',ii))
    xRange<-range(c(virus[[ii]][,'LTM_-12']-virus[[ii]][,'CHX_-12'],trainTest[[ii]][['test']][,'LTM_-12']-trainTest[[ii]][['test']][,'CHX_-12']))
    yRange<-range(c(virus[[ii]][,'sum3_DMSO']+virus[[ii]][,'sum3_CHX'],trainTest[[ii]][['test']][,'sum3_DMSO']+trainTest[[ii]][['test']][,'sum3_CHX']))
    plot(virus[[ii]][,'ltmChx'],virus[[ii]][,'sum3_DMSO']+virus[[ii]][,'sum3_CHX'],col=1+svms[[ii]][[3]],xlim=xRange,ylim=yRange)
    plot(trainTest[[ii]][['test']][,'ltmChx'],trainTest[[ii]][['test']][,'sum3_DMSO']+trainTest[[ii]][['test']][,'sum3_CHX'],col=1+svms[[ii]][[2]],xlim=xRange,ylim=yRange)
  dev.off()
  pdf(sprintf('out/svm/ltm_%s.pdf',ii))
    humanQuant<-quantile(trainTest[[ii]][['train']][,'ltmChx'],humanQuantileCut)
    virusBreak<-virus[[ii]][,'ltmChx']>humanQuant#&virus[[ii]][,'sum3_DMSO']+virus[[ii]][,'sum3_CHX']>.1
    pos<-as.numeric(sub('.*_pos([0-9]+)_.*$','\\1',rownames(virus[[ii]])))
    plot(virus[[ii]][,'ltmChx'],pos,col=1+virusBreak,xlim=xRange,ylab='Position',xlab='Virus LTM - CHX',main=ii)
    abline(v=.05,lty=2)
    humanBreak<-trainTest[[ii]][['test']][,'ltmChx']>humanQuant #& trainTest[[ii]][['test']][,'sum3_DMSO']+trainTest[[ii]][['test']][,'sum3_CHX']>.1
    plot(trainTest[[ii]][['test']][,'ltmChx'],trainTest[[ii]][['test']][,'sum3_DMSO']+trainTest[[ii]][['test']][,'sum3_CHX'],col=1+humanBreak,xlim=xRange,ylim=yRange,xlab='Host LTM - CHX',ylab='3-base codon sum',main=ii)
    abline(v=.05,lty=2)
  dev.off()
  message(ii)
  thisNames<-rownames(virus[[ii]])[virusBreak]
  print(thisNames)
}

#hivFiles<-list.files('work/virusCounts/','txt$',full.name=TRUE)
#hivData<-lapply(hivFiles,function(x)as.numeric(readLines(x)))
#names(hivData)<-sub('.txt$','',basename(hivFiles))
#viruses<-sub('_[^_]+$','',names(hivData))
#drugs<-sub('^.*_','',names(hivData))

drugs<-unique(sub('_[0-9-]+','',colnames(virus[[1]])[grepl('_[0-9-]+$',colnames(virus[[1]]))]))
bestFeatures<-lapply(features,function(dat){
  selector<-dat[,'total_Total']>1000&dat[,'total_LTM']>100&dat[,'total_CHX']>100&dat[,'total_PatA']>100 &dat[,'total_DMSO']>100&dat[,'LTM_-12']>.06 & dat[,'CHX_-12']>.06# &x[,'ltmChx']>.03
  dat[selector,]
})
drugProfiles<-lapply(bestFeatures,function(dat){
  meanDrug<-lapply(drugs,function(drug){
    apply(dat[,grepl(sprintf('%s_[0-9-]+',drug),colnames(dat))],2,mean)
  })
  names(meanDrug)<-drugs
  return(meanDrug)
})
jsDists<-lapply(names(bestFeatures),function(strain){
  dat<-bestFeatures[[strain]]
  js<-do.call(cbind,lapply(drugs,function(drug){
    apply(dat[,grepl(sprintf('%s_[0-9-]+',drug),colnames(dat))],1,jensenShannon,drugProfiles[[strain]][[drug]])
  }))
  colnames(js)<-drugs
  virusDat<-virus[[strain]]
  jsVirus<-do.call(cbind,lapply(drugs,function(drug){
    apply(virusDat[,grepl(sprintf('%s_[0-9-]+',drug),colnames(virusDat))],1,jensenShannon,drugProfiles[[strain]][[drug]])
  }))
  colnames(jsVirus)<-drugs
  return(list('host'=js,'virus'=jsVirus))
})
names(jsDists)<-names(bestFeatures)

pdf('out/jsDist.pdf')
for(ii in names(jsDists)){
  dat<-jsDists[[ii]]
  virusDat<-virus[[ii]]
  drugNames<-rep(colnames(dat[['host']]),each=nrow(dat[['host']]))
  virusDrugNames<-sprintf("%s_%s",rep(colnames(dat[['virus']]),each=nrow(dat[['virus']])),rep(sub('^.*_','',rownames(dat[['virus']])),ncol(dat[['virus']])))
  hostQuant<-apply(dat[['host']],2,quantile,.75)
  #goodVirus<-apply(dat[['virus']][,'LTM',drop=FALSE],1,function(x)all(x<hostQuant['LTM'])) 
  goodVirus<-apply(dat[['virus']],1,function(x)all(x<hostQuant))&virusDat[,'LTM_-12']-virusDat[,'CHX_-12']>0
  par(mar=c(7,4,1,.1))
  pos<-vpPlot(c(drugNames,virusDrugNames),c(as.vector(dat[['host']]),as.vector(dat[['virus']])),ylab='Jensen Shannon divergence',main=ii,las=2,cex=.5)
  abline(v=seq(4.5,20,4))
  vPos<-pos[(length(drugNames)+1):length(pos)]
  points(vPos[rep(goodVirus,ncol(dat[['virus']]))],as.vector(dat[['virus']][goodVirus,]),col='red',cex=1.1)
  message(ii)
  mtext(paste(sub('.*_pos([0-9]+).*','\\1',rownames(dat[['virus']])[goodVirus]),collapse=', '),3,line=-1,col='red',cex=.8)
}
dev.off()


pdf('out/virusLtmChx.pdf',width=20)
  for(ii in names(virus)){
    dat<-virus[[ii]]
    ltmChx<-dat[,'LTM_-12']-dat[,'CHX_-12']
    ltmChx<-ifelse(ltmChx<0,0,ltmChx)
    plot(ltmChx,col=ifelse(grepl('_start$',rownames(dat)),'red',ifelse(grepl('_oneOff$',rownames(dat)),'blue','black')),main=ii)
    abline(h=.05,lty=2)
    plot(ltmChx,dat[,'total_LTM']*dat[,'LTM_-12']+1,log='y',col=ifelse(grepl('_start$',rownames(dat)),'red',ifelse(grepl('_oneOff$',rownames(dat)),'blue','black')))
    plot(ltmChx,dat[,'total_CHX']*dat[,'CHX_-12']+1,log='y',col=ifelse(grepl('_start$',rownames(dat)),'red',ifelse(grepl('_oneOff$',rownames(dat)),'blue','black')))
  }
dev.off()


#output ltm-chx
ltmChxOut<-lapply(names(virus),function(ii){
  pos<-as.numeric(sub('.*_pos([0-9]+)_.*$','\\1',rownames(virus[[ii]])))
  start<-sub('^.*_([^_]+)','\\1',rownames(virus[[ii]]))
  out<-data.frame('pos'=pos,'ltmMinusChx'=virus[[ii]][,'ltmChx'],'start'=start)
  out$passCutoff<-out[,'ltmMinusChx']>.05
  write.csv(out,sprintf('out/ltmChx/%s.csv',ii))
  return(out)
})
names(ltmChxOut)<-names(virus)


if(FALSE){
  ltm1<-do.call(rbind,lapply(goodCounts,function(x)x['LTM1_HIV_CH0236',]))
  ltm1[is.na(ltm1)]<-0
  ltmProp<-t(apply(ltm1,1,function(x)x/ifelse(sum(x)==0,1,sum(x))))
  colnames(ltmProp)<-c(-100:99)
  prData<-ltmProp[,as.character(-20:30)]
  prData<-cbind(
    prData,
    'sum3'=apply(prData[,as.character(seq(-9,30,3))],1,sum),
    'readCount'=apply(ltm1,1,sum)
  )
  xx<-prcomp(prData)
  pdf('test.pdf');plot(xx);dev.off()
  png('test.png');biplot(xx);dev.off()
}

