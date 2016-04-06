library(parallel)
library(e1071)
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
  virusPredicts<-predict(model,virus[[datName]][,selectColumns])
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
    humanQuant<-quantile(trainTest[[ii]][['test']][,'ltmChx'],humanQuantileCut)
    virusBreak<-virus[[ii]][,'ltmChx']>humanQuant#&virus[[ii]][,'sum3_DMSO']+virus[[ii]][,'sum3_CHX']>.1
    plot(virus[[ii]][,'ltmChx'],virus[[ii]][,'sum3_DMSO']+virus[[ii]][,'sum3_CHX'],col=1+virusBreak,xlim=xRange,ylim=yRange)
    humanBreak<-trainTest[[ii]][['test']][,'ltmChx']>humanQuant #& trainTest[[ii]][['test']][,'sum3_DMSO']+trainTest[[ii]][['test']][,'sum3_CHX']>.1
    plot(trainTest[[ii]][['test']][,'ltmChx'],trainTest[[ii]][['test']][,'sum3_DMSO']+trainTest[[ii]][['test']][,'sum3_CHX'],col=1+humanBreak,xlim=xRange,ylim=yRange)
  dev.off()
  message(ii)
  thisNames<-rownames(virus[[ii]])[virusBreak]
  print(thisNames)
}

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

