library(parallel)
library(e1071)
if(!exists('features')){
  featureFiles<-list.files('out','data_.*.csv',full.names=TRUE)
  features<-lapply(featureFiles,read.csv,stringsAsFactors=FALSE,row.names=1)
  names(features)<-sub('out/data_([^.]+).csv','\\1',featureFiles)
}

set.seed(1234)
trainTest<-lapply(features,function(x){
  trainIds<-sample(nrow(x),round(sum(nrow(x)*.5)))
  list('train'=x[trainIds,],'test'=x[-trainIds,])
})

svms<-mclapply(names(trainTest),function(datName){
  message(datName)
  dat<-trainTest[[datName]]
  model<-svm(dat[['train']],as.factor(rep('TIS',nrow(dat[['train']]))),type='one-classification')
  print(model)
  summary(model)
  predicts<-predict(model,dat[['test']])
  print(summary(predicts))
  return(list(model,predicts))
},mc.cores=6)

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

