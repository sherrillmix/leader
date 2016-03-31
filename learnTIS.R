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

svms<-lapply(names(trainTest)[2],function(datName){
  message(datName)
  dat<-trainTest[[datName]]
  model<-svm(dat[['train']],as.factor(rep('TIS',nrow(dat[['train']]))),type='one-classification')
  print(model)
  summary(model)
  predicts<-predict(model,dat[['test']])
  print(summary(predicts))
  return(list(model,predicts))
})
