options(java.parameters = "-Xmx5g")
require(MASS)
require(AER)
require(tree)
require(mgcv)
require(earth)
require(BayesTree)
require(fBasics)
require(randomForest)
require(hydroGOF)
require(neuralnet)
setwd("~/Google Drive/Docs Kevin/National Water Model/NWM/Scripts")

df = read.csv("50hours.txt",sep='\t')
df2 = read.csv('IFIS_NWM_data_new.txt',sep='\t')
df = na.omit(df)
ind = df$gage_height>0
df = df[ind,]

df$ifis_id = factor(df$ifis_id,levels = levels(df2$ifis_id))
df$storm = as.logical(df$storm)
#test = as.logical(df$storm)
for (i in 1:length(df2$ifis_id)){
  a = which(df$ifis_id==df2$ifis_id[i])
  b = which(df$storm==TRUE)
  ind = intersect(a,b)
  if (df2$dynamics[i]<0.8){
    for (j in 1:length(ind)){
      df$storm[ind[j]] = FALSE
    }
  }
  print(i)
}


times = as.POSIXct(as.character(df$time),format='%Y-%m-%d %H:%M:%S')
#need to fix the input file so that hours show up in string. probably need to work with a datetime index in python
uniq_times = unique(times)
ids = unique(df$ifis_id)
all_ids = df$ifis_id
#df = df[,c(4:33)]

feb1 = as.POSIXct('02/01/2017 00:00:00',format='%m/%d/%Y %H:%M:%S')
a = which(times<feb1)
b = which(df$storm==TRUE)
train = intersect(a,b)
a = which(times>=feb1)
new_ids = df2$ifis_id[df2$dynamics>0.8]
n = length(new_ids)
b = which(df$ifis_id %in% new_ids)
bart_test = intersect(a,b)

NSE_6 = data.frame(linear = double(n),mars=double(n),bart=double(n),randforest = double(n))#,ann = double(n))
NSE_12 = data.frame(linear = double(n),mars=double(n),bart=double(n),randforest = double(n))#,ann = double(n))
NSE_18 = data.frame(linear = double(n),mars=double(n),bart=double(n),randforest = double(n))#,ann = double(n))
NSE_24 = data.frame(linear = double(n),mars=double(n),bart=double(n),randforest = double(n))#,ann = double(n))
NSE_30 = data.frame(linear = double(n),mars=double(n),bart=double(n),randforest = double(n))#,ann = double(n))
NSE_36 = data.frame(linear = double(n),mars=double(n),bart=double(n),randforest = double(n))#,ann = double(n))

formula1 = paste(colnames(df[,c(5,7,8,12,19:31,64:76)]),collapse='+')
formula1 = paste('gage_height~',formula1,sep='')
f6 = lm(formula1,data = df[train,])
mars = invisible(earth(df[train,c(5,7,8,12,19:31,64:76)],df$gage_height[train]))
bart1 = invisible(bart(df[train,c(5,7,8,12,19:31,64:76)],df$gage_height[train],sigest=1,nskip=100,ndpost=100,printevery=10,x.test = df[bart_test,c(5,7,8,12,19:31,64:76)]))
randforest = randomForest(df[train,c(5,7,8,12,19:31,64:76)],df$gage_height[train],importance=TRUE,keep.forest = TRUE,ntree = 200)
#ann = neuralnet(formula1,data=df[train,])
for (i in 1:n){
  b = which(all_ids==new_ids[i])
  test = intersect(a,b)
  ind = bart_test %in% test
  if (length(test)>0){
    observed = df$gage_height[test]
    # estimated = predict(f6,newdata=df[test,c(5:12,19:31,64:76)])
    # NSE_6$linear[i] = NSE(estimated,observed)
    
    estimated = predict(mars,newdata = df[test,],type = 'response')
    #NSE_6$mars[i] = NSE(estimated[,1],observed)
    
    # estimated = bart1$yhat.test.mean[ind]
    # NSE_6$bart[i] = NSE(estimated,observed)
    # 
    # estimated = predict(randforest,newdata = df[test,c(5,7,8,12,19:31,64:76)],type='response')
    # NSE_6$randforest[i] = NSE(estimated,observed)
    # 
    #estimated = compute(ann,df[test,])
    #NSE_6$ann = NSE(estimated,observed)
    
    timestamps = times[test]
    pdf(paste('~/Google Drive/Grad School/Fall 2016/IOE 691/Project/hydrographs/',new_ids[i],'_mars_6.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',main=paste('6 hour MARS',new_ids[i]))
    lines(timestamps,estimated[,1],col='red')
    dev.off()
  }
}
rm(f6)

formula1 = paste(colnames(df[,c(5,7,8,12,25:37,64:82)]),collapse='+')
formula1 = paste('gage_height~',formula1,sep='')
f12 = lm(formula1,data = df[train,])
mars = invisible(earth(df[train,c(5,7,8,12,25:37,64:82)],df$gage_height[train]))
bart1 = invisible(bart(df[train,c(5,7,8,12,25:37,64:82)],df$gage_height[train],sigest=1,nskip=100,ndpost=100,printevery=10,x.test = df[bart_test,c(5,7,8,12,25:37,64:82)]))
randforest = randomForest(df[train,c(5,7,8,12,25:37,64:82)],df$gage_height[train],importance=TRUE,keep.forest = TRUE,ntree = 200)
#ann = neuralnet(formula1,data=df[train,])
for (i in 1:n){
  b = which(all_ids==new_ids[i])
  test = intersect(a,b)
  ind = bart_test %in% test
  if (length(test)>0){
    observed = df$gage_height[test]
    estimated = predict(f12,newdata=df[test,c(5,7,8,12,25:37,64:82)])
    NSE_12$linear[i] = NSE(estimated,observed)
    
    estimated = predict(mars,newdata = df[test,],type = 'response')
    NSE_12$mars[i] = NSE(estimated[,],observed)
    
    estimated = bart1$yhat.test.mean[ind]
    NSE_12$bart[i] = NSE(estimated,observed)
    
    estimated = predict(randforest,newdata = df[test,c(5,7,8,12,25:37,64:82)],type='response')
    NSE_12$randforest[i] = NSE(estimated,observed)
    
    # estimated = compute(ann,df[test,])
    # NSE_12$ann = NSE(estimated,observed)
  }
}
rm(f12)

formula1 = paste(colnames(df[,c(5,7,8,12,31:43,64:88)]),collapse='+')
formula1 = paste('gage_height~',formula1,sep='')
f18 = lm(formula1,data = df[train,])
mars = invisible(earth(df[train,c(5,7,8,12,31:43,64:88)],df$gage_height[train]))
bart1 = invisible(bart(df[train,c(5,7,8,12,31:43,64:88)],df$gage_height[train],sigest=1,nskip=100,ndpost=100,printevery=10,x.test = df[bart_test,c(5,7,8,12,31:43,64:88)]))
randforest = randomForest(df[train,c(5,7,8,12,31:43,64:88)],df$gage_height[train],importance=TRUE,keep.forest = TRUE,ntree = 200)
#ann = neuralnet(formula1,data=df[train,])
for (i in 1:n){
  b = which(all_ids==new_ids[i])
  test = intersect(a,b)
  ind = bart_test %in% test
  if (length(test)>0){
    observed = df$gage_height[test]
    estimated = predict(f18,newdata=df[test,c(5,7,8,12,31:43,64:88)])
    NSE_18$linear[i] = NSE(estimated,observed)
    
    estimated = predict(mars,newdata = df[test,],type = 'response')
    NSE_18$mars[i] = NSE(estimated[,1],observed)
    
    estimated = bart1$yhat.test.mean[ind]
    NSE_18$bart[i] = NSE(estimated,observed)
    
    estimated = predict(randforest,newdata = df[test,c(5,7,8,12,31:43,64:88)],type='response')
    NSE_18$randforest[i] = NSE(estimated,observed)
    
    # estimated = compute(ann,df[test,])
    # NSE_18$ann = NSE(estimated,observed)
  }
}
rm(f18)

formula1 = paste(colnames(df[,c(5,7,8,12,37:49,64:94)]),collapse='+')
formula1 = paste('gage_height~',formula1,sep='')
f24 = lm(formula1,data = df[train,])
mars = invisible(earth(df[train,c(5,7,8,12,37:49,64:94)],df$gage_height[train]))
bart1 = invisible(bart(df[train,c(5,7,8,12,37:49,64:94)],df$gage_height[train],sigest=1,nskip=100,ndpost=100,printevery=10,x.test = df[bart_test,c(5,7,8,12,37:49,64:94)]))
randforest = randomForest(df[train,c(5,7,8,12,37:49,64:94)],df$gage_height[train],importance=TRUE,keep.forest = TRUE,ntree = 200)
# ann = neuralnet(formula1,data=df[train,])
for (i in 1:n){
  b = which(all_ids==new_ids[i])
  test = intersect(a,b)
  ind = bart_test %in% test
  if (length(test)>0){
    observed = df$gage_height[test]
    estimated = predict(f24,newdata=df[test,c(5,7,8,12,37:49,64:94)])
    NSE_24$linear[i] = NSE(estimated,observed)
    
    estimated = predict(mars,newdata = df[test,],type = 'response')
    NSE_24$mars[i] = NSE(estimated[,1],observed)
    
    estimated = bart1$yhat.test.mean[ind]
    NSE_24$bart[i] = NSE(estimated,observed)
    
    estimated = predict(randforest,newdata = df[test,c(5,7,8,12,37:49,64:94)],type='response')
    NSE_24$randforest[i] = NSE(estimated,observed)
    
    # estimated = compute(ann,df[test,])
    # NSE_24$ann = NSE(estimated,observed)
  }
}
rm(f24)

formula1 = paste(colnames(df[,c(5,7,8,12,43:55,64:100)]),collapse='+')
formula1 = paste('gage_height~',formula1,sep='')
f30 = lm(formula1,data = df[train,])
mars = invisible(earth(df[train,c(5,7,8,12,43:55,64:100)],df$gage_height[train]))
bart1 = invisible(bart(df[train,c(5,7,8,12,43:55,64:100)],df$gage_height[train],sigest=1,nskip=100,ndpost=100,printevery=10,x.test = df[bart_test,c(5,7,8,12,43:55,64:100)]))
randforest = randomForest(df[train,c(5,7,8,12,43:55,64:100)],df$gage_height[train],importance=TRUE,keep.forest = TRUE,ntree = 200)
# ann = neuralnet(formula1,data=df[train,])
for (i in 1:n){
  b = which(all_ids==new_ids[i])
  test = intersect(a,b)
  ind = bart_test %in% test
  if (length(test)>0){
    observed = df$gage_height[test]
    estimated = predict(f30,newdata=df[test,c(5,7,8,12,43:55,64:100)])
    NSE_30$linear[i] = NSE(estimated,observed)
    
    estimated = predict(mars,newdata = df[test,],type = 'response')
    NSE_30$mars[i] = NSE(estimated[,1],observed)
    
    estimated = bart1$yhat.test.mean[ind]
    NSE_30$bart[i] = NSE(estimated,observed)
    
    estimated = predict(randforest,newdata = df[test,c(5,7,8,12,43:55,64:100)],type='response')
    NSE_30$randforest[i] = NSE(estimated,observed)
    
    # estimated = compute(ann,df[test,])
    # NSE_30$ann = NSE(estimated,observed)
  }
}
rm(f30)

formula1 = paste(colnames(df[,c(5,7,8,12,49:61,64:106)]),collapse='+')
formula1 = paste('gage_height~',formula1,sep='')
f36 = lm(formula1,data = df[train,])
mars = invisible(earth(df[train,c(5,7,8,12,49:61,64:106)],df$gage_height[train]))
bart1 = invisible(bart(df[train,c(5,7,8,12,49:61,64:106)],df$gage_height[train],sigest=1,nskip=100,ndpost=100,printevery=10,x.test = df[bart_test,c(5,7,8,12,49:61,64:106)]))
randforest = randomForest(df[train,c(5,7,8,12,49:61,64:106)],df$gage_height[train],importance=TRUE,keep.forest = TRUE,ntree = 200)
#ann = neuralnet(formula1,data=df[train,])
for (i in 1:n){
  b = which(all_ids==new_ids[i])
  test = intersect(a,b)
  ind = bart_test %in% test
  if (length(test)>0){
    observed = df$gage_height[test]
    estimated = predict(f36,newdata=df[test,c(5,7,8,12,49:61,64:106)])
    NSE_36$linear[i] = NSE(estimated,observed)
    
    estimated = predict(mars,newdata = df[test,],type = 'response')
    NSE_36$mars[i] = NSE(estimated[,1],observed)
    
    estimated = bart1$yhat.test.mean[ind]
    NSE_36$bart[i] = NSE(estimated,observed)
    
    estimated = predict(randforest,newdata = df[test,c(5,7,8,12,49:61,64:106)],type='response')
    NSE_36$randforest[i] = NSE(estimated,observed)
    
    # estimated = compute(ann,df[test,])
    # NSE_36$ann = NSE(estimated,observed)
  }
}
rm(f36)

