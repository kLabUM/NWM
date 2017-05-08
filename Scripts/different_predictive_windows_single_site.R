require(quantregForest)
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
# might have to do this in case there are too many na's
# for (i in 1:length(df$gage_height)){
#   for (j in 5:16){
#     if (is.na(df[i,j])){
#       if (j==16){
#         df[i,j] = df[i,j-1]
#       }
#       else if (is.finite(df[i,j+1] && is.finite(df[i,j-1]))){
#         df[i,j] = (df[i,j-1]+df[i,j+1])/2
#       }
#       else if (is.na(df[i,j+1] && is.finite(df[i,j-1]))){
#         df[i,j] = df[i,j-1]
#       }
#     }
#   }
# }


times = as.POSIXct(as.character(df$time),format='%Y-%m-%d %H:%M:%S')
#need to fix the input file so that hours show up in string. probably need to work with a datetime index in python
uniq_times = unique(times)
ids = unique(df$ifis_id)
all_ids = df$ifis_id
#df = df[,c(4:33)]

feb1 = as.POSIXct('02/01/2017 00:00:00',format='%m/%d/%Y %H:%M:%S')

a = which(times>=feb1)
new_ids = df2$ifis_id[df2$dynamics>0.8]
for (i in 1:length(new_ids)){
  a = which(times<feb1)
  b = which(df$storm==TRUE)
  c = which(df$ifis_id==new_ids[i])
  train = intersect(a,b)
  train = intersect(train,c)
  a = which(times>=feb1)
  test = intersect(a,b)
  test = intersect(test,c)
  if (length(test)>0){
    f6_tree = quantregForest(df[train,c(19:31,64:76)],df$gage_height[train],nthreads=4,do.trace = TRUE,ntree=200,importance = TRUE,keep.forest = TRUE)
    f12_tree = quantregForest(df[train,c(25:37,64:82)],df$gage_height[train],nthreads=4,do.trace = TRUE,ntree=200,importance = TRUE,keep.forest = TRUE)
    f18_tree = quantregForest(df[train,c(31:43,64:88)],df$gage_height[train],nthreads=4,do.trace = TRUE,ntree=200,importance = TRUE,keep.forest = TRUE)
    f24_tree = quantregForest(df[train,c(37:49,64:94)],df$gage_height[train],nthreads=4,do.trace = TRUE,ntree=200,importance = TRUE,keep.forest = TRUE)
    f30_tree = quantregForest(df[train,c(43:55,64:100)],df$gage_height[train],nthreads=4,do.trace = TRUE,ntree=200,importance = TRUE,keep.forest = TRUE)
    f36_tree = quantregForest(df[train,c(49:61,64:106)],df$gage_height[train],nthreads=4,do.trace = TRUE,ntree=200,importance = TRUE,keep.forest = TRUE)
    #f50 = quantregForest(df[train,c(5:12,19:31,64:76)],df$gage_height[train],nthreads=4,do.trace = TRUE,ntree=200,importance = TRUE,keep.forest = TRUE)
    
    formula = paste(colnames(df[,c(19:31,64:76)]),collapse='+')
    formula = paste('gage_height~',formula,sep='')
    f6 = lm(formula,data = df[train,])
    
    formula = paste(colnames(df[,c(25:37,64:82)]),collapse='+')
    formula = paste('gage_height~',formula,sep='')
    f12 = lm(formula,data = df[train,])
    
    formula = paste(colnames(df[,c(31:43,64:88)]),collapse='+')
    formula = paste('gage_height~',formula,sep='')
    f18 = lm(formula,data = df[train,])
    
    formula = paste(colnames(df[,c(37:49,64:94)]),collapse='+')
    formula = paste('gage_height~',formula,sep='')
    f24 = lm(formula,data = df[train,])
    
    formula = paste(colnames(df[,c(43:55,64:100)]),collapse='+')
    formula = paste('gage_height~',formula,sep='')
    f30 = lm(formula,data = df[train,])
    
    formula = paste(colnames(df[,c(49:61,64:106)]),collapse='+')
    formula = paste('gage_height~',formula,sep='')
    f36 = lm(formula,data = df[train,])
    
    
    est6 = predict(f6,newdata=df[test,c(19:31,64:76)],what=.5)
    est12 = predict(f12,newdata=df[test,c(25:37,64:82)],what=.5)
    est18 = predict(f18,newdata=df[test,c(31:43,64:88)],what=.5)
    est24 = predict(f24,newdata=df[test,c(37:49,64:94)],what=.5)
    est30 = predict(f30,newdata=df[test,c(43:55,64:100)],what=.5)
    est36 = predict(f36,newdata=df[test,c(49:61,64:106)],what=.5)
    
    est6_tree = predict(f6_tree,newdata=df[test,c(19:31,64:76)],what=.5)
    est12_tree = predict(f12_tree,newdata=df[test,c(25:37,64:82)],what=.5)
    est18_tree = predict(f18_tree,newdata=df[test,c(31:43,64:88)],what=.5)
    est24_tree = predict(f24_tree,newdata=df[test,c(37:49,64:94)],what=.5)
    est30_tree = predict(f30_tree,newdata=df[test,c(43:55,64:100)],what=.5)
    est36_tree = predict(f36_tree,newdata=df[test,c(49:61,64:106)],what=.5)
    #est50
    observed = df$gage_height[test]
    
    timestamps = times[test]
    
    ymax = max(c(est6,observed))
    ymin = min(c(est6,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f6_linear.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('6 hour',ids[i]))
    lines(timestamps,est6,col='red')
    dev.off()
    ymax = max(c(est12,observed))
    ymin = min(c(est12,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f12_linear.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('12 hour',ids[i]))
    lines(timestamps,est12,col='red')
    dev.off()
    ymax = max(c(est18,observed))
    ymin = min(c(est18,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f18_linear.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('18 hour',ids[i]))
    lines(timestamps,est18,col='red')
    dev.off()
    ymax = max(c(est24,observed))
    ymin = min(c(est24,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f24_linear.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('24 hour',ids[i]))
    lines(timestamps,est24,col='red')
    dev.off()
    ymax = max(c(est30,observed))
    ymin = min(c(est30,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f30_linear.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('30 hour',ids[i]))
    lines(timestamps,est30,col='red')
    dev.off()
    ymax = max(c(est36,observed))
    ymin = min(c(est36,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f36_linear.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('36 hour',ids[i]))
    lines(timestamps,est36,col='red')
    dev.off()
    
    ymax = max(c(est6_tree,observed))
    ymin = min(c(est6_tree,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f6.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('6 hour',ids[i]))
    lines(timestamps,est6_tree,col='red')
    dev.off()
    ymax = max(c(est12_tree,observed))
    ymin = min(c(est12_tree,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f12.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('12 hour',ids[i]))
    lines(timestamps,est12_tree,col='red')
    dev.off()
    ymax = max(c(est18_tree,observed))
    ymin = min(c(est18_tree,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f18.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('18 hour',ids[i]))
    lines(timestamps,est18_tree,col='red')
    dev.off()
    ymax = max(c(est24_tree,observed))
    ymin = min(c(est24_tree,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f24.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('24 hour',ids[i]))
    lines(timestamps,est24_tree,col='red')
    dev.off()
    ymax = max(c(est30_tree,observed))
    ymin = min(c(est30_tree,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f30.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('30 hour',ids[i]))
    lines(timestamps,est30_tree,col='red')
    dev.off()
    ymax = max(c(est36_tree,observed))
    ymin = min(c(est36_tree,observed))
    pdf(paste('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs trained on single site/',ids[i],'_f36.pdf',sep=''))
    plot(timestamps,observed,col='black',type='l',ylim=c(ymin,ymax),main=paste('36 hour',ids[i]))
    lines(timestamps,est36_tree,col='red')
    dev.off()
    #lines(timestamps,est50,col='magenta')
    
    #legend('topright',c('observed','6hr','12hr','18hr','24hr','30hr','36hr','50hr'),col=c('black','red','blue','orange','green','purple','brown'))
  }
  
  print(paste(i,'/',length(new_ids),sep=''))
}