setwd("~/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/")
require(pspline)
require(signal)

ifis = read.csv('IFIS_site_id',header = FALSE)
ids = ifis$V1
for (i in 1:length(ids)){
  check = read.csv(paste('./IFIS_Data/',ids[i],'.txt',sep = ''),sep='\t',nrows = 2)
  if (length(check[1,])!=2){
    next
  }
  df = read.csv(paste('./IFIS_Data/',ids[i],'.txt',sep = ''),sep='\t',col.names = c('time','stage'))
  
  times = as.POSIXct(df$time,format='%Y-%m-%d %H:%M:%S')
  ind = times>as.POSIXct('10/10/2016 00:00:00',format='%m/%d/%Y %H:%M:%S')
  times = times[ind]
  df = df[ind,]
  ind = df$stage>0
  df = df[ind,]
  times = times[ind]
  
  # need to resample to hourly before applying filter to get rid of noise at beginning of timeseries
  
  
  bf = butter(2,.03)
  new_data = filtfilt(bf,df$stage)
  smooth = sm.spline(times,new_data)
  deriv = predict(smooth,times,1)
  smooth_deriv = sm.spline(times,deriv)
  pdf(paste('~/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs with derivatives/',ids[i],'.pdf',sep=''))
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(smooth$x,smooth$ysmth,type='l',col='red',xlab='time',ylab='smoothed stage')
  par(new=TRUE)
  plot(smooth_deriv$x,smooth_deriv$ysmth,type='l',col='blue',axes=FALSE,ylab='',xlab='')
  axis(side=4,at=pretty(range(deriv)))
  mtext('smoothed derivative',side=4)
  dev.off()
  print(paste(i,'/',length(ids),sep=''))
}
