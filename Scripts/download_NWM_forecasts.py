import urllib
import datetime
import zipfile
import os
import ftplib
from netCDF4 import Dataset
import pandas as pd
import numpy as np

#need to figure out how to store data in panels with items being stations, major axis being time, and minor axis being
#t+1 through t+15. then i append each item to a csv file over and over with index of time and columns of the t+1 to t+15

#def download_NWM_forecasts():
print(datetime.datetime.now())
f = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
ids = f['nwm_id']
g = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d.txt' % ids[0],sep='\t')
last_time = pd.to_datetime(g.iloc[len(g)-1,0])
now = datetime.datetime.now()
if last_time <= datetime.datetime(now.year,now.month,now.day-1): #if i missed a day of data or something
    start_time = datetime.datetime(now.year,now.month,now.day-1)
    end_time = datetime.datetime(year=start_time.year,month=start_time.month,day=start_time.day,hour=23)
elif last_time.hour == 23: #if the last time was the last of yesterday
    start_time = datetime.datetime(now.year,now.month,now.day)
    end_time = datetime.datetime.now() - datetime.timedelta(hours=2) #give it a couple hours to get files online
elif last_time >= datetime.datetime(now.year,now.month,now.day): #if the last time was sometime today
    start_time = last_time + datetime.timedelta(hours=1)
    end_time = datetime.datetime.now() - datetime.timedelta(hours=2)
else: #if the last time was sometime yesterday
    start_time = last_time + datetime.timedelta(hours=1)
    end_time = datetime.datetime(year=start_time.year,month=start_time.month,day=start_time.day+1)
ftp = ftplib.FTP('ftpprd.ncep.noaa.gov')
ftp.login()
#now = datetime.datetime.utcnow()
#yester = now - datetime.timedelta(days=1)
path = 'pub/data/nccf/com/nwm/prod/nwm.%d%02d%02d/short_range/' % (start_time.year,start_time.month,start_time.day)
ftp.cwd(path)
#all_flow = pd.Panel(items=ids,major_axis=range(0,24),minor_axis=range(1,16))
x = end_time-start_time
x = int(round(x.seconds/3600,0))
all_flow = np.ndarray(shape=(len(ids),x,15))
for i in range(0,x):
    for j in range(1,16):
        filename = 'nwm.t%02dz.short_range.channel_rt.f%03d.conus.nc.gz' % (start_time.hour+i,j)
        ftp.retrbinary('RETR '+filename,open(filename,'wb').write)
        os.system('gzip -d '+filename)
        filename = 'nwm.t%02dz.short_range.channel_rt.f%03d.conus.nc' % (start_time.hour+i,j)
        temp = Dataset(filename)
        os.remove(filename)

        stations = temp.variables['station_id']
#        stations = list(stations[:])
#        stations = pd.DataFrame(stations)
#        ind = stations[stations.isin(list(ids))]
#        ind = ind.dropna()
#        ind = list(ind.index)

        flow = temp.variables['streamflow']
        for k in range(0,len(ids)):
            all_flow[k,i,j-1] = flow[stations==ids[k]]

        print('%d,%d' % (start_time.hour+i,j))
times = pd.date_range('%d/%d/%d %02d' % (start_time.month,start_time.day,start_time.year,start_time.hour),periods=x,freq='H')
#all_flow = pd.Panel(data=all_flow,items=ids,major_axis=times,minor_axis=('t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12',
#                                                                         't13','t14','t15'))

for i in range(0,len(ids)):
    df = pd.DataFrame(data=all_flow[i,:,:],index=times,columns=('t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14','t15'))

    new_lines = df.to_csv(sep='\t',float_format='%.7f',header=False)

    with open('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d.txt' % ids[i],'ab') as old_file:
            old_file.write(new_lines)

print(datetime.datetime.now())

df = pd.DataFrame()
df.to_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/forecast_test.txt')
######
#need to figure out how to store data in panels
