import urllib
import datetime
import zipfile
import os
import ftplib
from netCDF4 import Dataset
import pandas as pd

def download_NWM():
    ftp = ftplib.FTP('ftpprd.ncep.noaa.gov')
    ftp.login()
    now = datetime.datetime.today()
    yester = now - datetime.timedelta(days=1)
    path = 'pub/data/nccf/com/nwm/prod/nwm.%d%02d%02d/analysis_assim/' % (yester.year,yester.month,yester.day)
    ftp.cwd(path)
    f = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/nwm_stations.txt',sep='\t',header=None)
    ids = f[0]
    all_flow = pd.DataFrame([])
    for i in range(0,24):
        #filename = 'nwm.t%02dz.analysis_assim.channel_rt.tm00.conus.nc.gz' % (i)
        #print(filename)
        try:
            filename = 'nwm.t%02dz.analysis_assim.channel_rt.tm00.conus.nc' % (i)
            print(filename)
            ftp.retrbinary('RETR '+filename,open(filename,'wb').write)
            #os.system('gzip -d '+filename)
            temp = Dataset(filename)
            os.remove(filename)
            stations = temp.variables['feature_id']
            stations = list(stations[:])
            stations = pd.DataFrame(stations)
            ind = stations[stations.isin(list(ids))]
            ind = ind.dropna()
            ind = list(ind.index)

            flow = temp.variables['streamflow']
            flow = list(flow[ind])
            #flow = flow*.01
            all_flow['%d%02d%02d%02d' % (yester.year,yester.month,yester.day,i)] = pd.Series(flow,index = stations.iloc[ind])
        except:
            pass

    all_flow.to_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (yester.year,yester.month,yester.day),sep='\t',float_format='%.7f')
