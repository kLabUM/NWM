import pandas as pd
import pytz
import datetime
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage,signal

constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
ids = constant_data.ifis_id
count = 1
for id in ids:
    #id = ids[10]
    ifis_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id,sep='\t',index_col = 0)
    ifis_data.index = pd.DatetimeIndex(ifis_data.index)
    ifis_data = ifis_data.resample('H').mean()
    ifis_data.index.tz = pytz.timezone('US/Central')
    ifis_data.index = ifis_data.index.tz_convert('UTC')
    ind = ifis_data.index[ifis_data.index>=datetime.datetime(2016,10,1)]
    ifis_data = ifis_data.loc[ind]
    if len(ifis_data.index)>0 and len(ifis_data.index)<100000:
        corr = signal.correlate(ifis_data.values.squeeze(),ifis_data.values.squeeze())
        corr = corr[round(len(corr)/2):]
        plt.plot(corr[0:336])
        plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/autocorrelation/%s' % (id),papertype='letter')
        plt.close()

    print('%i/%i' % (count,len(ids)))
    count = count+1

