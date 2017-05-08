import pandas as pd
import matplotlib.pyplot as plt
import datetime
import pytz
import glob
import time


df = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/50hours.txt',sep='\t',usecols=['time','storm','ifis_id','gage_height','t0'])
constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
ids = constant_data.ifis_id

ind = df.index[df.storm==True]
df = df.loc[ind]
df.time = pd.DatetimeIndex(df.time)

count = 1
for id in ids:
    ind = df.index[df['ifis_id'] == id]

    ifis = df['gage_height'][ind]
    time = df.time[ind]
    nwm = df['t0'][ind]

    ind = ifis>0
    ifis = ifis[ind]
    time = time[ind]
    nwm = nwm[ind]


    fig, ax1 = plt.subplots(figsize=(15,7))
    ax1.plot(time,ifis, 'b-')
    ax1.set_xlabel('time')
    # Make the y-axis label, ticks and tick labels match the line color.
    ax1.set_ylabel('height', color='b')
    ax1.tick_params('y', colors='b')
    #ax1.plot(ind,filt.loc[ind],'g-')

    ax2 = ax1.twinx()
    ax2.plot(time,nwm, 'r-')
    ax2.set_ylabel('flow', color='r')
    ax2.tick_params('y', colors='r')

    fig.tight_layout()
    plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs with flows/%s' % (id),papertype='letter')
    plt.close()
    print('%i/%i' % (count,len(ids)))
    count = count+1
