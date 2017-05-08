import pandas as pd
import matplotlib.pyplot as plt
import datetime
import pytz
import glob
import time

def changeIndex(x):
    x = x.replace("'","")
    x = x.replace("[","")
    x = x.replace("]","")
    return(x)
#data = pd.DataFrame()
ids = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_ids.txt',sep='\t')
files = glob.glob('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/201*.txt')

for i in range(0,len(ids)):
    #nwm = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
    data = pd.DataFrame()
    filename = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % (ids['ifis_id'][i])
    ifis_data = pd.read_csv(filename,sep='\t',index_col = 0)
    ifis_data.index = pd.to_datetime(ifis_data.index)
    ifis_data = ifis_data.resample('H').mean()
    ifis_data.index.tz = pytz.timezone('US/Central')
    ifis_data.index = ifis_data.index.tz_convert('UTC')
    for file in files:
        nwm = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
        nwm_data = nwm.loc[ids['NWM_id'][i]]
        nwm_data.index = pd.to_datetime(nwm_data.index,format='%Y%m%d%H')
        nwm_data.index.tz = pytz.timezone('UTC')
        try:
            x = ifis_data.loc[nwm_data.index]
            df = {'ifis': x['stage (ft)'],'nwm': nwm_data}
            df = pd.DataFrame(data = df,index = nwm_data.index)
            data = data.append(df)
        except:
            pass
        print(file,ids['ifis_id'][i])
    if len(data)>0:
        plt.scatter(data['nwm'],data['ifis'])
        time.sleep(.1)
        plt.xlim(max(0,data.nwm.quantile(.01)),data.nwm.quantile(.99))
        plt.ylim(max(0,data.ifis.quantile(.01)),data.ifis.quantile(.99))
        plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/figures/%s' % (ids['ifis_id'][i]))
        plt.close()

########################################################
#data = pd.DataFrame()
ids = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_ids.txt',sep='\t')
files = glob.glob('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/201*.txt')

for i in range(0,len(ids)):
    #nwm = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
    data = pd.DataFrame()
    filename = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts//IFIS_Data/%s.txt' % (ids['ifis_id'][i])
    ifis_data = pd.read_csv(filename,sep='\t',index_col = 0)
    ifis_data.index = pd.to_datetime(ifis_data.index)
    ifis_data = ifis_data.resample('H').mean()
    ifis_data.index.tz = pytz.timezone('US/Central')
    ifis_data.index = ifis_data.index.tz_convert('UTC')
    for file in files:
        nwm = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
        nwm_data = nwm.loc[ids['NWM_id'][i]]
        nwm_data.index = pd.to_datetime(nwm_data.index,format='%Y%m%d%H')
        nwm_data.index.tz = pytz.timezone('UTC')
        try:
            x = ifis_data.loc[nwm_data.index]
            df = {'ifis': x['stage (ft)'],'nwm': nwm_data}
            df = pd.DataFrame(data = df,index = nwm_data.index)
            data = data.append(df)
        except:
            pass
        print(file,ids['ifis_id'][i])

    if len(data)>0:
        fig, ax1 = plt.subplots(figsize=(15,7))
        ax1.plot(data.index,data['ifis'], 'b-')
        ax1.set_xlabel('time')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('height', color='b')
        ax1.tick_params('y', colors='b')
        #ax1.plot(ind,filt.loc[ind],'g-')

        ax2 = ax1.twinx()
        ax2.plot(data.index,data['nwm'], 'r-')
        ax2.set_ylabel('flow', color='r')
        ax2.tick_params('y', colors='r')

        fig.tight_layout()
        plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs with derivatives/%s' % (id),papertype='letter')
        plt.close()

# plt.scatter(data['ifis'],data['nwm'])
# plt.xlim(0,30)
# plt.show()
