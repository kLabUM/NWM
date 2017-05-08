import pandas as pd
from scipy import ndimage,signal
import pytz
import datetime
import matplotlib.pyplot as plt
import numpy as np

def changeIndex(x):
    x = x.replace("'","")
    x = x.replace("[","")
    x = x.replace("]","")
    return(x)

def get_times(ifis_data):
    #constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
    #ids = constant_data.ifis_id

    #for id in ids:
        #id = ids[0]
    #ifis_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id,sep='\t',index_col = 0)
    if len(ifis_data.index)>0:
        # ind = ifis_data.index[ifis_data['stage (ft)']>0]
        # ifis_data = ifis_data.loc[ind]
        #
        # ifis_data.index = pd.DatetimeIndex(ifis_data.index)
        # ifis_data = ifis_data.resample('H').mean()
        # ifis_data.index.tz = pytz.timezone('US/Central')
        # ifis_data.index = ifis_data.index.tz_convert('UTC')
        #
        # ind = ifis_data.index[ifis_data.index>=start_time and ifis_data.index<=end_time]
        # ifis_data = ifis_data.loc[ind]

        filt = ndimage.filters.gaussian_filter1d(ifis_data['stage (ft)'],10)

        deriv = np.diff(filt)
        deriv = np.append(deriv,0)

        filt = pd.DataFrame(filt,index=ifis_data.index)
        deriv = pd.DataFrame(deriv,index=ifis_data.index)

        ind=np.empty(1,dtype=int)
        i=0
        #std = np.nanstd(ifis_data.values)
        size = max(ifis_data.values) - min(ifis_data.values)
        while i < len(ifis_data.index):
            if deriv[0][i] > 0:
                j=i
                while deriv[0][j]>0:
                    j=j+1
                if (ifis_data.values[j]-ifis_data.values[i-1])>size*.15:
                    temp = list(range(i-(j-i),j+2*(j-i)))
                    ind = np.append(ind,temp)
            i = i+1

        ind = np.unique(ind)
        ind = ind[ind>0]
        ind = ind[ind<len(ifis_data.index)]
        ind = ifis_data.index[ind]

        return ind
    else:
        return []

        #plt.plot(ifis_data.index,filt)
        #plt.plot(ifis_data.index[0:len(ifis_data.index)-1],deriv)
        #plt.show()
        #
        # fig, ax1 = plt.subplots(figsize=(15,7))
        # ax1.plot(ifis_data.index,filt, 'b-')
        # ax1.set_xlabel('time')
        # # Make the y-axis label, ticks and tick labels match the line color.
        # ax1.set_ylabel('height', color='b')
        # ax1.tick_params('y', colors='b')
        # ax1.plot(ind,filt.loc[ind],'g-')
        #
        # ax2 = ax1.twinx()
        # ax2.plot(ifis_data.index,deriv, 'r-')
        # ax2.set_ylabel('derivative', color='r')
        # ax2.tick_params('y', colors='r')
        #
        # fig.tight_layout()
        # plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/hydrographs with derivatives/%s' % (id),papertype='letter')
        # plt.close()

nwm_start = 0
nwm_end = 50
ifis_start = 1
ifis_end = 50
start_time = datetime.datetime(2016,10,10)
end_time = datetime.datetime(2017,2,14)
output_file = '50hours'
constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
nwm_hours = range(nwm_start,nwm_end+1)
ifis_hours = range(ifis_start,ifis_end+1)
times = pd.DatetimeIndex(start=start_time,end=end_time,freq='H')
times.tz = pytz.timezone('UTC')




# need to edit this so that it yields two columns, one with the values, and one that's boolean to tell if it goes in
# the training data or not

ids = constant_data.ifis_id[0:220]
all_ifis = pd.DataFrame(index=times,columns=ids)
training = pd.DataFrame(index=times,columns=ids)
count = 0
for id in ids:
    ifis_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id,sep='\t',index_col = 0)
    ifis_data.index = pd.DatetimeIndex(ifis_data.index)
    ifis_data = ifis_data.resample('H').mean()
    ifis_data.index.tz = pytz.timezone('US/Central')
    ifis_data.index = ifis_data.index.tz_convert('UTC')

    ind = ifis_data.index.isin(times)
    ind = ifis_data.index[ind]
    all_ifis.loc[ind,id] = ifis_data.loc[ind,'stage (ft)']

    ifis_data = ifis_data.loc[ind]
    x = get_times(ifis_data)
    x = ifis_data.index.isin(x)
    training.loc[ind,id] = x



    #for time in times:
    #    if time in ifis_data.index:
    #        all_ifis.loc[time,id] = float(ifis_data.loc[time])
    print('generating ifis dataframe: %i/%i'%(count,len(ids)))
    count = count+1


cols=list(constant_data.columns)
cols.insert(0,'time')
cols.insert(1,'storm')
cols.append('gage_height')

for i in ifis_hours:
    cols.append('h%i'%i)

for i in nwm_hours:
    cols.append('t%i'%i)

dates = pd.DatetimeIndex(start=start_time,end=end_time,freq='d')
first_time = True
for date in dates:
    #pull in NWM files if they exist. if not, go to next loop
    #print(date)
    try:
        file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (date.year,date.month,date.day)
        nwm1 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
        x = date-datetime.timedelta(hours=nwm_end)
        a = date-datetime.timedelta(days=1)
        b = date-datetime.timedelta(days=2)
        c = date-datetime.timedelta(days=3)
        if x.day == date.day:
            pass
        elif (x.day == a.day):
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (date.year,date.month,date.day-1)
            nwm2 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
        elif (x.day == b.day):
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (date.year,date.month,date.day-1)
            nwm2 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (date.year,date.month,date.day-2)
            nwm3 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
        elif (x.day == c.day):
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (date.year,date.month,date.day-1)
            nwm2 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (date.year,date.month,date.day-2)
            nwm3 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (date.year,date.month,date.day-3)
            nwm4 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
    except:
        print('error, cannot read nwm files')
        continue

    times = pd.DatetimeIndex(start=date,end=datetime.datetime(date.year,date.month,date.day,hour=23),freq='H')
    all_data = pd.DataFrame(columns=cols)
    for time in times:
        for id in ids:
            #print(id)
            try:
                if time-datetime.timedelta(hours=ifis_end) in all_ifis.index:
                    x = np.ones(len(all_data.columns))
                    data = pd.DataFrame([x],columns=all_data.columns)
                    data.loc[0,'time'] = time
                    ind = constant_data.index[constant_data.ifis_id == id]
                    data.loc[0,constant_data.columns] = constant_data.loc[ind[0],constant_data.columns]
                    data.loc[0,'gage_height'] = float(all_ifis.loc[time,id])
                    data.loc[0,'storm'] = training.loc[time,id]
                    for i in ifis_hours:
                        data.loc[0,'h%i'%i] = float(all_ifis.loc[time-datetime.timedelta(hours=i),id])

                    for i in nwm_hours:
                        x = time - datetime.timedelta(hours=i)
                        a = date-datetime.timedelta(days=1)
                        b = date-datetime.timedelta(days=2)
                        c = date-datetime.timedelta(days=3)
                        if x.day == time.day:
                            data.loc[0,'t%i'%i] = float(nwm1['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[ind]])

                        elif (x.day == a.day):
                            data.loc[0,'t%i'%i] = float(nwm2['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[ind]])

                        elif (x.day == b.day):
                            data.loc[0,'t%i'%i] = float(nwm3['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[ind]])

                        elif (x.day == c.day):
                            data.loc[0,'t%i'%i] = float(nwm4['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[ind]])

                    if first_time:
                        data.to_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/%s.txt' % output_file,sep='\t',index=False)
                        first_time = False
                    else:
                        newlines = data.to_csv(sep='\t',index=False,header=False)
                        with open('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/%s.txt' % output_file,'a') as old_file:
                            old_file.write(newlines)
            except:
                pass

        print(time)

    try:
        del(nwm1)
        del(nwm2)
        del(nwm3)
        del(nwm4)
    except:
        pass

#create_data_file(0,5,1,5,datetime.datetime(2016,10,1),datetime.datetime(2016,11,1),'test')

