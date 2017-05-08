import pandas as pd
import generate_inputs_with_height
import datetime

#read old data
#find most recent time
#set that as start time
#make end time current time
#def append_full_inputs():
f = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/inputs_with_heights.txt',sep='\t',parse_dates=False)
g = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/inputs_full.txt',sep='\t',parse_dates=False)
constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
ids = constant_data.ifis_id
time1 = pd.DatetimeIndex(f.time)
time2 = pd.DatetimeIndex(g.time)
start_time = max(time1) + datetime.timedelta(hours=1)
end_time = max(time2)

numhours = end_time - start_time
numhours = int(numhours.days*24+numhours.seconds/3600)
times = [start_time+datetime.timedelta(hours=x) for x in range(0,numhours)]

for time in times:
    ind = g.index[time2==time]
    newdata = g.loc[ind]
    newdata = newdata.assign(h1=pd.Series(),h2=pd.Series(),h3=pd.Series(),h4=pd.Series(),h5=pd.Series(),h6=pd.Series(),h7=pd.Series(),h8=pd.Series(),h9=pd.Series(),h10=pd.Series(),h11=pd.Series(),h12=pd.Series())
    newdata.index = range(0,len(newdata))
    for j in range(1,13):
        x = g.loc[g.index[time2==time-datetime.timedelta(hours=j)]]
        if len(x.ifis_id)==len(newdata.ifis_id):
            y=pd.Series(index=x.index)
            ind1 = x.ifis_id==newdata.ifis_id
            y[ind1] = x.gage_height[ind1]
            newdata['h%i'%j] = list(y)
        else:
            for i in range(0,len(newdata)):
                a = x.index[x.ifis_id==newdata.ifis_id[i]]
                if len(a)==1:
                    newdata.loc[i,'h%i'%j] = float(x.gage_height[a])
                else:
                    newdata.loc[i,'h%i'%j] = float('NaN')

    newdata = newdata[f.columns]
    new_lines = newdata.to_csv(sep='\t',index = False,header=False)

    with open('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/inputs_with_heights.txt','ab') as old_file:
                old_file.write(new_lines)

    print(time)


