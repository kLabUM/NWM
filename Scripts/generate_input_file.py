import pandas as pd
import datetime
import numpy as np
import csv
import pytz
import os

# function with a start time, end time, and list of IFIS ids
# columns: gage height, previous 12 hours, stream characteristics (total of 23 columns)

def changeIndex(x):
    x = x.replace("'","")
    x = x.replace("[","")
    x = x.replace("]","")
    return(x)

def generate_inputs(start_time,end_time,ids = []): #start and end times must be datetime objects
#start_time = datetime.datetime(2016,10,10)
#end_time = datetime.datetime(2016,11,14)
#ids = []
    constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
    numhours = end_time - start_time
    numhours = int(numhours.days*24+numhours.seconds/3600)
    times = [start_time+datetime.timedelta(hours=x) for x in range(0,numhours)]
    times = pd.to_datetime(times)
    if len(ids)==0:
        x = constant_data.ifis_id
        ids=[]
        for id in x:
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id
            check = os.path.getmtime(file)
            check = datetime.datetime(1970,1,1)+datetime.timedelta(seconds=check)
            if check>start_time-datetime.timedelta(days=1):
                ids.append(id)

    ifis_id = []
    nwm_id = []
    gage_height = []
    t0 = []
    t1 = []
    t2 = []
    t3 = []
    t4 = []
    t5 = []
    t6 = []
    t7 = []
    t8 = []
    t9 = []
    t10 = []
    t11 = []
    t12 = []
    order = []
    elevation = []
    mannings = []
    slope = []
    timestamp = []
    for time in times:
        try:
            file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (time.year,time.month,time.day)
            nwm1 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
            x = time-datetime.timedelta(hours=12)
            if x.day == time.day:
                pass
            else:
                file = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (time.year,time.month,time.day-1)
                nwm2 = pd.read_csv(file,sep='\t',converters={0:changeIndex},index_col=0)
            for id in ids:
                try:
                    flows = []
                    ind = constant_data.ifis_id.index[constant_data.ifis_id == id]
                    ifis_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id,sep='\t',index_col = 0)
                    ifis_data.index = pd.to_datetime(ifis_data.index)
                    ifis_data = ifis_data.resample('H').mean()
                    ifis_data.index.tz = pytz.timezone('US/Central')
                    ifis_data.index = ifis_data.index.tz_convert('UTC')
                    for i in range(0,13):
                        x = time - datetime.timedelta(hours=i)
                        if x.day == time.day:
                            flows.append(nwm1['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[np.asscalar(ind)]])
                        else:
                            flows.append(nwm2['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[np.asscalar(ind)]])
                    if ifis_data['stage (ft)'][time] < -10:
                        continue
                    gage_height.append(ifis_data['stage (ft)'][time])
                    t0.append(flows[0])
                    t1.append(flows[1])
                    t2.append(flows[2])
                    t3.append(flows[3])
                    t4.append(flows[4])
                    t5.append(flows[5])
                    t6.append(flows[6])
                    t7.append(flows[7])
                    t8.append(flows[8])
                    t9.append(flows[9])
                    t10.append(flows[10])
                    t11.append(flows[11])
                    t12.append(flows[12])
                    order.append(constant_data.order[np.asscalar(ind)])
                    elevation.append(constant_data.elevation[np.asscalar(ind)])
                    mannings.append(constant_data.mannings[np.asscalar(ind)])
                    slope.append(constant_data.slope[np.asscalar(ind)])
                    timestamp.append(time)
                    ifis_id.append(constant_data.ifis_id[np.asscalar(ind)])
                    nwm_id.append(constant_data.nwm_id[np.asscalar(ind)])

                except:
                    pass
        except:
            pass

        print(time)
    all_data = {'time':timestamp,'ifis_id':ifis_id,'nwm_id':nwm_id,'gage_height':gage_height,'t0':t0,'t1':t1,'t2':t2,'t3':t3,'t4':t4,'t5':t5,'t6':t6,'t7':t7,'t8':t8,'t9':t9,'t10':t10,'t11':t11,'t12':t12,'elevation':elevation,'mannings':mannings,'order':order,'slope':slope}
    #all_data = np.array([timestamp,ifis_id,nwm_id,gage_height,hr1,hr2,hr3,hr4,hr5,hr6,hr7,hr8,hr9,hr10,hr11,hr12,bottom_width,channel_slope,elevation,mannings,musk_coeff,order,rout_time,slope])
    df = pd.DataFrame.from_dict(all_data,orient='columns')
    df = df[['time','ifis_id','nwm_id','gage_height','t0','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','elevation','mannings','order','slope']]
    return df

#df.to_csv('inputs.txt',sep='\t',float_format='%.7f',index=False)

#generate_inputs(datetime.datetime(2016,11,9),datetime.datetime(2016,11,10))
