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

#    ifis_id = []
#    nwm_id = []
#    gage_height = []
#    t0 = []
#    t1 = []
#    t2 = []
#    t3 = []
#    t4 = []
#    t5 = []
#    t6 = []
#    t7 = []
#    t8 = []
#    t9 = []
#    t10 = []
#    t11 = []
#    t12 = []
#    h1 = []
#    h2 = []
#    h3 = []
#    h4 = []
#    h5 = []
#    h6 = []
#    h7 = []
#    h8 = []
#    h9 = []
#    h10 = []
#    h11 = []
#    h12 = []
#    order = []
#    elevation = []
#    mannings = []
#    slope = []
#    timestamp = []
    df = pd.DataFrame(columns=['time','ifis_id','nwm_id','gage_height','t0','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10',
                               't11','t12','h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','order',
                               'elevation','mannings','slope'])
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
                    ifis_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id,sep='\t',index_col = 0,parse_dates=True)
                    #ifis_data.index = pd.DatetimeIndex(ifis_data.index)
                    ifis_data = ifis_data.resample('H').mean()
                    #ifis_data.index.tz = pytz.timezone('US/Central')
                    #ifis_data.index = ifis_data.index.tz_convert('UTC')
                    for i in range(0,13):
                        x = time - datetime.timedelta(hours=i)
                        if x.day == time.day:
                            flows.append(nwm1['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[np.asscalar(ind)]])
                        else:
                            flows.append(nwm2['%d%02d%02d%02d' % (x.year,x.month,x.day,x.hour)].loc[constant_data.nwm_id[np.asscalar(ind)]])
                    if ifis_data['stage (ft)'][time] < -10:
                        continue
                    #gage_height.append(ifis_data['stage (ft)'][time])
                    #h1.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=1)])
                    #h2.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=2)])
                    #h3.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=3)])
                    #h4.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=4)])
                    #h5.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=5)])
                    #h6.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=6)])
                    #h7.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=7)])
                    #h8.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=8)])
                    #h9.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=9)])
                    #h10.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=10)])
                    #h11.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=11)])
                    #h12.append(ifis_data['stage (ft)'][time-datetime.timedelta(hours=12)])
                    #t0.append(flows[0])
                    #t1.append(flows[1])
                    #t2.append(flows[2])
                    #t3.append(flows[3])
                    #t4.append(flows[4])
                    #t5.append(flows[5])
                    #t6.append(flows[6])
                    #t7.append(flows[7])
                    #t8.append(flows[8])
                    #t9.append(flows[9])
                    #t10.append(flows[10])
                    #t11.append(flows[11])
                    #t12.append(flows[12])
                    #order.append(constant_data.order[np.asscalar(ind)])
                    #elevation.append(constant_data.elevation[np.asscalar(ind)])
                    #mannings.append(constant_data.mannings[np.asscalar(ind)])
                    #slope.append(constant_data.slope[np.asscalar(ind)])
                    #timestamp.append(time)
                    #ifis_id.append(constant_data.ifis_id[np.asscalar(ind)])
                    #nwm_id.append(constant_data.nwm_id[np.asscalar(ind)])

                    newdata = pd.DataFrame([[time,constant_data.ifis_id[np.asscalar(ind)],constant_data.nwm_id[np.asscalar(ind)],ifis_data['stage (ft)'][time],
                                            flows[0],flows[1],flows[2],flows[3],flows[4],flows[5],flows[6],flows[7],flows[8],flows[9],flows[10],flows[11],flows[12],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=1)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=2)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=3)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=4)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=5)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=6)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=7)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=8)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=9)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=10)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=11)],
                                            ifis_data['stage (ft)'][time-datetime.timedelta(hours=12)],
                                            constant_data.order[np.asscalar(ind)],
                                            constant_data.elevation[np.asscalar(ind)],
                                            constant_data.mannings[np.asscalar(ind)],
                                            constant_data.slope[np.asscalar(ind)]]],
                                           columns=['time','ifis_id','nwm_id','gage_height','t0','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','order','elevation','mannings','slope'])

                    df = df.append(newdata)

                except:
                    pass
        except:
            pass
        print(time)
    #all_data = {'time':timestamp,'ifis_id':ifis_id,'nwm_id':nwm_id,'gage_height':gage_height,'t0':t0,'t1':t1,'t2':t2,'t3':t3,'t4':t4,'t5':t5,'t6':t6,'t7':t7,'t8':t8,'t9':t9,'t10':t10,'t11':t11,'t12':t12,'elevation':elevation,'mannings':mannings,'order':order,'slope':slope,'h1':h1,'h2':h2,'h3':h3,'h4':h4,'h5':h5,'h6':h6,'h7':h7,'h8':h8,'h9':h9,'h10':h10,'h11':h11,'h12':h12}
    #all_data = np.array([timestamp,ifis_id,nwm_id,gage_height,hr1,hr2,hr3,hr4,hr5,hr6,hr7,hr8,hr9,hr10,hr11,hr12,bottom_width,channel_slope,elevation,mannings,musk_coeff,order,rout_time,slope])
    #df = pd.DataFrame.from_dict(all_data,orient='columns')
    #df = df[['time','ifis_id','nwm_id','gage_height','h1','h2','h3','h4','h5','h6','h7','h8','h9','h10','h11','h12','t0','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','elevation','mannings','order','slope']]
    return df

#df.to_csv('inputs.txt',sep='\t',float_format='%.7f',index=False)

#generate_inputs(datetime.datetime(2016,11,9),datetime.datetime(2016,11,10))
