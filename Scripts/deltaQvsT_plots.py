import pandas as pd
from scipy import ndimage,signal,optimize
import pytz
import datetime
import matplotlib.pyplot as plt
import numpy as np

# def changeIndex(x):
#     x = x.replace("'","")
#     x = x.replace("[","")
#     x = x.replace("]","")
#     return(x)

#def get_times(ifis_data):
constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
ids = constant_data.ifis_id
count= 1
for id in ids:
    #id = ids[0]
    deltaH = np.empty(0)
    t = np.empty(0)
    ifis_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id,sep='\t',index_col = 0)
    if len(ifis_data.index)>0:
        ind = ifis_data.index[ifis_data['stage (ft)']>0]
        ifis_data = ifis_data.loc[ind]

        ifis_data.index = pd.DatetimeIndex(ifis_data.index)
        ifis_data = ifis_data.resample('H').mean()
        ifis_data.index.tz = pytz.timezone('US/Central')
        ifis_data.index = ifis_data.index.tz_convert('UTC')

        ind = ifis_data.index[ifis_data.index>=datetime.datetime(2016,10,1)]
        ifis_data = ifis_data.loc[ind]

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
                    k = i + 1
                    time = 1
                    sum = 0
                    while k <= j:
                        sum = sum+deriv.values[k]
                        deltaH = np.append(deltaH,sum)
                        t = np.append(t,time)
                        time = time + 1
                        k = k + 1
                i = j
            else:
                i = i + 1

        ind = deltaH>=0
        deltaH = deltaH[ind]
        t = t[ind]
        plt.scatter(t,deltaH)
        plt.xlabel('look forward time')
        plt.ylabel('change in height')
        #plt.show()

        plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/HvsT/%s' % (id),papertype='letter')
        plt.close()

    print('%i/%i' % (count,len(ids)))
    count = count + 1


#################################
constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
ids = constant_data.ifis_id
count= 1
for id in ids:
# id = ids[0]
    ifis_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % id,sep='\t',index_col = 0)
    storm = 1
    if len(ifis_data.index)>0:
        ind = ifis_data.index[ifis_data['stage (ft)']>0]
        ifis_data = ifis_data.loc[ind]

        ifis_data.index = pd.DatetimeIndex(ifis_data.index)
        ifis_data = ifis_data.resample('H').mean()
        ifis_data.index.tz = pytz.timezone('US/Central')
        ifis_data.index = ifis_data.index.tz_convert('UTC')

        ind = ifis_data.index[ifis_data.index>=datetime.datetime(2016,10,1)]
        ifis_data = ifis_data.loc[ind]

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
                deltaH = np.empty(0)
                t = np.empty(0)
                while deriv[0][j]>0:
                    j=j+1

                if (ifis_data.values[j]-ifis_data.values[i-1])>size*.15:
                    k = i
                    while k <= j:
                        m = k+1
                        time = 1
                        while m<=j:
                            deltaH = np.append(deltaH,float(filt.values[m] - filt.values[k]))
                            t = np.append(t,time)
                            time = time + 1
                            m = m + 1
                            #print(m)
                        k = k + 1
                    ind = deltaH>=0
                    deltaH = deltaH[ind]
                    t = t[ind]
                    plt.scatter(t,deltaH)
                    plt.xlabel('look forward time')
                    plt.ylabel('change in height')
                    #plt.show()

                    plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/HvsT/%s_density_%i' % (id,storm),papertype='letter')
                    plt.close()

                    storm = storm+1
                    #break
                i = j
            else:
                i = i + 1

        # ind = deltaH>=0
        # deltaH = deltaH[ind]
        # t = t[ind]
        # plt.scatter(t,deltaH)
        # plt.xlabel('look forward time')
        # plt.ylabel('change in height')
        # #plt.show()
        #
        # plt.savefig('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/HvsT/%s_density' % (id),papertype='letter')
        # plt.close()

    print('%i/%i' % (count,len(ids)))
    count = count + 1

##############################
def rating_fxn(Q,a,b,c):
    return (Q/c)**(1/b)+a

constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
ids = constant_data.ifis_id
all_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/test.txt',sep='\t',index_col=2)

#for id in ids:
id = ids[0]
heights = all_data['gage_height'].loc[id].values
Q = all_data['t0'].loc[id].values

ind = heights>0
heights = heights[ind]
Q = Q[ind]

fit = optimize.curve_fit(rating_fxn,Q,heights)

plt.scatter(np.log(Q),np.log(heights))
plt.close()
