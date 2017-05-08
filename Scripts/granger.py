import pandas as pd
import numpy as np
import statsmodels.tsa.stattools as stattools

df = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/dec5.txt',sep='\t')
constant_data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_NWM_data.txt',sep='\t')
granger = pd.Series(index=constant_data.ifis_id)
count = 1
for id in constant_data.ifis_id:
    #id = constant_data.ifis_id[34]
    try:
        ind = df.index[df.ifis_id == id]
        x=np.array([df.gage_height[ind],df.t0[ind]])
        x = x.transpose()
        tests = stattools.grangercausalitytests(x,24,verbose=False)
        for i in range(1,25):
            #print(i)
            if tests[i][0]['lrtest'][1]<=.001:
                granger.loc[id] = 1
                break
            if tests[i][0]['params_ftest'][1]<=.001:
                granger.loc[id] = 1
                break
            if tests[i][0]['ssr_chi2test'][1]<=.001:
                granger.loc[id] = 1
                break
            if tests[i][0]['ssr_ftest'][1]<=.001:
                granger.loc[id] = 1
                break
            if i==24:
                granger.loc[id] = 0
    except:
        pass
    print('%s/%s' % (count,len(constant_data.ifis_id)))
    count = count+1
