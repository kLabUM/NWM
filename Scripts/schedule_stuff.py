#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
from main_ifis import main_ifis
#from append_full_inputs import append_full_inputs
from download_NWM import download_NWM
import pandas as pd
import datetime
import os

#data = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/inputs_full.txt',sep='\t')
y = datetime.datetime.now()
if os.path.isfile('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/NWM_Data/%d%02d%02d.txt' % (y.year,y.month,y.day-1)):
    pass
else:
    main_ifis()
    download_NWM()

    #append_full_inputs()

df = pd.DataFrame()
df.to_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/schedule_test.txt')
