import pandas as pd
import generate_new_inputs
import datetime

#read old data
#find most recent time
#set that as start time
#make end time current time
#def append_full_inputs():
f = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/inputs_new.txt',sep='\t')
start_time = pd.to_datetime(max(f.time)) + datetime.timedelta(hours=1)
end_time = datetime.datetime.now()

new_data = generate_new_inputs.generate_new_inputs(start_time,end_time)
new_lines = new_data.to_csv(sep='\t',index = False,header=False)

with open('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/inputs_new.txt','ab') as old_file:
            old_file.write(new_lines)

