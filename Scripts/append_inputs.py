import pandas as pd
import generate_input_file
import datetime

#read old data
#find most recent time
#set that as start time
#make end time current time

f = pd.read_csv('inputs.txt',sep='\t')
start_time = pd.to_datetime(max(f.time)) + datetime.timedelta(hours=1)
end_time = datetime.datetime.now()

new_data = generate_input_file.generate_inputs(start_time,end_time)
new_lines = new_data.to_csv(sep='\t',index = False,header=False)

with open('inputs.txt','ab') as old_file:
            old_file.write(new_lines)

