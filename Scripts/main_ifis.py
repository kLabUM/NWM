# IFIS Data download
import pandas as pd
import csv

def main_ifis():
    f = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/IFIS_site_id.txt',sep='\t')
    ids = f.site_id
    ids = ids[1:len(ids)-1]

    for site in ids:
        filename = 'http://ifis.iowawis.org/ws/sites.php?site=%s&period=720&format=tab' % (site)
        f = pd.read_table(filename,sep='\t',skip_blank_lines=True,comment='#',names = ['timestamp (CST)','stage (ft)'])
        filename = '/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_Data/%s.txt' % (site)
        #f = f[len(f):None:-1]
        try:
            g = pd.read_table(filename,sep='\t')
            if len(g)==0:
                f = f[len(f):None:-1]
                new_lines = f.to_csv(sep='\t',index = False,header=False)
                with open(filename,'ab') as old_file:
                    old_file.write(new_lines)
            else:
                last_time = g.iloc[len(g)-1,0]
                for i in range(0,len(f)):
                    if last_time == f.iloc[i,0]:
                        x = i
                        break
                f = f[0:x]
                f = f[len(f):None:-1]
                new_lines = f.to_csv(sep='\t',index = False,header=False)
                with open(filename,'ab') as old_file:
                    old_file.write(new_lines)
        except:
            pass
            # old_file = open(filename,'a')
        # old_file.write(new_lines)
        # old_file.close()
