import pandas as pd
import csv

f = pd.read_csv('IFIS_site_id.txt',sep='\t')
ids = f.site_id
ids = ids[1:len(ids)-1]

filename = 'http://ifis.iowafloodcenter.org/ifis/ws/rgsm/sites.php'
f = pd.read_table(filename,sep=',',skip_blank_lines=True,comment='#',names = ['ifis_id','ifis_name','lat','long','time','temp1[F]','temp2[F]','temp3[F]','temp4[F]','water1[%]','water2[%]','water3[%]','water4[%]','rainfall[inch]'])

for site in ids:
