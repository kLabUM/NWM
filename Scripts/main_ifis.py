# IFIS Data download
import pandas as pd

f = pd.read_csv('IFIS_site_id.txt',sep='\t')
ids = f.site_id

for site in ids:
	filename = 'http://ifis.iowawis.org/ws/sites.php?site=%s&period=720&format=tab' % (site)
	f = pd.read_table(filename,sep='\t',skip_blank_lines=True,comment='#',names = ['timestamp (CST)','stage (ft)'])
	filename = '%s.txt' % (site)
	f.to_csv(filename,sep='\t')