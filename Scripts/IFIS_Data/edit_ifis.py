import pandas as pd

f = pd.read_csv('IFIS_site_id.txt',sep='\t')
ids = f.site_id
#ids = ids[1:len(ids)-1]

for site in ids:
    filename = '%s.txt' % (site)
    f = pd.read_table(filename,sep='\t',skip_blank_lines=True,comment='#',index_col=0)
    test = f[len(f):None:-1]
    test.to_csv(filename,sep='\t',index=False)