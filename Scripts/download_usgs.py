# USGS data download

import pandas as pd

f = pd.read_csv('USGS_site_id.txt',sep='\t')
ids = f.site_id

for site in ids:
	filename = 'http://waterservices.usgs.gov/nwis/iv/?sites=%s&period=PT720H&format=rdb' % (site)
	#need to put something in about how to change the names of the param codes to be something understandable
	f = pd.read_table(filename,sep='\t',skip_blank_lines=True,comment='#')
	names = ['site_no','datetime (CST)']
	cols = f.columns
	for name in f.columns:
		if(name.find('00045')!=-1):
			if(name.find('00045_')==-1)
				cols.replace(name,'Precipitation (in)')
				names = names.append('Precipitation (in)')
		elif(name.find('00055')!=-1):
			if(nname.find('00055_')==-1):
				cols.replace(name,'Stream Velocity (ft/s)')
				names = names.append('Stream Velocity (ft/s)')
		elif(name.find('00060')!=-1):
			if(nname.find('00060_')==-1):
				cols.replace(name,'Discharge (ft^3/s)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('00061')!=-1):
			if(nname.find('00061_')==-1):
				cols.replace(name,'Instantaneous Discharge (ft^3/s)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('00065')!=-1):
			if(nname.find('00065_')==-1):
				cols.replace(name,'Gage height (ft)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('00072')!=-1):
			if(nname.find('00072_')==-1):
				cols.replace(name,'Stage (m)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('30208')!=-1):
			if(nname.find('30208_')==-1):
				cols.replace(name,'Discharge (m^3/s)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('30209')!=-1):
			if(nname.find('30209_')==-1):
				cols.replace(name,'Instantaneous Discharge (m^3/s)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('30212')!=-1):
			if(nname.find('30212_')==-1):
				cols.replace(name,'Stage, sonic range (ft)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('30213')!=-1):
			if(nname.find('30213_')==-1):
				cols.replace(name,'Stage, sonic range (m)')
				names = names.append('Discharge (ft^3/s)')
		if(name.find('46529')!=-1):
			if(name.find('46529_')==-1)
				cols.replace(name,'Precipitation (in)')
				names = names.append('Discharge (ft^3/s)')
		elif(name.find('72149')!=-1):
			if(nname.find('72149_')==-1):
				cols.replace(name,'Stream Velocity (m/s)')
				names = names.append('Discharge (ft^3/s)')

	f.columns = cols


	filename = '%s.txt' % (site)
	f.to_csv(filename,sep='\t')