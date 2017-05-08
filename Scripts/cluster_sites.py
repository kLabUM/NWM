import pandas as pd
import matplotlib.pyplot as plt
import datetime
import pytz
import glob
import time
from sklearn import cluster,metrics
import numpy as np

df = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/cluster_data.txt',sep='\t')
ids = df.ifis_id
# latlon = pd.read_csv('/Users/kjfries/Google Drive/Docs Kevin/National Water Model/NWM/Scripts/IFIS_site_id.txt',sep='\t',index_col=0)
# lat = pd.Series()
# lon = pd.Series()
# constant_data = constant_data.assign(lat=lat,lon=lon)
# for id in ids:
#     x = latlon.lat.loc[id]
#     y = latlon.lon.loc[id]
#     ind = constant_data.index[constant_data.ifis_id == id]
#     constant_data.lat.loc[ind] = x
#     constant_data.lon.loc[ind] = y

# standardize shit
df.covs_train = df.covs_train.fillna(-1)
df.lags_train = (df.lags_train - df.lags_train.mean())/df.lags_train.std()
df.covs_test = df.covs_test.fillna(-1)
df.lags_test = (df.lags_test - df.lags_test.mean())/df.lags_test.std()
df.bottom_width = (df.bottom_width - df.bottom_width.mean())/df.bottom_width.std()
df.elevation = (df.elevation - df.elevation.mean())/df.elevation.std()
df.mannings = (df.mannings - df.mannings.mean())/df.mannings.std()
df.tf21 = df.tf21/100
df.tf21 = df.tf21.fillna(-1)
df.tf21[df.tf21<-1] = -1

df.tf31 = df.tf31/100
df.tf31 = df.tf31.fillna(-1)
df.tf31[df.tf31<-1] = -1

df.tf32 = df.tf32/100
df.tf32 = df.tf32.fillna(-1)
df.tf32[df.tf32<-1] = -1

df.ss1 = df.ss1/100
df.ss1 = df.ss1.fillna(-1)
df.ss1[df.ss1<-1] = -1

df.ss2 = df.ss2/100
df.ss2 = df.ss2.fillna(-1)
df.ss2[df.ss2<-1] = -1

df.ss3 = df.ss3/100
df.ss3 = df.ss3.fillna(-1)
df.ss3[df.ss3<-1] = -1

df.ss4 = df.ss4/100
df.ss4 = df.ss4.fillna(-1)
df.ss4[df.ss4<-1] = -1

df.oe221 = df.oe221/100
df.oe221 = df.oe221.fillna(-1)
df.oe221[df.oe221<-1] = -1

df.oe211 = df.oe211/100
df.oe211 = df.oe211.fillna(-1)
df.oe211[df.oe211<-1] = -1

df.usgs_dist = (df.usgs_dist - df.usgs_dist.mean())/df.usgs_dist.std()

X = df.as_matrix()
cols = np.array(range(1,19))
cols = np.append(cols,21)
X = X[:,cols]
# X_norm = (X - X.min(0))/X.ptp(0)

A = metrics.pairwise_distances(X)
[aff_centers,aff_labels] = cluster.affinity_propagation(A)
[k_centroids,k_labels,inertia] = cluster.k_means(A,n_clusters=2)

ind1 = aff_labels==0
ind2 = aff_labels == 1

plt.plot(constant_data.lon[ind[ind1]],constant_data.lat[ind[ind1]],'bo')
plt.plot(constant_data.lon[ind[ind2]],constant_data.lat[ind[ind2]],'rx')

plt.show()

ind1 = k_labels==0
ind2 = k_labels == 1

plt.plot(constant_data.lon[ind1],constant_data.lat[ind1],'bo')
plt.plot(constant_data.lon[ind2],constant_data.lat[ind2],'rx')

plt.show()


# cluster only on "good" sites

ind = constant_data.index[constant_data.dynamics>=0.8]
X = constant_data.loc[ind].as_matrix()

X = X[:,np.array([2,4,5,7,9])]
X_norm = (X - X.min(0))/X.ptp(0)

A = metrics.pairwise_distances(X_norm)
[aff_centers,aff_labels] = cluster.affinity_propagation(A)
[k_centroids,k_labels,inertia] = cluster.k_means(A,n_clusters=2)

ind1 = aff_labels==0
ind2 = aff_labels == 1

plt.figure(1)
plt.plot(constant_data.lon[ind[ind1]],constant_data.lat[ind[ind1]],'bo')
plt.plot(constant_data.lon[ind[ind2]],constant_data.lat[ind[ind2]],'rx')

plt.show()

ind1 = k_labels==0
ind2 = k_labels == 1

plt.figure(2)
plt.plot(constant_data.lon[ind[ind1]],constant_data.lat[ind[ind1]],'bo')
plt.plot(constant_data.lon[ind[ind2]],constant_data.lat[ind[ind2]],'rx')

plt.show()

# Plotting with dynamics colorbar

cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(constant_data.lon,constant_data.lat,cmap=cm,c=constant_data.dynamics)
plt.colorbar(sc)
plt.show()


plt.figure(1)
plt.subplot(421)
plt.scatter(constant_data.bottom_width,constant_data.dynamics)
plt.xlabel('bottom width')
# plt.show()

plt.subplot(422)
plt.scatter(constant_data.channel_slope,constant_data.dynamics)
plt.xlabel('channel slope')
# plt.show()

plt.subplot(423)
plt.scatter(constant_data.elevation,constant_data.dynamics)
plt.xlabel('elevation')
# plt.show()

plt.subplot(424)
plt.scatter(constant_data.mannings,constant_data.dynamics)
plt.xlabel('mannings')
# plt.show()

plt.subplot(425)
plt.scatter(constant_data.musk_coeff,constant_data.dynamics)
plt.xlabel('muskingum coeff')
# plt.show()

plt.subplot(426)
plt.scatter(constant_data.order,constant_data.dynamics)
plt.xlabel('order')
# plt.show()

plt.subplot(427)
plt.scatter(constant_data.rout_time,constant_data.dynamics)
plt.xlabel('rout_time')
# plt.show()

plt.subplot(428)
plt.scatter(constant_data.slope,constant_data.dynamics)
plt.xlabel('slope')
plt.show()

