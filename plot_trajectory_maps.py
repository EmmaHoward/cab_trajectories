import numpy as np
import pickle
from trajectory_grid import trajectory_grid
from matplotlib.collections import LineCollection
import iris
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
path = "/storage/silver/metstudent/msc/users_2021/yc810080/data/traj_data2/"
files = ["utraj-rf_cab2009101012", "utraj-rf_cab2010102512", "utraj-rf_cab2011100812", "utraj-rf_cab2012103112", "utraj-rf_cab2013102212", "utraj-rf_cab2014100712", "utraj-rf_cab2014113012", "utraj-rf_cab2016111312", "utraj-rf_cab2009101412", "utraj-rf_cab2010103012", "utraj-rf_cab2011102412", "utraj-rf_cab2012110512", "utraj-rf_cab2013102712", "utraj-rf_cab2014101712", "utraj-rf_cab2015102212", "utraj-rf_cab2016111812", "utraj-rf_cab2009102212", "utraj-rf_cab2010110512", "utraj-rf_cab2011103112", "utraj-rf_cab2012111912", "utraj-rf_cab2013102912", "utraj-rf_cab2014103112", "utraj-rf_cab2015103112", "utraj-rf_cab2016113012", "utraj-rf_cab2009110812", "utraj-rf_cab2010110812", "utraj-rf_cab2011111312", "utraj-rf_cab2012112112", "utraj-rf_cab2013110912", "utraj-rf_cab2014110612", "utraj-rf_cab2015112312", "utraj-rf_cab2009112212", "utraj-rf_cab2010111212", "utraj-rf_cab2011112512", "utraj-rf_cab2012112612", "utraj-rf_cab2013111712", "utraj-rf_cab2014110712", "utraj-rf_cab2015112612", "utraj-rf_cab2009112712", "utraj-rf_cab2010111812", "utraj-rf_cab2011113012", "utraj-rf_cab2012112912", "utraj-rf_cab2013111812", "utraj-rf_cab2014111812", "utraj-rf_cab2015113012", "utraj-rf_cab2009113012", "utraj-rf_cab2010113012", "utraj-rf_cab2012101912", "utraj-rf_cab2013101812", "utraj-rf_cab2013112712", "utraj-rf_cab2014112012", "utraj-rf_cab2016102312"]

files = ["utraj-rf_cab2009101012","utraj-rf_cab2009101412","utraj-rf_cab2009112712","utraj-rf_cab2010103012","utraj-rf_cab2010111212","utraj-rf_cab2011100812","utraj-rf_cab2011102412","utraj-rf_cab2011103112","utraj-rf_cab2011111312","utraj-rf_cab2012101912","utraj-rf_cab2012103112","utraj-rf_cab2012110512","utraj-rf_cab2012112112","utraj-rf_cab2012112912","utraj-rf_cab2013102912","utraj-rf_cab2013111712","utraj-rf_cab2013112712","utraj-rf_cab2014100712","utraj-rf_cab2014101712","utraj-rf_cab2014103112","utraj-rf_cab2014110612","utraj-rf_cab2014112012","utraj-rf_cab2015102212","utraj-rf_cab2015112612","utraj-rf_cab2016111312","utraj-rf_cab2016111812"]
print(len(files))

def modmean(array,mod,axis=0):
  a  = np.sin(2*np.pi*array/mod).mean(axis=axis)
  b  = np.cos(2*np.pi*array/mod).mean(axis=axis)
  r = mod*np.arctan2(a,b)/2/np.pi
  r[r>mod//2] -= mod
  return r

keep  = pickle.load(open("prop_trajectories.pickle","rb"))
summary = {}
norm = plt.Normalize(100,1000)
cubes=iris.cube.CubeList()
for j,f in enumerate(files[:]):
  print(f)
  traj=trajectory_grid(filename=path+f,shape=(20,34,41,61))[:,8,::,10:51:]
  traj.data["lon"][traj.data["lon"]>180] -= 360
  traj.data["q"] = np.maximum(traj.data["q"],1e-9)
  traj.data["e"]= traj.data["p"]*100*traj.data["q"]/0.622
  traj.data["Td"]=(1/273.15-1.844*0.0001*np.log(traj.data["e"]/611.3))**-1
  dd = traj.data["T"][0]-traj.data["Td"][0] 
  mask = (dd>=8) #* (traj.data["lon"] <50).all(axis=0)
  summary[f] = {"time":traj.time}
  for key in traj.data.keys():
    summary[f][key] = (mask[np.newaxis]*traj.data[key]).sum(axis=(1,2))/mask[np.newaxis].sum()
  iy,ix=np.where(mask)
  nt,ny,nx=traj.shape
#  for key in traj.data.keys():
#    traj.data[key] = np.ma.masked_array(traj.data[key],mask = np.ones(traj.shape)*mask)
  ax=plt.subplot(5,6,j+1,projection=ccrs.PlateCarree())
  ax.coastlines()
  mask2 = keep[f[-10:]][2]
  nmask2 = (1-keep[f[-10:]][2]).astype(bool)
  segments1 = np.array([[traj.data["lon"][:,i,j] for (i,j) in zip(iy[mask2],ix[mask2]) if i%2==0 and j%2==0],
                      [traj.data["lat"][:,i,j] for (i,j) in zip(iy[mask2],ix[mask2])  if i%2==0 and j%2==0]]).transpose(1,2,0)
  segments2 = np.array([[traj.data["lon"][:,i,j] for (i,j) in zip(iy[nmask2],ix[nmask2])  if i%2==0 and j%2==0],
                      [traj.data["lat"][:,i,j] for (i,j) in zip(iy[nmask2],ix[nmask2]) if i%2==0 and j%2==0 ]]).transpose(1,2,0)
#  segments = np.ma.masked_array([traj.data["lon"].reshape(nt,ny*nx),traj.data["lat"].reshape(nt,ny*nx)]).T#.reshape(-1,1,2)
  #ax.plot(traj.data["lon"][0].flatten(),traj.data["lat"][0].flatten(),"o")
  ax.plot(segments2[:,:,0].T,segments2[:,:,1].T,"r",lw=0.4)
  ax.plot(segments1[:,:,0].T,segments1[:,:,1].T,"b",lw=0.4)
#  lc = LineCollection(segments, cmap='viridis', norm=norm)
#  lc.set_array(np.array([traj.data["p"][:,i,j] for (i,j) in zip(iy,ix)]).flatten())
#  lc.set_linewidth(0.3)
#  line=ax.add_collection(lc)
  ax.set_xlim(traj.data["lon"].min(),traj.data["lon"].max())
  ax.set_ylim(traj.data["lat"].min(),traj.data["lat"].max())
  ax.set_title(f[12:-2])


plt.show()
#pickle.dump(summary,open("mean_trajectories.pickle","wb"))
