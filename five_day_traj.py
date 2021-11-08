import numpy as np
from trajectory_grid import trajectory_grid
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
path = "/storage/silver/metstudent/msc/users_2021/yc810080/data/traj_data2/"
files = ["utraj-rf_cab2009101012", "utraj-rf_cab2010102512", "utraj-rf_cab2011100812", "utraj-rf_cab2012103112", "utraj-rf_cab2013102212", "utraj-rf_cab2014100712", "utraj-rf_cab2014113012", "utraj-rf_cab2016111312", "utraj-rf_cab2009101412", "utraj-rf_cab2010103012", "utraj-rf_cab2011102412", "utraj-rf_cab2012110512", "utraj-rf_cab2013102712", "utraj-rf_cab2014101712", "utraj-rf_cab2015102212", "utraj-rf_cab2016111812", "utraj-rf_cab2009102212", "utraj-rf_cab2010110512", "utraj-rf_cab2011103112", "utraj-rf_cab2012111912", "utraj-rf_cab2013102912", "utraj-rf_cab2014103112", "utraj-rf_cab2015103112", "utraj-rf_cab2016113012", "utraj-rf_cab2009110812", "utraj-rf_cab2010110812", "utraj-rf_cab2011111312", "utraj-rf_cab2012112112", "utraj-rf_cab2013110912", "utraj-rf_cab2014110612", "utraj-rf_cab2015112312", "utraj-rf_cab2009112212", "utraj-rf_cab2010111212", "utraj-rf_cab2011112512", "utraj-rf_cab2012112612", "utraj-rf_cab2013111712", "utraj-rf_cab2014110712", "utraj-rf_cab2015112612", "utraj-rf_cab2009112712", "utraj-rf_cab2010111812", "utraj-rf_cab2011113012", "utraj-rf_cab2012112912", "utraj-rf_cab2013111812", "utraj-rf_cab2014111812", "utraj-rf_cab2015113012", "utraj-rf_cab2009113012", "utraj-rf_cab2010113012", "utraj-rf_cab2012101912", "utraj-rf_cab2013101812", "utraj-rf_cab2013112712", "utraj-rf_cab2014112012", "utraj-rf_cab2016102312"]

files = ["utraj-rf_cab2009101012","utraj-rf_cab2009101412","utraj-rf_cab2009112712","utraj-rf_cab2010103012","utraj-rf_cab2010111212","utraj-rf_cab2011100812","utraj-rf_cab2011102412","utraj-rf_cab2011103112","utraj-rf_cab2011111312","utraj-rf_cab2012101912","utraj-rf_cab2012103112","utraj-rf_cab2012110512","utraj-rf_cab2012112112","utraj-rf_cab2012112912","utraj-rf_cab2013102912","utraj-rf_cab2013111712","utraj-rf_cab2013112712","utraj-rf_cab2014100712","utraj-rf_cab2014101712","utraj-rf_cab2014103112","utraj-rf_cab2014110612","utraj-rf_cab2014112012","utraj-rf_cab2015102212","utraj-rf_cab2015112612","utraj-rf_cab2016111312","utraj-rf_cab2016111812"]
print(len(files))


files = [f for f in files if f[-6:] != "113012"]
print(len(files))
fig,ax=plt.subplots(2,2)
ax=ax.flatten()


def modmean(array,mod,axis=0):
  a  = np.sin(2*np.pi*array/mod).mean(axis=axis)
  b  = np.cos(2*np.pi*array/mod).mean(axis=axis)
  r = mod*np.arctan2(a,b)/2/np.pi
  r[r>mod//2] -= mod
  return r

"""
for j,f in enumerate(files[:30]):
  print(f)
  traj=trajectory_grid(filename=path+f,shape=(20,34,41,61))[:,8,10:,20:40]
  traj.data["lon"][traj.data["lon"]>180] -= 360
  day = (traj.traj_basetime - traj.traj_basetime.replace(day=1,month=10)).days
  mask1 = (traj.data["lon"]<20)
  mask2 = (traj.data["lat"]<-25)
  ax[0].plot([traj.data["q"][0].mean()] ,[(mask1+mask2).mean(axis=(1,2))[-1]],"o",label = f[12:-2],c=get_cmap("viridis")(day/60))
  for i,var in enumerate(["lon","lat","p"]):
    traj.data[var] = np.ma.masked_array(traj.data[var])
    traj.data[var].mask += 1-(mask1+mask2).any(axis=0)[np.newaxis,:,:]
    if var == "lon":
      ax[i+1].plot([traj.data["q"][0].mean()] ,modmean(traj.data[var][-1].reshape(1,31*20),360,axis=1),"o",c=get_cmap("viridis")(day/60))
    else:
      ax[i+1].plot([traj.data["q"][0].mean()],[traj.data[var][-1].reshape(31*20).mean(axis=0)],"o", c=get_cmap("viridis")(day/60))

ax[3].set_ylim(850,500)
fig.legend(loc="right")
plt.show()
"""

for j,f in enumerate(files[:30]):
  print(f)
  traj=trajectory_grid(filename=path+f,shape=(20,34,41,61))[:,8,:,10:51]
  traj.data["lon"][traj.data["lon"]>180] -= 360
  traj.data["q"] = np.maximum(traj.data["q"],1e-9)
  traj.data["e"]= traj.data["p"]*100*traj.data["q"]/0.622
  traj.data["Td"]=(1/273.15-1.844*0.0001*np.log(traj.data["e"]/611.3))**-1
  dd = traj.data["T"][0]-traj.data["Td"][0] 
  mask = dd<8
  iy,ix=np.where(dd>=8)
  nt,ny,nx=traj.shape
  for key in traj.data.keys():
    traj.data[key] = np.ma.masked_array(traj.data[key],mask = np.ones(traj.shape)*mask)
  day = (traj.data["q"][0].mean()-0.0045)/0.0055
#  day = (traj.traj_basetime - traj.traj_basetime.replace(day=1,month=10)).days/60
  mask1 = (traj.data["lon"]<20)
  mask2 = (traj.data["lat"]<-25)
  ax[0].plot(np.arange(0,-5,-0.25),(mask1+mask2).mean(axis=(1,2)),label = f[12:-2],c=get_cmap("viridis")(day))
  for i,var in enumerate(["lon","lat","p"]):
    traj.data[var] = np.ma.masked_array(traj.data[var])
    traj.data[var].mask = traj.data[var].mask +  1-(mask1+mask2).any(axis=0)[np.newaxis,:,:]
    if var == "lon":
      ax[i+1].plot(np.arange(0,-5,-0.25),modmean(traj.data[var].reshape(20,31*20),360,axis=1),c=get_cmap("viridis")(day))
    else:
      ax[i+1].plot(np.arange(0,-5,-0.25),traj.data[var].reshape(20,31*20).mean(axis=1), c=get_cmap("viridis")(day))

ax[3].set_ylim(850,500)
fig.legend(loc="right")
plt.show()
exit()
"""
fig=plt.figure()
norm = plt.Normalize(100,1000)
for j,f in enumerate(files):
  print(f)
  traj=trajectory_grid(filename=path+f,shape=(20,34,41,61))[:,14,10::2,20:40:2]
  traj.data["lon"][traj.data["lon"]>180] -= 360
  print(traj)
  ax=plt.subplot(5,6,j+1,projection=ccrs.PlateCarree())
  ax.coastlines()
  #points = np.array([traj.data["lon"][:,0,0],traj.data["lat"][:,0,0]]).T.reshape(-1,1,2)
  segments = np.array([traj.data["lon"].reshape(20,16*10),traj.data["lat"].reshape(20,16*10)]).T#.reshape(-1,1,2)
#  segments = np.concatenate([points[:,:-1], points[:,1:]], axis=1)
  lc = LineCollection(segments, cmap='viridis', norm=norm)
  lc.set_array(traj.data["p"].T.flatten())
  lc.set_linewidth(0.8)
  line=ax.add_collection(lc)
  ax.set_xlim(traj.data["lon"].min(),traj.data["lon"].max())
  ax.set_ylim(traj.data["lat"].min(),traj.data["lat"].max())
  ax.set_title(f[12:-2])

plt.show()
"""
