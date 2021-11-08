import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
from trajectory_grid import trajectory_grid
from matplotlib.collections import LineCollection
import datetime as dt
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import iris
import cartopy.crs as ccrs
from iris.analysis.calculus import differentiate
from iris.coord_categorisation import add_day_of_year
import pickle
norm = plt.Normalize(-50,0)
path = "/storage/silver/NCAS-Weather/oa921437/tmp/"

tpath = "/storage/silver/metstudent/msc/users_2021/yc810080/data/traj_data2/"
tpath3 = "/storage/silver/metstudent/msc/users_2021/yc810080/data/traj_data3/"

files = ["2009101012","2009101412","2009112712","2010103012","2010111212","2011100812","2011102412","2011103112","2011111312","2012101912","2012103112","2012110512","2012112112","2012112912","2013102912","2013111712","2013112712","2014100712","2014101712","2014103112","2014110612","2014112012","2015102212","2015112612","2016111312","2016111812"]

files3 = ["2012103112","2013111712","2014103112","2016111312"]
   
keep  = pickle.load(open("prop_trajectories.pickle","rb"))
traj  = [dt.datetime(int(f[:4]),int(f[4:6]),int(f[6:8])) for f in files]

array = np.array([keep[x][2].mean() for x in keep])
order = array.argsort()
ranks = order.argsort()

dewpoint = []
for i,f in enumerate(files):
      print(f)
      t = dt.datetime(int(f[:4]),int(f[4:6]),int(f[6:8]))
#      s = summaries["utraj-rf_cab"+f]
#      plt.plot(s["lon"], [(tt-dt.datetime(year,10,1)).total_seconds()/3600/24 for tt in s["time"]],"k")
      traj=trajectory_grid(filename=tpath+"utraj-rf_cab"+f,shape=(20,34,41,61))[:,8]
      nt,ny,nx = traj.shape
      traj.data["lon"][traj.data["lon"]>180] -= 360
      traj.data["q"] = np.maximum(traj.data["q"],1e-9)
      traj.data["e"]= traj.data["p"]*100*traj.data["q"]/0.622
      traj.data["Td"]=(1/273.15-1.844*0.0001*np.log(traj.data["e"]/611.3))**-1
      dd = traj.data["T"][0]-traj.data["Td"][0] 
      mask2 = (dd>=8) #* (traj.data["lon"] <50).all(axis=0)
      dewpoint.append((dd*mask2).sum()/mask2.sum())
      ax=plt.subplot(4,7,i+1,projection=ccrs.PlateCarree())
      #ax=plt.subplot(4,7,ranks[i]+1,projection=ccrs.PlateCarree())
      ax.coastlines()
      plt.pcolormesh(traj.lon,traj.lat,dd,vmin=0,vmax=30)
      plt.title("%s %0.3f"%(f,keep[f][2].mean()))
      x2,y2=np.meshgrid(traj.lon[10:51],traj.lat)
      plt.plot(np.ma.masked_array(x2.flatten()[mask2[:,10:51].flatten()],mask=keep[f][2]),    np.ma.masked_array(y2.flatten()[mask2[:,10:51].flatten()],mask=keep[f][2]),"k.",ms=0.5)
      plt.plot(np.ma.masked_array(x2.flatten()[mask2[:,10:51].flatten()],mask=(1-keep[f][2])),np.ma.masked_array(y2.flatten()[mask2[:,10:51].flatten()],mask=(1-keep[f][2])),"r.",ms=0.5)
      plt.colorbar()  
      plt.show()

#plt.plot([keep[x][2].mean() for x in keep],dewpoint,".")
