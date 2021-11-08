import numpy as np
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
from trajectory_grid import trajectory_grid
from matplotlib.collections import LineCollection
import datetime as dt
import iris.plot as iplt
import matplotlib.pyplot as plt
import iris
from iris.analysis.calculus import differentiate
from iris.coord_categorisation import add_day_of_year
import pickle
norm = plt.Normalize(-50,0)
path = "/storage/silver/NCAS-Weather/oa921437/tmp/"

tpath = "/storage/silver/metstudent/msc/users_2021/yc810080/data/traj_data2/"
tpath3 = "/storage/silver/metstudent/msc/users_2021/yc810080/data/traj_data3/"

files = ["2009101012","2009101412","2009112712","2010103012","2010111212","2011100812","2011102412","2011103112","2011111312","2012101912","2012103112","2012110512","2012112112","2012112912","2013102912","2013111712","2013112712","2014100712","2014101712","2014103112","2014110612","2014112012","2015102212","2015112612","2016111312","2016111812"]

files3 = ["2012103112","2013111712","2014103112","2016111312"]
 
keep = {}  
#keep  = pickle.load(open("prop_trajectories.pickle","rb"))
traj  = [dt.datetime(int(f[:4]),int(f[4:6]),int(f[6:8])) for f in files]

summaries = pickle.load(open("mean_trajectories.pickle","rb"))
cp = iris.Constraint(pressure_level=200)
fig=plt.figure(figsize=(12,8))
ct=iris.Constraint(time = lambda t: t.point.month in [10,11])
cy = iris.Constraint(latitude = lambda y: -50<=y<=-25)
for year in range(2009,2017):
  print(year)
  z=iris.load_cube(path+"era5_z_%d.nc"%year,ct)
  add_day_of_year(z,"time","doyr")
  #z=z.aggregated_by("doyr",iris.analysis.MEAN)#
  z=z.extract(cy&cp)
  z.coord("latitude").guess_bounds()
  z.coord("longitude").guess_bounds()
  weights = iris.analysis.cartography.area_weights(z)
  z=z.collapsed("latitude",iris.analysis.MEAN,weights=weights)
  z = z.rolling_window("time",iris.analysis.MEAN,24)
  z.coord("time").convert_units("days since %d-10-01"%year)
  dz = differentiate(z,"longitude")
  v=iris.load_cube(path+"era5_v_%d.nc"%year,ct)
  add_day_of_year(v,"time","doyr")
  #v=v.aggregated_by("doyr",iris.analysis.MEAN)#
  v=v.extract(cy)
  v.coord("latitude").guess_bounds()
  v.coord("longitude").guess_bounds()
  weights = iris.analysis.cartography.area_weights(v)
  v=v.collapsed("latitude",iris.analysis.MEAN,weights=weights)
  v.coord("time").convert_units("days since %d-10-01"%year)
  v = v.rolling_window("time",iris.analysis.MEAN,24)
  dv = differentiate(v,"longitude")
  mask = dz.copy(data = (dz.data < 0) * (dv.data > 0))
  w=iris.load_cube(path+"era5_w_%d.nc"%year,ct)
  add_day_of_year(w,"time","doyr")
  w=w.aggregated_by("doyr",iris.analysis.MEAN)#
  cy = iris.Constraint(latitude = lambda y: -40<=y<=-20)
  w=w.extract(cy)
  w.coord("latitude").guess_bounds()
  w.coord("longitude").guess_bounds()
  weights = iris.analysis.cartography.area_weights(w)
  w=w.collapsed("latitude",iris.analysis.MEAN,weights=weights)
  w.coord("time").convert_units("days since %d-10-01"%year)
  w.coord("doyr").points = w.coord("doyr").points - 273 - int(year%4==0)  +0.5
  mask.coord("doyr").points = [(t - dt.datetime(year,10,1)).total_seconds()/3600/24 for t in mask.coord("time").units.num2date(mask.coord("time").points)]
  ax=plt.subplot(1,8,year-2008)
  a=iplt.contourf(w,np.arange(-0.4,0.4,0.05),cmap="BrBG_r",coords=["longitude","doyr"],extend="both")
  b=iplt.contour(mask,[0.5],color="k",hatches=["XX"],coords=["longitude","doyr"])
  plt.title(year)
#  b=iplt.contourf(v,np.arange(-35,36,5),cmap="PiYG",coords=["longitude","doyr"])
  if year==2009:
    plt.yticks([1.5,10.5,20.5,31.5,41.5,51.5,61.5],["1 Oct","10 Oct","20 Oct","1 Nov","10 Nov","20 Nov","30 Nov"])
  else:
    plt.yticks([1.5,10.5,20.5,31.5,41.5,51.5,61.5],["","","","","","",""])
  plt.xlim(0,60)
  lon = mask.coord("longitude").points
  for f in files:
    if int(f[:4])==year:
      t = dt.datetime(int(f[:4]),int(f[4:6]),int(f[6:8]))
#      s = summaries["utraj-rf_cab"+f]
#      plt.plot(s["lon"], [(tt-dt.datetime(year,10,1)).total_seconds()/3600/24 for tt in s["time"]],"k")
      plt.plot([0,60],[(t-dt.datetime(year,10,1)).days+0.5,(t-dt.datetime(year,10,1)).days+0.5],"r")
      if f in files3:
       traj=trajectory_grid(filename=tpath3+"utraj-rf_cab"+f,shape=(40,34,41,61))[:,8,::,10:51:]
      else:
        traj=trajectory_grid(filename=tpath+"utraj-rf_cab"+f,shape=(20,34,41,61))[:,8,::,10:51:]
      nt,ny,nx = traj.shape
      traj.data["lon"][traj.data["lon"]>180] -= 360
      traj.data["q"] = np.maximum(traj.data["q"],1e-9)
      traj.data["e"]= traj.data["p"]*100*traj.data["q"]/0.622
      traj.data["Td"]=(1/273.15-1.844*0.0001*np.log(traj.data["e"]/611.3))**-1
      dd = traj.data["T"][0]-traj.data["Td"][0] 
      mask2 = (dd>=8) #* (traj.data["lon"] <50).all(axis=0)
      traj.data["lon"] = np.ma.masked_array(traj.data["lon"])
      traj.data["lat"] = np.ma.masked_array(traj.data["lat"])
      traj.data["lon"].mask += (1-mask2)
      phase = []
      for i in range(nt):
        ct_ = iris.Constraint(time = lambda t: (t.point.month == traj.time[i].month) & (t.point.day == traj.time[i].day) & (t.point.hour == traj.time[i].hour) )
        phase.append(mask.extract(ct_).data[[np.argmin(np.abs(lon - x)) for x in traj.data["lon"][i][mask2]]])
      phase = (np.array(phase[1:])*  (traj.data["lon"][1:][:,mask2] < traj.data["lon"][:-1][:,mask2]) * (traj.data["lat"][:-1][:,mask2]<-20 ) ).any(axis=(0))
      print(phase.shape)
      x = traj.data["lon"][:,mask2]
#      [time,x,phase] = keep[f]
      keep[f] = [traj.time,x,phase,mask2]
      for i in np.random.randint(0,mask2.sum(),50):
        plt.plot(x[:,i],np.array([(tt-dt.datetime(year,10,1)).total_seconds()/3600/24 for tt in traj.time])[:,np.newaxis],["r","b"][phase[i]],lw=0.2)
#      times = np.array([(tt-dt.datetime(year,10,1)).total_seconds()/3600/24 for tt in traj.time])
#      iy,ix=np.where(mask)
#      segments = np.array([[traj.data["lon"][:,i,j] for (i,j) in zip(iy,ix)],
#                      [times for (i,j) in zip(iy,ix)]]).transpose(1,2,0)
#      lc = LineCollection(segments, cmap='magma', norm=norm)
#      lc.set_array(np.array([traj.data["lat"][:,i,j] for (i,j) in zip(iy,ix)]).flatten())
#      lc.set_linewidth(0.2)
#      line=ax.add_collection(lc)
#cax1=fig.add_axes([0.1,0.05,0.35,0.05])
#fig.colorbar(a,cax=cax1,orientation="horizontal")

pickle.dump(keep,open("prop_trajectories2.pickle","wb"))

plt.subplots_adjust(bottom=0.15,left=0.08)
cax2=fig.add_axes([0.3,0.05,0.4,0.05])
fig.colorbar(a,cax=cax2,orientation="horizontal")
plt.show()
      

