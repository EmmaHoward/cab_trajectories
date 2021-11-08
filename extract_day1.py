import numpy as np
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

cubes=iris.cube.CubeList()
for j,f in enumerate(files):
  print(f)
  traj=trajectory_grid(filename=path+f,shape=(20,34,41,61))[:]
  traj.data["q"] = np.maximum(traj.data["q"],1e-9)
  traj.data["e"]= traj.data["p"]*100*traj.data["q"]/0.622
  traj.data["Td"]=(1/273.15-1.844*0.0001*np.log(traj.data["e"]/611.3))**-1
  dd = traj.data["T"][0]-traj.data["Td"][0] 
  sub = (traj.data["z"][0]-traj.data["z"][4])/24/3600
  pc = iris.coords.DimCoord(traj.p,standard_name="air_pressure",units="hPa")
  lonc = iris.coords.DimCoord(traj.lon,standard_name="longitude",units="degrees")
  latc = iris.coords.DimCoord(traj.lat[::-1],standard_name="latitude",units="degrees")
  t0 = dt.datetime(2009,1,1)
  timec=iris.coords.DimCoord([(t-t0).total_seconds()/24/3600 for t in traj.time],units="days since %04d-%02d-%02d"%(t0.year,t0.month,t0.day),standard_name="time")
  ddc =iris.cube.Cube(dd[:,::-1],standard_name="dew_point_depression", units="K",dim_coords_and_dims=[(pc,0),(latc,1),(lonc,2)])
  ddc.add_aux_coord(timec[0])
  subc =iris.cube.Cube(sub[:,::-1],long_name="subsidence", units="m/s",dim_coords_and_dims=[(pc,0),(latc,1),(lonc,2)])
  subc.add_aux_coord(timec[0])
  cubes += [ddc,subc]

cubes = cubes.merge()
iris.save(cubes,"dp_sub.nc")




