import iris
from iris.experimental.equalise_cubes import equalise_attributes
from scipy.ndimage import convolve
from iris.coord_categorisation import add_day_of_year,add_year
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import pandas as pd
path = "/storage/silver/metstudent/msc/users_2021/yc810080/data/"
dp = iris.load(path+"surface_td_data.nc")

equalise_attributes(dp)
dp=dp.concatenate_cube()

add_day_of_year(dp,"time","doyr")
add_year(dp,"time","year")

cab = pd.read_csv("~/josh_msc/era5_cab_0p25_1979-2018.csv",parse_dates=True,index_col=0)
dp = dp.aggregated_by(["doyr","year"],iris.analysis.MEAN)



def cmap_discretise(cmap, N,cyc=False):
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.cm import get_cmap
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.
    """
    if type(cmap) == str:
        cmap = get_cmap(cmap)
    if cyc:
      colors_i = np.concatenate((np.linspace(0, 1.-1./N, N), (0.,0.,0.,0.)))
    else:
      colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki, key in enumerate(('red','green','blue')):
        cdict[key] = [(indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1)]
    # Return colormap object.
    return LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)



files = ["2009101012","2009101412","2009112712","2010103012","2010111212","2011100812","2011102412","2011103112","2011111312","2012101912","2012103112","2012110512","2012112112","2012112912","2013102912","2013111712","2013112712","2014100712","2014101712","2014103112","2014110612","2014112012","2015102212","2015112612","2016111312","2016111812"]

traj  = [dt.datetime(int(f[:4]),int(f[4:6]),int(f[6:8])) for f in files] 

ct = iris.Constraint(time=lambda t: t.point.month in [10,11])
cx = iris.Constraint(longitude = lambda x: 15<=x<=25)
cy = iris.Constraint(latitude = lambda y: -25<=y<=-10)

new = iris.cube.CubeList()
for year in range(2009,2017):
  cyr = iris.Constraint(time = lambda t: t.point.year==year)
  new.append(dp.extract(cyr&ct))
  new[-1].remove_coord("time")
  new[-1].remove_coord("year")
  iris.util.promote_aux_coord_to_dim_coord(new[-1],"doyr")
  new[-1].add_aux_coord(iris.coords.DimCoord(year,long_name="year"))
  if year%4==0:
     new[-1].coord("doyr").points = new[-1].coord("doyr").points - 1

dp = new.merge_cube()

dp=dp.extract(cx&cy)#.collapsed(["longitude","latitude"],iris.analysis.MEAN)
dpm = convolve(dp.data.mean(axis=0),np.ones((7,1,1))/7)
#dp.data =np.ma.masked_array(dp.data,np.ones(dp.shape)*(dpm>273.15+18))
dp=dp.collapsed(["longitude","latitude"],iris.analysis.MEAN)
dp = dp - 273.15
lat=dp.coord("latitude").points
ny=len(lat)
plt.ion()
for year in range(2009,2017):
  plt.subplot(4,2,year-2008)
  #iplt.contourf(dp.extract(cyr)-dpm,cmap="bwr",vmin=-5,vmax=5)
  iplt.plot(dp[year-2009]/dp[year-2009].data.max())
  gc = cab["CAB Grid Cells"][dt.datetime(year,10,1):dt.datetime(year,11,30)]
  plt.plot(dp.coord("doyr").points,gc/gc.max(),".-")
  for t in traj:
    if t.year==year:
       tt = (t - dt.datetime(year,1,1)).days+int(year%4!=0)
       plt.plot([tt,tt],[0,1],'k')
       #plt.plot([tt,tt],[dp.data.min(),dp.data.max()],'k')


plt.show()
