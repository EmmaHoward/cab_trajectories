import iris
from iris.experimental.equalise_cubes import equalise_attributes
from iris.coord_categorisation import add_day_of_year,add_year
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
path = "/storage/silver/metstudent/msc/users_2021/yc810080/data/"
path2 = "/storage/silver/NCAS-Weather/oa921437/tmp/"
ct=iris.Constraint(time=lambda t:t.point.month in[10,11])
v=iris.load(path2+"era5_w_*.nc",ct)
equalise_attributes(v)
z=v.concatenate_cube()

#z = iris.load(path+"gp300.nc")[0]
dp = iris.load(path+"surface_td_data.nc")

equalise_attributes(dp)
dp=dp.concatenate_cube()

add_day_of_year(z,"time","doyr")
add_year(z,"time","year")
z = z.aggregated_by(["doyr","year"],iris.analysis.MEAN)

add_day_of_year(dp,"time","doyr")
add_year(dp,"time","year")
dp = dp.aggregated_by(["doyr","year"],iris.analysis.MEAN)


z = z - z.collapsed("longitude",iris.analysis.MEAN)


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

z.remove_coord("doyr")
z.remove_coord("year")
ct = iris.Constraint(time=lambda t: t.point.month in [10,11])
cx = iris.Constraint(longitude = lambda x: 18<=x<=30)
cy = iris.Constraint(latitude = lambda y: -34<=y<=-10)

zz=z.extract(cx&cy).collapsed("longitude",iris.analysis.MEAN)
dp=dp.extract(cx&cy&ct).collapsed("longitude",iris.analysis.MEAN)

lat=dp.coord("latitude").points
ny=len(lat)
time = zz.coord("time").units.num2date(zz.coord("time").points)
from scipy.ndimage import convolve
zz=zz.data.reshape(8,61,ny)
dp=dp[:488].data.reshape(8,61,ny)

zz2=zz#-convolve(zz.mean(axis=0),np.ones((7,1))/7)
dp2=dp-convolve(dp.mean(axis=0),np.ones((7,1))/7)
#zz2 = zz.copy()
#dp2 = dp.copy()
#zz2 = zz2.mean(axis=2)

zm = zz2.mean(axis=(0,1))
dm = dp2.mean(axis=(0,1))

zs = zz2.std(axis=(0,1),ddof=0)
ds = dp2.std(axis=(0,1),ddof=0)

corr = [0 for i in range(15)]

corr[7] = ((zz2-zm)[:,:]*(dp2-dm)).mean(axis=(0,1))/ds/zs

for i in range(1,8):
  corr[7-i] = ((zz2[:,i:]-zm)*(dp2[:,:-i]-dm)).mean(axis=(0,1))/ds/zs
  corr[7+i] = ((zz2[:,:-i]-zm)*(dp2[:,i:]-dm)).mean(axis=(0,1))/ds/zs

corr=np.array(corr)  
plt.clf()
plt.contourf(np.arange(-7,8),lat,corr.T,np.arange(-0.425,0.45,0.05),cmap="bwr")

plt.colorbar()
plt.ylabel("latitude")
plt.xlabel("dewpoint leads omega      omega leads dewpoint \n lag (days)")
plt.title("lead-lag correlation between upper \n level vertical velocity and dewpoint")
plt.subplots_adjust(bottom=0.15)
plt.show()


zz2 = zz2.mean(axis=2)


zm = zz2.mean(axis=(0,1))
dm = dp2.mean(axis=(0,1))

zs = zz2.std(axis=(0,1),ddof=0)
ds = dp2.std(axis=(0,1),ddof=0)

corr = [0 for i in range(15)]

corr[7] = ((zz2-zm)[:,:,np.newaxis]*(dp2-dm)).mean(axis=(0,1))/ds/zs

for i in range(1,8):
  corr[7-i] = ((zz2[:,i:]-zm)[:,:,np.newaxis]*(dp2[:,:-i]-dm)).mean(axis=(0,1))/ds/zs
  corr[7+i] = ((zz2[:,:-i]-zm)[:,:,np.newaxis]*(dp2[:,i:]-dm)).mean(axis=(0,1))/ds/zs

corr=np.array(corr)  
plt.clf()
plt.contourf(np.arange(-7,8),lat,corr.T,np.arange(-0.425,0.45,0.05),cmap="bwr")

plt.ylabel("latitude")
plt.xlabel("dewpoint leads omega      omega leads dewpoint \n lag (days)")
plt.title("lead-lag correlation between upper \n level vertical velocity and dewpoint")
plt.colorbar()
plt.subplots_adjust(bottom=0.15)
plt.show()


dp = dp-273.15

cmap = cmap_discretise("RdBu_r",17)

fig=plt.figure(figsize=(9,7))
plat = np.array([lat[0]-0.5]+list(lat+0.5))
for i in range(8):
  plt.subplot(4,2,i+1)
  a=plt.pcolormesh(np.arange(-0.5,61),plat,zz[i].T,vmin=-0.5,vmax=0.5,cmap=cmap)
  b=plt.contour(np.arange(61),lat,dp[i].T,np.arange(7,26,2),linewidths=1,cmap="viridis_r")
  plt.ylim(-33,-10)
  plt.yticks([-30,-20,-10])
  if i%2==0:
    plt.ylabel("latitude")
  if i>5:
    plt.xlabel("days since 1/10")
  plt.title(i+2009)
  for t in traj:
    if t.year == i+2009:
       days = (t-dt.datetime(i+2009,10,1)).days
       plt.plot([days,days],[-34,0],":k")

fig.subplots_adjust(bottom=0.23,top=0.92,left=0.08,right=0.97,wspace=0.15,hspace=0.5)
cax1=fig.add_axes([0.15,0.1,0.3,0.05])
cax2=fig.add_axes([0.62,0.1,0.3,0.05])
cb1=fig.colorbar(a,cax=cax1,orientation="horizontal")
cb2=fig.colorbar(b,cax=cax2,orientation="horizontal")
cb1.set_label("500 hPa vertical velocity")
cb2.set_label("2m dewpoint temperature (C)")
plt.show()







comp = np.gradient(lcab,axis=1)

z1 = z[comp.flatten()>0]


z1 = z[comp.flatten()>0]
z2 = z[comp.flatten()<0]
z3 = z[comp.flatten()==0]

ax1=plt.subplot(211,projection=ccrs.PlateCarree())
ax1.coastlines()
iplt.pcolormesh(z1.collapsed("time",iris.analysis.MEAN)-z.collapsed("time",iris.analysis.MEAN),vmin=-100,vmax=100,cmap="PiYG")
plt.colorbar()

ax2=plt.subplot(212,projection=ccrs.PlateCarree())
ax2.coastlines()
iplt.pcolormesh(z2.collapsed("time",iris.analysis.MEAN)-z.collapsed("time",iris.analysis.MEAN),vmin=-100,vmax=100,cmap="PiYG")
plt.colorbar()


ax3=plt.subplot(313,projection=ccrs.PlateCarree())
ax3.coastlines()
iplt.pcolormesh(z3.collapsed("time",iris.analysis.MEAN)-z.collapsed("time",iris.analysis.MEAN))

plt.show()




