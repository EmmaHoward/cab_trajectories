import numpy as np
import datetime as dt
attr_keys = {1:"T",3:"PV",4:"q",10:"z",159:"blh"}

def date(s):
  return dt.datetime(int(s[:4]),int(s[4:6]),int(s[6:8]),int(s[8:10]))

class trajectory_grid:
  def __init__(self,data={},filename="unset",shape = ()):
  
    self.filename = filename
    if len(data)==0:
      self.shape = shape
      self.read()
      self.init_coords()
    else:   
      self.data = data
      self.shape = data[list(data.keys())[0]].shape
       
  def read(self):
     nt,nz,ny,nx = self.shape
     with open(self.filename,"r") as f:
       a = f.readline()
       self.traj_basetime = date(a.split()[-1])
       a = f.readline()
       self.data_basetime = date(a.split()[-1])
       a = f.readline()       
       self.time_interval = int(a.split()[3])
       self.time_units = a.split()[4]
       self.interval_timesteps = int(a.split()[-2])
       a = f.readline()     
       self.ntraj = int(a.split()[-1])
       assert (nx*ny*nz == self.ntraj),   \
              "number of trajectories in file ({0}) doesn't match shape ({1},{2},{3})"\
              .format(self.ntraj,nz,ny,nx)
       a = f.readline()       
       self.n_attr = int(a.split()[-1])
       a = f.readline()     
       a = f.readline()              
       attr_indices = [int(i) for i in a.split()]
       keys = ["lat","lon","p"]+[attr_keys[i] for i in attr_indices]
       a = f.readline()              
       nz2 =  int(a.split()[-1])
       assert (nz == nz2), "number of clusters in file ({1}) doesn't match nz ({2})".format(nz2,nz)
       for i in np.arange(np.ceil(34/10)+1):
         a = f.readline()              
       a = f.readline()              
       a = f.readline()              
       a = f.readline()              
       self.forward = a.split()[-1] 
       data = {}
       for key in keys:
         data[key] = np.zeros((nt,nz,ny,nx))
       for iz in range(nz):
         for iy in range(ny):
           for ix in range(nx):
             a = f.readline()              
             a = f.readline()        
             nt_t = int(a.split()[-2])
             if nt_t < nt:             
               for key in keys:
                 data[key][nt_t:,iz,iy,ix]=np.nan
             a = f.readline()    
             a = f.readline()                           
             for it in range(nt_t):
                a = f.readline().split()              
                for j,key in enumerate(keys):
                  data[key][it,iz,iy,ix] = float(a[j+2])
             a = f.readline()              
       f.close()
       self.data=data 
   
  def init_coords(self):
      self.coords = ["time","p","lat","lon"]
      self.lat = self.data["lat"][0,0,:,0]
      self.lon = self.data["lon"][0,0,0]
      self.p = self.data["p"][0,:,:,:].max(axis=(1,2))
      if self.forward == "T":
        self.time = np.array([self.traj_basetime+dt.timedelta(self.time_interval*i/24) for i in range(self.shape[0])])
      elif self.forward == "F":
        self.time = np.array([self.traj_basetime-dt.timedelta(self.time_interval*i/24) for i in range(self.shape[0])])
    
  def __getitem__(self,items):
    data = {}
    for key in self.data.keys():
      data[key]=self.data[key][items]
    new = trajectory_grid(data,filename=self.filename,shape=data[key].shape)
    new_coords = []
    if type(items)!=tuple:
      items = (items,)
    for i,coord in enumerate(self.coords): 
        if i >= len(items):
           new_coords.append(coord)
           setattr(new,coord,getattr(self,coord))
        elif type(items[i])==slice: 
           new_coords.append(coord)
           setattr(new,coord,getattr(self,coord)[items[i]])
        elif type(items[i])==int:
           setattr(new,coord,getattr(self,coord)[items[i]])
        else:
           assert 1==0, "shouldn't be here"
    new.coords = new_coords
    for key in  ["traj_basetime","data_basetime","time_units","time_interval","interval_timesteps","n_attr","forward"]:
       setattr(new,key,getattr(self,key))
    for key in ["time","p","lat","lon"]:
      if key not in dir(new):
        setattr(new,key,getattr(self,key))
    return new

  def __repr__(self):
     t=self.traj_basetime
     s = "trajectory_grid(base=%04d-%02d-%02d_%02d, "%(t.year,t.month,t.day,t.hour)
     for key in ["p","lat","lon"]:
        if key in self.coords:
          s +="%s: %0.1f-%0.1f, "%(key, getattr(self,key)[0], getattr(self,key)[-1])
        else:
          s +="%s: %0.1f, "%(key, getattr(self,key))
     s += "%s)"%(self.shape,)
     return s

