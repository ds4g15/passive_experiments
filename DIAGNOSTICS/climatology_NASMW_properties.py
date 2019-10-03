import numpy as np
import netCDF4 as nc

PATH_TO_MESH_MASK             ='./'
PATH_TO_SUBBASINS             ='./'
PATH_TO_TRAJECTORY_CLIMATOLOGY='./'

MESH     =nc.Dataset(PATH_TO_MESH_MASK+'mesh_mask.nc')
SUBBASINS=nc.Dataset(PATH_TO_SUBBASINS+'subbasins.nc')
TRAJ_CLIM=nc.Dataset(PATH_TO_TRAJECTORY_CLIMATOLOGY+'TRAJ_CLIMATOLOGY_60y.nc')

# Climatological temperature and salinity
print('loading trajectory T')
tn=TRAJ_CLIM.variables['tn'][0:365,:]
print('loading trajectory S')
sn=TRAJ_CLIM.variables['sn'][0:365,:]

# Grid info and basin mask
atlmsk=SUBBASINS.variables['atlmsk'][:]
atlmsk[SUBBASINS.variables['navlat'][:]<0]=0 #North Atlantic only
e1t=MESH.variables['e1t'][:]
e2t=MESH.variables['e2t'][:]
e3t=MESH.variables['e3t'][:]
lon=MESH.variables['glamt'][:]

# Finding NASMW in climatology
############################
#Temperature constraint
NASMW_t=tn*atlmsk
NASMW_t[ (tn*atlmsk>19) | (tn*atlmsk<17) ]=100
NASMW_t[NASMW_t!=100]=1
NASMW_t[NASMW_t==100]=0

#Salinity constraint
NASMW_s=sn*atlmsk
NASMW_s[ (NASMW_s>36.6) | (NASMW_s<36.4) ]=100
NASMW_s[NASMW_s!=100]=1
NASMW_s[NASMW_s==100]=0

#Zonal boundary
NASMW_lon=lon*atlmsk
NASMW_lon[NASMW_lon<-35]=1000
NASMW_lon[NASMW_lon!=1000]=0
NASMW_lon[NASMW_lon==1000]=1

#All water satisfying above 3 constraints:
NASMW_TSlon=NASMW_t*NASMW_s*NASMW_lon

#Thickness constraint
NASMW_thck=np.sum(NASMW_TSlon*e3t,axis=1,keepdims=True) #Thickness of mode water
NASMW_thck[NASMW_thck<125]=0                           #Mask if <125m
NASMW_thck[NASMW_thck!=0 ]=1                           

# Result:
NASMW_clim=NASMW_TSlon*NASMW_thck                       #Apply mask to mode water

np.save('NASMW_climavg.npy',NASMW_clim)

# Time series of NASMW volume (m^3) and outcrop area (km^2)
NASMW_volume =np.sum( (NASMW_clim*e1t*e2t*e3t).reshape(365,-1) ,axis=1)  
NASMW_outcrop=np.sum( (NASMW_clim[:,0,:,:]*e1t*e2t).reshape(365,-1),axis=1)/1e6 
