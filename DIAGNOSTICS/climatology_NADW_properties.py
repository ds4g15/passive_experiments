import numpy as np
import netCDF4 as nc

PATH_TO_MESH_MASK             ='./'
PATH_TO_SUBBASINS             ='./'
PATH_TO_TRAJECTORY_CLIMATOLOGY='./'

MESH     =nc.Dataset(PATH_TO_MESH_MASK+'mesh_mask.nc')
SUBBASINS=nc.Dataset(PATH_TO_SUBBASINS+'subbasins.nc')
TRAJ_CLIM=nc.Dataset(PATH_TO_TRAJECTORY_CLIMATOLOGY+'TRAJ_CLIMATOLOGY_60y.nc')

# Climatological temperature and salinity

tn=TRAJ_CLIM.variables['tn'][0:365,:,:,:]
sn=TRAJ_CLIM.variables['sn'][0:365,:,:,:]

# Grid info and basin mask
atlmsk=SUBBASINS.variables['atlmsk'][:]
e1t=MESH.variables['e1t'][:]
e2t=MESH.variables['e2t'][:]
e3t=MESH.variables['e3t'][:]


# Finding NADW in climatology
############################
# Temperature constraint (2,4)*C
NADW_t=tn*atlmsk
NADW_t[ (NADW_t>4) | (NADW_t<2) ]=100
NADW_t[NADW_t!=100]=1
NADW_t[NADW_t==100]=0

# Salinity constraint (34.9,35.0) psu
NADW_s=sn*atlmsk
NADW_s[ (NADW_s>35.0) | (NADW_s<34.9) ]=100
NADW_s[NADW_s!=100]=1
NADW_s[NADW_s==100]=0

NADW_clim=NADW_t*NADW_s #Masks non-NADW locations in the climatology

np.save('NADW_climavg.npy',NADW_clim)

### Climatological NADW volume (m^3) and outcrop area (km^2) time series:
NADW_volume =np.sum( (NADW_clim*e1t*e2t*e3t).reshape(365,-1) ,axis=1)  
NADW_outcrop=np.sum( (NADW_clim[:,0,:,:]*e1t*e2t).reshape(365,-1),axis=1)/1e6 
