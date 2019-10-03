import sys
import numpy as np
import netCDF4 as nc

PATH_TO_MESH_MASK             ='./'
PATH_TO_SUBBASINS             ='./'
PATH_TO_TRAJECTORY_CLIMATOLOGY='./'

MESH     =nc.Dataset(PATH_TO_MESH_MASK+'mesh_mask.nc')
SUBBASINS=nc.Dataset(PATH_TO_SUBBASINS+'subbasins.nc')
TRAJ_CLIM=nc.Dataset(PATH_TO_TRAJECTORY_CLIMATOLOGY+'TRAJ_CLIMATOLOGY_60y.nc')

print('loading grid data')

e1v=MESH.variables['e1v'][0,:,:]
e3v=MESH.variables['e3v'][0,:,:,:]

print('loading basin mask')
atlmsk=SUBBASINS.variables['atlmsk'][:]
atlmsk[MESH.variables['nav_lat'][:]<0]=0 #North Atlantic only

########################################
### BAROTROPIC AND MERIDIONAL STREAM FUNCTION CALCULATIONS
print('getting mean velocity fields')
print('zonal')
u_bar=np.mean(TRAJ_CLIM.variables['un'][:],axis=0,keepdims=False) #Annual mean zonal velocity
print('meridional')
v_bar=np.mean(TRAJ_CLIM.variables['vn'][:],axis=0,keepdims=False) #Annual mean meridional velocity

print('calculating barotropic stream function')

# Calculate Atlantic barotropic stream function (Sv, integrated from W to E)
BSFv=np.cumsum(np.sum(e3v[:,:,:]*e1v[:,:]*v_bar[:,:,:]*atlmsk[:,:],axis=0),axis=1)*1e-6
print('calculating meridional overturning stream function')
#Calculate Atlantic meridional overturning stream function (Sv):
MSF=np.cumsum(np.sum(atlmsk*e1v*e3v*v_bar,axis=2),axis=0)*1e-6
