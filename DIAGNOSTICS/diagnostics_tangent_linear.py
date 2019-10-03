import sys
import numpy as np
import netCDF4 as nc
################################################################################
#                             DESCRIPTION
################################################################################
''' 
Script to calculate quantities related to water masses tracked in the 
tangent-linear mode of NEMOTAM using a passive tracer. Calculates:
tracer_volume                 
       (t,z,y,x) : time series of 3D tracer volume
tracer_depth_integrated_volume
       (t,  y,x) : time series of 2D depth-integrated tracer volume
tracer_initial_volume
       ()        : Scalar value of passive tracer volume at run start
tracer_total_volume
       (t,)      : Time series of total passive tracer volume through run
tracer_depth_integrated_prdens
       (t,  y,x) : Time-evolving probability density (/m^2) that tracer is
                   found at a given location
lat_bar         
       (t,)      : time series of tracer weighted centre of mass latitude
lon_bar         
       (t,)      : time series of tracer weighted centre of mass longitude
dep_bar         
       (t,)      : time series of tracer weighted centre of mass depth
tem_bins
       (T,)      : temperature values used to bin water tagged by tracer
sal_bins
       (S,)      : salinity values used to bin water tagged by tracer
tracer_TS_volume_histogram
       (t,S,T)   : Time-evolving distribution of passive tracer volume 
                   in TS space
'''

################################################################################
#                  LOAD TANGENT-LINEAR OUTPUT AND GRID DATA:
################################################################################

# NB: For each run in WATER_MASS_RUNS, PTTAM_output_????????_????.nc
# should be stitched across mpp tiles (using TOOLS/rebuild_nemo) & concatenated
# in time (using ncrcat) to produce a file "<cn_exp>_output.nc", where
# cn_exp is the name of the run in the namelist.

PATH_TO_OUTPUTS  ='WATER_MASS_OUTPUTS/'
PATH_TO_MESH_MASK='./'

## TAM output:
OUTPUT_NC=nc.Dataset(PATH_TO_OUTPUTS+'WATER_MASS_tan_output.nc')

tracer_conc=OUTPUT_NC.variables['pt_conc_tl'][:] # Passive tracer concentration
traj_tn   =OUTPUT_NC.variables['tn'        ][:] # Trajectory temperature
traj_sn    =OUTPUT_NC.variables['sn'        ][:] # Trajectory salinity

## Grid data:
MESH=nc.Dataset(PATH_TO_MESH_MASK+'mesh_mask.nc')
e1t=MESH.variables['e1t'  ][:]
e2t=MESH.variables['e2t'  ][:]
e3t=MESH.variables['e3t'  ][:]
lat=MESH.variables['gphit'][0,:]
lon=MESH.variables['glamt'][0,:]
dep=np.cumsum(e3t[0,:],axis=0)

################################################################################
#                              DIAGNOSTICS
################################################################################

# DEPTH-INTEGRATED PROBABILITY DENSITY
noutputs=np.shape(tracer_conc)[0] #number of outputs
tracer_volume                 =tracer_conc*(e1t*e2t*e3t)    #Volume in each grid cell
tracer_depth_integrated_volume=np.sum(tracer_volume,axis=1) #Depth-integrated volume
tracer_initial_volume         =np.sum(tracer_volume[0,:])   #Injected tracer volume
tracer_total_volume           =np.sum( (tracer_volume).reshape(noutputs,-1) ,axis=1)
tracer_depth_integrated_prdens=(tracer_depth_integrated_volume/tracer_initial_volume)\
                             /(e1t*e2t)

# TRACER CENTRE OF MASS:
## Get mean lat and lon:
### Project grid to Cartesian coordinates using inverse spherical projection:
R_earth=6371e3
X=R_earth*np.sin(np.deg2rad(lat+90))*np.cos(np.deg2rad(lon+180))
Y=R_earth*np.sin(np.deg2rad(lat+90))*np.sin(np.deg2rad(lon+180))
Z=R_earth*np.cos(np.deg2rad(lat+90))

### Mean position in cartesian coordinates = sum(volume*{X,Y,Z})/sum(volume)
X_bar=np.sum( (X*tracer_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       /tracer_total_volume
Y_bar=np.sum( (Y*tracer_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       /tracer_total_volume
Z_bar=np.sum( (Z*tracer_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       /tracer_total_volume

### Project mean position back to spherical coordinates and get lat,lon:
lat_bar  =(np.rad2deg(np.arccos(Z_bar/np.sqrt(X_bar**2 + Y_bar**2 + Z_bar**2)))-90)
lon_bar  =(np.rad2deg(np.arctan2(Y_bar,X_bar))-180)

## Get mean depth:
### Mean depth = sum(volume * depth)/sum(volume)
dep_bar = np.sum( (tracer_volume*dep).reshape(noutputs,-1),axis=1 )\
          /tracer_total_volume


# TS PROPERTIES OF WATER OCCUPIED BY TRACER:
## Initialise histogram
tem_bins=np.linspace(- 2  , 5  ,29)
sal_bins=np.linspace( 34.5,35.5,21)
tracer_TS_volume_histogram=np.zeros((noutputs,len(tem_bins)-1,len(sal_bins)-1))

## Populate histogram
for ii in np.arange(noutputs):
    tracer_TS_volume_histogram[ii,:,:],_,_=\
                    np.histogram2d(traj_tn[ii,:].flatten(),\
                                   traj_sn[ii,:].flatten(),\
                                   [tem_bins,sal_bins],\
                                   weights=tracer_volume[ii,:].flatten())
                          
