import sys
import numpy as np
import netCDF4 as nc
################################################################################
#                             DESCRIPTION
################################################################################
''' 
Script to calculate quantities related to water masses tracked in the adjoint
mode of NEMOTAM using a passive tracer. In this script, outputs are flipped along
the time axis, so the time dimension here can be intepreted as "age".

tracer_initial_volume
       ()        : Scalar value of passive tracer volume at run start
tracer_ventilation_prdens
       (t,  y,x) : probability density that tracer of age "t" has reached the 
                   surface at location (x,y)
tem_bins
       (T,)      : temperature values used to bin surface water tagged by tracer
sal_bins
       (S,)      : salinity values used to bin surface water tagged by tracer
tracer_TS_volume_histogram
       (t,S,T)   : cumulative volume of tracer (age "t" or younger) which has
                   ventilated in surface water of salinity S and temperature T
tracer_age_probability
       (t,)      : cumulative probability that tracer of  age "t" or younger 
                   has reached the surface 
'''

################################################################################
#                  LOAD ADJOINT OUTPUT AND GRID DATA:
################################################################################

# NB: For each run in WATER_MASS_RUNS, PTTAM_output_????????_????.nc
# should be stitched across mpp tiles (using TOOLS/rebuild_nemo) & concatenated
# in time (using ncrcat) to produce a file "<cn_exp>_output.nc", where
# cn_exp is the name of the run in the namelist.

PATH_TO_OUTPUTS  ='WATER_MASS_OUTPUTS/'
PATH_TO_MESH_MASK='./'

## TAM output:
OUTPUT_NC=nc.Dataset(PATH_TO_OUTPUTS+'WATER_MASS_adj_output.nc') # output file 

tracer_vol =np.flip(OUTPUT_NC.variables['pt_vol_ad' ][:],0) # Passive tracer concentration
tracer_vent=np.flip(OUTPUT_NC.variables['pt_vent_ad'][:],0) # Tracer ventilation volume
traj_sst   =np.flip(OUTPUT_NC.variables['tn'        ][:,0,:,:],0) # Trajectory SST
traj_sss   =np.flip(OUTPUT_NC.variables['sn'        ][:,0,:,:],0) # Trajectory SSS

## Grid data:
MESH=nc.Dataset(PATH_TO_MESH_MASK+'mesh_mask.nc')
e1t=MESH.variables['e1t'  ][:]
e2t=MESH.variables['e2t'  ][:]

################################################################################
#                            DIAGNOSTICS
################################################################################
noutputs=tracer_vol.shape[0]

# VENTILATION LOCATION PROBABILITY DENSITY:
tracer_initial_volume     = np.sum(tracer_vol[0,:])
tracer_ventilation_prdens = tracer_vent/(tracer_initial_volume*(e1t*e2t))

# VENTILATION TS PROBABILITY DENSITY:
## Initialise histogram
tem_bins=np.linspace(- 2  , 5  ,29)
sal_bins=np.linspace( 34.5,35.5,21)
tracer_TS_volume_histogram=np.zeros((noutputs,len(tem_bins)-1,len(sal_bins)-1))
## Populate histogram
for ii in np.arange(noutputs):
    tracer_TS_volume_histogram[ii,:,:],_,_=\
                    np.histogram2d(traj_sst[ii,:].flatten(),\
                                   traj_sss[ii,:].flatten(),\
                                   [tem_bins,sal_bins],\
                                   weights=tracer_vent[ii,:].flatten())

# TRACER AGE PROBABILITY DISTRIBUTION
tracer_age_probability=np.sum(tracer_vent.reshape(noutputs,-1),axis=1)\
                        /tracer_initial_volume


