import sys
import numpy as np
import netCDF4 as nc
################################################################################
#                             DESCRIPTION
################################################################################
'''
Script to calculate quantities allowing the comparison of different advection 
schemes in NEMOTAM. Schemes are abbreviated as:
TVD: Total variation diminishing scheme (NEMO default, nonlinear)
CE : Centred scheme (NEMOTAM default, forward in time, centred in space)
UW : Trajectory-upstream scheme
WM : Weighted-mean scheme (c.f. Fiadeiro and Veronis, 1977)

The script calculates the following quantities (where "X" is any of TVD,CE,UW,WM):
X_volume
       (t,z,y,x) : time series of 3D tracer volume with advection scheme X
X_total_volume
       (t,)      : time series of total passive tracer volume through run
X_depth_integrated_volume
       (t,  y,x) : Lateral distribution of total passive tracer volume
X_horiz_integrated_volume
       (t,z)     : Vertical distribution of total passive tracer volume
X_total_positive_volume
       (t,)      : The volume of tracer with concentration of 0 or higher
X_total_negative_volume
       (t,)      : The volume of tracer with concentration less than 0
lat_bar_X
       (t,)      : Time series of tracer weighted centre of mass latitude
lon_bar_X
       (t,)      : Time series of tracer weighted centre of mass longitude
dep_bar_X
       (t,)      : Time series of tracer weighted centre of mass depth
lateral_STD_X
       (t,)      : The lateral spread of passive tracer, as standard deviation
                   about its mean horizontal location
vertical_STD_X
       (t,)      : The vertical spread of passive tracer, as standard deviation
                   about its mean location in depth
'''
################################################################################
#                  LOAD TANGENT-LINEAR OUTPUT AND GRID DATA:
################################################################################
# NB: For each run in ADV_SCHEME_RUNS, tl_trajectory_????????_????.nc
# should be stitched across mpp tiles (using TOOLS/rebuild_nemo) & concatenated
# in time (using ncrcat) to produce a file "<cn_exp>_output.nc", where
# cn_exp is the name of the run in the namelist.

PATH_TO_OUTPUTS  ='ADV_OUTPUTS/'
PATH_TO_MESH_MASK='./'

# Load tangent-linear outputs of demo run with each advection scheme

TVD_NC=nc.Dataset(PATH_TO_OUTPUTS+'adv_TVD_output.nc')
WM_NC =nc.Dataset(PATH_TO_OUTPUTS+'adv_weighted_mean_output.nc')
UW_NC =nc.Dataset(PATH_TO_OUTPUTS+'adv_upwind_output.nc')
CE_NC =nc.Dataset(PATH_TO_OUTPUTS+'adv_centred_output.nc')

TVD_conc=TVD_NC.variables['pt_conc_tl'][:]
WM_conc = WM_NC.variables['pt_conc_tl'][:]
UW_conc = UW_NC.variables['pt_conc_tl'][:]
CE_conc = CE_NC.variables['pt_conc_tl'][:]

MESH=nc.Dataset(PATH_TO_MESH_MASK+'mesh_mask.nc')
e1t =MESH.variables['e1t'  ][:]
e2t =MESH.variables['e2t'  ][:]
e3t =MESH.variables['e3t'  ][:]
lat =MESH.variables['gphit'][0,:]
lon =MESH.variables['glamt'][0,:]
dep0=MESH.variables['gdept_0'][0,:]
dep =np.cumsum(e3t[0,:],axis=0)



################################################################################
#                              DIAGNOSTICS
################################################################################
# FUNCTIONS:

def haversine(lon1, lat1, lon2, lat2):
    """Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)"""
    # convert decimal degrees to radians
    
    lon1, lat1, lon2, lat2 = map(np.deg2rad, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return c * r


#################################################################################
# TRACER VOLUMES:
noutputs=TVD_conc.shape[0]

TVD_volume                  = TVD_conc*(e1t*e2t*e3t)
TVD_total_volume            = np.sum(                  (TVD_volume ).reshape(noutputs,-1),axis=1)
TVD_depth_integrated_volume = np.sum(   TVD_volume                                       ,axis=1)
TVD_horiz_integrated_volume = np.sum(   TVD_volume.reshape(noutputs,31,-1)               ,axis=2)
TVD_total_positive_volume   = np.sum( ((TVD_volume>=0)*(TVD_volume)).reshape(noutputs,-1),axis=1)
TVD_total_negative_volume   = np.abs(\
                              np.sum( ((TVD_volume< 0)*(TVD_volume)).reshape(noutputs,-1),axis=1) )

CE_volume                   = CE_conc *(e1t*e2t*e3t)
CE_total_volume             = np.sum(                  (CE_volume  ).reshape(noutputs,-1),axis=1)
CE_depth_integrated_volume  = np.sum(   CE_volume                                        ,axis=1)
CE_horiz_integrated_volume  = np.sum(   CE_volume.reshape(noutputs,31,-1)                ,axis=2)
CE_total_positive_volume    = np.sum( ((CE_volume >=0)*(CE_volume )).reshape(noutputs,-1),axis=1)
CE_total_negative_volume    = np.abs(\
                              np.sum( ((CE_volume < 0)*(CE_volume )).reshape(noutputs,-1),axis=1) )

UW_volume                   = UW_conc *(e1t*e2t*e3t)
UW_total_volume             = np.sum(                  (UW_volume  ).reshape(noutputs,-1),axis=1)
UW_depth_integrated_volume  = np.sum(   UW_volume                                        ,axis=1)
UW_horiz_integrated_volume  = np.sum(   UW_volume.reshape(noutputs,31,-1)                ,axis=2)
UW_total_positive_volume    = np.sum( ((UW_volume >=0)*(UW_volume )).reshape(noutputs,-1),axis=1)
UW_total_negative_volume    = np.abs(\
                              np.sum( ((UW_volume < 0)*(UW_volume )).reshape(noutputs,-1),axis=1) )

WM_volume                   = WM_conc *(e1t*e2t*e3t)
WM_total_volume             = np.sum(                  (WM_volume  ).reshape(noutputs,-1),axis=1)
WM_depth_integrated_volume  = np.sum(   WM_volume                                        ,axis=1)
WM_horiz_integrated_volume  = np.sum(   WM_volume.reshape(noutputs,31,-1)                ,axis=2)
WM_total_positive_volume    = np.sum( ((WM_volume >=0)*(WM_volume )).reshape(noutputs,-1),axis=1)
WM_total_negative_volume    = np.abs(\
                              np.sum( ((WM_volume < 0)*(WM_volume )).reshape(noutputs,-1),axis=1) )


# TRACER CENTRE OF MASS:
## Get mean lat and lon:
### Project grid to Cartesian coordinates using inverse spherical projection:
R_earth=6371e3
X=R_earth*np.sin(np.deg2rad(lat+90))*np.cos(np.deg2rad(lon+180))
Y=R_earth*np.sin(np.deg2rad(lat+90))*np.sin(np.deg2rad(lon+180))
Z=R_earth*np.cos(np.deg2rad(lat+90))

### Mean position in cartesian coordinates = sum(volume*{X,Y,Z})/sum(volume)
#### TVD 
X_bar_TVD=np.sum( (X*TVD_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       /TVD_total_volume
Y_bar_TVD=np.sum( (Y*TVD_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       /TVD_total_volume
Z_bar_TVD=np.sum( (Z*TVD_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       /TVD_total_volume
#### Centred
X_bar_CE =np.sum( (X* CE_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / CE_total_volume
Y_bar_CE =np.sum( (Y* CE_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / CE_total_volume
Z_bar_CE =np.sum( (Z* CE_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / CE_total_volume
#### Upwind
X_bar_UW =np.sum( (X* UW_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / UW_total_volume
Y_bar_UW =np.sum( (Y* UW_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / UW_total_volume
Z_bar_UW =np.sum( (Z* UW_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / UW_total_volume
#### Weighted-mean
X_bar_WM =np.sum( (X* WM_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / WM_total_volume
Y_bar_WM =np.sum( (Y* WM_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / WM_total_volume
Z_bar_WM =np.sum( (Z* WM_depth_integrated_volume).reshape(noutputs,-1),axis=1)\
       / WM_total_volume
### Project mean position back to spherical coordinates and get lat,lon:
lat_bar_TVD =(np.rad2deg(np.arccos(Z_bar_TVD/np.sqrt(X_bar_TVD**2 + Y_bar_TVD**2 + Z_bar_TVD**2)))-90)
lat_bar_CE  =(np.rad2deg(np.arccos(Z_bar_CE /np.sqrt(X_bar_CE **2 + Y_bar_CE **2 + Z_bar_CE **2)))-90)
lat_bar_UW  =(np.rad2deg(np.arccos(Z_bar_UW /np.sqrt(X_bar_UW **2 + Y_bar_UW **2 + Z_bar_UW **2)))-90)
lat_bar_WM  =(np.rad2deg(np.arccos(Z_bar_WM /np.sqrt(X_bar_WM **2 + Y_bar_WM **2 + Z_bar_WM **2)))-90)

lon_bar_TVD =(np.rad2deg(np.arctan2(Y_bar_TVD,X_bar_TVD))-180)
lon_bar_CE  =(np.rad2deg(np.arctan2(Y_bar_CE ,X_bar_CE ))-180)
lon_bar_UW  =(np.rad2deg(np.arctan2(Y_bar_UW ,X_bar_UW ))-180)
lon_bar_WM  =(np.rad2deg(np.arctan2(Y_bar_WM ,X_bar_WM ))-180)

## Get mean depth:
### Mean depth = sum(volume * depth)/sum(volume)
dep_bar_TVD = np.sum( (TVD_volume*dep).reshape(noutputs,-1),axis=1 )/ TVD_total_volume
dep_bar_CE  = np.sum( (CE_volume *dep).reshape(noutputs,-1),axis=1 )/  CE_total_volume
dep_bar_UW  = np.sum( (UW_volume *dep).reshape(noutputs,-1),axis=1 )/  UW_total_volume
dep_bar_WM  = np.sum( (WM_volume *dep).reshape(noutputs,-1),axis=1 )/  WM_total_volume

# TRACER LATERAL STANDARD DEVIATION
## Distance of all grid points from centre of mass of tracer patch
lateral_distance_from_centre_TVD = haversine(lon                              ,lat,\
                                             lon_bar_TVD.reshape(noutputs,1,1),lat_bar_TVD.reshape(noutputs,1,1))
lateral_distance_from_centre_CE  = haversine(lon                              ,lat,\
                                             lon_bar_CE. reshape(noutputs,1,1),lat_bar_CE .reshape(noutputs,1,1) )
lateral_distance_from_centre_UW  = haversine(lon                              ,lat,\
                                             lon_bar_UW. reshape(noutputs,1,1),lat_bar_UW .reshape(noutputs,1,1) )
lateral_distance_from_centre_WM  = haversine(lon                              ,lat,\
                                             lon_bar_WM. reshape(noutputs,1,1),lat_bar_WM .reshape(noutputs,1,1) )

## Weight distance from centre of mass by proportion of tracer volume and calculate total STD
lateral_STD_TVD = np.sqrt(np.sum( ( (lateral_distance_from_centre_TVD**2)*(TVD_depth_integrated_volume) \
                                    /( TVD_total_volume.reshape(noutputs,1,1)) ).reshape(noutputs,-1) ,axis=1 ))
lateral_STD_CE  = np.sqrt(np.sum( ( (lateral_distance_from_centre_CE **2)*( CE_depth_integrated_volume) \
                                    /(  CE_total_volume.reshape(noutputs,1,1)) ).reshape(noutputs,-1) ,axis=1 ))
lateral_STD_UW  = np.sqrt(np.sum( ( (lateral_distance_from_centre_UW **2)*( UW_depth_integrated_volume) \
                                    /(  UW_total_volume.reshape(noutputs,1,1)) ).reshape(noutputs,-1) ,axis=1 ))
lateral_STD_WM  = np.sqrt(np.sum( ( (lateral_distance_from_centre_WM **2)*( WM_depth_integrated_volume) \
                                    /(  WM_total_volume.reshape(noutputs,1,1)) ).reshape(noutputs,-1) ,axis=1 ))
# TRACER VERTICAL STANDARD DEVIATION
## Distance of vertical grid points from mean depth of tracer 
vertical_distance_from_centre_TVD= np.abs(dep0-dep_bar_TVD.reshape(noutputs,1))
vertical_distance_from_centre_CE = np.abs(dep0-dep_bar_CE .reshape(noutputs,1))
vertical_distance_from_centre_UW = np.abs(dep0-dep_bar_UW .reshape(noutputs,1))
vertical_distance_from_centre_WM = np.abs(dep0-dep_bar_WM .reshape(noutputs,1))

## Weight distance from centre of mass by proportion of tracer volume and calculate total STD
vertical_STD_TVD=np.sqrt(np.sum( ((vertical_distance_from_centre_TVD**2)*(TVD_horiz_integrated_volume) \
                                  /( TVD_total_volume.reshape(noutputs,1  ))  ).reshape(noutputs,-1),axis=1  ))
vertical_STD_CE =np.sqrt(np.sum( ((vertical_distance_from_centre_CE **2)*( CE_horiz_integrated_volume) \
                                  /(  CE_total_volume.reshape(noutputs,1  ))  ).reshape(noutputs,-1),axis=1  ))
vertical_STD_UW =np.sqrt(np.sum( ((vertical_distance_from_centre_UW **2)*(UW_horiz_integrated_volume) \
                                  /(  UW_total_volume.reshape(noutputs,1  ))  ).reshape(noutputs,-1),axis=1  ))
vertical_STD_WM =np.sqrt(np.sum( ((vertical_distance_from_centre_WM**2)*(WM_horiz_integrated_volume) \
                                  /(  WM_total_volume.reshape(noutputs,1  ))  ).reshape(noutputs,-1),axis=1  ))

################################################################################
