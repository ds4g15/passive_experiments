#!/bin/bash

################################################################################
# Create a climatology from NEMOTAM trajectory files
# Takes files TRAJECTORY_DIRECTORY/t_????????_????.nc and produces TRAJ_CLIMATOLOGY_Ny.nc

# First averages in time along each processor tile, then ensures (using python script) that
# output is compatible with TOOLS/rebuild_nemo and stitches in space, before
# concatenating in time.
################################################################################

###
REBUILD_NEMO_DIRECTORY=[PATH TO "TOOLS"]
TRAJECTORY_DIRECTORY=[PATH TO TRAJECTORY]
NSTEPS_PER_YEAR=5475
NSTEPS_TOTAL=328500
NPROCESSORS=64 # Number of processors trajectory was run on

for D in {000..002}; #{000..365}
do
    for T in `seq -f %04g 0 ${NPROCESSORS}`; 	        
    do
    echo "Tile ${T}, day ${D}"
    ncra `seq -f "${TRAJECTORY_DIRECTORY}/t_%08g_${T}.nc" $(((0+$((10#$D)))*15)) 5475 328500` "./CLIMATOLOGY_${D}_${T}.nc"
    python2.7 rearrange_climatology_for_rebuild_nemo.py "CLIMATOLOGY_${D}_${T}.nc"
    rm "./CLIMATOLOGY_${D}_${T}.nc"
    done;    
    ${REBUILD_NEMO_DIRECTORY}/rebuild_nemo ./TRAJ_CLIMATOLOGY_${D} 64;
    rm ./TRAJ_CLIMATOLOGY_${D}_????.nc
done

ncrcat -d t,0,,1 TRAJ_CLIMATOLOGY_???.nc TRAJ_CLIMATOLOGY_$(( 10#${NSTEPS_TOTAL}/10#${NSTEPS_PER_YEAR} )).nc
