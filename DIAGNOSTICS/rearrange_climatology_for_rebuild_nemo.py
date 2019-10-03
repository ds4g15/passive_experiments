import netCDF4 as nc
import numpy as np
import sys
''' 
Rearrange netCDF dimensions produced by trajectory_climatology.sh
so the files can be stitched together by TOOLS/rebuild_nemo
'''
fname=sys.argv[1]
        
infile=fname
outfile=('./TRAJ_'+fname)
print('writing to '+outfile)        
IN=nc.Dataset(infile)    
#Copy global attributes:
OUT=nc.Dataset(outfile,'w')
OUT.setncatts(IN.__dict__)
        
# Create dimensions in (I believe) the correct order:
OUT.createDimension('x', IN.dimensions['x'].size)
OUT.createDimension('y', IN.dimensions['y'].size)
OUT.createDimension('z', IN.dimensions['z'].size)
OUT.createDimension('t',None)
        
# Copy variables from input file
for name, variable in IN.variables.items():
    x = OUT.createVariable(name, variable.datatype, variable.dimensions)
    OUT[name][:] = IN[name][:]
    # copy variable attributes all at once via dictionary
    OUT[name].setncatts(IN[name].__dict__)
OUT.close()
    
