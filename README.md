# North Atlantic water mass tracking experiments
## Overview
This is an archive of scripts necessary to reproduce the experiments detailed in our manuscript, submitted to _Geoscientific Model Development_ under the title  "Tracking water masses using passive-tracer transport in NEMO v3.4 with NEMOTAM: application to North Atlantic Deep Water and North Atlantic Subtropical Mode Water"

Our configuration is exactly as preserved at [DOI: 10.5281/zenodo.3459314](https://doi.org/10.5281/zenodo.3459314). It was compiled with the following keys:

`key_trabbl` `key_orca_r2` `key_lim2` `key_dynspg_flt` `key_diaeiv` `key_ldfslp` `key_traldf_c2d` `key_traldf_eiv` `key_dynldf_c3d` `key_zdftke` `key_zdftmx`  `key_mpp_mpi`  `key_mpp_rep` `key_nosignedzero` `key_tam` `key_diainstant`

## Directory structure
There are four subdirectories, relating to different model runs:

### `SPINUP_RUNS`
Contains scripts to produce the exact namelist files as used in our experiments to spin up the model for 950 years. 

### `TRAJECTORY_RUNS`
Contains scripts to produce a trajectory using the nonlinear model. Our experiments used a 400-year trajectory run in 50-year increments on 64 CPUs.

### `ADV_SCHEME_RUNS`
Contains the input file and namelists used to compare advection schemes (as in Section 2.3 of our submitted manuscript). 

### `WATER_MASS_RUNS`
Contains input files and namelists used to run water-mass-tracking experiments, as follows:

- `NASMW_tan.{nc,namelist}`: corresponding to tangent-linear run (described in Section 3.2)
- `NASMW_adj_above_MLD.{nc,namelist}`: corresponding to adjoint run above the mixed layer (described in Section 3.3)
- `NASMW_adj_below_MLD.{nc,namelist}`: corresponding to adjoint run below the mixed layer (described in Section 3.3)

- `SPNADW_tan.{nc,namelist}`: corresponding to tangent-linear run of subpolar outcropping NADW (described in Section 4.2.1)
- `ANADW_tan.{nc,namelist}`: corresponding to tangent-linear run of Arctic outcropping NADW (described in Section 4.2.2)
- `NADW_adj.{nc,namelist}`: corresponding to adjoint run of NADW (described in Section 4.3)
