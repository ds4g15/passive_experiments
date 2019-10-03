# North Atlantic water mass tracking experiments
## Overview
This is an archive of scripts necessary to reproduce the experiments detailed in our manuscript, submitted to _Geoscientific Model Development_ under the title  "Tracking water masses using passive-tracer transport in NEMO v3.4 with NEMOTAM: application to North Atlantic Deep Water and North Atlantic Subtropical Mode Water"

Our configuration is exactly as preserved at [DOI: 10.5281/zenodo.3459314](https://doi.org/10.5281/zenodo.3459314). The full model, along with the forcing and input files used in our experiments is available at [DOI: 10.5281/zenodo.1471702](https://doi.org/10.5281/zenodo.1471702) (file ORCA2_LIM_v3.4.tar). It was compiled with the following keys:

`key_trabbl` `key_orca_r2` `key_lim2` `key_dynspg_flt` `key_diaeiv` `key_ldfslp` `key_traldf_c2d` `key_traldf_eiv` `key_dynldf_c3d` `key_zdftke` `key_zdftmx`  `key_mpp_mpi`  `key_mpp_rep` `key_nosignedzero` `key_tam` `key_diainstant`

## Directory structure
There are five subdirectories, relating to different model runs:

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


### `DIAGNOSTICS`
Contains python and bash scripts used to produce the diagnostics in our manuscript, as follows:

- `produce_trajectory_climatology.sh` : a bash script which takes the raw NEMOTAM trajectory and produces a single netCDF file corresponding to its average year
- `rearrange_climatology_for_rebuild_nemo.py` : a python script called within `produce_trajectory_climatology.sh` which corrects for `ncra` re-arranging dimensions when time-averaging individual NEMO output tiles. Uncorrected, the tile averages cannot be stitched together with `rebuild_nemo`
- `climatology_stream_functions.py` : calculates the time-averaged barotropic and meridional overturning stream functions of the North Atlantic (as shown in Figs. 1 & 2)
- `climatology_NADW_properties` : calculates the location, volume and outcrop area of NADW over the climatology (as shown in Figs. 1, 2 & 5)
- `climatology_NASMW_properties` : calculates the location, volume and outcrop area of NASMW over the climatology (as shown in Figs. 1, 2 & 5)
- `compare_advection_schemes.py` : calculates the lateral and vertical spread of tracer when the same passive-tracer injection is propagated using different advection schemes. Also calculates the total volume of tracer with positive-valued and negative-valued concentration in these runs (as shown in Fig. 4).
- `diagnostics_tangent_linear.py` :  calculates the probability density that a water mass can be found at a given location or in a given TS class at a given time (as shown in Figs. 6 & 7 for NASMW, Figs. 11, 12, 13 & 14 for SPNADW and Figs. 15 & 16 for ANADW). Also calculates the average location and depth of a water mass based on its volume (as shown in Fig. 6 for NASMW, Figs. 11 & 12 for SPNADW and Fig. 15 for ANADW)
- `diagnostics_adjoint.py` : calculates the probability density that a water mass of a given age has originated from a given location or TS class (as shown in Figs. 8 & 9 for NASMW and Figs. 17 and 18 for NADW). Also calculates the probability distribution that a water mass is of a certain age (as in Fig. 10 for NASMW and Fig. 19 for NADW)