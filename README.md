# sourcecode_GAexport

This toolbox contains the programs used to reproduce results and figures from Messié et al. (submitted).   
The programs are written for Matlab. The Statistics toolbox is required for multilinear regression in `ga_validate_gridded_export_vs_SedTraps` (l. 93).

* * *

### Get started ###

See script:  

	start_sourcecode_GAexport

The code generates figures, saved in `figures/`, reproducing results from Fig. 1 and 2 from the paper. The figures are simplified versions of the paper figures, with Matlab native colorbars and functions. Please contact me if you want to be able to reproduce the exact same figures as in the paper. The programs also output statistics from Table 1.
  
### Description of functions ### 

`ga_analyze_SedTraps_trajectories`: reproduces results using backward Lagrangian trajectories originating from the CCE-LTER sediment trap positions. Trajectories are available in `data/CCE_SedTrap_trajectories.mat`.  
`ga_validate_gridded_export_vs_SedTraps`: reproduces the validation of the GA gridded product published on Zenodo, against the CCE-LTER dataset. The program also displays results using the EF-OC dataset (available online at https://spg-satdata.ucsd.edu/wc_productivity/wc_productivity.htm). 

### Description of data ###

**Data to download and save in the `data/` folder:** 

`Nsupply_California_CCMP_REP.nc`: download from https://doi.org/10.5281/zenodo.14641977 (nutrient supply used to force the plankton model)   
`export_CCMP3km_GlobCurrent_daily_REP.nc`: download from https://doi.org/10.5281/zenodo.14084208 (output of the GA model, forced by Nsupply and GlobCurrent).

**CCE-LTER datasets:**  

`CCE_SedTrap_Data.csv`:  
Columns 1-14 contains the CCE-LTER sediment trap dataset at the base of the euphotic zone (available at https://doi.org/10.6073/pasta/cdee03ef7b17c2a4027a4a8b33c5b09b). Carbon_flux_mgm2day is the POC flux measured at Depth_m (closest to Zeu); Carbon_flux_corr_mgm2day is the POC flux extrapolated to the base of the euphotic zone (Zeu) following Stukel et al. (2004). Carbon_flux_standard_error_mgm2day and Carbon_flux_standard_error_corr_mgm2day are the corresponding standard errors.  
Columns 15-18 contains the along-trajectory outputs (water age, coastal Nsupply, modeled export) reproduced by `ga_analyze_SedTraps_trajectories`.  
Column 19 contains the EF-OC export derived from ocean color corresponding to each sediment trap calculated from 5-day data (https://spg-satdata.ucsd.edu/wc_productivity/wc_productivity.htm) interpolated daily and shifted by 2 days (so that each day is based on primary production during the previous 5 days).  
Columns 20-21 contains the gridded GA output corresponding to each sediment trap, reproduced by `ga_validate_gridded_export_vs_SedTraps`  

`CCE_SedTrap_trajectories.mat`:  
backward Lagrangian trajectories originating from the CCE-LTER sediment trap positions, computed using a 2D custom version of Ariane. (Please contact me if you want to reproduce these trajectories).


* * *

### Reference ###

Please refer this paper when using these scripts:  

Messié, M., C.L. Huffard, M. Stukel, and H.A. Ruhl (2024). **Spatial and temporal interplay between oceanic circulation and biological production in shaping carbon export off the California coast**.  *Geophysical Research Letters*, submitted. Preprint available at https://doi.org/10.22541/essoar.173272956.60470361/v1

* * *

### Contact ###

monique@mbari.org

Do not hesitate to contact me if you cannot run the code in `start_sourcecode_GAexport`, if you notice bugs, or if you need help implementing the code for a custom application. Note that the GA model is separately available at https://github.com/messiem/toolbox_GrowthAdvection.
