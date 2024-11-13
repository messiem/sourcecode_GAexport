%% START_SOURCECODE_GAEXPORT





%% --------------------------------------------- Set up GA toolbox ----------------------------------------------- %%

% Download the GA toolbox from https://github.com/messiem/toolbox_GrowthAdvection into directory toolbox_GrowthAdvection.
% The toolbox is only needed for using the plankton model (function ga_model_2P2Z_fromNsupply.m).
% Trajectories were calculated using ga_advection_ariane and can be reproduced by downloading GlobCurrent 
% from https://doi.org/10.48670/mds-00327 and following the toolbox instructions.

% Download data and save in data/:
% Nsupply_California_CCMP.nc from 
% export_CCMP3km_GlobCurrent_daily_REP.nc from 

addpath('toolbox_GrowthAdvection/')
addpath('toolbox_GrowthAdvection/utils/')



%% --------------------------------- Reproducing CCE-LTER dataset results ------------------------------------ %%

% The main dataset is data/CCE_SedTrap_Data.csv. 
% Columns 1-12 contains the CCE-LTER sediment trap dataset at the base of the euphotic zone 
%    (available at https://doi.org/10.6073/pasta/cdee03ef7b17c2a4027a4a8b33c5b09b). 
%    Carbon_flux_mgm2day is the POC flux measured at Depth_m (closest to Zeu); 
%    Carbon_flux_corr_mgm2day is the POC flux extrapolated to the base of the euphotic zone (Zeu) following Stukel et al. (2004).
% Columns 13-16 contains the along-trajectory outputs (water age, coastal Nsupply, modeled export) reproduced by ga_analyze_SedTraps_trajectories
% Column 17 contains the EF-OC export derived from ocean color corresponding to each sediment trap calculated from 5-day data 
%    (https://spg-satdata.ucsd.edu/wc_productivity/wc_productivity.htm) interpolated daily and shifted by 2 days 
%    (so that each day is based on primary production during the previous 5 days)
% Columns 18-19 contains the gridded GA output corresponding to each sediment trap, reproduced by ga_validate_gridded_export_vs_SedTraps


% Reproduce along-trajectory outputs (water age, coastal Nsupply, modeled export)
ga_analyze_SedTraps_trajectories

% Reproduce correlations with gridded GA outputs
ga_validate_gridded_export_vs_SedTraps

