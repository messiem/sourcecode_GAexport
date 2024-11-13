function ga_analyze_SedTraps_trajectories


%% GA_ANALYZE_SEDTRAPS_TRAJECTORIES: look at sediment trap source water characteristics 
% This function recalculates water_age, coastal Nsupply, and along-trajectory modeled export.
% Outputs are compared against CCE_SedTrap_Data.csv (differences are due to precision in the csv file).
%
% ga_analyze_SedTraps_trajectories


distcoast_lim=10;		% distance to coast used as "source data". This is a trade-off between too far from the coast, and unreliable/no current data near shore.



%% ------------------------- (1) Load datasets

% Sediment trap dataset
SedTrap=readtable('data/CCE_SedTrap_Data.csv');
SedTrap_ini=SedTrap;		% keeping the original structure to double-check that outputs are correctly reproduced (see section 7)

% Backward trajectories starting at each sediment trap deployment and recovery positions (start position is the last one ie time2D(:,end))
load('data/CCE_SedTrap_trajectories.mat','deployment','recovery')
[nb_pts,nb_time] = size(deployment.lon2D);			% nb_pts is the number of sediment traps (identical to size to SedTrap), nb_time the number of advection time steps

% Coastline
load('data/coastline_California.mat','coast_x','coast_y')	

% Coastal nitrate supply: download file https://data.mbari.org/products/satellite-derived/nitrate-supply/Nsupply_California_CCMP.nc and save in inputs/
ncfile='data/Nsupply_California_CCMP.nc';
Nsupply=struct(); Nsupply.unit=struct();
for varname={'latitude','time','Nsupply_coastal_pervolume'}, varname=varname{:};
	Nsupply.(varname)=ncread(ncfile,varname);
	Nsupply.unit.(varname)=ncreadatt(ncfile,varname,'units');
end
Nsupply.Nsupply_coastal_pervolume=Nsupply.Nsupply_coastal_pervolume';			% matlab's ncread transposes dimensions
Nsupply.time=days(Nsupply.time)+datetime(1970,1,1); Nsupply.unit.time='';		% moving to datetime format




%% ------------------------- (2) Remove trajectory data within the coast


% Sort coast data by latitude
sort_coast=sortrows([coast_y,coast_x]);
sort_coast=sort_coast(~isnan(sort_coast(:,1)),:);
[~,iok]=unique(sort_coast(:,1)); 
sort_coast=sort_coast(iok,:);

% Loop on all trajectories, deployment and recovery
for ii=1:nb_pts

	% find the coast longitude for each trajectory latitude (at the coast point closest to that latitude)
	coast_lon=interp1(sort_coast(:,1),sort_coast(:,2),deployment.lat2D(ii,:));
	% identify the first time the coastline is crossed
	iscoast_first=find(deployment.lon2D(ii,:)>coast_lon,1,'last');
	% set all latitudes/longitudes past that point at NaN, keeping the first coastal point for future interpolation
	deployment.lat2D(ii,1:iscoast_first-1)=NaN;
	deployment.lon2D(ii,1:iscoast_first-1)=NaN;
	deployment.distcoast2D(ii,1:iscoast_first-1)=NaN;
	% for the coastline distance, change the distance to a negative number within the coast
	deployment.distcoast2D(ii,iscoast_first)=-deployment.distcoast2D(ii,iscoast_first);

	% same calculations for recovery
	coast_lon=interp1(sort_coast(:,1),sort_coast(:,2),recovery.lat2D(ii,:));
	iscoast_first=find(recovery.lon2D(ii,:)>coast_lon,1,'last');
	recovery.lat2D(ii,1:iscoast_first-1)=NaN;
	recovery.lon2D(ii,1:iscoast_first-1)=NaN;
	recovery.distcoast2D(ii,1:iscoast_first-1)=NaN;
	recovery.distcoast2D(ii,iscoast_first)=-recovery.distcoast2D(ii,iscoast_first);

end



%% ------------------------- (3) Get trajectory origin and corresponding Nsupply and water_age



% Set up variables
deployment.source_time=NaT(nb_pts,1); 
recovery.source_time=NaT(nb_pts,1);
deployment.water_age=NaT(nb_pts,1)-datetime(2000,1,1); 
recovery.water_age=NaT(nb_pts,1)-datetime(2000,1,1);
deployment.unit=struct(); recovery.unit=struct();
for varname={'source_lon','source_lat','Nsupply_mmolCm3d'}, varname=varname{:};
	deployment.(varname)=NaN(nb_pts,1); recovery.(varname)=NaN(nb_pts,1);
end
deployment.unit.Nsupply=Nsupply.unit.Nsupply_coastal_pervolume;
recovery.unit.Nsupply=Nsupply.unit.Nsupply_coastal_pervolume;


% Loop on all trajectories, deployment and recovery
for ii=1:nb_pts

	% find the decimal position of the target distcoast_lim, identifying source waters
	icoast=max(ga_find_index(deployment.distcoast2D(ii,:),distcoast_lim));
	if ~isnan(icoast)
		% interpolate lon, lat and time trajectories to the source water location
		deployment.source_lon(ii)=interp1(1:nb_time,deployment.lon2D(ii,:),icoast);
		deployment.source_lat(ii)=interp1(1:nb_time,deployment.lat2D(ii,:),icoast);
		deployment.source_time(ii)=interp1(1:nb_time,deployment.time2D(ii,:),icoast);
		% compute water_age as deployment time minus source time
		deployment.water_age(ii)=deployment.time_trap(ii)-deployment.source_time(ii);
		% find the closest Nsupply time step
		dt=abs(Nsupply.time-deployment.source_time(ii));
		itimeN=find(dt==min(dt),1);
		% interpolate to the source latitude
		deployment.Nsupply_mmolCm3d(ii)=interp1(Nsupply.latitude,Nsupply.Nsupply_coastal_pervolume(:,itimeN),deployment.source_lat(ii));
	end

	% do the same for recovery
	icoast=max(ga_find_index(recovery.distcoast2D(ii,:),distcoast_lim));
	if ~isnan(icoast)
		recovery.source_lon(ii)=interp1(1:nb_time,recovery.lon2D(ii,:),icoast);
		recovery.source_lat(ii)=interp1(1:nb_time,recovery.lat2D(ii,:),icoast);
		recovery.source_time(ii)=interp1(1:nb_time,recovery.time2D(ii,:),icoast);
		recovery.water_age(ii)=recovery.time_trap(ii)-recovery.source_time(ii);
		dt=abs(Nsupply.time-recovery.source_time(ii));
		itimeN=find(dt==min(dt),1);
		recovery.Nsupply_mmolCm3d(ii)=interp1(Nsupply.latitude,Nsupply.Nsupply_coastal_pervolume(:,itimeN),recovery.source_lat(ii));
	end

end



%% ------------------------- (4) Average deployment and recovery information and match with sediment trap data


% Define SedTrap lon,lat,time as the mean of deployment and recovery positions
SedTrap.lon=mean([SedTrap.Deployment_Longitude_deg,SedTrap.Recovery_Longitude_deg],2,'omitnan');
SedTrap.lat=mean([SedTrap.Deployment_Latitude_deg,SedTrap.Recovery_Latitude_deg],2,'omitnan');
SedTrap.time=mean([SedTrap.Deployment_Datetime,SedTrap.Recovery_Datetime],2,'omitnan');

% Define SedTrap source_lat, water_age, Nsupply as the mean of deployment and recovery trajectories
% Set up new variables
for varname={'source_lat','source_time','water_age','Nsupply_mmolCm3d'}, varname=varname{:};
	SedTrap.(varname)=mean([deployment.(varname),recovery.(varname)],2,'omitnan');
end

% Retain information on how far apart deployment and recovery trajectories were
SedTrap.source_dlat=abs(deployment.source_lat-recovery.source_lat);
SedTrap.source_dtime=abs(deployment.water_age-recovery.water_age);



%% ------------------------- (5) Get the modeled export for each sediment trap as a function of Nsupply and water age


% Use the krill configuration used in Messié et al. (2022)
suff_krill={'gmax_big',0.6*0.6,'eZ',0.1*0.6,'mZ',0.05*16/106*0.6};					

% Run the model for each sediment trap
for ipts=1:nb_pts
	% plankton model as a function of Nsupply
	output=ga_model_2P2Z_fromNsupply(SedTrap.Nsupply_mmolCm3d(ipts),suff_krill{:},'nbdays_advec',90);
	% interpolate the modeled export to water_age
	SedTrap.GAalongtraj_Cproduction_mgCm3day_unscaled(ipts)=interp1(output.time,output.Cproduction,SedTrap.water_age(ipts));
end

% Compute GAalongtraj_export at Zeu, assuming it is representative of a 30-m surface layer
% The attenuation coefficient (0.0072) is based on Stukel et al. (2023) and used to correct in situ trap data in SedTrap
SedTrap.GAalongtraj_CZeu_mgCm2day_unscaled=SedTrap.GAalongtraj_Cproduction_mgCm3day_unscaled*30.*exp((30-SedTrap.Zeu_m)*0.0072);



%%  ------------------------- (6) Results


%% Statistics (Table 1)

% Exclude data points for which no coastal origin was identified
% or deployment and recovery backward trajectories are too different 
% (source latitude separated by more than 1°, or water age different by more than 45 days)
inocoast = isnan(SedTrap.source_lat); 
idiverge = SedTrap.source_dlat>1 | SedTrap.source_dtime>days(45);
iok = ~idiverge & ~inocoast; 
% For export estimates, also exclude points above 30m (the model is expected to be representative of the top ~ 30m)
ibelow30m = iok & SedTrap.Zeu_m>=30;

% Projection such that Carbon_flux_corr_mgm2day = a * GAalongtraj_CZeu_mgCm2day
% in order to calculate the scaling factor
scaling_factor = SedTrap.Carbon_flux_corr_mgm2day(ibelow30m)'/SedTrap.GAalongtraj_CZeu_mgCm2day_unscaled(ibelow30m)';
disp(['Scaling factor between in situ export and Zeu GA export: ',num2str(scaling_factor,'%.3f')])
disp(['Corresponding scaling factor between Cproduction and GA Zbig: ',...
	num2str(scaling_factor*(1-output.attributs.arg.epsilon)*output.attributs.arg.eZ*12,'%.3f'),' gC molC^{-1} d^{-1}'])
SedTrap.GAalongtraj_Cproduction_mgCm3day=SedTrap.GAalongtraj_Cproduction_mgCm3day_unscaled*scaling_factor;
SedTrap.GAalongtraj_CZeu_mgCm2day=SedTrap.GAalongtraj_CZeu_mgCm2day_unscaled*scaling_factor;

% Table 1 results
disp(' ')
disp('Results:')
disp('--------')

% Correlations
corr=corrcoef(SedTrap.Nsupply_mmolCm3d(iok),SedTrap.Carbon_flux_corr_mgm2day(iok));
disp(['R² in situ export / Nsupply = ',num2str(corr(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok)),')'])
corr=corrcoef(SedTrap.GAalongtraj_CZeu_mgCm2day(iok),SedTrap.Carbon_flux_corr_mgm2day(iok));
corr10=corrcoef(log10(SedTrap.Carbon_flux_corr_mgm2day(iok)),log10(SedTrap.GAalongtraj_CZeu_mgCm2day(iok)));
disp(['R² in situ export / modeled Zeu export = ',num2str(corr(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok)),'), ',...
		'log10 R² = ',num2str(corr10(1,2)^2,'%.2f'),', ',...
		'RMSE = ',num2str(round(rmse(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.GAalongtraj_CZeu_mgCm2day(iok))))])
disp(' ')
disp('Excluding points where Zeu is shallower than 30m:')
corr=corrcoef(SedTrap.GAalongtraj_CZeu_mgCm2day(ibelow30m),SedTrap.Carbon_flux_corr_mgm2day(ibelow30m));
corr10=corrcoef(log10(SedTrap.Carbon_flux_corr_mgm2day(ibelow30m)),log10(SedTrap.GAalongtraj_CZeu_mgCm2day(ibelow30m)));
disp(['R² in situ export / modeled Zeu export = ',num2str(corr(1,2)^2,'%.2f'),' (N = ',num2str(sum(ibelow30m)),'), ',...
		'log10 R² = ',num2str(corr10(1,2)^2,'%.2f'),', ',...
		'RMSE = ',num2str(round(rmse(SedTrap.Carbon_flux_corr_mgm2day(ibelow30m),SedTrap.GAalongtraj_CZeu_mgCm2day(ibelow30m))))])
	

%% Figure

% Compute distance between deployment and recovery using Haversine formula (for plot)
xlat=(SedTrap.Deployment_Latitude_deg-SedTrap.Recovery_Latitude_deg)/180*pi;
xlon=(SedTrap.Deployment_Longitude_deg-SedTrap.Recovery_Longitude_deg)/180*pi;
a = (sin(xlat/2)).^2 + cos(SedTrap.Deployment_Latitude_deg/180*pi).*cos(SedTrap.Recovery_Latitude_deg/180*pi).*(sin(xlon/2)).^2;
c = 2*atan2(sqrt(a),sqrt(1-a));
SedTrap.distdrift_km = 6378*c;

% Set up figure
figure, set(gcf,'Units','Centimeters'); pos=get(gcf,'Position');
pos(3)=18; pos(4)=6; set(gcf,'Position',pos)
set(gcf,'PaperPositionMode','Auto','PaperUnits','Centimeters','PaperSize',[pos(3), pos(4)])
fontsize=8;

% Map (Fig. 1a)
lonlim=[-128 -120]; latlim=[32 40.5];
axes('Position',[0.02 0.06 0.28 0.93],'FontSize',fontsize-1), hold on
for ii=1:nb_pts
	if ~isnan(deployment.water_age(ii)), plot(deployment.lon2D(ii,:),deployment.lat2D(ii,:),'Color',0.5*[1 1 1]), end
	if ~isnan(recovery.water_age(ii)), plot(recovery.lon2D(ii,:),recovery.lat2D(ii,:),'Color',0.7*[1 1 1]), end
end
scatter(SedTrap.lon(iok),SedTrap.lat(iok),150-SedTrap.Zeu_m(iok),SedTrap.Carbon_flux_corr_mgm2day(iok),'filled','MarkerEdgeColor','k')
plot(coast_x,coast_y,'k')
xlim(lonlim), ylim(latlim), clim([0 600])
hbar1=colorbar; hbar1.Title.String={'{\it In situ} Zeu export','[mgC m^{-2} d^{-1}]'};
hbar1.Position=[0.24 0.58 0.01 0.3];

% Export vs water age and Nsupply (Fig. 1b)
axes('Position',[0.37 0.6 0.25 0.34],'FontSize',fontsize-1), hold on
scatter(days(SedTrap.water_age(iok)),SedTrap.Carbon_flux_corr_mgm2day(iok),150-SedTrap.Zeu_m(iok),...
	SedTrap.Nsupply_mmolCm3d(iok),'filled','MarkerEdgeColor','k')
clim([0 50])
xlabel('Water age (days)'), ylabel({'{\it In situ} export','[mg m^{-2} d^{-1}]'})
hbar2=colorbar; hbar2.Title.String='mmolC m^{-3} d^{-1}'; ylabel(hbar2,'Nsupply','Rotation',-90,'FontSize',fontsize)
hbar2.Position=[0.63 0.6 0.01 0.34];

% Model output (Fig. 1c)
axes('Position',[0.37 0.12 0.25 0.34],'FontSize',fontsize-1), hold on
output=ga_model_2P2Z_fromNsupply(20,suff_krill{:},'nbdays_advec',90);
h1=plot(output.time,output.P_big,'k--','LineWidth',3);
h2=plot(output.time,output.Z_big,'k','LineWidth',3);
list_Nsupply=0:5:50;
couleurs=parula(length(list_Nsupply)+1);
for iNsupply=1:length(list_Nsupply)
	output=ga_model_2P2Z_fromNsupply(list_Nsupply(iNsupply),suff_krill{:},'nbdays_advec',100);
	plot(output.time,output.Z_big,'Linewidth',1.5,'Color',couleurs(iNsupply+1,:))
end
legend([h1,h2],{'P_{big}','Z_{big}'},'FontSize',fontsize-1)
xlabel('Elapsed time (days)'), ylabel({'Model output',['[',output.units.P_big,']']})
clim([min(list_Nsupply) max(list_Nsupply)])
hbar3=colorbar; hbar3.Title.String='mmolC m^{-3} d^{-1}'; ylabel(hbar3,'Nsupply','Rotation',-90,'FontSize',fontsize)
hbar3.Position=[0.63 0.12 0.01 0.34];

% Scatter plot (Fig. 2a)
axes('Position',[0.75 0.27 0.24 0.7],'FontSize',fontsize-1), hold on
scatter(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.GAalongtraj_CZeu_mgCm2day(iok),150-SedTrap.Zeu_m(iok),SedTrap.distdrift_km(iok),'filled','MarkerEdgeColor','w')
scatter(SedTrap.Carbon_flux_corr_mgm2day(iok & ~ibelow30m),SedTrap.GAalongtraj_CZeu_mgCm2day(iok & ~ibelow30m),150,'k')
plot(xlim,xlim,'k')
xlabel({'{\it In situ} Zeu export [mg m^{-2} d^{-1}]'}), ylabel('Along-trajectory C_{Zeu} [mg m^{-2} d^{-1}]')
hbar4=colorbar('SouthOutside'); hbar4.Title.String='Drift [km]';
hbar4.Position=[0.75 0.06 0.19 0.03];

% Save the figure
print('-djpeg','-r300','figures/results_alongtraj.jpg')


%%  ------------------------- (7) Reproducibility

disp(' ')
disp('Reproducibility:')
disp('----------------')
disp(['Max diff reproducing water_age: ',num2str(max(abs(days(SedTrap.water_age)-SedTrap_ini.water_age)))])
disp(['Max diff reproducing Nsupply: ',num2str(max(abs(SedTrap.Nsupply_mmolCm3d-SedTrap_ini.Nsupply_mmolCm3d)))])
disp(['Max diff reproducing alongtraj GA surface export: ',num2str(max(abs(SedTrap.GAalongtraj_Cproduction_mgCm3day-SedTrap_ini.GAalongtraj_Cproduction_mgCm3day)))])
disp(['Max diff reproducing alongtraj GA Zeu export: ',num2str(max(abs(SedTrap.GAalongtraj_CZeu_mgCm2day-SedTrap_ini.GAalongtraj_CZeu_mgCm2day)))])
disp(' ')
