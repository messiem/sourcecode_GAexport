function ga_validate_gridded_export_vs_SedTraps


%% GA_VALIDATE_GRIDDED_EXPORT_VS_SEDTRAPS: match gridded export product with the CCE-LTER sediment trap dataset
% This function matches gridded export products (EF-OC and GA) with the CCE-LTER dataset, 
% and reproduces Fig. 2b & c in Messié et al. (2025) as well as results from Table 1.
% Outputs are compared against CCE_SedTrap_Data.csv.
%
% ga_validate_gridded_export_vs_SedTraps
% 
% Monique Messié, 2024
	
	

%% ------------------------- (1) Load datasets


% Load sediment trap dataset
SedTrap = readtable('data/CCE_SedTrap_Data.csv');
SedTrap_ini=SedTrap;		% keeping the original structure to double-check that outputs are correctly reproduced (see section 6)


% Compute deployment date and prepare GAgridded_Cproduction_mgm3day field
SedTrap.Deployment_Date = dateshift(SedTrap.Deployment_Datetime,'Start','Day');
SedTrap.Recovery_Date = dateshift(SedTrap.Recovery_Datetime,'Start','Day');
SedTrap.GAgridded_Cproduction_mgm3day = NaN(size(SedTrap,1),1);


% Load 3D daily GA export (downloaded from Zenodo, see start_sourcecode_GAexport.m)
% Cproduction is already calibrated in this file.
ncfile='data/export_CCMP3km_GlobCurrent_daily_REP.nc';
GA=struct(); GA.unit=struct();
for varname={'longitude','latitude','time','Cproduction'}, varname=varname{:};
	GA.(varname)=ncread(ncfile,varname);
	GA.unit.(varname)=ncreadatt(ncfile,varname,'units');
end
GA.Cproduction=permute(GA.Cproduction,[3 2 1]);			% matlab's ncread transposes dimensions
GA.time=days(GA.time)+datetime(1970,1,1); GA.unit.time='';		% moving to datetime format



%% ------------------------- (2) Match 3D GA product with sediment traps


for ipts=1:size(SedTrap,1)

	% get daily times from deployment to recovery
	itime = GA.time>=dateshift(SedTrap.Deployment_Date(ipts),'Start','Day') ...
		& GA.time<=dateshift(SedTrap.Recovery_Date(ipts),'Start','Day');

	% get longitude and latitude deployment/recovery boundaries
	lat_trap = sort( [SedTrap.Deployment_Latitude_deg(ipts),SedTrap.Recovery_Latitude_deg(ipts)] );
	lon_trap = sort( [SedTrap.Deployment_Longitude_deg(ipts),SedTrap.Recovery_Longitude_deg(ipts)] );

	% get corresponding indices in gridded GA product as data points within boundaries (or closest to the mean, if none)
	ilat = GA.latitude>=lat_trap(1) & GA.latitude<=lat_trap(2);
	if sum(ilat)==0
		diff_lat = abs(GA.latitude-mean(lat_trap));
		ilat = diff_lat==min(diff_lat); 
	end
	ilon = GA.longitude>=lon_trap(1) & GA.longitude<=lon_trap(2);
	if sum(ilon)==0
		diff_lon = abs(GA.longitude-mean(lon_trap));
		ilon = diff_lon==min(diff_lon); 
	end

	% average the corresponding GA data points
	SedTrap.GAgridded_Cproduction_mgm3day(ipts) = mean(GA.Cproduction(ilat,ilon,itime),'all','omitnan');

end



%% ------------------------- (3) Estimate Zeu GA export from surface export


% Compute modeled_export at Zeu, assuming it is representative of a 30-m surface layer
% The attenuation coefficient (0.0072) is based on Stukel et al. (2023) and used to correct in situ trap data in SedTrap
SedTrap.GAgridded_CZeu_mgm2day = SedTrap.GAgridded_Cproduction_mgm3day*30.*exp(-0.0072*(SedTrap.Zeu_m-30));



%% ------------------------- (4) Multilinear regression to combine both products


% Compute distance between deployment and recovery using Haversine formula
xlat=(SedTrap.Deployment_Latitude_deg-SedTrap.Recovery_Latitude_deg)/180*pi;
xlon=(SedTrap.Deployment_Longitude_deg-SedTrap.Recovery_Longitude_deg)/180*pi;
a = (sin(xlat/2)).^2 + cos(SedTrap.Deployment_Latitude_deg/180*pi).*cos(SedTrap.Recovery_Latitude_deg/180*pi).*(sin(xlon/2)).^2;
c = 2*atan2(sqrt(a),sqrt(1-a));
SedTrap.distdrift_km = 6378*c;

% Multilinear regression excluding data points where the drift is <2km (not representative of 12 km GA product)
% Note that regress.m requires the Statistics toolbox.
iok = ~isnan(SedTrap.Carbon_flux_corr_mgm2day) & ~isnan(SedTrap.export_EFOC_mgm2day) & ~isnan(SedTrap.GAgridded_CZeu_mgm2day) & SedTrap.distdrift_km>=2;
b=regress(SedTrap.Carbon_flux_corr_mgm2day(iok),[SedTrap.export_EFOC_mgm2day(iok),SedTrap.GAgridded_CZeu_mgm2day(iok)]);
disp(['Multilinear regression: combined export = ',num2str(b(1)),'*EF + ',num2str(b(2)),'*Czeu'])
SedTrap.export_combined_mgm2day=b(1)*SedTrap.export_EFOC_mgm2day+b(2)*SedTrap.GAgridded_CZeu_mgm2day;



%% ------------------------- (5) Results


disp(' ')
disp('Results (Table 1):')
disp('----------------')
iok = SedTrap.distdrift_km>=2; 	% exclude data points where the drift is <2km

% Statistics GA
corr=corrcoef(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.GAgridded_CZeu_mgm2day(iok));
inot0=SedTrap.GAgridded_CZeu_mgm2day>0;
corr10=corrcoef(log10(SedTrap.Carbon_flux_corr_mgm2day(iok & inot0)),log10(SedTrap.GAgridded_CZeu_mgm2day(iok & inot0)));
disp(['R² in situ export / modeled GA Zeu export = ',num2str(corr(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok)),'), ',...
	'log10 R² = ',num2str(corr10(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok & inot0)),'), ',...
	'RMSE = ',num2str(round(rmse(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.GAgridded_CZeu_mgm2day(iok))))])

% Statistics EF-OC
corr=corrcoef(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.export_EFOC_mgm2day(iok));
corr10=corrcoef(log10(SedTrap.Carbon_flux_corr_mgm2day(iok)),log10(SedTrap.export_EFOC_mgm2day(iok)));
disp(['R² in situ export / EF-OC export = ',num2str(corr(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok)),'), ',...
	'log10 R² = ',num2str(corr10(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok)),'), ',...
	'RMSE = ',num2str(round(rmse(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.export_EFOC_mgm2day(iok))))])

% Statistics combined
corr=corrcoef(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.export_combined_mgm2day(iok));
corr10=corrcoef(log10(SedTrap.Carbon_flux_corr_mgm2day(iok)),log10(SedTrap.export_combined_mgm2day(iok)));
disp(['R² in situ export / combined export = ',num2str(corr(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok)),'), ',...
	'log10 R² = ',num2str(corr10(1,2)^2,'%.2f'),' (N = ',num2str(sum(iok)),'), ',...
	'RMSE = ',num2str(round(rmse(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.export_combined_mgm2day(iok))))])

% Reproduce Fig. 2b & 2c
figure, set(gcf,'Units','Centimeters'); pos=get(gcf,'Position');
	pos(3)=12; pos(4)=7; set(gcf,'Position',pos)
	set(gcf,'PaperPositionMode','Auto','PaperUnits','Centimeters','PaperSize',[pos(3), pos(4)])
	fontsize=7;
axes('Position',[0.08 0.27 0.4 0.7],'FontSize',fontsize), hold on
	for ipts=find(iok)'
		plot([SedTrap.Carbon_flux_corr_mgm2day(ipts) SedTrap.Carbon_flux_corr_mgm2day(ipts)],...
		[SedTrap.export_EFOC_mgm2day(ipts) SedTrap.GAgridded_CZeu_mgm2day(ipts)],'Color',[0.7 0.7 0.7])
	end
	h1=scatter(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.export_EFOC_mgm2day(iok),50,'r','filled','MarkerEdgeColor','w');
	h2=scatter(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.GAgridded_CZeu_mgm2day(iok),50,'b','filled','MarkerEdgeColor','w');
	set(gca,'Xlim',[0 900],'YLim',[0 900],'XTick',0:200:800,'YTick',0:200:800)
	plot(xlim,xlim,'k')
	xlabel({'{\it In situ} Zeu export [mg m^{-2} d^{-1}]'}), ylabel('Gridded export products [mg m^{-2} d^{-1}]')
	legend([h1,h2],{'EF-OC','GA C_{Zeu}'},'FontSize',fontsize,'Location','northwest')
axes('Position',[0.58 0.27 0.4 0.7],'FontSize',fontsize), hold on
	ibelow30m=SedTrap.Zeu_m>=30;
	errSat=abs(SedTrap.export_EFOC_mgm2day*b(1)-SedTrap.GAgridded_CZeu_mgm2day*b(2))*sqrt(2);
	for ipts=find(iok)'
		plot([SedTrap.Carbon_flux_corr_mgm2day(ipts)-SedTrap.Carbon_flux_standard_error_corr_mgm2day(ipts)/2 ...
			SedTrap.Carbon_flux_corr_mgm2day(ipts)+SedTrap.Carbon_flux_standard_error_corr_mgm2day(ipts)/2],...
			[SedTrap.export_combined_mgm2day(ipts) SedTrap.export_combined_mgm2day(ipts)],...
			'Color',[0.7 0.7 0.7])
		plot([SedTrap.Carbon_flux_corr_mgm2day(ipts) SedTrap.Carbon_flux_corr_mgm2day(ipts)],...
			[SedTrap.export_combined_mgm2day(ipts)-errSat(ipts)/2 SedTrap.export_combined_mgm2day(ipts)+errSat(ipts)/2],...
			'Color',[0.7 0.7 0.7])
	end
	scatter(SedTrap.Carbon_flux_corr_mgm2day(iok),SedTrap.export_combined_mgm2day(iok),150-SedTrap.Zeu_m(iok),SedTrap.distdrift_km(iok),'filled','MarkerEdgeColor','w')
	scatter(SedTrap.Carbon_flux_corr_mgm2day(iok & ~ibelow30m),SedTrap.export_combined_mgm2day(iok & ~ibelow30m),150,'k')
	set(gca,'Xlim',[0 900],'YLim',[0 900],'XTick',0:200:800,'YTick',0:200:800)
	plot(xlim,xlim,'k')
	xlabel('{\it In situ} Zeu export [mgC m^{-2} d^{-1}]'), ylabel('Combined satellite export [mgC m^{-2} d^{-1}]')
	hbar4=colorbar('SouthOutside'); hbar4.Title.String='Drift [km]';
	hbar4.Position=[0.58 0.06 0.4 0.03];
print('-djpeg','-r300','figures/Fig2bc.jpg')



%% ------------------------- (6) Reproducibility


disp(' ')
disp('Reproducibility:')
disp('----------------')
disp('(differences below 1e-5 are due to precision differences)')
disp(['Max diff reproducing GA POC production: ',num2str(max(abs(SedTrap.GAgridded_Cproduction_mgm3day-SedTrap_ini.GAgridded_Cproduction_mgm3day)))])
disp(['Max diff reproducing GA Zeu export: ',num2str(max(abs(SedTrap.GAgridded_CZeu_mgm2day-SedTrap_ini.GAgridded_CZeu_mgm2day)))])

