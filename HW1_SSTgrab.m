

% Load in data
addpath('/Users/gyl/Downloads/02_Research/01_Data')
load('Hist_FluxAll_Masked.mat')
load('Historical_SST_masked_Avg.mat')

% Load in MLD Climatology and get point Lat/Lon
% Select Point to model at
lonf = 330; % [0 360]
latf = 50; % [-90 90]
[loi,lai] = findcoords(lonf,latf,2,{LON,LAT});

% Grab FSNS
fsns = tmean_FSNS(loi,lai)

% Grab mean temp
t = SST(loi,lai,:);
t = squeeze(t);
t = reshape(t,12,length(t)/12);
tm = nanmean(t,2)-273.15;
