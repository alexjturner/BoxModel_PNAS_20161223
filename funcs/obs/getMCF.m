%%% =======================================================================
%%% = getMCF.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the d13C observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getMCF( dataDir )

%%% Diagnostic
fprintf('   * MCF\n');

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;


%%% =======================================================================
%%% NOAA
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/mcf/NOAA/combined/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fName = 'HATS_global_MC.txt';
nHDR  = 85;
% Specify the site names and latitudes
sNames = {'alt','brw','nwr','kum','mlo','smo','cgo','psa','spo','mhd','thd','ush','sum'};
sLat   = [ 82.5, 71.3,40.05, 19.5, 19.5, 14.3,-40.7,-64.6,  -90,  53,   41, -54.8, 72.6];
sCol   = [    9,   13,   19,   21,   23,   25,   27,   31,   33,  15,   17,    29,   11];

%%% Load the data
dat   = importdata(sprintf('%s%s',dataDirU,fName),' ',nHDR);
dat   = dat.data;
tDatO = datenum(dat(:,1),dat(:,2),ones(size(dat(:,1))));

%%% Read the data
for i = 1:length(sNames)
    % Get the column for this site the data
    yDat = dat(:,sCol(i));
    tDat = tDatO;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_NOAA',sNames{i})) = yDat;
        out.tim.(sprintf('%s_NOAA',sNames{i})) = tDat;
        out.lat.(sprintf('%s_NOAA',sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% GAGE
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/mcf/GAGE/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = '%s-gage.mon';
sNames = { 'CGO', 'MHD', 'RPB', 'SMO', 'ORG'};
sLat   = [-40.68, 53.33, 13.17,-14.23, 45.00];
nHDR   = [     6,     6,     6,     6,     6];

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Load the data
    dat   = importdata(fName,' ',nHDR(i));
    dat   = dat.data;
    tDat  = datenum(dat(:,3),dat(:,2),ones(size(dat(:,1))));
    yDat  = dat(:,14);
    yDat(yDat == 0) = NaN;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_GAGE',sNames{i})) = yDat;
        out.tim.(sprintf('%s_GAGE',sNames{i})) = tDat;
        out.lat.(sprintf('%s_GAGE',sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% AGAGE
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/mcf/AGAGE/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = '%s-gcmd.mon';
sNames = { 'CGO', 'MHD', 'RPB', 'SMO', 'THD'};
sLat   = [-40.68, 53.33, 13.17,-14.23, 41.05];
nHDR   = [    16,    16,    16,    16,    16];

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Load the data
    dat   = importdata(fName,' ',nHDR(i));
    dat   = dat.data;
    tDat  = datenum(dat(:,3),dat(:,2),ones(size(dat(:,1))));
    yDat  = dat(:,11);
    yDat(yDat == 0) = NaN;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_AGAGE',sNames{i})) = yDat;
        out.tim.(sprintf('%s_AGAGE',sNames{i})) = tDat;
        out.lat.(sprintf('%s_AGAGE',sNames{i})) = sLat(i);
    end
end


end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================