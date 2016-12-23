%%% =======================================================================
%%% = getEdObs.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads Ed's observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getEdObs( dataDir, ajt_obs, St, tAvg, reread )

%%% Diagnostic
fprintf('\n *** UPDATE THE OBSERVATION STRUCTURE *** \n');

%%% Create the output structure
out_Ed = struct;


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/EJD_InputData_%4i-%4i_%s-%s.mat',...
                  reread.dir,reread.sYear,reread.eYear,reread.tRes,reread.tAvg);

%%% Load pre-existing data file
if ~reread.flag
    % Check if a file exists
    if exist(OutName, 'file') == 2
        fprintf('   * LOADING OLD OBS STRUCTURE\n');
        load(OutName);
        return % Don't need to read the data
    end
end


%%% =======================================================================
%%% CH4
%%% =======================================================================

%%% Diagnostic
fprintf('   * CH4\n');

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/ch4/dlugokencky/',dataDir);

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = 'zone_%s.mbl.ch4';
sNames = {'nh','sh'};
sLat   = [45,-45];
nHDR   = 0;

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Load the data
    dat  = importdata(fName,' ',nHDR);
    tDat = fracDateConv(dat(:,1));
    yDat = dat(:,2);
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_Ed',sNames{i})) = yDat;
        out.tim.(sprintf('%s_Ed',sNames{i})) = tDat;
        out.lat.(sprintf('%s_Ed',sNames{i})) = sLat(i);
    end
end

%%% Save the methane data
ch4_out = out;


%%% =======================================================================
%%% CH4
%%% =======================================================================

%%% Diagnostic
fprintf('   * CH4C13\n');

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/ch4c13/dlugokencky/',dataDir);

%%% Conversion from VPDB to NIWA scale (d13C(NIWA) = d13C + 0.169 permil)
NIWAconv = 0.169; % permil

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = 'zone_%s.mbl.ch4c13';
sNames = {'nh','sh'};
sLat   = [45,-45];
nHDR   = 0;

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Load the data
    dat  = importdata(fName,' ',nHDR);
    tDat = fracDateConv(dat(:,1));
    yDat = dat(:,2) - NIWAconv;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_Ed',sNames{i})) = yDat;
        out.tim.(sprintf('%s_Ed',sNames{i})) = tDat;
        out.lat.(sprintf('%s_Ed',sNames{i})) = sLat(i);
    end
end

%%% Save the d13C data
ch4c13_out = out;


%%% =======================================================================
%%% MCF
%%% =======================================================================

%%% Diagnostic
fprintf('   * MCF\n');

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/mcf/NOAA/combined/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fName = 'HATS_global_MC.txt';
nHDR  = 85;
% Specify the site names and latitudes
sNames = {'nh','sh'};
sLat   = [ 45, -45];
sCol   = [  3,   5];

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
        out.obs.(sprintf('%s_Ed',sNames{i})) = yDat;
        out.tim.(sprintf('%s_Ed',sNames{i})) = tDat;
        out.lat.(sprintf('%s_Ed',sNames{i})) = sLat(i);
    end
end

%%% Save the MCF data
mcf_out = out;


%%% =======================================================================
%%% Create the observation structure
%%% =======================================================================

%%% Fill the CH4 data
% Get obs
[oDat_NH, oDat_SH] = ReadObsStruct(St,tAvg,ch4_out);
eDat_NH            = ajt_obs.nh_ch4_err;
eDat_SH            = ajt_obs.sh_ch4_err;
% Throw out NaNs
ind_NH          = isnan(oDat_NH) | isnan(eDat_NH);
ind_SH          = isnan(oDat_SH) | isnan(eDat_SH);
oDat_NH(ind_NH) = NaN;
oDat_SH(ind_SH) = NaN;
eDat_NH(ind_NH) = NaN;
eDat_SH(ind_SH) = NaN;
% Put the data in the output structure
out_Ed.nh_ch4     = oDat_NH;
out_Ed.sh_ch4     = oDat_SH;
out_Ed.nh_ch4_err = eDat_NH;
out_Ed.sh_ch4_err = eDat_SH;

%%% Fill the CH4C13 data
% Get obs
[oDat_NH, oDat_SH] = ReadObsStruct(St,tAvg,ch4c13_out);
eDat_NH            = ajt_obs.nh_ch4c13_err;
eDat_SH            = ajt_obs.sh_ch4c13_err;
% Throw out NaNs
ind_NH          = isnan(oDat_NH) | isnan(eDat_NH);
ind_SH          = isnan(oDat_SH) | isnan(eDat_SH);
oDat_NH(ind_NH) = NaN;
oDat_SH(ind_SH) = NaN;
eDat_NH(ind_NH) = NaN;
eDat_SH(ind_SH) = NaN;
% Put the data in the output structure
out_Ed.nh_ch4c13     = oDat_NH;
out_Ed.sh_ch4c13     = oDat_SH;
out_Ed.nh_ch4c13_err = eDat_NH;
out_Ed.sh_ch4c13_err = eDat_SH;

%%% Fill the MCF data
% Get obs
[oDat_NH, oDat_SH] = ReadObsStruct(St,tAvg,mcf_out);
eDat_NH            = ajt_obs.nh_mcf_err;
eDat_SH            = ajt_obs.sh_mcf_err;
% Throw out NaNs
ind_NH          = isnan(oDat_NH) | isnan(eDat_NH);
ind_SH          = isnan(oDat_SH) | isnan(eDat_SH);
oDat_NH(ind_NH) = NaN;
oDat_SH(ind_SH) = NaN;
eDat_NH(ind_NH) = NaN;
eDat_SH(ind_SH) = NaN;
% Put the data in the output structure
out_Ed.nh_mcf     = oDat_NH;
out_Ed.sh_mcf     = oDat_SH;
out_Ed.nh_mcf_err = eDat_NH;
out_Ed.sh_mcf_err = eDat_SH;


%%% =======================================================================
%%% SAVE THIS OBSERVATION FILE
%%% =======================================================================

%%% Rename the output structure
clear out;
out = out_Ed;

%%% Save the structure
fprintf('   * SAVING OBS STRUCTURE\n');
if exist(OutName, 'file') == 2
    delete(OutName);
end
save(OutName,'out');

end

%%% Read the observation structure
function [ oDat_NH, oDat_SH ] = ReadObsStruct( t, tAvg, obs )

%%% Define parameters for the block averaging
fDays = 365.25; % Number of days in the block average
if strcmp(tAvg,'month') || strcmp(tAvg,'MONTH') || strcmp(tAvg,'monthly')
    fDays = fDays / 12;
end

%%% Get the site names
sNames = fieldnames(obs.obs);

%%% Initialize a vector for the observation times and 
tDat_NH = [];
yDat_NH = [];
tDat_SH = [];
yDat_SH = [];

%%% Get the data from the structure and sort it
% Get the data
for i = 1:length(sNames);
    % Get data for this site
    lat  = obs.lat.(sNames{i});
    tDat = obs.tim.(sNames{i});
    yDat = obs.obs.(sNames{i});
    % Remove the seasonal cycle
    yDat_noSeas = DeseasonalizeData(tDat,yDat,fDays);
    % Throw out NaNs
    ind = ~isnan(tDat) & ~isnan(yDat_noSeas);
    % Require at least a 5-year record
    tDat        = tDat(ind);
    yDat_noSeas = yDat_noSeas(ind);
    if ( abs(tDat(end) - tDat(1)) > 365.25*5 )
        % Sort by latitude
        if lat > 0 % NH
            tDat_NH = [tDat_NH;tDat];
            yDat_NH = [yDat_NH;yDat_noSeas];
        end
        if lat < 0 % SH
            tDat_SH = [tDat_SH;tDat];
            yDat_SH = [yDat_SH;yDat_noSeas];
        end
    end
end
% Sort
[~,ind] = sort(tDat_NH);
tDat_NH = tDat_NH(ind);
yDat_NH = yDat_NH(ind);
[~,ind] = sort(tDat_SH);
tDat_SH = tDat_SH(ind);
yDat_SH = yDat_SH(ind);

%%% Put the data on our temporal grid
oDat_NH = interp1(tDat_NH,yDat_NH,t);
oDat_SH = interp1(tDat_SH,yDat_SH,t);

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================
