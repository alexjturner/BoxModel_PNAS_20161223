%%% =======================================================================
%%% = makeObs.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Aggregates the observations into Northern and Southern
%%% =        hemispheric averages and puts them onto my temporal grid.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St         -- Our time vector.
%%% =  ( 2): tAvg       -- String defining our averaging period.
%%% =  ( 3): ch4_obs    -- Structure containing the methane obs.
%%% =  ( 4): ch4c13_obs -- Structure containing the d13C obs.
%%% =  ( 5): mcf_obs    -- Structure containing the methylchloroform obs.
%%% =  ( 6): dataDir    -- Directory containing the data.
%%% =  ( 7): reread     -- Structure the says if we'll re-read the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = makeObs( St, tAvg, ch4_obs, ch4c13_obs, mcf_obs, dataDir, reread )

%%% Diagnostic
fprintf('\n *** MAKING THE OBSERVATION STRUCTURE *** \n');

%%% Create the output structure
out = struct;

%%% Get the years
yrs = datevec(St);
yrs = yrs(:,1);


%%% =======================================================================
%%% HAVE WE READ THIS DATA BEFORE?
%%% =======================================================================

%%% Build the filename
OutName = sprintf('%sobs/StoredData/InputData_%4i-%4i_%s-%s.mat',...
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

%%% Read the observation structure
fprintf('   * CH4\n');
[oDat_NH,oDat_SH,eDat_NH,eDat_SH] = ReadObsStruct(St,tAvg,ch4_obs);

%%% Make sure the errors aren't overly optimistic
min_err = [9999, 2;... % [year and ppb]
           1999, 3];
for i = 1:size(min_err,1);
    ind          = yrs < min_err(i,1) & eDat_NH < min_err(i,2);
    eDat_NH(ind) = min_err(i,2);
    ind          = yrs < min_err(i,1) & eDat_SH < min_err(i,2);
    eDat_SH(ind) = min_err(i,2);
end
% Make sure older obs are always less certain than new obs
for i = length(St)-1:-1:1
    if ~isnan(eDat_NH(i)) && eDat_NH(i) < nanmax(eDat_NH(i+1:end))
        eDat_NH(i) = nanmax(eDat_NH(i+1:end));
    end
    if ~isnan(eDat_SH(i)) && eDat_SH(i) < nanmax(eDat_SH(i+1:end))
        eDat_SH(i) = nanmax(eDat_SH(i+1:end));
    end
end

%%% Put the data in the output structure
out.nh_ch4     = oDat_NH;
out.sh_ch4     = oDat_SH;
out.nh_ch4_err = eDat_NH;
out.sh_ch4_err = eDat_SH;


%%% =======================================================================
%%% CH4C13
%%% =======================================================================

%%% Read the observation structure
fprintf('   * CH4C13\n');
[oDat_NH,oDat_SH,eDat_NH,eDat_SH] = ReadObsStruct(St,tAvg,ch4c13_obs);

%%% Make sure the errors aren't overly optimistic
min_err = [9999, 0.03;... % [year and permil]
           1998, 0.03;...
           1988, 1.10];
for i = 1:size(min_err,1);
    ind          = yrs < min_err(i,1) & eDat_NH < min_err(i,2);
    eDat_NH(ind) = min_err(i,2);
    ind          = yrs < min_err(i,1) & eDat_SH < min_err(i,2);
    eDat_SH(ind) = min_err(i,2);
end
% Make sure older obs are always less certain than new obs
for i = length(St)-1:-1:1
    if ~isnan(eDat_NH(i)) && eDat_NH(i) < nanmax(eDat_NH(i+1:end))
        eDat_NH(i) = nanmax(eDat_NH(i+1:end));
    end
    if ~isnan(eDat_SH(i)) && eDat_SH(i) < nanmax(eDat_SH(i+1:end))
        eDat_SH(i) = nanmax(eDat_SH(i+1:end));
    end
end

%%% Put the data in the output structure
out.nh_ch4c13     = oDat_NH;
out.sh_ch4c13     = oDat_SH;
out.nh_ch4c13_err = eDat_NH;
out.sh_ch4c13_err = eDat_SH;


%%% =======================================================================
%%% MCF
%%% =======================================================================

%%% Read the observation structure
fprintf('   * MCF\n');
[oDat_NH,oDat_SH,eDat_NH,eDat_SH] = ReadObsStruct(St,tAvg,mcf_obs);

%%% Make sure the errors aren't overly optimistic
[eDat_NH,eDat_SH] = ReadMCFerr(St,tAvg,eDat_NH,eDat_SH,dataDir);

%%% Put the data in the output structure
out.nh_mcf     = oDat_NH;
out.sh_mcf     = oDat_SH;
out.nh_mcf_err = eDat_NH;
out.sh_mcf_err = eDat_SH;


%%% =======================================================================
%%% SAVE THIS OBSERVATION FILE
%%% =======================================================================

%%% Save the structure
fprintf('   * SAVING OBS STRUCTURE\n');
if exist(OutName, 'file') == 2
    delete(OutName);
end
save(OutName,'out');


end

%%% Read the observation structure
function [ oDat_NH, oDat_SH, eDat_NH, eDat_SH ] = ReadObsStruct( t, tAvg, obs )

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
iDat_NH = [];
tDat_SH = [];
yDat_SH = [];
iDat_SH = [];

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
            iDat_NH = [iDat_NH;repmat(i,size(tDat))];
        end
        if lat < 0 % SH
            tDat_SH = [tDat_SH;tDat];
            yDat_SH = [yDat_SH;yDat_noSeas];
            iDat_SH = [iDat_SH;repmat(i,size(tDat))];
        end
    end
end
% Sort
[~,ind] = sort(tDat_NH);
tDat_NH = tDat_NH(ind);
yDat_NH = yDat_NH(ind);
iDat_NH = iDat_NH(ind);
[~,ind] = sort(tDat_SH);
tDat_SH = tDat_SH(ind);
yDat_SH = yDat_SH(ind);
iDat_SH = iDat_SH(ind);

%%% Bootstrap the mean and error
nBoot    = 50;
tDatO_NH = [];
yDatO_NH = [];
tDatO_SH = [];
yDatO_SH = [];
for i = 1:nBoot
    
    %%% NH
    uniqID = unique(iDat_NH);   % Get the unique IDs
    nDraw  = length(uniqID);    % How many different timeseries do we have?
    % Randomly sample from our sites
    tDat = [];
    yDat = [];
    for j = 1:nDraw
        ind  = randsample(uniqID,1) == iDat_NH;
        tDat = [tDat;tDat_NH(ind)];
        yDat = [yDat;yDat_NH(ind)];
    end
    % Block average
    [tDat, yDat] = BlockAverage(tDat,yDat,ones(size(tDat)),fDays);
    % Store the data
    tDatO_NH = [tDatO_NH;tDat];
    yDatO_NH = [yDatO_NH;yDat];
    
    %%% SH
    uniqID = unique(iDat_SH);   % Get the unique IDs
    nDraw  = length(uniqID);    % How many different timeseries do we have?
    % Randomly sample from our sites
    tDat = [];
    yDat = [];
    for j = 1:nDraw
        ind  = randsample(uniqID,1) == iDat_SH;
        tDat = [tDat;tDat_SH(ind)];
        yDat = [yDat;yDat_SH(ind)];
    end
    % Block average
    [tDat, yDat] = BlockAverage(tDat,yDat,ones(size(tDat)),fDays);
    % Store the data
    tDatO_SH = [tDatO_SH;tDat];
    yDatO_SH = [yDatO_SH;yDat];
    
end

%%% Block average to get the mean and error
% Remove NaNs
ind      = ~isnan(tDatO_NH) & ~isnan(yDatO_NH);
tDatO_NH = tDatO_NH(ind);
yDatO_NH = yDatO_NH(ind);
ind      = ~isnan(tDatO_SH) & ~isnan(yDatO_SH);
tDatO_SH = tDatO_SH(ind);
yDatO_SH = yDatO_SH(ind);
% Sort
[~,ind]  = sort(tDatO_NH);
tDatO_NH = tDatO_NH(ind);
yDatO_NH = yDatO_NH(ind);
[~,ind]  = sort(tDatO_SH);
tDatO_SH = tDatO_SH(ind);
yDatO_SH = yDatO_SH(ind);
% Block average
[tDat_NH, yDat_NH, eDat_NH] = BlockAverage_AltError(tDatO_NH,yDatO_NH,ones(size(tDatO_NH)),fDays);
[tDat_SH, yDat_SH, eDat_SH] = BlockAverage_AltError(tDatO_SH,yDatO_SH,ones(size(tDatO_SH)),fDays);

%%% Put the data on our temporal grid
oDat_NH = interp1(tDat_NH,yDat_NH,t);
eDat_NH = interp1(tDat_NH,eDat_NH,t);
oDat_SH = interp1(tDat_SH,yDat_SH,t);
eDat_SH = interp1(tDat_SH,eDat_SH,t);

end

%%% Read the NH/SH errors from the NOAA network
function [ eOut_NH, eOut_SH ] = ReadMCFerr( t, tAvg, eDat_NH, eDat_SH, dataDir )

%%% Define the minimum error
min_err = 0.1; % ppt

%%% Append the directory onto the dataDir
dataDir = sprintf('%sobs/MCF/NOAA/combined/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fName = 'HATS_global_MC.txt';
nHDR  = 85;

%%% Load the data
dat     = importdata(sprintf('%s%s',dataDir,fName),' ',nHDR);
dat     = dat.data;
tDatO   = datenum(dat(:,1),dat(:,2),ones(size(dat(:,1))));
eMCF_NH = dat(:,4); % NH uncertainty
eMCF_SH = dat(:,6); % SH uncertainty

%%% Get the observations onto my temporal grid
tDat_NH = tDatO;
tDat_SH = tDatO;
% Remove NaNs
ind_NH  = ~isnan(tDat_NH) & ~isnan(eMCF_NH);
ind_SH  = ~isnan(tDat_SH) & ~isnan(eMCF_SH);
tDat_NH = tDat_NH(ind_NH);
eMCF_NH = eMCF_NH(ind_NH);
tDat_SH = tDat_SH(ind_SH);
eMCF_SH = eMCF_SH(ind_SH);
% Put it on my grid (Smooth it first)
fDays = 365; % Number of days in the block average
if strcmp(tAvg,'month') || strcmp(tAvg,'MONTH') || strcmp(tAvg,'monthly')
    fDays = fDays / 12;
end
[tDat_NH, eMCF_NH] = BlockAverage(tDat_NH,eMCF_NH,ones(size(tDat_NH)),fDays);
[tDat_SH, eMCF_SH] = BlockAverage(tDat_SH,eMCF_SH,ones(size(tDat_SH)),fDays);
eOut_NH = interp1(tDat_NH,eMCF_NH,t);
eOut_SH = interp1(tDat_SH,eMCF_SH,t);

%%% Assume errors before the first obs are twice the max
eOut_NH(isnan(eOut_NH)) = 2*nanmax(eOut_NH);
eOut_SH(isnan(eOut_SH)) = 2*nanmax(eOut_SH);

%%% Compare this to our previously estimated errors
eOut_NH = nanmax([eOut_NH,eDat_NH],[],2);
eOut_SH = nanmax([eOut_SH,eDat_SH],[],2);

%%% Make sure this is more than the absolute min value
eOut_NH(eOut_NH < min_err) = min_err;
eOut_SH(eOut_SH < min_err) = min_err;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================