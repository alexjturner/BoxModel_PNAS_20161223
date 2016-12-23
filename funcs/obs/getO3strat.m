%%% =======================================================================
%%% = getO3strat.m
%%% = Alex Turner
%%% = 06/10/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the stratospheric ozone observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getO3strat( dataDir )

%%% Diagnostic
fprintf('   * STRATOSPHERIC OZONE\n');

%%% Append the directory onto the dataDir
dataDir = sprintf('%sobs/o3_strat/NOAA/',dataDir);

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = '%sTOTo3.txt';
% Specify the site names and latitudes
sNames = { 'ams',   'bdr',  'bis',  'bna', 'brw',  'car', 'fbk','hnx',    'ldr', 'mlo', 'ohp',   'pth',    'smo',  'wai'};
sLat   = [-89.99,40.01667,46.7667,36.2469,71.323,46.5667,64.859,36.31,-45.04278,19.533,43.931,-31.9222,-14.25022,37.8594];
nHDR   = 2;
fspec  = '%s %s %s %s %f %s %s';

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDir,sprintf(fNameS,sNames{i}));
    % Load the data
    fid  = fopen(fName);
    dat  = textscan(fid,fspec,'HeaderLines',nHDR,'Delimiter',' ','MultipleDelimsAsOne',1);
    fclose(fid);
    tDat = datenum(dat{1},'yymmdd');
    yDat = dat{5};
    % Make sure year is never greater than 2020
    for j = 1:length(tDat)
        [YYYY,MM,DD,HH,MIN,SEC] = datevec(tDat(j));
        if YYYY > 2020
            tDat(j) = datenum(YYYY-100,MM,DD,HH,MIN,SEC);
        end
    end
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat) & 0 < yDat & yDat < 999;
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_NOAA',sNames{i})) = yDat;
        out.tim.(sprintf('%s_NOAA',sNames{i})) = tDat;
        out.lat.(sprintf('%s_NOAA',sNames{i})) = sLat(i);
    end
end

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================