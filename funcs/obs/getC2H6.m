%%% =======================================================================
%%% = getC2H6.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the ethane observations.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- A structure containing the observation information.
%%% =======================================================================

function [ out ] = getC2H6( dataDir )

%%% Diagnostic
fprintf('   * C2H6\n');

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;


%%% =======================================================================
%%% WDCGG
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/c2h6/WDCGG/event/',dataDir);

%%% Get the filenames
files  = dir(sprintf('%s*.dat',dataDirU));
nFiles = length(files);
sNames = cell(nFiles,2);
sLat   = nan(nFiles,1);
for i = 1:nFiles
    tmp         = strsplit(files(i).name,'.');
    sNames{i,1} = files(i).name(1:3);
    sNames{i,2} = tmp{2};
    sNames{i,3} = tmp{4};
end

%%% Read the data
for i = 1:nFiles
    % Current filename
    fName = sprintf('%s%s',dataDirU,files(i).name);
    % Get the number of header lines and the latitude
    [~,tDat] = grep('-s','LATITUDE',fName);
    tDat     = strsplit(tDat.match{1});
    sLat(i)  = str2double(tDat{3});
    [~,tDat] = grep('-s','HEADER LINES',fName);
    tDat     = strsplit(tDat.match{1});
    nHDR     = str2double(tDat{4});
    % Load the data
    fspec = '%s %s %s %s %f %f %f %s %s %s';
    fid = fopen(fName);
    dat = textscan(fid,fspec,'HeaderLines',nHDR,'Delimiter',' ','MultipleDelimsAsOne',true);
    fclose(fid);
    yDat = dat{5};
    eDat = dat{8};
    % Convert dates to a usable form
    tDat = zeros(size(yDat));
    for j = 1:length(tDat)
        tDat(j) = datenum(dat{1}{j},'yyyy-mm-dd');
    end
    % Check flags
    ind = ones(size(yDat));
    for j = 1:length(yDat);
        ind(j) = strcmp(eDat{j},'...') || strcmp(eDat{j},'0') || strcmp(eDat{j},'1') || strcmp(eDat{j},'2');
    end
    yDat(~ind) = NaN;
    tDat(~ind) = NaN;
    % Check for negatives (non-physical)
    yDat(yDat <= 0) = NaN;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_%s_%s',sNames{i,1},sNames{i,2},sNames{i,3})) = yDat;
        out.tim.(sprintf('%s_%s_%s',sNames{i,1},sNames{i,2},sNames{i,3})) = tDat;
        out.lat.(sprintf('%s_%s_%s',sNames{i,1},sNames{i,2},sNames{i,3})) = sLat(i);
    end
end

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================