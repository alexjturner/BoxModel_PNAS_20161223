%%% =======================================================================
%%% = getCH4C13.m
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

function [ out ] = getCH4C13( dataDir )

%%% Diagnostic
fprintf('   * CH4C13\n');

%%% Create the output structure
out = struct;
out.obs = struct;
out.tim = struct;
out.lat = struct;

%%% Conversion from VPDB to NIWA scale (d13C(NIWA) = d13C + 0.169 permil)
NIWAconv = 0.169; % permi


%%% =======================================================================
%%% NOAA (INSTAAR)
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/ch4c13/INSTAAR/month/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = 'ch4c13_%s_surface-flask_1_sil_month.txt';
% Make the site list with: "ls month/ch4c13_*_month.txt | cut -d'_' -f2 > site_list.csv"
fid = fopen(sprintf('%s/../site_list.csv',dataDirU));
dat = textscan(fid,'%s %f','delimiter',',');
fclose(fid);
sNames = dat{1};
sLat   = dat{2};

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Get the number of header lines
    fid  = fopen(fName);
    nHDR = textscan(fid, '%s', 1,'delimiter','\n');
    nHDR = strsplit(char(nHDR{1}));
    nHDR = str2double(nHDR(end));
    fclose(fid);
    % Load the data
    dat  = importdata(fName,' ',nHDR);
    dat  = dat.data;
    tDat = datenum(dat(:,1),dat(:,2),ones(size(dat(:,1))));
    yDat = dat(:,3) - NIWAconv;
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_INSTAAR',sNames{i})) = yDat;
        out.tim.(sprintf('%s_INSTAAR',sNames{i})) = tDat;
        out.lat.(sprintf('%s_INSTAAR',sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% Heidelberg
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/ch4c13/Heidelberg/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = 'd13Cdata_%s.txt';
sNames = {'Alert','Izana','Neumayer'};
sLat   = [ 82.45,  28.30,     -70.65];
nHDR   = [    66,     66,         68];

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Load the data
    dat = importdata(fName,'\t',nHDR(i));
    dat = dat.data;
    if (nHDR(i) == 68) % Neumayer file format is slightly different
        tDat = datenum(dat(:,1),dat(:,2),dat(:,3)); % Gives year, month, day
        yDat = dat(:,6);
    else
        tDat = mean([dat(:,1),dat(:,2)],2); % Gives a range of fractional dates
        tDat = fracDateConv(tDat);
        yDat = dat(:,5);
    end
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_Heidelberg',sNames{i})) = yDat;
        out.tim.(sprintf('%s_Heidelberg',sNames{i})) = tDat;
        out.lat.(sprintf('%s_Heidelberg',sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% UCI
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/ch4c13/UCI/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = 'use_%sch4.dat';
sNames = {'nwr','mdo'};
sLat   = [ 41.0,  35.0];
nHDR   = [   29,    24];

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Load the data
    if strcmp(sNames{i},'nwr')
        fspec   = '%s %s %s %f %s %s %s %s %s';
        yDatInd = 4;
    else
        fspec   = '%s %s %s %s %f %s %s %s %s %s %s';
        yDatInd = 5;
    end
    fid = fopen(fName);
    dat = textscan(fid,fspec,'HeaderLines',nHDR(i),'Delimiter','\t');
    fclose(fid);
    dDat = dat{1};
    yDat = dat{yDatInd} - NIWAconv;
    % Convert dates to a usable form
    tDat = zeros(size(dDat));
    for j = 1:length(tDat)
        tDat(j) = datenum(dDat{j},'mm/dd/yyyy');
    end
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_UCI',sNames{i})) = yDat;
        out.tim.(sprintf('%s_UCI',sNames{i})) = tDat;
        out.lat.(sprintf('%s_UCI',sNames{i})) = sLat(i);
    end
end


%%% =======================================================================
%%% UW
%%% =======================================================================

%%% Append the directory onto the dataDir
dataDirU = sprintf('%sobs/ch4c13/UW/',dataDir);

%%% Define the site names, header lengths, and latitudes
% Filename structure
fNameS = 'use_%sExtract.txt';
sNames = {'CG','CP','FD','MI','ML','NZ','PB','SM'};
sLat   = [ -41,  48,  50,   7,  19, -41,  71, -14];
nHDR   = [  18,  18,  16,  18,  18,  18,  18,  18];

%%% Read the data
for i = 1:length(sNames)
    % Current filename
    fName = sprintf('%s%s',dataDirU,sprintf(fNameS,sNames{i}));
    % Load the data
    fspec   = '%s %s %s %s %f';
    yDatInd = 5;
    if strcmp(sNames{i},'FD');
        fspec   = '%s %s %s %f';
        yDatInd = 4;
    end
    fid   = fopen(fName);
    dat = textscan(fid,fspec,'HeaderLines',nHDR(i),'Delimiter',' ','MultipleDelimsAsOne',1);
    fclose(fid);
    dDat = dat{2};
    yDat = dat{yDatInd};
    % Convert dates to a usable form
    tDat = zeros(size(dDat));
    for j = 1:length(tDat)
        date_info = str2double(strsplit(dDat{j},'-'));
        tDat(j)   = datenum(1900+date_info(3),date_info(2),date_info(1));
    end
    % Remove NaNs
    ind  = ~isnan(yDat) & ~isnan(tDat);
    if sum(ind) > 0
        tDat = tDat(ind);
        yDat = yDat(ind);
        % Put the data in a structure
        out.obs.(sprintf('%s_UW',sNames{i})) = yDat;
        out.tim.(sprintf('%s_UW',sNames{i})) = tDat;
        out.lat.(sprintf('%s_UW',sNames{i})) = sLat(i);
    end
end

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================
