%%% =======================================================================
%%% = getMCF.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads in the methylchloroform emissions and puts them onto our
%%% =        temporal grid.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the methylchloroform emissions.
%%% =======================================================================

function [ out ] = getMCFems( St, tRes, dataDir )

%%% Diagnostic
fprintf('   * MCF\n');

%%% Smooth the MCF emissions with a 5-year filter? (sensitivity test)
global smooth_MCF set_MCF_EMS MCF_EMS_val

%%% Read the Prinn/McCulloch emissions
% Append the directory onto the dataDir
dataDir = sprintf('%s/ems/mcf/',dataDir);
% Filenames and header lengths
fnames = {'mcf_Prinn.dat','mcf_McCulloch.dat'};
% Construct the output matrix
outS = nan(length(St),length(fnames));
% Read the data
for i = 1:length(fnames)
    % Load the data
    dat   = load(sprintf('%s%s',dataDir,fnames{i}));
    tDat  = datenum(dat(:,1),ones(size(dat(:,1))),ones(size(dat(:,1))));
    yDat  = dat(:,2);
    % Interpolate it to my temporal grid
    if strcmp(tRes,'year') || strcmp(tRes,'YEAR') || strcmp(tRes,'yearly')
        fDays        = 365.25; % Number of days in the block average
        % Smooth the emissions with a 5-year filter?
        if smooth_MCF
            fDays = fDays * 5;
        end
        [tDat, yDat] = BlockAverage(tDat,yDat,ones(size(tDat)),fDays);
        % Set post-2000 emissions to fixed values
        if set_MCF_EMS
            ind       = datenum(2000,1,1) <= tDat;
            yDat(ind) = MCF_EMS_val;
        end
    end
    ind  = ~isnan(tDat) & ~isnan(yDat);
    oDat = interp1(tDat(ind),yDat(ind),St,'spline');
    % Save the data
    outS(:,i) = oDat;
end
% Make the structure
out.prinn     = outS(:,1);
out.mcculloch = outS(:,2);

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================