%%% =======================================================================
%%% = getCH4ems.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Create the methane emissions and puts them onto our temporal 
%%% =        grid.  Can either use constant emissions, constant with a 30
%%% =        Tg/yr jump in 2007, or EDGAR v4.2FT2010.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the methane emissions.
%%% =======================================================================

function [ out ] = getCH4ems( St, tRes, dataDir )

%%% Diagnostic
fprintf('   * CH4\n');

%%% Which emissions do we want to use?
cfranken  = false;
edgar     = false;
const_ems = true;
if cfranken
    % Total methane source
    tot_ch4 = 550;  % Tg/yr
    frac_nh = 0.7;  % Fraction of emissions in the NH
    % CH4 Emissions
    ems_nh = zeros(size(St)) + tot_ch4 * frac_nh;
    ems_sh = zeros(size(St)) + tot_ch4 * (1 - frac_nh);
    % Add a step change of emissions after 2006
    ems_nh(St>datenum(2006,1,1)) = ems_nh(St>datenum(2006,1,1)) + 15/2;
    ems_sh(St>datenum(2006,1,1)) = ems_sh(St>datenum(2006,1,1)) + 15/2;
end
if edgar
    % Read in EDGAR
    fname       = sprintf('%s/ems/ch4/%s',dataDir,'edgarv42ft2010.xls');
    [dat, ~, ~] = xlsread(fname,'CH4_timeseries');
    tDat        = dat(8,:);     % Year
    tDat        = datenum(tDat,ones(size(tDat)),ones(size(tDat)));
    yDat        = dat(246,:);   % Global emissions (Gg/yr)
    % Interpolate it to my temporal grid
    if strcmp(tRes,'year') || strcmp(tRes,'YEAR') || strcmp(tRes,'yearly')
        fDays        = 365; % Number of days in the block average
        [tDat, yDat] = BlockAverage(tDat,yDat,ones(size(tDat)),fDays);
    end
    ind  = ~isnan(tDat) & ~isnan(yDat);
    oDat = interp1(tDat(ind),yDat(ind),St,'spline');
    oDat(isnan(oDat)) = nanmax(oDat);
    % Convert it to the units we need and add a natural source
    nat_source   = 215*1d3; % Tg/yr
    frac_nh_anth = 0.9;     % Fraction of anthropogenic emissions in the NH
    frac_nh_nat  = 0.4;     % Fraction of natural emissions in the NH
    ems_nh       = oDat*frac_nh_anth     + nat_source*frac_nh_nat;
    ems_sh       = oDat*(1-frac_nh_anth) + nat_source*(1-frac_nh_nat);
    ems_nh       = ems_nh * 1d-3;   % Tg/yr
    ems_sh       = ems_sh * 1d-3;   % Tg/yr
end
if const_ems
    % Total methane source (constant)
    tot_ch4 = 550;  % Tg/yr
    frac_nh = 0.75; % Fraction of emissions in the NH
    % CH4 Emissions
    ems_nh = zeros(size(St)) + tot_ch4 * frac_nh;
    ems_sh = zeros(size(St)) + tot_ch4 * (1 - frac_nh);
end
% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================