%%% =======================================================================
%%% = readData.m
%%% = Alex Turner
%%% = 06/06/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Reads a csv file that we made.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): baseName -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): obs      -- Structure containing the observations.
%%% =  ( 3): ems      -- Emission sources (and OH) for the box model.
%%% =  ( 4): mod      -- Structure containing the model results.
%%% =  ( 5): IC       -- Initial conditions.
%%% =======================================================================

function [ St, obs, mod, ems, IC ] = readData( baseName )

%%% Read the data
dat = importdata(sprintf(baseName,'Data'),',',2);
dat = dat.data;
% Fill the years
St = dat(:,1);
% Fill obs
obs.nh_ch4    = dat(:,2);
obs.nh_ch4c13 = dat(:,3);
obs.nh_mcf    = dat(:,4);
obs.nh_ch4    = dat(:,5);
obs.nh_ch4c13 = dat(:,6);
obs.nh_mcf    = dat(:,7);
% Fill obs error
obs.nh_ch4_err    = dat(:,8);
obs.nh_ch4c13_err = dat(:,9);
obs.nh_mcf_err    = dat(:,10);
obs.nh_ch4_err    = dat(:,11);
obs.nh_ch4c13_err = dat(:,12);
obs.nh_mcf_err    = dat(:,13);
% Fill obs
mod.nh_ch4    = dat(:,14);
mod.nh_ch4c13 = dat(:,15);
mod.nh_mcf    = dat(:,16);
mod.nh_ch4    = dat(:,17);
mod.nh_ch4c13 = dat(:,18);
mod.nh_mcf    = dat(:,19);
% Fill emissions
ems(:,1) = dat(:,20);       % NH CH4
ems(:,2) = dat(:,21);       % NH CH4C13
ems(:,3) = dat(:,22);       % NH MCF
ems(:,4) = dat(:,23);       % SH CH4
ems(:,5) = dat(:,24);       % SH CH4C13
ems(:,6) = dat(:,25);       % SH MCF
ems(:,7) = dat(:,26)/100+1; % NH OH
ems(:,8) = dat(:,27)/100+1; % SH OH

%%% Read the ICs
dat = importdata( sprintf(baseName,'ICs'),',',2);
dat = dat.data;
IC  = dat;

end


%%% =======================================================================
%%% = END
%%% =======================================================================
