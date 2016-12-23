%%% =======================================================================
%%% = getCH4C13ems.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Creates the isotopic composition for the Northern and Southern
%%% =        hemispheric methane emissions.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St      -- Our time vector.
%%% =  ( 2): tRes    -- String containing the temporal resolution.
%%% =  ( 3): dataDir -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure containing the d13C composition.
%%% =======================================================================

function [ out ] = getCH4C13ems( St, tRes, dataDir )

%%% Diagnostic
fprintf('   * CH4C13\n');

%%% Define the delta13C composition
% delta13C source
nh_ch4c13 = -52.15; % permil
sh_ch4c13 = -51.25; % permil
% delta13C composition
ems_nh = zeros(size(St)) + nh_ch4c13;
ems_sh = zeros(size(St)) + sh_ch4c13;
% Make the structure
out.nh = ems_nh;
out.sh = ems_sh;

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================
