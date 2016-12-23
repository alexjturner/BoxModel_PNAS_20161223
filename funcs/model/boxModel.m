%%% =======================================================================
%%% = boxModel.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): 2-box model for methane, delta13C, and methylchloroform.
%%% =  ( 2): Adapted from C. Frankenberg's box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): t      -- Box model time.
%%% =  ( 2): y      -- Concentrations at time t.
%%% =  ( 3): St     -- Out time vector.
%%% =  ( 4): S      -- Emission sources (and OH) for the box model.
%%% =  ( 5): params -- Structure with the parameters for the box model.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): dy -- Changes in the box model concentrations.
%%% =======================================================================

function [ dy ] = boxModel(t,y,St,S,params)

%%% Initialize dy
% - Column 1: 12CH4 in the Northern Hemisphere
% - Column 2: 13CH4 in the Northern Hemisphere
% - Column 3:   MCF in the Northern Hemisphere
% - Column 4: 12CH4 in the Southern Hemisphere
% - Column 5: 13CH4 in the Southern Hemisphere
% - Column 6:   MCF in the Southern Hemisphere
dy = zeros(6,1);

%%% Interpolate sources in time
S = interp1(St,S,t);

%%% Scale the rate constants (day^-1) by the OH factors
k_12ch4_NH = (S(7) * params.k_12ch4);   % NH
k_12ch4_SH = (S(8) * params.k_12ch4);   % SH
k_13ch4_NH = (S(7) * params.k_13ch4);   % NH
k_13ch4_SH = (S(8) * params.k_13ch4);   % SH
k_mcf_NH   = (S(7) * params.k_mcf  );   % NH
k_mcf_SH   = (S(8) * params.k_mcf  );   % SH

%%% Compute dy
% 12CH4
dy(1) = S(1) + (y(4) - y(1))/params.tau_NS - k_12ch4_NH*y(1);	% NH
dy(4) = S(4) + (y(1) - y(4))/params.tau_NS - k_12ch4_SH*y(4);	% SH
% 13CH4
dy(2) = S(2) + (y(5) - y(2))/params.tau_NS - k_13ch4_NH*y(2);   % NH
dy(5) = S(5) + (y(2) - y(5))/params.tau_NS - k_13ch4_SH*y(5);   % SH
% MCF
dy(3) = S(3) + (y(6) - y(3))/params.tau_NS - k_mcf_NH*y(3);     % NH
dy(6) = S(6) + (y(3) - y(6))/params.tau_NS - k_mcf_SH*y(6);     % SH

end


%%% =======================================================================
%%% = END
%%% =======================================================================