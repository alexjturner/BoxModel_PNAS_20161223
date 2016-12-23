%%% =======================================================================
%%% = boxModel_wrapper.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): A wrapper for the box model.
%%% =  ( 2): Precomputes terms for the box model for speed.
%%% =  ( 3): Runs the box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St     -- Our time vector.
%%% =  ( 2): S      -- Emission sources (and OH) for the box model.
%%% =  ( 3): IC     -- Initial conditions for the box model.
%%% =  ( 4): params -- Structure with parameters for the box model.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): out -- Structure with the simulated concentrations.
%%% =======================================================================

function [ out ] = boxModel_wrapper(St,S,IC,params)

%%% Set up the emissions for the box model
% Convert CH4 & MCF emissions to units of per day
S(:,[1,3,4,6]) = S(:,[1,3,4,6]) / params.YrToDay;
% 12CH4
S(:,1) = 2/params.mm_ch4*S(:,1);    % NH
S(:,4) = 2/params.mm_ch4*S(:,4);    % SH
% 13CH4
S(:,2) = S(:,1).*(1 + S(:,2)/1000); % NH (precompute S1*S2)
S(:,5) = S(:,4).*(1 + S(:,5)/1000); % SH (precompute S4*S5)
% MCF
S(:,3) = 2/params.mm_mcf*S(:,3);    % NH
S(:,6) = 2/params.mm_mcf*S(:,6);    % SH

%%% Run the box model with ode45
[T, F] = ode45(@(t,y) boxModel(t,y,St,S,params),params.Tspan,IC,params.odeOpts);

%%% Make the output structure
out.nh_ch4    = interp1(T,F(:,1),St,'spline');
out.sh_ch4    = interp1(T,F(:,4),St,'spline');
out.nh_ch4c13 = interp1(T,(F(:,2)./F(:,1)-1)*1000,St,'spline');
out.sh_ch4c13 = interp1(T,(F(:,5)./F(:,4)-1)*1000,St,'spline');
out.nh_mcf    = interp1(T,F(:,3),St,'spline');
out.sh_mcf    = interp1(T,F(:,6),St,'spline');

end


%%% =======================================================================
%%% = END
%%% =======================================================================
