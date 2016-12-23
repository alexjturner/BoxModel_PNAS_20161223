%%% =======================================================================
%%% = getParameters.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Computes the parameters that are needed for the box model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St -- Our time vector.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): params -- Structure with the parameters for the box model.
%%% =======================================================================

function [ params ] = getParameters(St)

%%% Box model parameters and unit conversions
m       = 5.15e21;              % Rough: Total mass of atmosphere in g
m_air   = 28.8;                 % Average molar mass of the atmosphere in g/mol (mostly N2 and O2)
m_ch4   = 16;                   % Molar mass of CH4
m_mcf   = 133.4;                % Molar mass of methyl-chloroform
mConv   = m/m_air/1e12*1e-9;	% Factor to convert mixing ratios to Tg
mm_ch4  = m_ch4 * mConv;        % Convert CH4 mixing ratios to Tg
mm_mcf  = m_mcf * mConv;        % Convert MCF mixing ratios to Tg
KIE     = 1.005;                % Kinetic isotope effect (k12/k13) from Burkholder et al. (2015)
iKIE    = round(1/KIE,3);       % Inverse of KIE (for later)
OH      = 1.0e6;                % Average OH concentration in molec/cm3
k_12ch4 = 3.395e-15;            % reaction rate of OH with 12CH4 (cm3/molec/s)
k_13ch4 = k_12ch4 * iKIE;       % reaction rate of OH with 13CH4 (cm3/molec/s)
DaysToS = 60 * 60 * 24;         % Days to Seconds
YrToDay = 365.25;               % Years to Days
k_mcf_A = 5.66e-15;             % reaction rate with MCF (computed such that lifetime is about 5.5 years, Talukdar et al 1992, taken at 273K)
k_mcf_B = 6.05e-15;             % reaction rate determined by AJT
tau_NS  = 1.0;                  % Interhemispheric exchange rate (years)
% Convert so everything has units of days (working with julian days)
tau_NS  = tau_NS * YrToDay;
OHfac   = OH * DaysToS;
% Which methyl chloroform reaction rate are we using?
global k_mcf_flag
k_mcf = k_mcf_A;
if k_mcf_flag
    k_mcf = k_mcf_B;
end

%%% Guess for the initial conditions for the box model
nh_12ch4 = 1575;                        % ppb
sh_12ch4 = 1520;                        % ppb
nh_13ch4 = nh_12ch4 * (1 - 47.50/1000); % ppb
sh_13ch4 = sh_12ch4 * (1 - 47.25/1000); % ppb
nh_mcf   = 85;                          % ppt
sh_mcf   = 75;                          % ppt
% Assemble the ICs into a vector
IC = [nh_12ch4, nh_13ch4, nh_mcf, sh_12ch4, sh_13ch4, sh_mcf];

%%% ODE45 parameters
Tspan = St;
opts  = odeset('MaxStep',YrToDay/12);   % Make sure the max timestep is 1 month

%%% Make a structure with the parameters
% Unit conversions
params.mm_ch4  = mm_ch4;
params.mm_mcf  = mm_mcf;
params.YrToDay = YrToDay;
% Rate constants/lifetimes
params.k_12ch4 = OHfac * k_12ch4;   
params.k_13ch4 = OHfac * k_13ch4;
params.k_mcf   = OHfac * k_mcf;
params.tau_NS  = tau_NS;
% ODE parameters
params.Tspan   = Tspan;
params.IC      = IC;
params.odeOpts = opts;

end


%%% =======================================================================
%%% = END
%%% =======================================================================
