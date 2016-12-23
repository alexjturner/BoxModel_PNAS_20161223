%%% =======================================================================
%%% = define_prior.m
%%% = Alex Turner
%%% = 03/30/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Defines the prior distribution.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): ems         -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC          -- Vector with the ICs.
%%% =  ( 3): input_param -- A structure containing inputs to the inversion.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): p_prior -- Probability of the prior distribution.
%%% =======================================================================

function [ p_prior ] = define_prior( ems, IC, input_param )

%%% Read the structure
St         = input_param.St;
ch4_ems    = input_param.ch4_ems;
ch4c13_ems = input_param.ch4c13_ems;
mcf_ems    = input_param.mcf_ems;
oh_scale   = input_param.oh_scale;
nT         = input_param.nT;
use_log    = input_param.use_log;


%%% Define the prior distributions for each source
% Temporal correlation
tau        = 1*365.25;  % 1-year temporal correlation
tau_ch4    = 5*tau;     % Temporal correlation for methane emissions
tau_ch4c13 = 5*tau;     % Temporal correlation for ch4c13 composition
tau_mcf    = 3*tau;     % Temporal correlation for mcf emissions
tau_oh     = 3*tau;     % Temporal correlation for the OH anomaly
% NH CH4 (bounded normal)
mu       = ch4_ems.nh;
sig      = 20.^2*eye(nT);
sig      = fillDiagonals(sig,tau_ch4,St);
xL       = 300*ones(nT,1);
xU       = 500*ones(nT,1);
p_nh_ch4 = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
% SH CH4 (uniform)
mu       = ch4_ems.sh;
sig      = 20.^2*eye(nT);
sig      = fillDiagonals(sig,tau_ch4,St);
xL       = 100*ones(nT,1);
xU       = 250*ones(nT,1);
p_sh_ch4 = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
% NH CH4C13 (uniform)
xL          = -60*ones(nT,1);
xU          = -45*ones(nT,1);
p_nh_ch4c13 = @(x) p_uniform(x,xL,xU,use_log);
% SH CH4C13 (bounded normal)
xL          = -60*ones(nT,1);
xU          = -45*ones(nT,1);
p_sh_ch4c13 = @(x) p_uniform(x,xL,xU,use_log);
% NH MCF (bounded normal)
mu       = mcf_ems.prinn;
sig      = max([1.5*ones(size(mu)),.2*mu],[],2).^2; sig = diag(sig);
sig      = fillDiagonals(sig,tau_mcf,St);
xL       = -0.5*ones(nT,1);
xU       = max([1.0*ones(size(mu)),2.0*mu],[],2);
p_nh_mcf = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
% SH MCF (uniform)
xL       = -1.0*ones(nT,1);
xU       =  1.0*ones(nT,1);
p_sh_mcf = @(x) p_uniform(x,xL,xU,use_log);
% NH OH rescaling (bounded normal)
mu      = oh_scale.nh;
sig     = 0.10^2*eye(nT);
sig     = fillDiagonals(sig,tau_oh,St);
xL      = 0.80*ones(nT,1);
xU      = 1.20*ones(nT,1);
p_nh_oh = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
% SH OH rescaling (bounded normal)
mu      = oh_scale.sh;
sig     = 0.10^2*eye(nT);
sig     = fillDiagonals(sig,tau_oh,St);
xL      = 0.80*ones(nT,1);
xU      = 1.20*ones(nT,1);
p_sh_oh = @(x) p_normalB(x,mu,sig,xL,xU,use_log);
% Initial conditions (uniform: nh_12ch4, nh_13ch4, nh_mcf, sh_12ch4, sh_13ch4, sh_mcf)
nh_ch4_B    = [ 1540, 1620]; sh_ch4_B    = [ 1480, 1580];
nh_ch4c13_B = [-48.2,-46.6]; sh_ch4c13_B = [-48.2,-46.6];
nh_mcf_B    = [   15,  135]; sh_mcf_B    = [   15,  135];
xL          = [nh_ch4_B(1), nh_ch4_B(1)*(1+nh_ch4c13_B(1)/1000), nh_mcf_B(1),...
               sh_ch4_B(1), sh_ch4_B(1)*(1+sh_ch4c13_B(1)/1000), sh_mcf_B(1)];
xU          = [nh_ch4_B(2), nh_ch4_B(2)*(1+nh_ch4c13_B(2)/1000), nh_mcf_B(2),...
               sh_ch4_B(2), sh_ch4_B(2)*(1+sh_ch4c13_B(2)/1000), sh_mcf_B(2)];
p_IC        = @(x) p_uniform(x,xL,xU,use_log);
%%% Construct the full prior distribution
priors = [p_nh_ch4(ems(:,1)),    p_sh_ch4(ems(:,4)),    ...
          p_nh_ch4c13(ems(:,2)), p_sh_ch4c13(ems(:,5)), ...
          p_nh_mcf(ems(:,3)),    p_sh_mcf(ems(:,6)),    ...
          p_nh_oh(ems(:,7)),     p_sh_oh(ems(:,8)),     ...
          p_IC(IC)];
% Diagnostic
if any(isnan(priors))
    fprintf('NH ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f | SH ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f | NH oh = %4.0f, SH oh = %4.0f, IC = %4.0f \n',...
            priors(1),priors(3),priors(5),priors(2),priors(4),priors(6),priors(7),priors(8),priors(9));
end
if use_log
    p_prior = sum(priors);
else
    p_prior = prod(priors);
end

end

%%% Fill diagonals
function [ sig ] = fillDiagonals(sig,tau,St)

%%% Populate the diagonals of the covariance matrices (temporal correlations)
n = length(sig);
for i = 1:n
for j = 1:n
    if i ~= j
      sig(i,j) = sqrt(sig(i,i))*sqrt(sig(j,j)) * exp(-abs(St(i) - St(j))/tau);
    end
end
end

end

%%% Distributioins
function [ p ] = p_normalB(x,mu,sig,xL,xU,use_log)
if any(x < xL) || any(xU < x)
    p = NaN;
else
	if use_log
        p = logmvnpdf(x,mu,sig);
    else
        p = mvnpdf(x,mu,sig);
    end
end
end
function [ p ] = p_normal(x,mu,sig,use_log)
if use_log
    p = logmvnpdf(x,mu,sig);
else
    p = mvnpdf(x,mu,sig);
end
end
function [ p ] = p_uniform(x,xL,xU,use_log)
if any(x < xL) || any(xU < x)
    p = NaN;
else
    p = 1;
	if use_log
        p = log(p);
    end
end
end


%%% =======================================================================
%%% = END
%%% =======================================================================
