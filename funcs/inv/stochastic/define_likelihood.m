%%% =======================================================================
%%% = define_likelihood.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Defines the likelihood distribution.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): ems         -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC          -- Vector with the ICs.
%%% =  ( 3): input_param -- A structure containing inputs to the inversion.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): p_like -- Probability of the likelihood distribution.
%%% =======================================================================

function [ p_like ] = define_likelihood( ems, IC, input_param )

%%% Read the structure
obs     = input_param.obs;
St      = input_param.St;
params  = input_param.params;
use_log = input_param.use_log;

%%% Run the box model
out = boxModel_wrapper(St,ems,IC,params);

%%% Define different types of distributions to use
if use_log
    p_normal  = @(x,mu,sig) logmvnpdf(x,mu,sig);
else
    p_normal  = @(x,mu,sig) mvnpdf(x,mu,sig);
end

%%% Define the likelihood distributions for each source
% NH CH4 (normal)
y        = out.nh_ch4;
mu       = obs.nh_ch4;
sig      = obs.nh_ch4_err;
ind      = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_nh_ch4 = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% SH CH4 (normal)
y        = out.sh_ch4;
mu       = obs.sh_ch4;
sig      = obs.sh_ch4_err;
ind      = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_sh_ch4 = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% NH CH4C13 (normal)
y           = out.nh_ch4c13;
mu          = obs.nh_ch4c13;
sig         = obs.nh_ch4c13_err;
ind         = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_nh_ch4c13 = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% SH CH4C13 (normal)
y           = out.sh_ch4c13;
mu          = obs.sh_ch4c13;
sig         = obs.sh_ch4c13_err;
ind         = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_sh_ch4c13 = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% NH MCF (normal)
y        = out.nh_mcf;
mu       = obs.nh_mcf;
sig      = obs.nh_mcf_err;
ind      = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_nh_mcf = p_normal(y(ind),mu(ind),diag(sig(ind).^2));
% SH MCF (normal)
y        = out.sh_mcf;
mu       = obs.sh_mcf;
sig      = obs.sh_mcf_err;
ind      = ~isnan(mu) & ~isnan(y) & ~isnan(sig);
p_sh_mcf = p_normal(y(ind),mu(ind),diag(sig(ind).^2));

%%% Construct the full likelihood distribution
likeli = [p_nh_ch4,    p_sh_ch4,    ...
          p_nh_ch4c13, p_sh_ch4c13, ...
          p_nh_mcf,    p_sh_mcf];
% Diagnostic
if any(isnan(likeli))
    fprintf('NH ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f | ch4 = %4.0f, ch4c13 = %4.0f, mcf = %4.0f SH\n',...
             likeli(1),likeli(3),likeli(5),likeli(2),likeli(4),likeli(6));
end
if use_log
    p_like = sum(likeli);
else
    p_like = prod(likeli);
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================
