%%% =======================================================================
%%% = define_Jacobian.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Constructs the Jacobian for a deterministic inversion.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St           -- Our time vector.
%%% =  ( 2): ems          -- Matrix with the emission sources (and OH).
%%% =  ( 3): IC           -- Vector with the initial conditions.
%%% =  ( 4): params       -- Structure with parameters for the box model.
%%% =  ( 5): run_parallel -- Are we running in parallel?  (True/False).
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): jacobian_ems -- Jacobian for the emission sources (and OH).
%%% =  ( 2): jacobian_IC  -- Jacobian for the ICs.
%%% =======================================================================

function [ jacobian_ems, jacobian_IC ] = define_Jacobian( St, ems, IC, params, run_parallel )

%%% Get the unperturbed output
out_base = assembleObs(boxModel_wrapper(St,ems,IC,params));

%%% Initialize the jacobians
nT = length(St);    % Number of timesteps
nE = size(ems,2);   % Number of emission sources
nI = length(IC);    % Number of initial conditions
nY = length(out_base);
% Jacobians
jacobian_ems = zeros(nY,nT,nE);
jacobian_IC  = zeros(nY,nI);

%%% Perturbation for the Jacobian
delta.ems = [5,1,5,5,1,5,0.01,0.01];
delta.IC  = [1,1,1,1,1,1];

%%% Jacobian for ICs
for i = 1:nI
    IC_pert          = IC;
    IC_pert(i)       = IC(i) + delta.IC(i);
    out_plus         = assembleObs(boxModel_wrapper(St,ems,IC_pert,params));
    jacobian_IC(:,i) = (out_plus - out_base)./delta.IC(i);
end

%%% Jacobian for sources (compute it in parallel?)
if run_parallel % Parallel
    parfor i = 1:nT
        for j = 1:nE
            ems_pert            = ems;
            ems_pert(i,j)       = ems(i,j) + delta.ems(j);
            out_plus            = assembleObs(boxModel_wrapper(St,ems_pert,IC,params));
            jacobian_ems(:,i,j) = (out_plus - out_base)./delta.ems(j);
        end
    end
else % Serial
    for i = 1:nT
        for j = 1:nE
            ems_pert            = ems;
            ems_pert(i,j)       = ems(i,j) + delta.ems(j);
            out_plus            = assembleObs(boxModel_wrapper(St,ems_pert,IC,params));
            jacobian_ems(:,i,j) = (out_plus - out_base)./delta.ems(j);
        end
    end
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================
