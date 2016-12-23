%%% =======================================================================
%%% = cmaes_fun_eval.m
%%% = Alex Turner
%%% = 04/22/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Function that can be evaluated by CMA-ES.
%%% =  ( 2): Allows for multiple box model evaluations so we can run CMA-ES
%%% =        in parallel.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): stateVector -- State vector(s) for the inversion.
%%% =  ( 2): fun_param   -- Structure with inputs for the box model.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): p_post -- Posterior probability (negative so we can minimize).
%%% =======================================================================

function [ p_post ] = cmaes_fun_eval( stateVector, fun_param )

%%%% How many state vectors are we looking at?
nSV = size(stateVector,2);

%%% Evaluating in parallel?
if (nSV > 1)
    p_post = nan(1,nSV);
    if fun_param.run_parallel
        parfor k = 1:size(stateVector,2)
            % Break up the state vector
            [ems,IC] = disassembleStateVector(stateVector(:,k),fun_param.nT,fun_param.nE,fun_param.nI);
            % Evaluate the function
            if fun_param.use_log
                p_post(1,k) = fun_param.p_prior(ems,IC,fun_param) + fun_param.p_like(ems,IC,fun_param);
            else
                p_post(1,k) = fun_param.p_prior(ems,IC,fun_param) * fun_param.p_like(ems,IC,fun_param);
            end
        end
    else
        for k = 1:size(stateVector,2)
            % Break up the state vector
            [ems,IC] = disassembleStateVector(stateVector(:,k),fun_param.nT,fun_param.nE,fun_param.nI);
            % Evaluate the function
            if fun_param.use_log
                p_post(1,k) = fun_param.p_prior(ems,IC,fun_param) + fun_param.p_like(ems,IC,fun_param);
            else
                p_post(1,k) = fun_param.p_prior(ems,IC,fun_param) * fun_param.p_like(ems,IC,fun_param);
            end
        end
    end
else
    % Break up the state vector
    [ems,IC] = disassembleStateVector(stateVector,fun_param.nT,fun_param.nE,fun_param.nI);
    % Evaluate the function
    if fun_param.use_log
        p_post = fun_param.p_prior(ems,IC,fun_param) + fun_param.p_like(ems,IC,fun_param);
    else
        p_post = fun_param.p_prior(ems,IC,fun_param) * fun_param.p_like(ems,IC,fun_param);
    end
end

%%% Need to find a minimum with CMAES, make this negative
p_post = -1 * p_post;

end


%%% =======================================================================
%%% = END
%%% =======================================================================
