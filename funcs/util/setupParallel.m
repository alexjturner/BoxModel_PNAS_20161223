%%% =======================================================================
%%% = setupParallel.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Sets up a parallel environment for the simulation.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): run_parallel -- Are we running in parallel?  (True/False).
%%% =  ( 2): nWorkers     -- Number of workers to use.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = setupParallel( run_parallel, nWorkers )

fprintf('\n *** SETTING UP PARALLEL SIMULATION ***\n')
if run_parallel
    loop = true;
    while (nWorkers > 1 && loop)
        try
            initializePool(gcp('nocreate'),nWorkers);
            loop = false;
        catch
            nWorkers = nWorkers - 1;
            fprintf('\n   * ERROR: Too many workers requested!  Reducing to %i workers\n',nWorkers);
        end
    end
    if nWorkers == 1
        delete(gcp('nocreate'));
    end
    fprintf('   * Using %i workers!\n',nWorkers);
else
    delete(gcp('nocreate'));
    fprintf('   * Running in serial!\n');
end

end

function [ ] = initializePool( p, nWorkers )

if length(p) < 1
    parpool(nWorkers);
elseif (p.NumWorkers ~= nWorkers)
    delete(p);
    parpool(nWorkers);
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================
