%%% =======================================================================
%%% = DriverScript.m
%%% = Alex Turner
%%% = Created on 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES:
%%% = 
%%% = This is the driver script for the 2-box model methane inversion.
%%% = There are currently two different inversions implemented: (1) a
%%% = linear or non-linear deterministic inversion following Rodgers (2000)
%%% = and (2) an inversion using the non-linear Covariance Matrix Adaptation 
%%% = Evolution Strategy (CMA-ES).  Case (1) requires us to compute gradients 
%%% = and only allows for Gaussian errors.  Case (2) is  a stochastic method 
%%% = that automatically tunes the proposal distribution to improve the sampling,
%%% = however it does not provide error statistics that are consistent with
%%% = the distributions.  Case (2) also allows one to specify non-analytic 
%%% = distributions (e.g., bounded Gaussians or uniform distributions).
%%% =======================================================================


%%
%%% =======================================================================
%%% 1. Initialize
%%% =======================================================================

%%% Clear the MatLab space
clf
clear all
close all
clc

%%% Header
fprintf('\n ***********************************\n')
fprintf(' *** STARTING GLOBAL 2-BOX MODEL ***\n')
fprintf(' ***********************************\n')

%%% Define the directories
baseDir = pwd;
utilDir = sprintf('%s/funcs/', baseDir);
dataDir = sprintf('%s/data/',  baseDir);
outDir  = sprintf('%s/output/',baseDir);

%%% Add the utility functions
addpath(utilDir);
addpath(sprintf('%s/obs',               utilDir));
addpath(sprintf('%s/ems',               utilDir));
addpath(sprintf('%s/model',             utilDir));
addpath(sprintf('%s/util',              utilDir));
addpath(sprintf('%s/plot',              utilDir));
addpath(sprintf('%s/inv',               utilDir));
addpath(sprintf('%s/inv/deterministic', utilDir));
addpath(sprintf('%s/inv/stochastic',    utilDir));

%%% Define the time period
sYear = 1980;
eYear = 2016;
tRes  = 'year';     % Can be 'year' or 'month' (year preferred)
tAvg  = 'year';     % Smooth the observations
St    = getTime(sYear,eYear,tRes); % Time vector
nT    = length(St);

%%% Execute in parallel?
run_parallel = true;
nWorkers     = 4;
setupParallel(run_parallel,nWorkers);

%%% What kind of inversions do we want to do?
do_deterministic = true;    % Rodgers (2000)
do_cmaes         = true;    % Covariance Matrix Adaptation Evolution Strategy

%%% For reading the observations
% Do we want to reread the raw data?
reread.flag  = false;
% Other flags for re-reading
reread.sYear = sYear;   
reread.eYear = eYear;
reread.tRes  = tAvg;
reread.tAvg  = tAvg;
reread.dir   = dataDir;

%%% Other options and flags
% Use globals for some flags
global fixedCH4 fixedOH onlyCH4 onlyMCF schaefer        % Linear inversion
global k_mcf_flag smooth_MCF set_MCF_EMS MCF_EMS_val    % Methyl Chloroform
% Plotting flags
ftype           = 'eps';    % Type of plots to make? (eps, pdf, tif, or png)
plot_prior      = true;     % Plot the prior?
plot_raw        = false;    % Plot the raw observations?
plot_old_cmaes  = true;     % Plot an old CMA-ES solution (false means run a new one)
% Linear inversion flags
det_linear      = true;     % Use a linear deterministic inversion?
fixedCH4        = false;    % Use fixed methane emissions
fixedOH         = false;    % Use fixed OH anomalies
onlyCH4         = false;    % Only invert for methane emissions
onlyMCF         = false;    % Only invert for MCF emissions
schaefer        = false;    % Case that is most similar to Schaefer et al.
% MCF sensitivity test flags
k_mcf_flag      = false;    % Use k_MCF that AJT derived
smooth_MCF      = false;    % Smooth the MCF emissions with a 5-year filter?
set_MCF_EMS     = false;    % Set post-2000 emissions to a fixed value?
MCF_EMS_val     = 0.0;      % Fixed post-2000 MCF emissions value (Gg/yr)
reduce_MCFerr   = false;    % Reduce the errors in MCF observations?
MCF_ERR_val     = 2.0;      % Error in MCF observations (ppt)
% Flags for other tests to run
use_OH_stratMLO = false;    % Use the OH derived from MLO strat ozone?
use_Ed          = false;    % Use Ed Dlugokencky's hemispheric averages?

%%% Set the seed for repeatability
rng('default');


%%
%%% =======================================================================
%%% 2. Load the obs
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** LOADING THE OBSERVATIONS *** \n');

%%% Load the observations
% Structures with with three fields:
% - "obs":  Observations from each NOAA site (ppb)
% - "tim":  Julian date for the observation
% - "lat":  Latitude of the NOAA site
try % Add a try-catch statement in case the user hasn't downloaded the data
    ch4_obs     = getCH4(dataDir);      % CH4 observations (ppb)
    ch4c13_obs  = getCH4C13(dataDir);   % delta13C observations (permil)
    mcf_obs     = getMCF(dataDir);      % Methylchloroform observations (ppt)
    c2h6_obs    = getC2H6(dataDir);     % Ethane observations (ppt)
    o3strat_obs = getO3strat(dataDir);  % Stratospheric ozone observations (DU)
catch % Some data is missing
    try % See if ethane is the only problem
        ch4_obs     = getCH4(dataDir);
        ch4c13_obs  = getCH4C13(dataDir);
        mcf_obs     = getMCF(dataDir);
        o3strat_obs = getO3strat(dataDir);
        c2h6_obs    = NaN;
    catch % Otherwise, set the observation structures to NaN
        ch4_obs     = NaN;
        ch4c13_obs  = NaN;
        mcf_obs     = NaN;
        c2h6_obs    = NaN;
        o3strat_obs = NaN;
    end
end

%%% Make the observation structure
% Structure with 12 fields:
% - NH/SH CH4    obs & err (ppb)
% - NH/SH CH4C13 obs & err (permil)
% - NH/SH MCF    obs & err (ppt)
obs = makeObs(St,tAvg,ch4_obs,ch4c13_obs,mcf_obs,dataDir,reread);

%%% Use Ed Dlugokencky's obs? (sensitivity test)
if use_Ed
    ajt_obs = obs;
    ed_obs  = getEdObs(dataDir,ajt_obs,St,tAvg,reread);
    obs     = ed_obs;
    plotEdObs(St,ajt_obs,ed_obs,sprintf('%s/%s/raw_EdObs.%s',outDir,tRes,ftype))
end

%%% Reduce MCF errors?  (sensitivity test)
if reduce_MCFerr
    obs.nh_mcf_err = min([obs.nh_mcf_err,MCF_ERR_val*ones(size(obs.nh_mcf_err))],[],2);
    obs.sh_mcf_err = min([obs.sh_mcf_err,MCF_ERR_val*ones(size(obs.sh_mcf_err))],[],2);
end

%%% Diagnostics (check the raw data)
if plot_raw
    deseasonalize  = true;
    plot_all_sites = false;
    plotAllObs(St,obs,ch4_obs,   tAvg, 'ch4' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,ch4c13_obs,tAvg, 'd13C',sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,mcf_obs,   tAvg, 'mcf' ,sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
    plotAllObs(St,obs,c2h6_obs,  tAvg, 'c2h6',sprintf('%s/%s/raw_%%s_%%s.%s',outDir,tRes,ftype),deseasonalize,plot_all_sites);
end


%%
%%% =======================================================================
%%% 3. Load the emissions (all will be arrays with a length matching "St")
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** LOADING THE EMISSIONS *** \n');

%%% Get the CH4 emissions
% Stucture with two fields
% - "nh": CH4 emissions from the Northern hemisphere (Tg/yr)
% - "sh": CH4 emissions from the Southern hemisphere (Tg/yr)
ch4_ems = getCH4ems(St,tRes,dataDir);

%%% Get the delta13C composition for NH/SH CH4 emissions
% Stucture with two fields
% - "nh": delta13C composition from the Northern hemisphere (permil)
% - "sh": delta13C composition from the Southern hemisphere (permil)
ch4c13_ems = getCH4C13ems(St,tRes,dataDir);

%%% Get the MCF emissions (assumed to be in NH only)
% Stucture with two fields
% - "prinn":     MCF emissions from Prinn (Gg/yr)
% - "mcculloch": MCF emissions from McCulloch (Gg/yr)
mcf_ems = getMCFems(St,tRes,dataDir);


%%
%%% =======================================================================
%%% 4. Initialize the 2-box model
%%% =======================================================================

%%% Diagnostic
fprintf('\n *** RUN THE 2-BOX MODEL WITH PRIOR FLUXES *** \n');

%%% OH scaling factor
oh_scale.nh = ones(nT,1);
oh_scale.sh = ones(nT,1);
% Derive OH from the stratospheric ozone?
if use_OH_stratMLO
    OH_sensitivity  = 4.2;          % a 1% increase in strat O3 leads to a 4.2% decrease in OH (Murray et al., 2013)
    O3_site         = 'mlo_NOAA';   % which site to use?
    fDays           = 365.25*2;     % How long of a smoothing?
    tO3             = o3strat_obs.tim.(O3_site);
    yO3             = o3strat_obs.obs.(O3_site);
    yO3             = DeseasonalizeData(tO3,yO3,fDays);
    [tO3, yO3, ~]   = BlockAverage_AltError(tO3,yO3,ones(size(tO3)),365.25);
    oh_change       = yO3 / nanmean(yO3); % Convert strat O3 to OH change
    oh_change       = 1 ./ ((oh_change - 1) * OH_sensitivity + 1);
    yOH             = interp1(tO3,oh_change,St);
    yOH(isnan(yOH)) = 1;
    % Store this OH
    oh_scale.nh = yOH;
    oh_scale.sh = yOH;
end

%%% Construct a matrix of sources
% Column 1: NH CH4 emissions
% Column 2: NH CH4C13 composition
% Column 3: NH MCF emissions
% Column 4: SH CH4 emissions
% Column 5: SH CH4C13 composition
% Column 6: SH MCF emissions
% Column 7: NH OH scaling factor
% Column 8: SH OH scaling factor
ems = [ch4_ems.nh, ch4c13_ems.nh, mcf_ems.prinn, ch4_ems.sh, ch4c13_ems.sh, zeros(nT,1), oh_scale.nh, oh_scale.sh];

%%% Run the box model
params = getParameters(St); % Only need to do this once
IC     = params.IC;         % Guess for the inital conditions
out    = boxModel_wrapper(St,ems,IC,params);
if plot_prior
    writeData(St,obs,out,ems,IC,sprintf('%s/%s/prior_%%s.csv',outDir,tRes));
    plotObs(St,out,obs,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype));
    plotDrivers(St,ems,NaN*ems,sprintf('%s/%s/prior_%%s.%s',outDir,tRes,ftype),dataDir);
end


%%
%%% =======================================================================
%%% 5. Deterministic inversion (Rodgers, 2000)
%%% =======================================================================

if do_deterministic
        
    %%% Diagnostic
    fprintf('\n *** DETERMINISTIC INVERSION *** \n');

    %%% Invert
    [anal_soln,jacobian_ems,jacobian_IC,reltol,abstol] = invert_methane(St,obs,ems,IC,params,det_linear,run_parallel);

    %%% Plot the Jacobians
    plotJacobian(St,jacobian_ems,tRes,sprintf('%s/%s/jacobian_%%s.%s',outDir,tRes,ftype));

    %%% Try plotting the solution
    ems_anal = anal_soln{1};
    IC_anal  = anal_soln{2};
    out_anal = boxModel_wrapper(St,ems_anal,IC_anal,params);
    writeData(St,obs,out_anal,ems_anal,IC_anal,sprintf('%s/%s/anal_%%s.csv',outDir,tRes));
    plotObs(St,out_anal,obs,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype));
    plotDrivers(St,ems_anal,ems,sprintf('%s/%s/anal_%%s.%s',outDir,tRes,ftype),dataDir);
    
end


%%
%%% =======================================================================
%%% 6. Invert with CMAES
%%% =======================================================================

if do_cmaes
    
    %%% Diagnostics
    fprintf('\n *** STARTING CMA-ES INVERSION *** \n');

    %%% Use a log-likelihood and get matrix dimensions
    use_log = true;
    nE      = size(ems,2);
    nI      = length(IC);

    %%% Set the function parameters for the box model
    fun_param.run_parallel = run_parallel;
    fun_param.p_prior      = @(ems,IC,input_param) define_prior(ems,IC,input_param);
    fun_param.p_like       = @(ems,IC,input_param) define_likelihood(ems,IC,input_param);
    fun_param.use_log      = use_log;
    fun_param.params       = params;
    fun_param.IC           = params.IC;
    fun_param.St           = St;
    fun_param.ch4_ems      = ch4_ems;
    fun_param.ch4c13_ems   = ch4c13_ems;
    fun_param.mcf_ems      = mcf_ems;
    fun_param.oh_scale     = oh_scale;
    fun_param.nT           = nT;
    fun_param.nE           = nE;
    fun_param.nI           = nI;
    fun_param.obs          = obs;

    %%% Set the options for CMAES
    CMAES_opts                   = cmaes;
    CMAES_opts.EvalParallel      = 'yes';
    CMAES_opts.SaveFilename      = sprintf('output/%s/cmaes_dat/variablescmaes.mat',tRes);
    CMAES_opts.LogFilenamePrefix = sprintf('output/%s/cmaes_dat/outcmaes',tRes);
    CMAES_opts.Resume            = 'no'; % Default is 'no'
    CMAES_opts.CMA.active        = 2;
    CMAES_opts.DiagonalOnly      = 100;
    % Full
    CMAES_opts.StopFunEvals      = 5000000;
    CMAES_opts.MaxIter           = 20000;
    CMAES_opts.Restarts          = 10;
    % Medium
    CMAES_opts.StopFunEvals      = 1000000;
    CMAES_opts.MaxIter           = 5000;
    CMAES_opts.Restarts          = 6;
    % Small
    CMAES_opts.StopFunEvals      = 10000;
    CMAES_opts.MaxIter           = 200;
    CMAES_opts.Restarts          = 2;
    
    %%% Get the starting point and standard deviations
    % Starting point
    if do_deterministic % Use the deterministic inversion as a starting point
        xstart = ems_anal;
        xstart(xstart(:,1)<0,1) = 0+eps; % Enforce positive CH4 emissions
        xstart(xstart(:,3)<0,3) = 0+eps; % Enforce positive CH4 emissions
        xstart = assembleStateVector(xstart,IC_anal);
    else
        xstart = [ems(:,1),ems(:,2),ems(:,3),ems(:,4),ems(:,5),0.5+zeros(nT,1),ems(:,7),ems(:,8)];
        xstart = assembleStateVector(xstart,IC);
    end
    % Standard deviations
    xsigma.ems      = ones(nT,nE);
    sig.ch4         =  1.00;         % Tg/yr
    sig.ch4c13      =  0.05;         % permil
    sig.mcf_nh      =  0.05;         % Tg/yr
    sig.mcf_sh      =  0.05;         % Tg/yr
    sig.oh          =  0.005;        % OH scale factor
    xsigma.IC       = [0.05,0.05,0.01,0.05,0.05,0.01]; % Initial conditions
    xsigma.ems(:,1) = xsigma.ems(:,1) * sig.ch4;
    xsigma.ems(:,2) = xsigma.ems(:,2) * sig.ch4c13;
    xsigma.ems(:,3) = xsigma.ems(:,3) * sig.mcf_nh;
    xsigma.ems(:,4) = xsigma.ems(:,4) * sig.ch4;
    xsigma.ems(:,5) = xsigma.ems(:,5) * sig.ch4c13;
    xsigma.ems(:,6) = xsigma.ems(:,6) * sig.mcf_sh;
    xsigma.ems(:,7) = xsigma.ems(:,7) * sig.oh;
    xsigma.ems(:,8) = xsigma.ems(:,8) * sig.oh;
    xsigma          = assembleStateVector(xsigma.ems,xsigma.IC);
    
    %%% Use an old estimate?
    use_old = true;
    if use_old
        fprintf('   * PLOTTING THE OLD CMA-ES RESULTS\n')
        fname = sprintf('./%s',CMAES_opts.SaveFilename);
        if (exist(fname,'file') == 2)
            load(fname);
            xstart = bestever.x;
        end
    end

    % Run the CMA-ES inversion or just plot an old one?
    if plot_old_cmaes
        [cmaes_res.ems,cmaes_res.IC] = disassembleStateVector(bestever.x,nT,nE,nI);
    else
        %%% Invert with CMAES
        % INPUTS
        %  - fitfun:  name of objective/fitness function
        %  - xstart:  objective variables initial point, determines N
        %  - insigma: initial coordinate wise standard deviation(s)
        %  - inopts:  options struct
        % OUTPUTS
        %  - xmin:      minimum search point of last iteration
        %  - fmin:      function value of xmin
        %  - counteval: number of function evaluations done
        %  - stopflag:  stop criterion reached
        %  - out:       struct with various histories and solutions
        %  - bestever:  struct containing overall best solution (for convenience)
        [xmin,fmin,counteval,stopflag,out,bestever] = cmaes('cmaes_fun_eval', xstart, xsigma, CMAES_opts, fun_param);
        [cmaes_res.ems,cmaes_res.IC]                = disassembleStateVector(bestever.x,nT,nE,nI);
    end

    %%% Plot the best one
    ems_best = cmaes_res.ems;
    IC_best  = cmaes_res.IC;
    out_best = boxModel_wrapper(St,ems_best,IC_best,params);
    writeData(St,obs,out_best,ems_best,IC_best,sprintf('%s/%s/cmaes_%%s.csv',outDir,tRes));
    plotObs(St,out_best,obs,sprintf('%s/%s/cmaes_%%s.%s',outDir,tRes,ftype));
    plotDrivers(St,ems_best,ems,sprintf('%s/%s/cmaes_%%s.%s',outDir,tRes,ftype),dataDir);
    
end

%%
%%% Finished simulation
fprintf('\n ***********************************\n')
fprintf(' ***            DONE!            ***\n')
fprintf(' ***********************************\n\n')


%%% =======================================================================
%%% =                              E N D                                  =
%%% =======================================================================
