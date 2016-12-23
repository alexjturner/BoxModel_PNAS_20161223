function [xmin, ...      % minimum search point of last iteration
	  fmin, ...      % function value of xmin
	  counteval, ... % number of function evaluations done
	  stopflag, ...  % stop criterion reached
	  out, ...     % struct with various histories and solutions
	  bestever ... % struct containing overall best solution (for convenience)
	 ] = cmaes( ...
    fitfun, ...    % name of objective/fitness function
    xstart, ...    % objective variables initial point, determines N
    insigma, ...   % initial coordinate wise standard deviation(s)
    inopts, ...    % options struct, see defopts below
    varargin )     % arguments passed to objective function 
% cmaes.m, Version 2.55, last change: July, 2007 
% CMAES implements an Evolution Strategy with Covariance Matrix
% Adaptation (CMA-ES) for nonlinear function minimization.  For
% introductory comments and copyright see end of file (type 'type
% cmaes'). It should run with MATLAB (Windows, Linux) and Octave
% (Linux, package octave-forge is needed).
%
% OPTS = CMAES returns default options. 
% OPTS = CMAES('defaults') returns default options quietly.
% OPTS = CMAES('displayoptions') displays options. 
% OPTS = CMAES('defaults', OPTS) supplements options OPTS with default 
% options. 
%
% XMIN = CMAES(FUN, X0, SIGMA[, OPTS]) locates the minimum XMIN of
% function FUN starting from column vector X0 with the initial
% coordinate wise search standard deviation SIGMA.
%
% Input arguments: 
%
%  FUN is a string function name like 'myfun'. FUN takes as argument a
%     column vector of size of X0 and returns a scalar. An easy way to
%     implement a hard non-linear constraint is to return NaN. Then,
%     this function evaluation is not counted and a newly sampled
%     point is tried immediately.
%
%   X0 is a column vector, or a matrix, or a string. If X0 is a matrix,
%     mean(X0, 2) is taken as initial point. If X0 is a string like
%     '2*rand(10,1)-1', the string is evaluated first.
%
%   SIGMA is a scalar, or a column vector of size(X0,1), or a string
%     that can be evaluated into one of these. SIGMA determines the
%     initial coordinate wise standard deviations for the search.
%     Setting SIGMA one third of the initial search region is
%     appropriate, e.g., the initial point in [0, 6]^10 and SIGMA=2
%     means cmaes('myfun', 3*rand(10,1), 2).  If SIGMA is missing and
%     size(X0,2) > 1, SIGMA is set to sqrt(var(X0')'). That is, X0 is
%     used as a sample for estimating initial mean and variance of the
%     search distribution.
%
%   OPTS (an optional argument) is a struct holding additional input
%     options. Valid field names and a short documentation can be
%     discovered by looking at the default options (type 'cmaes'
%     without arguments, see above). Empty or missing fields in OPTS
%     invoke the default value, i.e. OPTS needs not to have all valid
%     field names.  Capitalization does not matter and unambiguous
%     abbreviations can be used for the field names. If a string is
%     given where a numerical value is needed, the string is evaluated
%     by eval, where 'N' expands to the problem dimension
%     (==size(X0,1)) and 'popsize' to the population size. 
%
% [XMIN, FMIN, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = ...
%    CMAES(FITFUN, X0, SIGMA)
% returns the best (minimal) point XMIN (found in the last
% generation); function value FMIN of XMIN; the number of needed
% function evaluations COUNTEVAL; a STOPFLAG value as cell array,
% where possible entries are 'fitness', 'tolx', 'tolupx', 'tolfun',
% 'maxfunevals', 'maxiter', 'stoptoresume', 'manual',
% 'warnconditioncov', 'warnnoeffectcoord', 'warnnoeffectaxis',
% 'warnequalfunvals', 'warnequalfunvalhist', 'bug' (use
% e.g. any(strcmp(STOPFLAG, 'tolx')) or findstr(strcat(STOPFLAG,
% 'tolx')) for further processing); a record struct OUT with various
% output, where the array HIST.MEAN.X contains the evolution of the
% mean value and the struct SOLUTIONS.BESTEVER contains the overall
% best evaluated point X with function value F evaluated at evaluation
% count EVALS. BESTEVER equals OUT.SOLUTIONS.BESTEVER. 
%
% A regular manual stop can be achieved via the file signals.par. The
% program is terminated if the first two non-white sequences in any
% line of this file are 'stop' and the value of the SaveFileName
% option (by default 'variablescmaes.mat'). Also a run can be skipped.
% Given, for example, 'skip SaveFileName run 2', skips the second run
% (if option Restarts is at least 2) and another run will be started.
% 
% To run the code completely quietly set Display, VerboseModulo, and
% Plotting options to 0.  When OPTS.Saving==1 (default) everything is
% saved in file OPTS.SaveFileName (default 'variablescmaes.mat')
% permitting to investigate the recent result (e.g. plot with function
% plotcmaes) even while CMAES is still running (which can be quite
% useful on expensive objective functions) and to resume the search
% afterwards by using the resume option.
%
% To find the best ever evaluated point load the variables typing
% "es=load('variablescmaes')" and investigate the variable
% es.out.solutions.bestever. To further control data sampling behavior
% use SaveModulo option.
%
% The primary strategy parameter to play with is OPTS.PopSize, which
% can be increased from its default value.  Increasing the population
% size (by default together with the parent number OPTS.ParentNumber)
% improves global search properties in exchange to speed. Speed
% decreases, as a rule, at most linearely with increasing population
% size. It is advisable to begin with the default small population
% size. The options Restarts and IncPopSize can be used for an
% automated multistart where the population size is increased by the
% factor IncPopSize (two by default) before each restart. X0 (given as
% string) is reevaluated for each restart. Stopping options
% StopFunEvals, StopIter, MaxFunEvals, and Fitness terminate the
% program, all others including MaxIter invoke another restart, where
% the iteration counter is reset to zero.
%
% Examples: 
%
%   XMIN = cmaes('myfun', 5*ones(10,1), 1.5); starts the search at
%   10D-point 5 and initially searches mainly between 5-3 and 5+3
%   (+- two standard deviations), but this is not a strict bound.
%   'myfun' is a name of a function that returns a scalar from a 10D
%   column vector.
%
%   opts.LBounds = 0; opts.UBounds = 10; 
%   X=cmaes('myfun', 10*rand(10,1), 5, opts);
%   search within lower bound of 0 and upper bound of 10. Bounds can
%   also be given as column vectors. If the optimum is not located
%   on the boundary, use rather a penalty approach to handle bounds. 
%
%   opts=cmaes; opts.StopFitness=1e-10;
%   X=cmaes('myfun', rand(5,1), 0.5, opts); stops the search, if
%   the function value is smaller than 1e-10.
%   
%   [X, F, E, STOP, OUT] = cmaes('myfun2', 'rand(5,1)', 1, [], P1, P2); 
%   passes two additional parameters to the function MYFUN2.
%
% See also FMINSEARCH, FMINUNC, FMINBND.

cmaVersion = '2.54'; 

% ----------- Set Defaults for Input Parameters and Options -------------
% These defaults may be edited for convenience

% Input Defaults (obsolete, these are obligatory now)
definput.fitfun = 'felli'; % frosen; fcigar; see end of file for more
definput.xstart = rand(10,1); % 0.50*ones(10,1);
definput.sigma = 0.3;

% Options defaults: Stopping criteria % (value of stop flag)
defopts.StopFitness  = '-Inf % stop if f(xmin) < stopfitness, minimization';
defopts.MaxFunEvals  = 'Inf  % maximal number of fevals';
defopts.MaxIter      = '1e3*(N+5)^2/sqrt(popsize) % maximal number of iterations';
defopts.StopFunEvals = 'Inf  % stop after resp. evaluation to resume later';
defopts.StopIter     = 'Inf  % stop after resp. iteration to resume later';
defopts.TolX         = '1e-11*max(insigma) % stop if x-change smaller TolX';
defopts.TolUpX       = '1e3*max(insigma) % stop if x-changes larger TolUpX';
defopts.TolFun       = '1e-12 % stop if fun-changes smaller TolFun';
defopts.TolHistFun   = '1e-13 % stop if back fun-changes smaller TolHistFun';
defopts.StopOnWarnings = 'yes  % ''no''==''off''==0, ''on''==''yes''==1 ';

% Options defaults: Other
defopts.DiffMaxChange = 'Inf  % maximal variable change(s), can be Nx1-vector';
defopts.DiffMinChange = '0    % minimal variable change(s), can be Nx1-vector';
defopts.WarnOnEqualFunctionValues = ...
    'yes  % ''no''==''off''==0, ''on''==''yes''==1 ';
defopts.LBounds = '-Inf % lower bounds, scalar or Nx1-vector'; 
defopts.UBounds = 'Inf  % upper bounds, scalar or Nx1-vector'; 
defopts.EvalParallel = 'no   % objective function FUN accepts NxM matrix, with M>1?';
defopts.EvalInitialX = 'yes  % evaluation of initial solution';
defopts.Restarts     = '0    % number of restarts ';
defopts.IncPopSize   = '2    % multiplier for population size before each restart';

defopts.PopSize      = '(4 + floor(3*log(N)))  % population size, lambda'; 
defopts.ParentNumber = 'floor(popsize/2)     % popsize equals lambda';
defopts.RecombinationWeights = 'superlinear decrease % or linear, or equal';
defopts.Display  = 'on   % display messages like initial and final message';
defopts.Plotting = 'on   % plot while running';
defopts.VerboseModulo = '100  % >=0, messaging after every i-th iteration';
defopts.Resume  = 'no   % resume former run from SaveFile';  
defopts.Science = 'off  % off==do some additional (minor) problem capturing';
defopts.Saving  =    'on   % [on|final|off][-v6] save data to file';
defopts.SaveModulo = '1    % if >1 record data less frequently after gen=100';
defopts.SaveTime   = '25   % max. percentage of time for recording data';
defopts.SaveFileName = 'variablescmaes.mat'; % file name for saving
defopts.Seed = 'sum(100*clock)  % evaluated if it is a string';

%qqqkkk 
%defopts.varopt1 = ''; % 'for temporary and hacking purposes'; 
%defopts.varopt2 = ''; % 'for temporary and hacking purposes'; 
defopts.UserData = 'for saving data/comments associated with the run';
defopts.UserDat2 = ''; 'for saving data/comments associated with the run';

% ---------------------- Handling Input Parameters ----------------------

if nargin < 1 || isequal(fitfun, 'defaults') % pass default options
  if nargin < 1
    disp('Default options returned (type "help cmaes" for help).');
  end
  xmin = defopts;
  if nargin > 1 % supplement second argument with default options
    xmin = getoptions(xstart, defopts);
  end
  return;
end

if isequal(fitfun, 'displayoptions')
 names = fieldnames(defopts); 
 for name = names'
   disp([name{:} repmat(' ', 1, 20-length(name{:})) ': ''' defopts.(name{:}) '''']); 
 end
 return; 
end

input.fitfun = fitfun; % record used input
if isempty(fitfun)
  % fitfun = definput.fitfun; 
  % warning(['Objective function not determined, ''' fitfun ''' used']);
  error(['Objective function not determined']);
end
if ~ischar(fitfun)
  error('first argument FUN must be a string');
end


if nargin < 2 
  xstart = [];
end

input.xstart = xstart;
if isempty(xstart)
  % xstart = definput.xstart;  % objective variables initial point
  % warning('Initial search point, and problem dimension, not determined');
  error('Initial search point, and problem dimension, not determined');
end

if nargin < 3 
  insigma = [];
end
if isa(insigma, 'struct')
  error(['Third argument SIGMA must be (or eval to) a scalar '...
	   'or a column vector of size(X0,1)']);
end
input.sigma = insigma;
if isempty(insigma)
  if size(myeval(xstart),2) > 1
    insigma = std(xstart, 0, 2); 
    if any(insigma == 0)
      error(['Initial search volume is zero, choose SIGMA or X0 appropriate']);
    end
  else
    % will be captured later
    % error(['Initial step sizes (SIGMA) not determined']);
  end
end

% Compose options opts
if nargin < 4 || isempty(inopts) % no input options available
  inopts = []; 
  opts = defopts;
else
  opts = getoptions(inopts, defopts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
counteval = 0; countevalNaN = 0; 
irun = 0;
while irun <= myeval(opts.Restarts) % for-loop does not work with resume
  irun = irun + 1; 

% ------------------------ Initialization -------------------------------

% Handle resuming of old run
flgresume = myevalbool(opts.Resume);
if ~flgresume % not resuming a former run
  % Assign settings from input parameters and options for myeval...
  xmean = mean(myeval(xstart), 2); % in case of xstart is a population 
  N = size(xmean, 1); numberofvariables = N; 
  popsize = floor(myeval(opts.PopSize) * myeval(opts.IncPopSize)^(irun-1)); 
  lambda = popsize;
  insigma = myeval(insigma);
  if all(size(insigma) == [N 2]) 
    insigma = 0.5 * (insigma(:,2) - insigma(:,1));
  end
else % flgresume is true, do resume former run
  tmp = whos('-file', opts.SaveFileName);
  for i = 1:length(tmp)
    if strcmp(tmp(i).name, 'localopts');
      error('Saved variables include variable "localopts", please remove');
    end
  end
  local.opts = opts; % keep stopping and display options
  local.varargin = varargin;
  load(opts.SaveFileName); 
  varargin = local.varargin;
  flgresume = 1;

  % Overwrite old stopping and display options
  opts.StopFitness = local.opts.StopFitness; 
  %%opts.MaxFunEvals = local.opts.MaxFunEvals;
  %%opts.MaxIter = local.opts.MaxIter; 
  opts.StopFunEvals = local.opts.StopFunEvals; 
  opts.StopIter = local.opts.StopIter;  
  opts.TolX = local.opts.TolX;
  opts.TolUpX = local.opts.TolUpX;
  opts.TolFun = local.opts.TolFun;
  opts.TolHistFun = local.opts.TolHistFun;
  opts.StopOnWarnings = local.opts.StopOnWarnings; 
  opts.Display = local.opts.Display;
  opts.Plotting = local.opts.Plotting;
  opts.VerboseModulo = local.opts.VerboseModulo;
  opts.Saving = local.opts.Saving;
  opts.SaveModulo = local.opts.SaveModulo;
  opts.SaveTime = local.opts.SaveTime;
  clear local; % otherwise local would be overwritten during load
end
  
% Evaluate options
stopFitness = myeval(opts.StopFitness); 
stopMaxFunEvals = myeval(opts.MaxFunEvals);  
stopMaxIter = myeval(opts.MaxIter);  
stopFunEvals = myeval(opts.StopFunEvals);  
stopIter = myeval(opts.StopIter);  
stopTolX = myeval(opts.TolX);
stopTolUpX = myeval(opts.TolUpX);
stopTolFun = myeval(opts.TolFun);
stopTolHistFun = myeval(opts.TolHistFun);
stopOnWarnings = myevalbool(opts.StopOnWarnings); 
flgWarnOnEqualFunctionValues = myevalbool(opts.WarnOnEqualFunctionValues);
flgEvalParallel = myevalbool(opts.EvalParallel);
flgdisplay = myevalbool(opts.Display);
flgplotting = myevalbool(opts.Plotting);
verbosemodulo = myeval(opts.VerboseModulo);
flgscience = myevalbool(opts.Science);
flgsaving = [];
strsaving = [];
if strfind(opts.Saving, '-v6') 
  i = strfind(opts.Saving, '%');
  if isempty(i) || i == 0 || strfind(opts.Saving, '-v6') < i
    strsaving = '-v6';
    flgsaving = 1;
    flgsavingfinal = 1;
  end
end
if strncmp('final', opts.Saving, 5)
  flgsaving = 0;
  flgsavingfinal = 1;
end
if isempty(flgsaving)
  flgsaving = myevalbool(opts.Saving);
  flgsavingfinal = flgsaving;
end
savemodulo = myeval(opts.SaveModulo);
savetime = myeval(opts.SaveTime);

if (isfinite(stopFunEvals) || isfinite(stopIter)) && ~flgsaving
  warning('To resume later the saving option needs to be set');
end

% Do more checking and initialization 
if flgresume % resume is on
  time.t0 = clock;
  if flgdisplay
    disp(['  resumed from ' opts.SaveFileName ]); 
  end
  if counteval >= stopMaxFunEvals 
    error(['MaxFunEvals exceeded, use StopFunEvals as stopping ' ...
	  'criterion before resume']);
  end
  if countiter >= stopMaxIter 
    error(['MaxIter exceeded, use StopIter as stopping criterion ' ...
	  'before resume']);
  end
  
else
  xmean = mean(myeval(xstart), 2); % evaluate xstart again, because of irun
  maxdx = myeval(opts.DiffMaxChange); % maximal sensible variable change
  mindx = myeval(opts.DiffMinChange); % minimal sensible variable change 
				      % can both also be defined as Nx1 vectors
  lbounds = myeval(opts.LBounds);		     
  ubounds = myeval(opts.UBounds);
  if length(lbounds) == 1
    lbounds = repmat(lbounds, N, 1);
  end
  if length(ubounds) == 1
    ubounds = repmat(ubounds, N, 1);
  end
  if isempty(insigma) % last chance to set insigma
    if all(lbounds > -Inf) && all(ubounds < Inf)
      if any(lbounds>=ubounds)
	error('upper bound must be greater than lower bound');
      end
      insigma = 0.3*(ubounds-lbounds);
      stopTolX = myeval(opts.TolX);  % reevaluate these
      stopTolUpX = myeval(opts.TolUpX);
    else
      error(['Initial step sizes (SIGMA) not determined']);
    end
  end

  % Check all vector sizes
  if size(xmean, 2) > 1 || size(xmean,1) ~= N
    error(['intial search point should be a column vector of size ' ...
	   num2str(N)]);
  elseif ~(all(size(insigma) == [1 1]) || all(size(insigma) == [N 1]))
    error(['input parameter SIGMA should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(stopTolX, 2) > 1 || ~ismember(size(stopTolX, 1), [1 N])
    error(['option TolX should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(stopTolUpX, 2) > 1 || ~ismember(size(stopTolUpX, 1), [1 N])
    error(['option TolUpX should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(maxdx, 2) > 1 || ~ismember(size(maxdx, 1), [1 N])
    error(['option DiffMaxChange should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(mindx, 2) > 1 || ~ismember(size(mindx, 1), [1 N])
    error(['option DiffMinChange should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(lbounds, 2) > 1 || ~ismember(size(lbounds, 1), [1 N])
    error(['option lbounds should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  elseif size(ubounds, 2) > 1 || ~ismember(size(ubounds, 1), [1 N])
    error(['option ubounds should be (or eval to) a scalar '...
	   'or a column vector of size ' num2str(N)] );
  end
  
  % Strategy internal parameter setting: Selection  
  mu = myeval(opts.ParentNumber); % number of parents/points for recombination
  if strncmp(lower(opts.RecombinationWeights), 'equal', 3)
    weights = ones(mu,1); % (mu_I,lambda)-CMA-ES
  elseif strncmp(lower(opts.RecombinationWeights), 'linear', 3)
    weights = mu+1-(1:mu)';
  elseif strncmp(lower(opts.RecombinationWeights), 'superlinear', 3)
    weights = log(mu+1)-log(1:mu)'; % muXone array for weighted recombination
  else
    error(['Recombination weights to be "' opts.RecombinationWeights ...
	   '" is not implemented']);
  end
  mueff=sum(weights)^2/sum(weights.^2); % variance-effective size of mu
  if mueff == lambda
    error(['Combination of values for PopSize, ParentNumber and ' ...
	  ' and RecombinationWeights is not reasonable']);
  end
  
  % Strategy internal parameter setting: Adaptation
  cc = 4/(N+4);         % time constant for cumulation for covariance matrix
  cs = (mueff+2)/(N+mueff+3); % t-const for cumulation for step size control
  mucov = mueff;   % size of mu used for calculating learning rate ccov
  ccov = (1/mucov) * 2/(N+1.41)^2 ... % learning rate for covariance matrix
	 + (1-1/mucov) * min(1,(2*mueff-1)/((N+2)^2+mueff)); 
  % ||ps|| is close to sqrt(mueff/N) for mueff large on linear fitness
  damps = ... % damping for step size control, usually close to one 
      (1 + 2*max(0,sqrt((mueff-1)/(N+1))-1)) ... % limit sigma increase
      * max(0.3, ... % reduce damps, if max. iteration number is small
	  1 - N/min(stopMaxIter,stopMaxFunEvals/lambda)) + cs; 

  %qqq hacking of a different parameter setting, e.g. for ccov or damps,
  % can be done here. 
  % ccov = 0.0*ccov; disp(['CAVE: ccov=' num2str(ccov)]);

  % Initialize dynamic internal state parameters
  if any(insigma <= 0) 
    error(['Initial search volume (SIGMA) must be greater than zero']);
  end
  if max(insigma)/min(insigma) > 1e6
    error(['Initial search volume (SIGMA) badly conditioned']);
  end
  sigma = max(insigma);              % overall standard deviation
  pc = zeros(N,1); ps = zeros(N,1);  % evolution paths for C and sigma
  if length(insigma) == 1
    insigma = insigma * ones(N,1) ;
  end
  B = eye(N,N);                      % B defines the coordinate system
  D = diag(insigma/max(insigma));    % diagonal matrix D defines the scaling
  BD = B*D;                          % for speed up only
  C = BD*(BD)';                      % covariance matrix
  fitness.hist=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values
  fitness.histsel=NaN*ones(1,10+ceil(3*10*N/lambda)); % history of fitness values

  % Initialize boundary handling
  bnd.isactive = any(lbounds > -Inf) || any(ubounds < Inf); 
  if bnd.isactive
    if any(lbounds>ubounds)
      error('lower bound found to be greater than upper bound');
    end
    [xmean ti] = xintobounds(xmean, lbounds, ubounds); % just in case
    if any(ti)
      warning('Initial point was out of bounds, corrected');
    end
    bnd.weights = zeros(N,1);         % weights for bound penalty
    % scaling is better in axis-parallel case, worse in rotated
    bnd.flgscale = 0; % scaling will be omitted if zero 
    if bnd.flgscale ~= 0 
      bnd.scale = diag(C)/mean(diag(C));
    else
      bnd.scale = ones(N,1);
    end
    
    idx = (lbounds > -Inf) | (ubounds < Inf);
    if length(idx) == 1
      idx = idx * ones(N,1);
    end
    bnd.isbounded = zeros(N,1);
    bnd.isbounded(find(idx)) = 1; 
    maxdx = min(maxdx, (ubounds - lbounds)/2);
    if any(sigma*sqrt(diag(C)) > maxdx)
      fac = min(maxdx ./ sqrt(diag(C)))/sigma;
      sigma = min(maxdx ./ sqrt(diag(C)));
      warning(['Initial SIGMA multiplied by the factor ' num2str(fac) ...
	       ', because it was larger than half' ...
	       ' of one of the boundary intervals']);
    end
    idx = (lbounds > -Inf) & (ubounds < Inf);
    dd = diag(C);
    if any(5*sigma*sqrt(dd(idx)) < ubounds(idx) - lbounds(idx))
      warning(['Initial SIGMA is, in at least one coordinate, ' ...
	       'much smaller than the '...
	       'given boundary intervals. For reasonable ' ...
	       'global search performance SIGMA should be ' ...
	       'between 0.2 and 0.5 of the bounded interval in ' ...
	       'each coordinate. If all coordinates have ' ... 
	       'lower and upper bounds SIGMA can be empty']);
    end
    bnd.dfithist = 1;              % delta fit for setting weights
    bnd.aridxpoints = [];          % remember complete outside points
    bnd.arfitness = [];            % and their fitness
    bnd.validfitval = 0;
    bnd.iniphase = 1;
  end

  % ooo initial feval, for output only
  if irun == 1 
    out.solutions.bestever.x = xmean;
    out.solutions.bestever.f = Inf;  % for simpler comparison below
    out.solutions.bestever.evals = counteval;
    bestever = out.solutions.bestever;
  end
  if myevalbool(opts.EvalInitialX)
    fitness.hist(1)=feval(fitfun, xmean, varargin{:}); 
    fitness.histsel(1)=fitness.hist(1);
    counteval = counteval + 1;
    if fitness.hist(1) < out.solutions.bestever.f 
	out.solutions.bestever.x = xmean;
	out.solutions.bestever.f = fitness.hist(1);
	out.solutions.bestever.evals = counteval;
	bestever = out.solutions.bestever;
    end
  else
    fitness.hist(1)=NaN; 
    fitness.histsel(1)=NaN; 
  end
    
  % initialize random number generator
  if ischar(opts.Seed)
    randn('state', eval(opts.Seed));     % random number generator state
  else
    randn('state', opts.Seed);
  end
  %qqq
%  load(opts.SaveFileName, 'startseed');
%  randn('state', startseed);
%  disp(['SEED RELOADED FROM ' opts.SaveFileName]);
  startseed = randn('state');         % for retrieving in saved variables

  % Initialize further constants
  chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of 
				      %   ||N(0,I)|| == norm(randn(N,1))
  weights = weights/sum(weights);     % normalize recombination weights array
  
  countiter = 0;
  % Initialize records and output
  if irun == 1
    time.t0 = clock;
    
    % per evaluation
    %  out.hist.AllSolutions.eval = 1; 
    %  out.hist.AllSolutions.f = NaN;
    %  out.hist.AllSolutions.x = xmean;     % most memory critical
    
    % TODO: keep also median solution? 
    out.evals = counteval;  % should be first entry
    out.stopflag = {};
    out.hist.evals = counteval;
    out.hist.iterations = countiter;
    out.hist.mean.x = xmean;
    out.hist.mean.f = fitness.hist(1);
    out.hist.mean.evals = counteval;
    out.hist.recentbest.x = xmean;
    out.hist.recentbest.f = fitness.hist(1); % reevaluations make an array
    out.hist.recentbest.evals = counteval;
    out.hist.recentworst.x = xmean;
    out.hist.recentworst.f = fitness.hist(1); % reevaluations make an array
    out.hist.recentworst.evals = counteval;
    
    % Single Parameters
    out.hist.param.evals = counteval;
    out.hist.param.iterations = countiter;
    out.hist.param.sigma = sigma;
    out.hist.param.maxstd = sigma * sqrt(max(diag(C)));
    out.hist.param.minstd = sigma * sqrt(min(diag(C)));
    [muell out.hist.param.maxstdidx] = max(diag(C));
    [muell out.hist.param.minstdidx] = min(diag(C));
    out.hist.param.maxD = max(diag(D));
    out.hist.param.minD = min(diag(D));
    out.hist.param.comment = ['maxD=sqrt(max(EV)) where EV are' ...
		    ' the eigenvalues of C and sigma^2*C is the' ...
		    ' covariance matrix of the search distribution'];
    
    % Parameter Arrays
    out.histParamArr.evals = counteval;
    out.histParamArr.iterations = countiter;
    out.histParamArr.sigma = sigma; % for convenience and completeness
    out.histParamArr.diagD = diag(D); 
    out.histParamArr.stds = sigma * sqrt(diag(C));
    out.histParamArr.Bmax = B(:,out.hist.param.maxstdidx); 
    out.histParamArr.Bmin = B(:,out.hist.param.minstdidx); 
    out.histParamArr.comment = ...
	['diagD = sort(sqrt(EV)), Bmax = eigenvector of largest ' ...
		    ' eigenvalue (EV) of C'];

    %  out.x = 0;                   
    %  out.y1=[fitness.hist(1) sigma max(diag(D))/min(diag(D)) ...
    %	  sigma*[max(diag(D)) min(diag(D))] fitness.hist(1)];
    %  out.y2=xmean'; out.y2a=xmean';
    %  out.y3=sigma*sqrt(diag(C))';
    %  out.y4=sort(diag(D))'; 
    outiter = 0;
  end
  
end % else flgresume

% Display initial message
if flgdisplay
  if mu == 1
    strw = '100';
  elseif mu < 8
    strw = [num2str(100*weights(1:end-1)','%.0f ') ... 
	    num2str(100*weights(end)','%.0f') ];
  else
    strw = [num2str(100*weights(1:2)','%.2g ') ...
	    num2str(100*weights(3)','%.2g') '...' ...
	    num2str(100*weights(end-1:end)',' %.2g') ']%, '];
    
  end
  if irun > 1
    strrun = [', run ' num2str(irun)];
  else
    strrun = '';
  end
  disp(['  n=' num2str(N) ': (' num2str(mu) ',' ...
	num2str(lambda) ')-CMA-ES(w=[' ...
	strw ']%, ' ...
	'mu_eff=' num2str(mueff,'%.1f') ...
	') on function ' ...
	(fitfun) strrun]);
end

% -------------------- Generation Loop --------------------------------
stopflag = {};
while isempty(stopflag)
  countiter = countiter + 1; 
  flush;
  % Generate and evaluate lambda offspring
 
  fitness.raw = repmat(NaN, 1, lambda);

  % parallel evaluation
  if flgEvalParallel
      arz = randn(N,lambda);
      arx = repmat(xmean, 1, lambda) + sigma * (BD * arz);                % Eq. (1)

      % You may handle constraints here. You may either resample
      % arz(:,k) and/or multiply it with a factor between -1 and 1
      % (the latter will decrease the overall step size) and
      % recalculate arx accordingly. Do not change arx or arz in any
      % other way.
 
      if ~bnd.isactive
        arxvalid = arx;
      else
        arxvalid = xintobounds(arx, lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx and arxvalid should not be changed.
      fitness.raw = feval(fitfun, arxvalid, varargin{:}); 
      countevalNaN = countevalNaN + sum(isnan(fitness.raw));
      counteval = counteval + sum(~isnan(fitness.raw)); 
  end

  % non-parallel evaluation and remaining NaN-values
  for k=find(isnan(fitness.raw)), 
    % fitness.raw(k) = NaN; 
    tries = 0;
    % Resample, until fitness is not NaN
    while isnan(fitness.raw(k))
      arz(:,k) = randn(N,1);
      arx(:,k) = xmean + sigma * (BD * arz(:,k));                % Eq. (1)

      % You may handle constraints here. You may either resample
      % arz(:,k) and/or multiply it with a factor between -1 and 1
      % (the latter will decrease the overall step size) and
      % recalculate arx accordingly. Do not change arx or arz in any
      % other way.
 
      if ~bnd.isactive
        arxvalid(:,k) = arx(:,k);
      else
        arxvalid(:,k) = xintobounds(arx(:,k), lbounds, ubounds);
      end
      % You may handle constraints here.  You may copy and alter
      % (columns of) arxvalid(:,k) only for the evaluation of the
      % fitness function. arx should not be changed.
      fitness.raw(k) = feval(fitfun, arxvalid(:,k), varargin{:}); 
      tries = tries + 1;
      if isnan(fitness.raw(k))
	countevalNaN = countevalNaN + 1;
      end
      if mod(tries, 100) == 0
	warning([num2str(tries) ...
                 ' NaN objective function values at evaluation ' ...
                 num2str(counteval)]);
      end
    end
    counteval = counteval + 1; % retries due to NaN are not counted
  end

  fitness.sel = fitness.raw; 

  % ----- handle boundaries -----
  if 1 < 3 && bnd.isactive
    % Get delta fitness values
    val = myprctile(fitness.raw, [25 75]);
    % more precise would be exp(mean(log(diag(C))))
    val = (val(2) - val(1)) / N / mean(diag(C)) / sigma^2;
    %val = (myprctile(fitness.raw, 75) - myprctile(fitness.raw, 25)) ...
    %    / N / mean(diag(C)) / sigma^2;
    % Catch non-sensible values 
    if ~isfinite(val)
      warning('Non-finite fitness range');
      val = max(bnd.dfithist);  
    elseif val == 0 % happens if all points are out of bounds
      val = min(bnd.dfithist(bnd.dfithist>0)); 
    elseif bnd.validfitval == 0 % first sensible val
      bnd.dfithist = [];
      bnd.validfitval = 1;
    end

    % Store delta fitness values
    if length(bnd.dfithist) < 20+(3*N)/lambda
      bnd.dfithist = [bnd.dfithist val];
    else
      bnd.dfithist = [bnd.dfithist(2:end) val];
    end

    [tx ti]  = xintobounds(xmean, lbounds, ubounds);

    % Set initial weights
    if bnd.iniphase 
      if any(ti) 
        bnd.weights(find(bnd.isbounded)) = 2.0002 * median(bnd.dfithist);
	if bnd.flgscale == 0 % scale only initial weights then
	  dd = diag(C); 
	  idx = find(bnd.isbounded); 
	  dd = dd(idx) / mean(dd); %  remove mean scaling
	  bnd.weights(idx) = bnd.weights(idx) ./ dd; 
	end
	if bnd.validfitval && countiter > 2
          bnd.iniphase = 0;
	end
      end
    end

    % Increase weights
    if  1 < 3 && any(ti) % any coordinate of xmean out of bounds
      % judge distance of xmean to boundary
      tx = xmean - tx;
      idx = (ti ~= 0 & abs(tx) > 3*max(1,sqrt(N)/mueff) ... 
	     * sigma*sqrt(diag(C))) ;
      % only increase if xmean is moving away
      idx = idx & (sign(tx) == sign(xmean - xold));
      if ~isempty(idx) % increase
        % the factor became 1.2 instead of 1.1, because
	bnd.weights(idx) = 1.2^(max(1, mueff/10/N)) * bnd.weights(idx); 
      end
    end

    % Calculate scaling biased to unity, product is one
    if bnd.flgscale ~= 0 
      bnd.scale = exp(0.9*(log(diag(C))-mean(log(diag(C))))); 
    end

    % Assigned penalized fitness
    bnd.arpenalty = (bnd.weights ./ bnd.scale)' * (arxvalid - arx).^2; 

    fitness.sel = fitness.raw + bnd.arpenalty;

  end % handle boundaries
  % ----- end handle boundaries -----
  
  % Sort by fitness 
  [fitness.raw, fitness.idx] = sort(fitness.raw);  
  [fitness.sel, fitness.idxsel] = sort(fitness.sel);   % minimization
  fitness.hist(2:end) = fitness.hist(1:end-1);    % record short history of
  fitness.hist(1) = fitness.raw(1);               % best fitness values
  fitness.histsel(2:end) = fitness.histsel(1:end-1);    % record short history of
  fitness.histsel(1) = fitness.sel(1);               % best fitness values

  % Calculate new xmean, this is selection and recombination 
  xold = xmean; % for speed up of Eq. (2) and (3)
  xmean = arx(:,fitness.idxsel(1:mu))*weights; 
  zmean = arz(:,fitness.idxsel(1:mu))*weights;%==D^-1*B'*(xmean-xold)/sigma
  if mu == 1
    fmean = fitness.sel(1);
  else
    fmean = NaN; % [] does not work in the latter assignment
    % fmean = feval(fitfun, xintobounds(xmean, lbounds, ubounds), varargin{:});
    % counteval = counteval + 1;
  end
  
  % Cumulation: update evolution paths
  ps = (1-cs)*ps + (sqrt(cs*(2-cs)*mueff)) * (B*zmean);          % Eq. (4)
  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.4 + 2/(N+1);
%  hsig = norm(ps)/sqrt(1-(1-cs)^(2*countiter))/chiN < 1.5 + 1/(N-0.5);
%  hsig = norm(ps) < 1.5 * sqrt(N);
%  hsig = 1;
  pc = (1-cc)*pc ...
        + hsig*(sqrt(cc*(2-cc)*mueff)/sigma) * (xmean-xold);     % Eq. (2)
  if hsig == 0
    %disp([num2str(countiter) ' ' num2str(counteval) ' pc update stalled']);
  end

  % Adapt covariance matrix
  if ccov > 0                                                    % Eq. (3)
    C = (1-ccov+(1-hsig)*ccov*cc*(2-cc)/mucov) * C ... % regard old matrix 
        + ccov * (1/mucov) * pc*pc' ...                % plus rank one update
        + ccov * (1-1/mucov) ...                       % plus rank mu update
          * sigma^-2 * (arx(:,fitness.idxsel(1:mu))-repmat(xold,1,mu)) ...
          * diag(weights) * (arx(:,fitness.idxsel(1:mu))-repmat(xold,1,mu))';
  end
  
  if 1 < 2 && ~flgscience 
    % remove momentum in ps, if ps is large and fitness is getting worse.
    % this should rarely happen. 
    % this is questionable in dynamic environments
    if sum(ps.^2)/N > 1.5 + 10*(2/N)^.5 && ...
        fitness.histsel(1) > max(fitness.histsel(2:3))
      ps = ps * sqrt(N*(1+max(0,log(sum(ps.^2)/N))) / sum(ps.^2));
      if flgdisplay
        disp(['Momentum in ps removed at [niter neval]=' ...
              num2str([countiter counteval]) ']']);
      end
    end
  end

  % Adapt sigma
  sigma = sigma * exp((norm(ps)/chiN - 1)*cs/damps);             % Eq. (5)

  % Update B and D from C
  if ccov > 0 && mod(countiter, 1/ccov/N/10) < 1
    C=triu(C)+triu(C,1)'; % enforce symmetry
    [B,D] = eig(C);       % eigen decomposition, B==normalized eigenvectors

    if any(~isfinite(diag(D)))
      clear idx; % prevents error under octave
      save(['tmp' opts.SaveFileName]);
      error(['function eig returned non-finited eigenvalues, cond(C)=' ...
	     num2str(cond(C)) ]);
    end
    if any(any(~isfinite(diag(B))))
      clear idx; % prevents error under octave
      save(['tmp' opts.SaveFileName]);
      error(['function eig returned non-finited eigenvectors, cond(C)=' ...
	     num2str(cond(C)) ]);
    end

    % limit condition of C to 1e14 + 1
    if min(diag(D)) <= 0
	if stopOnWarnings
	  stopflag(end+1) = {'warnconditioncov'};
	else
	  warning(['Iteration ' num2str(countiter) ...
		   ': Eigenvalue (smaller) zero']);
	  D(D<0) = 0;
	  tmp = max(diag(D))/1e14;
	  C = C + tmp*eye(N,N); D = D + tmp*eye(N,N); 
	end
    end
    if max(diag(D)) > 1e14*min(diag(D)) 
	if stopOnWarnings
	  stopflag(end+1) = {'warnconditioncov'};
	else
	  warning(['Iteration ' num2str(countiter) ': condition of C ' ...
		   'at upper limit' ]);
	  tmp = max(diag(D))/1e14 - min(diag(D));
	  C = C + tmp*eye(N,N); D = D + tmp*eye(N,N); 
	end
    end

    % Align order of magnitude of scales of sigma and C for nicer output
    % needs to be carefully reviewed and tested yet
    if 11 < 2 && sigma > 1e10*sqrt(max(diag(D)))
      fac = sigma / sqrt(median(diag(D)));
      sigma = sigma/fac;
      pc = fac * pc;
      C = fac^2 * C;
      D = fac^2 * D;
    end

    D = diag(sqrt(diag(D))); % D contains standard deviations now
    % D = D / prod(diag(D))^(1/N);  C = C / prod(diag(D))^(2/N);
    BD = B*D; % for speed up only

  end % if mod

  % ----- numerical error management -----
  % Adjust maximal coordinate axis deviations
  if any(sigma*sqrt(diag(C)) > maxdx)
    sigma = min(maxdx ./ sqrt(diag(C)));
    %warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
    %         'deviation at upper limit of ' num2str(maxdx)]);
    % stopflag(end+1) = {'maxcoorddev'};
  end
  % Adjust minimal coordinate axis deviations
  if any(sigma*sqrt(diag(C)) < mindx)
    sigma = max(mindx ./ sqrt(diag(C))) * exp(0.05+cs/damps); 
    %warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
    %         'deviation at lower limit of ' num2str(mindx)]);
    % stopflag(end+1) = {'mincoorddev'};;
  end
  % Adjust too low coordinate axis deviations
  if any(xmean == xmean + 0.2*sigma*sqrt(diag(C))) 
    if stopOnWarnings
	stopflag(end+1) = {'warnnoeffectcoord'};
    else
      warning(['Iteration ' num2str(countiter) ': coordinate axis std ' ...
	       'deviation too low' ]);
	C = C + ccov * diag(diag(C) .* ...
			    (xmean == xmean + 0.2*sigma*sqrt(diag(C))));
	sigma = sigma * exp(0.05+cs/damps); 
    end
  end
  % Adjust step size in case of (numerical) precision problem 
  if all(xmean == xmean ...
	    + 0.1*sigma*BD(:,1+floor(mod(countiter,N))))
    i = 1+floor(mod(countiter,N));
    if stopOnWarnings
	stopflag(end+1) = {'warnnoeffectaxis'};
    else
      warning(['Iteration ' num2str(countiter) ...
	       ': main axis standard deviation ' ...
	       num2str(sigma*D(i,i)) ' has no effect' ]);
	sigma = sigma * exp(0.2+cs/damps); 
    end
  end
  % Adjust step size in case of equal function values (flat fitness)
  if fitness.sel(1) == fitness.sel(1+ceil(0.1+lambda/4))
    if flgWarnOnEqualFunctionValues && stopOnWarnings
	stopflag(end+1) = {'warnequalfunvals'};
    else
      if flgWarnOnEqualFunctionValues
	warning(['Iteration ' num2str(countiter) ...
		 ': equal function values f=' num2str(fitness.sel(1)) ...
		 ' at maximal main axis sigma ' ...
		 num2str(sigma*max(diag(D)))]);
      end
      sigma = sigma * exp(0.2+cs/damps); 
    end
  end
  % Adjust step size in case of equal function values
  if countiter > 2 && myrange([fitness.hist fitness.sel(1)]) == 0  
    if stopOnWarnings
	stopflag(end+1) = {'warnequalfunvalhist'};
    else
      warning(['Iteration ' num2str(countiter) ...
	       ': equal function values in history at maximal main ' ...
	       'axis sigma ' num2str(sigma*max(diag(D)))]);
	sigma = sigma * exp(0.2+cs/damps); 
    end
  end
    
  % ----- end numerical error management -----
  
  % Keep overall best solution
  out.evals = counteval;
  out.solutions.evals = counteval;
  out.solutions.mean.x = xmean;
  out.solutions.mean.f = fmean;
  out.solutions.mean.evals = counteval;
  out.solutions.recentbest.x = arxvalid(:, fitness.idx(1));
  out.solutions.recentbest.f = fitness.raw(1);
  out.solutions.recentbest.evals = counteval + fitness.idx(1) - lambda;
  out.solutions.recentworst.x = arxvalid(:, fitness.idx(end));
  out.solutions.recentworst.f = fitness.raw(end);
  out.solutions.recentworst.evals = counteval + fitness.idx(end) - lambda;
  if fitness.hist(1) < out.solutions.bestever.f
    out.solutions.bestever.x = arxvalid(:, fitness.idx(1));
    out.solutions.bestever.f = fitness.hist(1);
    out.solutions.bestever.evals = counteval + fitness.idx(1) - lambda;
    bestever = out.solutions.bestever;
  end

  % Set stop flag
  if fitness.raw(1) <= stopFitness, stopflag(end+1) = {'fitness'}; end
  if counteval >= stopMaxFunEvals, stopflag(end+1) = {'maxfunevals'}; end
  if countiter >= stopMaxIter, stopflag(end+1) = {'maxiter'}; end
  if all(sigma*(max(abs(pc), sqrt(diag(C)))) < stopTolX) 
    stopflag(end+1) = {'tolx'};
  end
  if any(sigma*sqrt(diag(C)) > stopTolUpX) 
    stopflag(end+1) = {'tolupx'};
  end
  if sigma*max(diag(D)) == 0  % should never happen
    stopflag(end+1) = {'bug'};
  end
  if countiter > 2 && myrange([fitness.sel fitness.hist]) <= stopTolFun 
    stopflag(end+1) = {'tolfun'};
  end
  if countiter >= length(fitness.hist) && myrange(fitness.hist) <= stopTolHistFun 
    stopflag(end+1) = {'tolhistfun'};
  end
  if counteval >= stopFunEvals || countiter >= stopIter
    stopflag(end+1) = {'stoptoresume'};
    if length(stopflag) == 1 && flgsaving == 0
      error('To resume later the saving option needs to be set');
    end
  end
  % read stopping message from file signals.par 
  fid = fopen('signals.par', 'rt');
  while fid > 0
    strline = fgets(fid, 300);
    if strline < 0 % fgets returns -1 at end of file
      break;
    end
    % 'stop filename' sets stopflag to manual
    str = sscanf(strline, ' %s %s', 2);
    if strcmp(str, ['stop' opts.SaveFileName])
      stopflag(end+1) = {'manual'};
      break;
    end
    % 'skip filename run 3' skips a run, but not the last
    str = sscanf(strline, ' %s %s %s', 3);
    if strcmp(str, ['skip' opts.SaveFileName 'run'])
      i = strfind(strline, 'run');
      if irun == sscanf(strline(i+3:end), ' %d ', 1) && irun <= myeval(opts.Restarts)
	stopflag(end+1) = {'skipped'};
      end	
    end
  end % while, break 
  if fid > 0
    fclose(fid);
    clear fid; % prevents strange error under octave
  end
  
  out.stopflag = stopflag;

  % ----- output generation -----
  if verbosemodulo > 0 && isfinite(verbosemodulo)
    if countiter == 1 || mod(countiter, 10*verbosemodulo) < 1 
      disp(['Iterat, #Fevals:   Function Value    (median,worst) ' ...
	    '|Axis Ratio|' ...
	    'idx:Min SD idx:Max SD']); 
    end
    if mod(countiter, verbosemodulo) < 1 ...
	  || (verbosemodulo > 0 && isfinite(verbosemodulo) && ...
	      (countiter < 3 || ~isempty(stopflag)))
      [minstd minstdidx] = min(sigma*sqrt(diag(C)));
      [maxstd maxstdidx] = max(sigma*sqrt(diag(C)));
      % format display nicely
      disp([repmat(' ',1,4-floor(log10(countiter))) ...
	    num2str(countiter) ' , ' ...
	    repmat(' ',1,5-floor(log10(counteval))) ...
	    num2str(counteval) ' : ' ...
            num2str(fitness.hist(1), '%.13e') ...
	    ' +(' num2str(median(fitness.raw)-fitness.hist(1), '%.0e ') ...
	    ',' num2str(max(fitness.raw)-fitness.hist(1), '%.0e ') ...
	    ') | ' ...
	    num2str(max(diag(D))/min(diag(D)), '%4.2e') ' | ' ...
	    repmat(' ',1,1-floor(log10(minstdidx))) num2str(minstdidx) ':' ...
	    num2str(minstd, ' %.1e') ' ' ...
	    repmat(' ',1,1-floor(log10(maxstdidx))) num2str(maxstdidx) ':' ...
	    num2str(maxstd, ' %.1e')]);
    end
  end

  % measure time for recording data
  if countiter < 3 
    time.c = 0.5;
    time.nonoutput = 0;
    time.recording = 0;
    time.saving  = 0.5; % first saving after 10 seconds
    time.plotting = 0;
  else
    time.c = min(1, time.nonoutput/3 + 1e-9); % set backward horizon
    time.c = max(1e-5, 1/countiter); % mean over all or 1e-5
  end
  % get average time per iteration
  time.t1 = clock;
  time.act = max(0,etime(time.t1, time.t0));
  time.nonoutput = (1-time.c) * time.nonoutput ...
      + time.c * time.act; 

  time.recording = (1-time.c) * time.recording;
  time.saving  = (1-time.c) * time.saving;
  time.plotting = (1-time.c) * time.plotting;
  
  % record output data, concerning time issues
  if countiter < 1e2 || ~isempty(stopflag) || ...
	countiter >= outiter + savemodulo
    outiter = countiter; 
      out.hist.evals(end+1) = counteval;
      out.hist.iterations(end+1) = countiter;
      out.hist.mean.x(:,end+1) = xmean;
      out.hist.mean.f(end+1) = fmean;
      out.hist.mean.evals(end+1) = counteval;
      out.hist.recentbest.x(:,end+1) = arx(:,fitness.idx(1));
      out.hist.recentbest.f(1,end+1) = fitness.raw(1); % reevaluations would make an array
      out.hist.recentbest.evals(end+1) = counteval + fitness.idx(1) - lambda;
      out.hist.recentworst.x(:,end+1) = arx(:,fitness.idx(end));
      out.hist.recentworst.f(1,end+1) = fitness.raw(end); 
      out.hist.recentworst.evals(end+1) = counteval + fitness.idx(end) - lambda;

      % Single Parameters
      out.hist.param.evals(end+1) = counteval;
      out.hist.param.iterations(end+1) = countiter;
      out.hist.param.sigma(end+1) = sigma;
      out.hist.param.maxstd(end+1) = sigma * sqrt(max(diag(C)));
      out.hist.param.minstd(end+1) = sigma * sqrt(min(diag(C)));
      [muell out.hist.param.maxstdidx(end+1)] = max(diag(C));
      [muell out.hist.param.minstdidx(end+1)] = min(diag(C));
      out.hist.param.maxD(end+1) = max(diag(D));
      out.hist.param.minD(end+1) = min(diag(D));
      out.hist.param.comment = '';
      
      % Parameter Arrays
      out.histParamArr.evals(end+1) = counteval;
      out.histParamArr.iterations(end+1) = countiter;
      out.histParamArr.sigma(end+1) = sigma; % for convenience and completeness
      out.histParamArr.diagD(:,end+1) = diag(D); 
      out.histParamArr.stds(:,end+1) = sigma * sqrt(diag(C));
      out.histParamArr.Bmax(:,end+1) = B(:,out.hist.param.maxstdidx(end)); 
      out.histParamArr.Bmin(:,end+1) = B(:,out.hist.param.minstdidx(end)); 
      out.histParamArr.comment = '';
      
    % get average time for recording data
    time.t2 = clock;
    time.recording = time.recording + time.c * max(0,etime(time.t2, time.t1)); 
    
    if flgplotting && countiter > 1
      if ~isempty(stopflag) || ...
	  time.plotting < 0.2 * time.nonoutput
	outplot(out); % outplot defined below
	if countiter > 3
	  time.plotting = time.plotting + time.c * max(0,etime(clock, time.t2)); 
	end
      end
    end
    if countiter > 100 && ...
	  time.recording > 1 && ...
	  time.recording > savetime * (time.nonoutput+time.recording) / 100 
      savemodulo = savemodulo + 1;
      % disp('++savemodulo'); %qqq
    end
  end % if output

  % save everything
  time.t3 = clock;
  if ~isempty(stopflag) || time.saving < 0.05 * time.nonoutput 
    xmin = arxvalid(:, fitness.idx(1));
    fmin = fitness.raw(1);
    if flgsaving && countiter > 2
      clear idx; % prevents error under octave
      % -v6 : non-compressed non-unicode for version 6 and earlier
      if ~isempty(strsaving) && ~isoctave
	save('-mat', strsaving, opts.SaveFileName); % for inspection and possible restart	
      else 
	save('-mat', opts.SaveFileName); % for inspection and possible restart
      end
      time.saving = time.saving + time.c * max(0,etime(clock, time.t3)); 
    end
  end
  time.t0 = clock;

  % ----- end output generation -----

end % while, end generation loop

% -------------------- Final Procedures -------------------------------

% Evaluate xmean and return best recent point in xmin
fmin = fitness.raw(1);
xmin = arxvalid(:, fitness.idx(1)); % Return best point of last generation.
if length(stopflag) > sum(strcmp(stopflag, 'stoptoresume')) % final stopping
  out.solutions.mean.f = ...
      feval(fitfun, xintobounds(xmean, lbounds, ubounds), varargin{:});
  counteval = counteval + 1;
  out.solutions.mean.evals = counteval;
  if out.solutions.mean.f < fitness.raw(1)
    fmin = out.solutions.mean.f;
    xmin = xintobounds(xmean, lbounds, ubounds); % Return xmean as best point
  end
  if out.solutions.mean.f < out.solutions.bestever.f
    out.solutions.bestever = out.solutions.mean; % Return xmean as bestever point
    out.solutions.bestever.x = xintobounds(xmean, lbounds, ubounds); 
    bestever = out.solutions.bestever;
  end
end

% Save everything and display final message
if flgsavingfinal
  clear idx; % prevents error under octave
  if ~isempty(strsaving) && ~isoctave
    save('-mat', strsaving, opts.SaveFileName); % for inspection and possible restart	
  else 
    save('-mat', opts.SaveFileName);    % for inspection and possible restart
  end
  message = [' (saved to ' opts.SaveFileName ')'];
else
  message = [];
end

if flgdisplay
  disp(['#Fevals:   f(returned x)   |    bestever.f     | stopflag' ...
        message]);
  if isoctave
    strstop = stopflag(:); 
  else
      strcat(stopflag(:), '.');
  end
  strstop = stopflag(:); %strcat(stopflag(:), '.');
  disp([repmat(' ',1,6-floor(log10(counteval))) ...
        num2str(counteval, '%6.0f') ': ' num2str(fmin, '%.11e') ' | ' ...
        num2str(out.solutions.bestever.f, '%.11e') ' | ' ...
	strstop{1:end}]);
  if exist('sfile', 'var') 
    disp(['Results saved in ' sfile]); 
  end
end

  out.arstopflags{irun} = stopflag;
  if any(strcmp(stopflag, 'fitness')) ...
	|| any(strcmp(stopflag, 'maxfunevals')) ...
	|| any(strcmp(stopflag, 'stoptoresume')) ...
	|| any(strcmp(stopflag, 'manual'))
    break; 
  end
end % while irun <= Restarts


% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function [x, idx] = xintobounds(x, lbounds, ubounds)
%
% x can be a column vector or a matrix consisting of column vectors
%
  if ~isempty(lbounds)
    if length(lbounds) == 1
      idx = x < lbounds;
      x(idx) = lbounds;
    else
      arbounds = repmat(lbounds, 1, size(x,2));
      idx = x < arbounds;
      x(idx) = arbounds(idx);
    end
  end
  if ~isempty(ubounds)
    if length(ubounds) == 1
      idx2 = x > ubounds;
      x(idx2) = ubounds;
    else
      arbounds = repmat(ubounds, 1, size(x,2));
      idx2 = x > arbounds;
      x(idx2) = arbounds(idx2);
    end
  end
  idx = idx2-idx;
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function opts=getoptions(inopts, defopts)
% OPTS = GETOPTIONS(INOPTS, DEFOPTS) handles an arbitrary number of
% optional arguments to a function. The given arguments are collected
% in the struct INOPTS.  GETOPTIONS matches INOPTS with a default
% options struct DEFOPTS and returns the merge OPTS.  Empty or missing
% fields in INOPTS invoke the default value.  Fieldnames in INOPTS can
% be abbreviated.
%
% The returned struct OPTS is first assigned to DEFOPTS. Then any
% field value in OPTS is replaced by the respective field value of
% INOPTS if (1) the field unambiguously (case-insensitive) matches
% with the fieldname in INOPTS (cut down to the length of the INOPTS
% fieldname) and (2) the field is not empty.
%
% Example:
%   In the source-code of the function that needs optional
%   arguments, the last argument is the struct of optional
%   arguments:
%
%   function results = myfunction(mandatory_arg, inopts)
%     % Define four default options
%     defopts.PopulationSize = 200;
%     defopts.ParentNumber = 50;
%     defopts.MaxIterations = 1e6;
%     defopts.MaxSigma = 1;
%  
%     % merge default options with input options
%     opts = getoptions(inopts, defopts);
%
%     % Thats it! From now on the values in opts can be used
%     for i = 1:opts.PopulationSize
%       % do whatever
%       if sigma > opts.MaxSigma
%         % do whatever
%       end
%     end
%   
%   For calling the function myfunction with default options:
%   myfunction(argument1, []);
%   For calling the function myfunction with modified options:
%   opt.pop = 100; % redefine PopulationSize option
%   opt.PAR = 10;  % redefine ParentNumber option
%   opt.maxiter = 2; % opt.max=2 is ambiguous and would result in an error 
%   myfunction(argument1, opt);

%
% 04/07/19: Entries can be structs itself leading to a recursive
%           call to getoptions. 
%

if nargin < 2 || isempty(defopts) % no default options available
  opts=inopts;
  return;
elseif isempty(inopts) % empty inopts invoke default options
  opts = defopts;
  return;
elseif ~isstruct(defopts) % handle a single option value
  if isempty(inopts) 
    opts = defopts;
  elseif ~isstruct(inopts)
    opts = inopts;
  else
    error('Input options are a struct, while default options are not');
  end
  return;
elseif ~isstruct(inopts) % no valid input options
  error('The options need to be a struct or empty');
end

  opts = defopts; % start from defopts 
  % if necessary overwrite opts fields by inopts values
  defnames = fieldnames(defopts);
  idxmatched = []; % indices of defopts that already matched
  for name = fieldnames(inopts)'
    name = name{1}; % name of i-th inopts-field
    if isoctave
      for i = 1:size(defnames, 1)
	idx(i) = strncmp(lower(defnames(i)), lower(name), length(name));
      end
    else
	idx = strncmpi(defnames, name, length(name));
    end
    if sum(idx) > 1
      error(['option "' name '" is not an unambigous abbreviation. ' ...
	     'Use opts=RMFIELD(opts, ''' name, ...
	     ''') to remove the field from the struct.']);
    end
    if sum(idx) == 1
      defname  = defnames{find(idx)}; 
      if ismember(find(idx), idxmatched)
	error(['input options match more than ones with "' ...
	       defname '". ' ...
	       'Use opts=RMFIELD(opts, ''' name, ...
	       ''') to remove the field from the struct.']);
      end
      idxmatched = [idxmatched find(idx)];
      val = getfield(inopts, name);
      % next line can replace previous line from MATLAB version 6.5.0 on and in octave
      % val = inopts.(name);
      if isstruct(val) % valid syntax only from version 6.5.0
	opts = setfield(opts, defname, ...
	    getoptions(val, getfield(defopts, defname))); 
      elseif isstruct(getfield(defopts, defname)) 
      % next three lines can replace previous three lines from MATLAB 
      % version 6.5.0 on
      %   opts.(defname) = ...
      %      getoptions(val, defopts.(defname)); 
      % elseif isstruct(defopts.(defname)) 
	warning(['option "' name '" disregarded (must be struct)']); 
      elseif ~isempty(val) % empty value: do nothing, i.e. stick to default
	opts = setfield(opts, defnames{find(idx)}, val);
	% next line can replace previous line from MATLAB version 6.5.0 on
	% opts.(defname) = inopts.(name); 
      end
    else
      warning(['option "' name '" disregarded (unknown field name)']);
    end
  end

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myeval(s)
  if ischar(s)
    res = evalin('caller', s);
  else
    res = s;
  end
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myevalbool(s)
  if ~ischar(s) % s may not and cannot be empty
    res = s;
  else % evaluation string s
    if strncmp(lower(s), 'yes', 3) || strncmp(lower(s), 'on', 2) ...
	  || strncmp(lower(s), 'true', 4) || strncmp(s, '1 ', 2)
      res = 1;
    elseif strncmp(lower(s), 'no', 2) || strncmp(lower(s), 'off', 3) ...
	  || strncmp(lower(s), 'false', 5) || strncmp(s, '0 ', 2)
      res = 0;
    else
      try res = evalin('caller', s); catch
	error(['String value "' s '" cannot be evaluated']);
      end
      try res ~= 0; catch
	error(['String value "' s '" cannot be evaluated reasonably']);
      end
    end
  end
  

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res = isoctave
% any hack to find out whether we are running octave
  s = version;
  res = 0;
  if exist('fflush', 'builtin') && eval(s(1)) < 7
    res = 1;
  end

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function flush
  if isoctave
    feval('fflush', stdout);
  end

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function outplot(o)
  persistent lasteval
  acteval = evalin('caller', 'counteval');
  foffset = 1e-99;

  if ~isfield(o, 'y1') % due to historical reason
    o.x = o.hist.param.evals;
    o.y1 = [o.hist.recentbest.f' ...
	    o.hist.param.sigma' ...
	    (o.hist.param.maxD ./ o.hist.param.minD)' ...
	    (o.hist.param.sigma.*o.hist.param.maxD)' ...
	    (o.hist.param.sigma.*o.hist.param.minD)' ...
	    o.hist.recentworst.f'];
    o.y2 = o.hist.mean.x'; 
    o.y2a = o.hist.recentbest.x';
    o.y3 = o.histParamArr.stds'; % cave: o.x comes from hist
    o.y4 = o.histParamArr.diagD'; 
  end
  minfitness = min(o.y1(:,1)); 
  dfit = o.y1(:,1)-minfitness; 
  dfit(find(dfit<1e-98)) = NaN;

%%%% OCTAVE %%%%
  if isoctave
    if exist('figure') % multiple windows require X11
      figure(327);
      text(); % in older Matlab versions: replace with text('')
      if isempty(lasteval) || acteval < lasteval 
	hold off;
      end
      semilogy(o.x, o.y4(:,1), ';principle axes lengths;'); hold on;
      semilogy(o.x, o.y4(:,2:end), ';;'); 

      figure(326);
      text(); % in older Matlab versions: replace with text('')
      if isempty(lasteval) || acteval < lasteval 
	hold off;
      end
      for i = 1:size(o.y3,2)
	text(o.x(end), o.y3(end, i), [' ' num2str(i)]);
      end

      semilogy(o.x, o.y3(:,1), ...
	       ';standard deviations of all variables;'); hold on;
      semilogy(o.x, o.y3(:,2:end), ';;'); 
      
      figure(325); 
      text(); % in older Matlab versions: replace with text('')
      if isempty(lasteval) || acteval < lasteval 
	hold off;
      end
      for i = 1:size(o.y2,2)
	text(o.x(end), o.y2(end, i), [' ' num2str(i)]);
      end
      plot(o.x, o.y2(:,1), ';object variables;'); hold on;
      plot(o.x, o.y2(:,2:end), ';;');
      
      figure(324); 
    end % exist figure

    text(); % in older Matlab versions: replace with text('')
    if isempty(lasteval) || acteval < lasteval 
      hold off;
    end
    
    semilogy(o.x, [o.y1(:,1)], ...
	     [';fitness=' num2str(o.y1(end, 1), '%.15g') ';']); hold on
    semilogy(o.x, [o.y1(:,2)], ';sigma;'); 
    semilogy(o.x, [dfit], [';fitness-[min(fitness)=' num2str(minfitness, '%.15g') '];']); 
    semilogy(o.x, o.y1(:,3), ';axis ratio;'); 
    semilogy(o.x, o.y1(:,4), 'k;max std in coordinate;'); 
    semilogy(o.x, o.y1(:,5), 'k;min std in coordinate;'); 
    semilogy(o.x, o.y1(:,6), ';worst fitness;'); 
    if size(o.y1,2) > 6
      semilogy(o.x, o.y1(:,7:end), ';more data;'); 
    end
    grid on;
    return;
  end

%%%% MATLAB %%%%
  if isempty(lasteval) || acteval < lasteval
    figure(324);
  elseif gcf ~= 324 % prevents repeated raise of figure
    if  ismember(324, findobj('Type', 'figure'))
      set(0, 'CurrentFigure', 324); % geht nur, wenn figure schon exisitiert
    else
      figure(324);
    end
  end
  lasteval = acteval;

  foffset = 1e-99;
  dfit = o.y1(:,1)-minfitness;
  dfit(find(dfit<1e-98)) = NaN;
  subplot(2,2,1); hold off; semilogy(o.x,dfit,'-c'); hold on;
  idx = find(o.y1(:,1)>1e-98);  % positive values
  subplot(2,2,1);semilogy(o.x(idx), o.y1(idx,1)+foffset, '.b'); hold on; 
  idx = find(o.y1(:,1) < -1e-98);  % negative values
  subplot(2,2,1);semilogy(o.x(idx), abs(o.y1(idx,1))+foffset,'.r');hold on; 
  subplot(2,2,1);semilogy(o.x,abs(o.y1(:,1))+foffset,'-b'); hold on;
  if size(o.y1, 2) > 3
    %qqq
    subplot(2,2,1);semilogy(o.x,abs(o.y1(:,4:5)),'-m'); hold on;
    subplot(2,2,1);semilogy(o.x,abs(o.y1(:,6:end))+foffset,'-k'); hold on;
    % subplot(2,2,1);semilogy(o.x,(o.y1(:,4:end)),'-k'); hold on;
    % subplot(2,2,1);semilogy(o.x,(o.y1(:,4:end)),'-m'); hold on;
    % subplot(2,2,1);semilogy(o.x,(o.y1(:,4:end))); hold on;
  end
  subplot(2,2,1);semilogy(o.x,(o.y1(:,3)),'-r'); hold on;
  subplot(2,2,1);semilogy(o.x,(o.y1(:,2)),'-g'); 
  ax = axis;
%  text(ax(1), 10^(log10(ax(3))+0.05*(log10(ax(4))-log10(ax(3)))), ...
%       [ ' f=' num2str(o.y1(end, 1), '%.15g')]);
  text(ax(1), 10^(log10(ax(3))+0.05*(log10(ax(4))-log10(ax(3)))), ...
       [ ' f=' num2str(o.y1(end, 1), ' %.15g') ...
	   ' ' num2str(minfitness, '(%.15g)')]);
  title('abs(f) (blue), f-min(f) (cyan), Sigma (green), Axis Ratio (red)');
  grid on; 

  subplot(2,2,2); hold off; plot(o.x, o.y2,'-'); 
  if size(o.y2, 2) < 100
    ax = axis;
    ax(2) = max(1.05*o.x(end), ax(2));
    axis(ax);
    yy = linspace(ax(3), ax(4), size(o.y2,2))';
    [yyl idx] = sort(o.y2(end,:));
    [muell idx2] = sort(idx);
    hold on;
    plot([o.x(end) ax(2)]', [o.y2(end,:)' yy(idx2)]', '-');
    plot(repmat(o.x(end),2), [ax(3) ax(4)], 'k-');
    for i = 1:length(idx)
      text(ax(2), yy(i), [' ' num2str(idx(i))]);
    end
  end
  title(['Object Variables (' num2str(size(o.y2, 2)) 'D)']);grid on;

  subplot(2,2,3); hold off; semilogy(o.x, o.y3, '-'); 
  if size(o.y2, 2) < 100
    ax = axis; 
    ax(2) = max(1.05*o.x(end), ax(2));
    axis(ax);
    yy = logspace(log10(ax(3)), log10(ax(4)), size(o.y3,2))';
    [yyl idx] = sort(o.y3(end,:));
    [muell idx2] = sort(idx);
    hold on;
    plot([o.x(end) ax(2)]', [o.y3(end,:)' yy(idx2)]', '-');
    plot(repmat(o.x(end),2), [ax(3) ax(4)], 'k-');
    for i = 1:length(idx)
      text(ax(2), yy(i), [' ' num2str(idx(i))]);
    end
  end
  title('Standard Deviations of All Variables'); grid on;
  xlabel('function evaluations'); 

  subplot(2,2,4); semilogy(o.x, o.y4, '-');
  title('Scaling (All Main Axes)'); grid on;
  xlabel('function evaluations'); 
  % zoom on; 
  drawnow;
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
% ----- replacements for statistic toolbox functions ------------
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res=myrange(x)
  res = max(x) - min(x);
  
% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res = myprctile(inar, perc, idx)
%
% Computes the percentiles in vector perc from vector inar
% returns vector with length(res)==length(perc)
% idx: optional index-array indicating sorted order
%

N = length(inar);
flgtranspose = 0;

% sizes 
if size(perc,1) > 1
  perc = perc';
  flgtranspose = 1;
  if size(perc,1) > 1
    error('perc must not be a matrix');
  end
end
if size(inar, 1) > 1 && size(inar,2) > 1
  error('data inar must not be a matrix');
end
 
% sort inar
if nargin < 3 || isempty(idx)
  [sar idx] = sort(inar);
else
  sar = inar(idx);
end

res = [];
for p = perc
  if p <= 100*(0.5/N)
    res(end+1) = sar(1);
  elseif p >= 100*((N-0.5)/N)
    res(end+1) = sar(N);
  else
    % find largest index smaller than required percentile
    availablepercentiles = 100*((1:N)-0.5)/N;
    i = max(find(p > availablepercentiles));
    % interpolate linearly
    res(end+1) = sar(i) ...
	+ (sar(i+1)-sar(i))*(p - availablepercentiles(i)) ...
	/ (availablepercentiles(i+1) - availablepercentiles(i));

  end
end

if flgtranspose
  res = res';
end


% ---------------------------------------------------------------  
% --------------- OBJECTIVE TEST FUNCTIONS ----------------------  
% ---------------------------------------------------------------  

%%% Unimodal functions

function f=fjens1(x)
%
% use population size about 2*N
%
  f = sum((x>0) .* x.^1, 1);
  if any(any(x<0))
    idx = sum(x < 0, 1) > 0;
    f(idx) = 1e3;
%    f = f + 1e3 * sum(x<0, 1);
%    f = f + 10 * sum((x<0) .* x.^2, 1);
    f(idx) = f(idx) + 1e-3*abs(randn(1,sum(idx)));
%    f(idx) = NaN;
  end

function f=fsphere(x)

  f=sum(x.^2, 1);

function f = fsphereoneax(x)
  f = x(1)^2;
  f = mean(x)^2;
  
function f=frandsphere(x)
  N = size(x,1);
  idx = ceil(N*rand(7,1));
  f=sum(x(idx).^2);

function f=fspherelb0(x, M) % lbound at zero for 1:M needed
  if nargin < 2 M = 0; end
  N = size(x,1);
  % M active bounds, f_i = 1 for x = 0
  f = -M + sum((x(1:M) + 1).^2);
  f = f + sum(x(M+1:N).^2);
  
function f=fspherehull(x)
  % Patton, Dexter, Goodman, Punch
  % in -500..500
  % spherical ridge through zeros(N,1)
  % worst case start point seems x = 2*100*sqrt(N)
  % and small step size
  N = size(x,1);
  f = norm(x) + (norm(x-100*sqrt(N)) - 100*N)^2;
  
function f=fellilb0(x, idxM, scal) % lbound at zero for 1:M needed
  N = size(x,1);
  if nargin < 3 || isempty(scal)
    scal = 100;
  end
  scale=scal.^((0:N-1)/(N-1));
  if nargin < 2 || isempty(idxM)
    idxM = 1:N;
  end
  %scale(N) = 1e0;
  % M active bounds
  xopt = 0.1;
  x(idxM) = x(idxM) + xopt;
  f = scale.^2*x.^2;
  f = f - sum((xopt*scale(idxM)).^2); 
%  f = exp(f) - 1;
%  f = log10(f+1e-19) + 19;

  f = f + 1e-19;
  
function f=fcornersphere(x)
  w = ones(size(x,1));
  w(1) = 2.5; w(2)=2.5;
  idx = x < 0;
  f = sum(x(idx).^2);
  idx = x > 0;
  f = f + 2^2*sum(w(idx).*x(idx).^2);
  
function f=fsectorsphere(x, scal)
%
% This is deceptive for cumulative sigma control in large dimension:
% The strategy (initially) diverges for N=50 and popsize = 150.  (Even
% for cs==1 this can be observed for larger settings of N and
% popsize.) The reason is obvious from the function topology. 
% Divergence can be avoided by setting boundaries or adding a
% penalty for large ||x||. Then, convergence can be observed again. 
% Conclusion: for popsize>N cumulative sigma control is not completely
% reasonable, but I do not know better alternatives.
%
  if nargin < 2 || isempty (scal)
    scal = 1e3;
  end
  f=sum(x.^2);
  idx = find(x<0);
  f = f + (scal-1)^2 * sum(x(idx).^2);
  
function f=fstepsphere(x, scal)
  if nargin < 2 || isempty (scal)
    scal = 1e0;
  end
  N = size(x,1);
  f=1e-11+sum(scal.^((0:N-1)/(N-1))*floor(x+0.5).^2);
  f=1e-11+sum(floor(scal.^((0:N-1)/(N-1))'.*x+0.5).^2);
%  f=1e-11+sum(floor(x+0.5).^2);

function f=fstep(x)
  % in -5.12..5.12 (bounded)
  N = size(x,1);
  f=1e-11+6*N+sum(floor(x));

function f=flnorm(x, scal, e)
if nargin < 2 || isempty(scal)
  scal = 1;
end
if nargin < 3 || isempty(e)
  e = 1;
end
if e==inf
  f = max(abs(x));
else
  N = size(x,1);
  scale = scal.^((0:N-1)/(N-1))';
  f=sum(abs(scale.*x).^e);
end

function f=fneumaier3(x) 
  % in -n^2..n^2
  % x^*-i = i(n+1-i)
  N = size(x,1);
%  f = N*(N+4)*(N-1)/6 + sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  f = sum((x-1).^2) - sum(x(1:N-1).*x(2:N));
  
function f=fchangingsphere(x)
  N = size(x,1);
  global scale_G; global count_G; if isempty(count_G) count_G=-1; end
  count_G = count_G+1;
  if mod(count_G,10) == 0
    scale_G = 10.^(2*rand(1,N));
  end
  %disp(scale(1));
  f = scale_G*x.^2;
  
function f= flogsphere(x)
 f = 1-exp(-sum(x.^2));
  
function f= fexpsphere(x)
 f = exp(sum(x.^2)) - 1;
  
function f=fbaluja(x)
  % in [-0.16 0.16]
  y = x(1);
  for i = 2:length(x)
    y(i) = x(i) + y(i-1);
  end
  f = 1e5 - 1/(1e-5 + sum(abs(y)));

function f=fschwefel(x)
  f = 0;
  for i = 1:size(x,1),
    f = f+sum(x(1:i))^2;
  end

function f=fcigar(x, ar)
  if nargin < 2 || isempty(ar)
    ar = 1e3;
  end
  f = x(1)^2 + ar^2*sum(x(2:end).^2);
  
function f=fcigtab(x)
  f = x(1,:).^2 + 1e8*x(end,:).^2 + 1e4*sum(x(2:(end-1),:).^2, 1);
  
function f=ftablet(x)
  f = 1e6*x(1,:).^2 + sum(x(2:end,:).^2, 1);
  
function f=felli(x, lgscal, expon, expon2)
  % lgscal: log10(axis ratio)
  % expon: x_i^expon, sphere==2
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  if nargin < 2 || isempty(lgscal), lgscal = 3; end
  if nargin < 3 || isempty(expon), expon = 2; end
  if nargin < 4 || isempty(expon2), expon2 = 1; end

  f=((10^(lgscal*expon)).^((0:N-1)/(N-1)) * abs(x).^expon).^(1/expon2);
%  if rand(1,1) > 0.015
%    f = NaN;
%  end


function f=fellii(x, scal)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  if nargin < 2
    scal = 1;
  end
  f= (scal*(1:N)).^2 * (x).^2;

function f=fellirot(x)
  N = size(x,1);
  global ORTHOGONALCOORSYSTEM_G
  if isempty(ORTHOGONALCOORSYSTEM_G) ...
	|| length(ORTHOGONALCOORSYSTEM_G) < N ...
	|| isempty(ORTHOGONALCOORSYSTEM_G{N})
    coordinatesystem(N);
  end
  f = felli(ORTHOGONALCOORSYSTEM_G{N}*x);
  
function coordinatesystem(N)
  if nargin < 1 || isempty(N)
    arN = 2:30;
  else
    arN = N;
  end
  global ORTHOGONALCOORSYSTEM_G
  ORTHOGONALCOORSYSTEM_G{1} = 1; 
  for N = arN
    ar = randn(N,N);
    for i = 1:N 
      for j = 1:i-1
	ar(:,i) = ar(:,i) - ar(:,i)'*ar(:,j) * ar(:,j);
      end
      ar(:,i) = ar(:,i) / norm(ar(:,i));
    end
    ORTHOGONALCOORSYSTEM_G{N} = ar; 
  end

function f=fplane(x)
  f=x(1);

function f=ftwoaxes(x)
  f = sum(x(1:floor(end/2)).^2) + 1e6*sum(x(floor(1+end/2):end).^2);

function f=fparabR(x)
  f = -x(1) + 100*sum(x(2:end).^2);

function f=fsharpR(x)
  f = abs(-x(1)) + 30*norm(x(2:end));
  
function f=frosen(x)
  if size(x,1) < 2 error('dimension must be greater one'); end
  f = 1e2*sum((x(1:end-1).^2 - x(2:end)).^2) + sum((x(1:end-1)-1).^2);
  
function f=frosenmodif(x)
  f = 74 + 100*(x(2)-x(1)^2)^2 + (1-x(1))^2 ...
      - 400*exp(-sum((x+1).^2)/2/0.05);
  
function f=fschwefelrosen1(x)
  % in [-10 10] 
  f=sum((x.^2-x(1)).^2 + (x-1).^2);
  
function f=fschwefelrosen2(x)
  % in [-10 10] 
  f=sum((x(2:end).^2-x(1)).^2 + (x(2:end)-1).^2);

function f=fdiffpow(x)
  N = size(x,1); if N < 2 error('dimension must be greater one'); end
  f=sum(abs(x).^(2+10*(0:N-1)'/(N-1)));

%%% Multimodal functions 

function f=fackley(x)
  % -32.768..32.768
  % Adding a penalty outside the interval is recommended,  
  % because for large step sizes, fackley imposes like frand
  % 
  N = size(x,1); 
  f = 20-20*exp(-0.2*sqrt(sum(x.^2)/N)); 
  f = f + (exp(1) - exp(sum(cos(2*pi*x))/N));
  % add penalty outside the search interval
  f = f + sum((x(x>32.768)-32.768).^2) + sum((x(x<-32.768)+32.768).^2);
  
function f = fbohachevsky(x)
 % -15..15
  f = sum(x(1:end-1).^2 + 2 * x(2:end).^2 - 0.3 * cos(3*pi*x(1:end-1)) ...
	  - 0.4 * cos(4*pi*x(2:end)) + 0.7);
  
function f=fconcentric(x)
  % in  +-600
  s = sum(x.^2);
  f = s^0.25 * (sin(50*s^0.1)^2 + 1);

function f=fgriewank(x)
  % in [-600 600]
  N = size(x,1);
  f = 1 - prod(cos(x'./sqrt(1:N))) + sum(x.^2)/4e3;
  % f = f + 1e4*sum(x(abs(x)>5).^2);
  % if sum(x(abs(x)>5).^2) > 0
  %   f = 1e4 * sum(x(abs(x)>5).^2) + 1e8 * sum(x(x>5)).^2;
  % end
  
function f=frastrigin(x, scal, skewfac, skewstart, amplitude)
% by default multi-modal about between -30 and 30
  if nargin < 5 || isempty(amplitude)
    amplitude = 10;
  end
  if nargin < 4 || isempty(skewstart)
    skewstart = 0;
  end
  if nargin < 3 || isempty(skewfac)
    skewfac = 1;
  end
  if nargin < 2 || isempty(scal)
    scal = 1;
  end
  N = size(x,1); 
  scale = 1;
  if N > 1
    scale=scal.^((0:N-1)'/(N-1));
  end
  % simple version: 
  % f = amplitude*(N - sum(cos(2*pi*(scale.*x)))) + sum((scale.*x).^2);

  % skew version: 
  y = scale.*x;
  idx = find(x > skewstart);
  if ~isempty(idx)
    y(idx) =  skewfac*x(idx);
  end
  f = amplitude * (N-sum(cos(2*pi*y))) + sum(y.^2);
  
function f = fschaffer(x)
 % -100..100
  N = size(x,1);
  s = x(1:N-1).^2 + x(2:N).^2;
  f = sum(s.^0.25 .* (sin(50*s.^0.1).^2+1));

function f=fschwefelmult(x)
  % -500..500
  N = size(x,1); 
  %   f = - sum(x.*sin(sqrt(abs(x))));
  f = 418.9829*N - 1.27275661e-5*N - sum(x.*sin(sqrt(abs(x))));
  % penalty term 
  f = f + sum(x(abs(x)>500).^2);
  
function f=ftwomax(x)
  % Boundaries at +/-5
  N = size(x,1); 
  f = -abs(sum(x)) + 5*N;

function f=ftwomaxtwo(x)
  % Boundaries at +/-10
  N = size(x,1); 
  f = abs(sum(x));
  if f > 30
    f = f - 30;
  end
  f = -f;

function f=frand(x)
  f=1/(1-rand) - 1;

% Changes: 
% 07/09: tolhistfun as termination criterion added, "<" changed to
%        "<=" also for TolFun to allow for stopping on zero difference. 
%        Name tolfunhist clashes with option tolfun. 
% 07/07: hsig threshold made slighly smaller for large dimension, 
%        useful for lambda < lambda_default. 
% 07/06: boundary handling: scaling in the boundary handling
%        is omitted now, see bnd.flgscale. This seems not to
%        have a big impact. Using the scaling is worse on rotated
%        functions, but better on separable ones. 
% 07/05: boundary handling: weight i is not incremented anymore
%        if xmean(i) moves towards the feasible space. Increment
%        factor changed to 1.2 instead of 1.1. 
% 07/05: boundary handling code simplified not changing the algorithm
% 07/04: bug removed for saving in octave
% 06/11/10: more testing of outcome of eig, fixed max(D) to max(diag(D))
% 06/10/21: conclusive final bestever assignment in the end 
% 06/10/21: restart and incpopsize option implemented for restarts
%        with increasing population size, version 2.50. 
% 06/09/16: output argument bestever inserted again for convenience and
%        backward compatibility
% 06/08: output argument out and struct out reorganized. 
% 06/01: Possible parallel evaluation included as option EvalParallel
% 05/11: Compatibility to octave implemented, package octave-forge
%   is needed. 
% 05/09: Raise of figure and waiting for first plots improved
% 05/01: Function coordinatesystem cleaned up. 
% 05/01: Function prctile, which requires the statistics toolbox,
%        replaced by myprctile. 
% 05/01: Option warnonequalfunctionvalues included. 
% 04/12: Decrease of sigma removed. Problems on fsectorsphere can 
%        be addressed better by adding search space boundaries. 
% 04/12: Boundary handling simpyfied. 
% 04/12: Bug when stopping criteria tolx or tolupx are vectors. 
% 04/11: Three input parameters are obligatory now. 
% 04/11: Bug in boundary handling removed: Boundary weights can decrease now. 
% 04/11: Normalization for boundary weights scale changed. 
% 04/11: VerboseModulo option bug removed. Documentation improved. 
% 04/11: Condition for increasing boundary weights changed.
% 04/10: Decrease of sigma when fitness is getting consistenly
%        worse. Addresses the problems appearing on fsectorsphere for
%        large population size.
% 04/10: VerboseModulo option included. 
% 04/10: Bug for condition for increasing boundary weights removed.
% 04/07: tolx depends on initial sigma to achieve scale invariance
%        for this stopping criterion. 
% 04/06: Objective function value NaN is not counted as function
%        evaluation and invokes resampling of the search point. 
% 04/06: Error handling for eigenvalue beeing zero (never happens
%        with default parameter setting)
% 04/05: damps further tuned for large mueff 
%      o Details for stall of pc-adaptation added (variable hsig 
%        introduced). 
% 04/05: Bug in boundary handling removed: A large initial SIGMA was
%        corrected not until *after* the first iteration, which could
%        lead to a complete failure.
% 04/05: Call of function range (works with stats toolbox only) 
%        changed to myrange. 
% 04/04: Parameter cs depends on mueff now and damps \propto sqrt(mueff)
%        instead of \propto mueff. 
%      o Initial stall to adapt C (flginiphase) is removed and
%        adaptation of pc is stalled for large norm(ps) instead.
%      o Returned default options include documentation. 
%      o Resume part reorganized.
% 04/03: Stopflag becomes cell-array. 

% ---------------------------------------------------------------
% CMA-ES: Evolution Strategy with Covariance Matrix Adaptation for
% nonlinear function minimization. To be used under the terms of the
% GNU General Public License (http://www.gnu.org/copyleft/gpl.html).
% Author: Nikolaus Hansen, 2001/3. e-mail: hansen@bionik.tu-berlin.de
% URL:http://www.bionik.tu-berlin.de/user/niko
% References: See below. 
% ---------------------------------------------------------------
%
% GENERAL PURPOSE: The CMA-ES (Evolution Strategy with Covariance
% Matrix Adaptation) is a robust search method which should be
% applied, if derivative based methods, e.g. quasi-Newton BFGS or
% conjucate gradient, (supposably) fail due to a rugged search
% landscape (e.g. noise, local optima, outlier, etc.). On smooth
% landscapes CMA-ES is roughly ten times slower than BFGS. For up to
% N=10 variables even the simplex direct search method (Nelder & Mead)
% is often faster, but far less robust than CMA-ES.  To see the
% advantage of the CMA, it will usually take at least 30*N and up to
% 300*N function evaluations, where N is the search problem dimension.
% On considerably hard problems the complete search (a single run) is
% expected to take at least 30*N^2 and up to 300*N^2 function
% evaluations.
%
% SOME MORE COMMENTS: 
% The adaptation of the covariance matrix (e.g. by the CMA) is
% equivalent to a general linear transformation of the problem
% coding. Nevertheless every problem specific knowlegde about the best
% linear transformation should be exploited before starting the
% search. That is, an appropriate a priori transformation should be
% applied to the problem. This also makes the identity matrix as
% initial covariance matrix the best choice.
%
% The strategy parameter lambda (population size, opts.PopSize) is the
% preferred strategy parameter to play with.  If results with the
% default strategy are not satisfactory, increase the population
% size. (Remark that the crucial parameter mu (opts.ParentNumber) is
% increased proportionally to lambda). This will improve the
% strategies capability of handling noise and local minima. We
% recomment successively increasing lambda by a factor of about three,
% starting with initial values between 5 and 20. Casually, population
% sizes even beyond 1000+100*N can be sensible.
%
%
% ---------------------------------------------------------------
%%% REFERENCES
%
% The equation numbers refer to 
% Hansen, N. and S. Kern (2004). Evaluating the CMA Evolution
% Strategy on Multimodal Test Functions.  Eighth International
% Conference on Parallel Problem Solving from Nature PPSN VIII,
% Proceedings, pp. 282-291, Berlin: Springer. 
% (http://www.bionik.tu-berlin.de/user/niko/ppsn2004hansenkern.pdf)
% 
% Further references:
% Hansen, N. and A. Ostermeier (2001). Completely Derandomized
% Self-Adaptation in Evolution Strategies. Evolutionary Computation,
% 9(2), pp. 159-195.
% (http://www.bionik.tu-berlin.de/user/niko/cmaartic.pdf).
%
% Hansen, N., S.D. Mueller and P. Koumoutsakos (2003). Reducing the
% Time Complexity of the Derandomized Evolution Strategy with
% Covariance Matrix Adaptation (CMA-ES). Evolutionary Computation,
% 11(1).  (http://mitpress.mit.edu/journals/pdf/evco_11_1_1_0.pdf).
%

