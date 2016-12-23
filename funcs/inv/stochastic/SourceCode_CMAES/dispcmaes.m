function data = dispcmaes3(filename, areval, arfunval)
%
% dispcmaes3(filename) displays (prints) time-data
%   from saved variables from cmaes.m.
% 
% dispcmaes3(cma) displays (prints) time-data from cma-variables,
%   previously loadede, e.g. 
%   cma = load(filename); dispcmaes3(cma);
%
% dispcmaes3(filename, areval, arfunval) displays (prints) time-data
%
% So far still a hack. 
% Not yet implemented: AREVAL gives the function evaluation numbers
% (times) of the data that should be plotted. Alternatively ARFUNEVAL
% asks for the respective first entry where the function value was
% lower (i.e. better).
%

  if nargin < 1 || isempty(filename)
    filename = 'variablescmaes.mat';
  end
  if ~ischar(filename)
    cma = filename; % input is a struct rather a filename
  else
    cma = load(filename);
  end
  
  comment = ['#Fevals:   Function Value   (d-worst) ' ...
	'|Axis Ratio|' ...
	'idx:Min SD idx:Max SD']; 

  o = cma.out;
  if ~isfield(o, 'y1')  % convert out to old format
    out = o;
    if ~isfield(out, 'solutions')
      out = setfield(out, 'solutions', out.solution)
    end
    cma.bestever = out.solutions.bestever;
    cma.bestever.counteval = out.solutions.bestever.evals;
    clear o;
    o.x = out.hist.evals;
    o.y1 = [out.hist.recentbest.f' ...
	    out.hist.param.sigma' ...
	    (out.hist.param.maxD ./ out.hist.param.minD)' ...
	    (out.hist.param.sigma.*out.hist.param.maxD)' ...
	    (out.hist.param.sigma.*out.hist.param.minD)' ...
	    out.hist.recentworst.f'];
    o.y2 = out.hist.mean.x'; 
    o.y2a = out.hist.recentbest.x';
    o.y3 = out.histParamArr.stds'; 
    o.y4 = out.histParamArr.diagD'; 
  end
  
  
%   history record OUTHIST with columns (1) function evaluation count,
%  (2) function value, (3) axis ratio of search distribution, (4)
%  maximal coordinate wise standard deviation
%  (sigma*sqrt(max(diag(C)))), (5) minimal coordinate wise standard
%  deviation, (6) maximal standard deviation in covariance matrix C;
 
  %   Iterat, #Fevals:   Function Value    (median,worst) |Axis Ratio|idx:Min SD idx:Max SD

  % process input argument areval as max or min/max #Fevals
  if nargin > 1 && ~isempty(areval) && length(areval) < 3
    if length(areval) == 1
      ii = find(o.x <= areval);
    elseif length(areval) == 2
      ii = intersect(find(o.x >= areval(1)), ...
		     find(o.x <= areval(2)));
    end
    o.x = o.x(ii);
    o.y1 = o.y1(ii,:);
    o.y2 = o.y2(ii,:);
    o.y3 = o.y3(ii,:);
    o.y4 = o.y4(ii,:);
  end

  % find indices 
  if nargin > 1 && isequal(areval, 'all')
    idx = 1:length(o.x);
  elseif nargin > 1 && length(areval) > 2
    error('not yet implemented');
  else
    l = length(o.x);
    idx = [2 floor(3:l/20:l-2) l-1 l];
    if ~isnan(o.y1(1,1))
      idx = [1 idx];
    end
  end
  
  if o.x(end) < cma.bestever.counteval
    [muell ibestever] = min(o.y1(:,1));
  else
    ibestever = find(o.x == cma.bestever.counteval);
  end
  if isempty(ibestever) % bestever is not in the displayable indices
    flgaddbestever = 1;
    ibestever = 0; % this is a hack
  else
    idx = unique([ibestever idx]);
    flgaddbestever = 0;
  end

  % Construct the string
  data = [comment repmat(' ', 1, 91-length(comment))]; % this makes sure it becomes a string

  for i = idx
    counteval = o.x(i);
    fitness = o.y1(i, 1);
    fitmax = [];
    if size(o.y1, 2) > 5
      fitmax = o.y1(i, 6);
    end
    axisratio = o.y1(i, 3);
    [minstd minstdidx] = min(o.y3(i,:));
    [maxstd maxstdidx] = max(o.y3(i,:));

    if flgaddbestever ~= 0 && (i == ibestever || counteval > cma.bestever.counteval)
      flgaddbestever = 0;
      dd = [ ...
	    '*' repmat(' ',1,4-floor(log10(max(1,cma.bestever.counteval)))) ...
	    num2str(cma.bestever.counteval) ' : ' ...
            num2str(cma.bestever.f, '%+.12e') ...
	  ];
      data(end+1,:) = [dd repmat(' ', 1, 91-length(dd))];
    end
    if i == ibestever
      str = ['*' repmat(' ',1,4-floor(log10(max(1,counteval))))]; 
    else
      str = repmat(' ',1,5-floor(log10(max(1,counteval)))); 
    end
    
    dd = ...
	[ ...
	    str ... 
	    num2str(counteval) ' : ' ...
            num2str(fitness, '%+.12e') ...
	    ' +(' num2str(fitmax-fitness, '%.0e ') ...
	    ') | ' ...
	    num2str(axisratio, '%4.2e') ' | ' ...
	    repmat(' ',1,1-floor(log10(minstdidx))) num2str(minstdidx) ':' ...
	    num2str(minstd, ' %.1e') ' ' ...
	    repmat(' ',1,1-floor(log10(maxstdidx))) num2str(maxstdidx) ':' ...
	    num2str(maxstd, ' %.1e') ...
	];
    
    data(end+1,:) = [dd repmat(' ', 1, 91-length(dd))];
  end

  if counteval < cma.bestever.counteval % real bestever was not in the displayed part
    dd = [ ...
	'*' repmat(' ',1,4-floor(log10(max(1,cma.bestever.counteval)))) ...
	num2str(cma.bestever.counteval) ' : ' ...
	num2str(cma.bestever.f, '%+.12e') ...
	 ];
    data(end+1,:) = [dd repmat(' ', 1, 91-length(dd))];
  end
  
  disp(data(:,1:85));
