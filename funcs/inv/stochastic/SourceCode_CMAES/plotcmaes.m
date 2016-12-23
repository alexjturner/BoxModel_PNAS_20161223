function plotcmaes(figNb, datei, out, ii)
% PLOTCMAES;
% PLOTCMAES(FIGURENUMBER, FILENAME, OUTSTRUCT, X_MIN_MAX_IDX);
%   plots data from runs of function CMAES. 
%
% PLOTCMAES is in particular practical in plotting the most recent
% data while CMAES is (still) running (in a different Matlab/Octave
% shell) and the plotting option was turned off.
%
% All input arguments are optional and can be omitted or empty:
%   FIGURENUMBER (default is 325): Number of figure where to plot.
%   FILENAME (default is 'variablescmaes.mat'): Filename where to
%     get the data. The function CMAES writes per default into
%     'variablescmaes.mat'.
%   OUTSTRUCT output struct as written by cmaes. Giving this 
%     argument, FILENAME is ignored. 
%   X_MIN_MAX_IDX allows to plot only part of the data. If it is 
%     a scalar it defines the maximum x-value (#FEs). If it has
%     two values it defines min and max x-value. Otherwise it 
%     interpreted as index array for choosing the data from
%     x-vector. 

  stfitfun = [];

  if nargin < 1 || isempty(figNb)
    figNb = 325;
  end
  if nargin < 2 || isempty(datei)
    datei = 'variablescmaes.mat';
  else
    stfitfun = load (datei, 'fitfun'); 
  end
  if nargin < 3 || isempty(out)
    out = load (datei, 'out');
    out = out.out;
    stfitfun = load (datei, 'fitfun'); 
  end

  if isfield(out, 'y1')
    o = out;
  else
    if isfield(out.hist.param, 'evals') % for historical reasons
      o.x = out.hist.param.evals;
    else
      o.x = out.hist.evals;
    end      
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
  
  % choose max x-value
  if nargin > 3 && ~isempty(ii)
    if length(ii) == 1
      ii = find(o.x <= ii);
    elseif length(ii) == 2
      ii = intersect(find(o.x >= ii(1)), ...
		     find(o.x <= ii(2)));
    end
    o.x = o.x(ii);
    o.y1 = o.y1(ii,:);
    o.y2 = o.y2(ii,:);
    o.y3 = o.y3(ii,:);
    o.y4 = o.y4(ii,:);
  end

%%%% OCTAVE %%%%
  if isoctave
    foffset = 1e-99;
    [minfitness minfitnessidx] = min(o.y1(:,1)); 
    dfit = o.y1(:,1)-minfitness; 
    dfit(find(dfit<1e-98)) = NaN;
    if exist('figure') % multiple windows require X11
      figure(427); 
      text(''); % in older Matlab versions: replace with text('')
      hold off;
      semilogy(o.x, o.y4(:,1), ';principle axes lengths;'); hold on;
      semilogy(o.x, o.y4(:,2:end), ';;'); 
      grid('on'); replot;

      figure(426);
      text(''); % in older Matlab versions: replace with text('')
      hold off;
      for i = 1:size(o.y3,2)
	text(o.x(end), o.y3(end, i), [' ' num2str(i)]);
      end
      semilogy(o.x, o.y3(:,1), ';sigma*sqrt(diag(C));'); hold on;
      semilogy(o.x, o.y3(:,2:end), ';;'); 
      
      figure(425); 
      text(''); % in older Matlab versions: replace with text('')
      hold off;
      for i = 1:size(o.y2,2)
	text(o.x(end), o.y2(end, i), [' ' num2str(i)]);
      end
      plot(o.x, o.y2(:,1), ';object variables;'); hold on;
      plot(o.x, o.y2(:,2:end), ';;');
      grid('on'); replot;
      
      figure(424); 
    end % exist figure

    text(''); % in older Matlab versions: replace with text('')
    hold off;
    
    semilogy(o.x, [o.y1(:,1)], ...
	     [';fitness=' num2str(o.y1(end, 1), '%.15g') ';']); hold on
    semilogy(o.x, [o.y1(:,2)], ';sigma;'); 
    semilogy(o.x, [dfit], [';fitness-[min(fitness)=' num2str(minfitness, '%.15g') '];']); 
    % plot best fitness point
    semilogy(o.x(minfitnessidx),abs(o.y1(minfitnessidx,1)),'*c;'); 
    semilogy(o.x, o.y1(:,3), ';axis ratio;'); 
    semilogy(o.x, o.y1(:,4), 'k;max std in coordinate;'); 
    semilogy(o.x, o.y1(:,5), 'k;min std in coordinate;'); 
    semilogy(o.x, o.y1(:,6), ';worst fitness;'); 
    if size(o.y1,2) > 6
      semilogy(o.x, o.y1(:,7:end), ';more data;'); 
    end
    grid('on'); replot;
    return;
  end


%%%% MATLAB %%%%

  % plottrajectory(o.y2, figNb+1000);

  figure(figNb); 
  if size(o.y2, 2) < 100
    minxend = 1.03*o.x(end);
  else
    minxend = 0;
  end

  foffset = 1e-99;
  [minfitness minfitnessidx] = min(o.y1(:,1)); 
  dfit = o.y1(:,1)-minfitness; 
  dfit(dfit<1e-98) = NaN;

  subplot(2,2,1); hold off; 
  if size(o.y1, 2) > 4
    subplot(2,2,1);semilogy(o.x,abs(o.y1(:,4:5)),'-m'); hold on;
  end
  if size(o.y1, 2) > 6
    subplot(2,2,1);semilogy(o.x,abs(o.y1(:,6:7)-min(o.y1(:,1))),'-k'); hold on;
  end
  if size(o.y1, 2) > 5
      subplot(2,2,1);semilogy(o.x,abs(o.y1(:,6:end))+foffset,'-k'); hold on;
  end

  semilogy(o.x,dfit,'-c');hold on;   % fitness difference to best
  idx = find(o.y1(:,1)>1e-98);  % positive values
  subplot(2,2,1);semilogy(o.x(idx), o.y1(idx,1)+foffset, '.b'); hold on; 
  idx = find(o.y1(:,1) < -1e-98);  % negative values
  subplot(2,2,1);semilogy(o.x(idx), abs(o.y1(idx,1))+foffset,'.r'); hold on; 
  subplot(2,2,1);semilogy(o.x,abs(o.y1(:,1))+foffset,'-b'); hold on;
  % plot best fitness point
  subplot(2,2,1);semilogy(o.x(minfitnessidx),abs(o.y1(minfitnessidx,1))+foffset,'*c'); hold on;
  subplot(2,2,1);semilogy(o.x,(o.y1(:,3)),'-r'); hold on;
  subplot(2,2,1);semilogy(o.x,(o.y1(:,2)),'-g'); 
  ax = axis;
  ax(2) = max(minxend, ax(2)); 
  axis(ax);
  if ~isempty(stfitfun)
    if ~ischar(stfitfun.fitfun)
      stfitfun.fitfun = func2str(stfitfun.fitfun);
    end
    text(ax(1), 10^(log10(ax(3))+0.05*(log10(ax(4))-log10(ax(3)))), ...
	 [ stfitfun.fitfun '=' num2str(o.y1(end, 1), ' %.15g') ...
	   ' ' num2str(min(o.y1(:,1)), '(%.15g)')]);
  end
  title('abs(f) (blue), f-min(f) (cyan), Sigma (green), Axis Ratio (red)');
  grid on; 

  subplot(2,2,2); hold off; plot(o.x, o.y2,'-'); 
  ax = axis;
  ax(2) = max(minxend, ax(2)); 
  axis(ax);
  if size(o.y2, 2) < 100
    yy = linspace(ax(3), ax(4), size(o.y2,2))';
    [yyl idx] = sort(o.y2(end,:));
    [muell idx2] = sort(idx);
    hold on;
    plot([o.x(end) ax(2)]', [o.y2(end,:)' yy(idx2)]', '-');
    plot(repmat(o.x(end),2), [ax(3) ax(4)], 'k-');
    %plot(repmat(o.x(end),2), [yyl(end) ax(4)], 'k-');
    %plot(repmat(o.x(end),2), [ax(3) yyl(1)], 'k-');
    for i = 1:length(idx)
      text(ax(2), yy(i), [' ' num2str(idx(i))]);
    end
  end
  title(['Object Variables (' num2str(size(o.y2, 2)) 'D)']);grid on;

  subplot(2,2,3); hold off; semilogy(o.x, o.y3, '-'); 
  ax = axis;
  ax(2) = max(minxend, ax(2)); 
  axis(ax);
  if size(o.y2, 2) < 100
    yy = logspace(log10(ax(3)), log10(ax(4)), size(o.y3,2))';
    [yyl idx] = sort(o.y3(end,:));
    [muell idx2] = sort(idx);
    hold on;
    plot([o.x(end) ax(2)]', [o.y3(end,:)' yy(idx2)]', '-');
    plot(repmat(o.x(end),2), [ax(3) ax(4)], 'k-');
    %plot(repmat(o.x(end),2), [yyl(end) ax(4)], 'k-');
    %plot(repmat(o.x(end),2), [ax(3) yyl(1)], 'k-');
    for i = 1:length(idx)
      text(ax(2), yy(i), [' ' num2str(idx(i))]);
    end
  end
  title('Standard Deviations in All Coordinates');grid on;
  xlabel('function evaluations'); 

  subplot(2,2,4); semilogy(o.x, o.y4, '-');
  ax = axis;
  ax(2) = max(minxend, ax(2)); 
  axis(ax);
  title('Scaling (All Main Axes)');grid on;
  xlabel('function evaluations'); 
  zoom('on'); drawnow;

  if isfield(out.hist, 'user')
    if gcf ~= 1000+figNb % prevents repeated raise of figure
      if  ismember(1000+figNb, findobj('Type', 'figure'))
        set(0, 'CurrentFigure', 1000+figNb); % geht nur, wenn figure schon exisitiert
      else
        figure(1000+figNb);
      end
    end
    semilogy(out.hist.evals(1:end), (out.hist.user')');
    %dat = diff(out.hist.user')';
    %plot(out.hist.evals(2:end), filter(ones(1,100)/100,1,[dat;mean(dat,1)]));
    %plot(out.hist.evals(2:end), filter(ones(1,100)/100,1,[dat]')');
    %plot(out.hist.evals(2:end), (diff(out.hist.user')'));
    grid; 
    drawnow; 
  end
  
  
function plottrajectory(x, figNb)
disp(figNb)
N = size(x, 2); 
if N < 5
  figure(figNb);
  for i = 1:N-1
    for j = i+1:N
      subplot(N-1,N-1,(i-0)+(N-1)*(j-2)); hold off;
      plot(x(:,i), x(:,j), '.-'); hold on;
      plot(x(:,i)-0.2*diff([x(1,i); x(:,i)]), x(:,j)-0.2*diff([x(1,j); x(:,j)]), '.g'); hold on;
      plot(x(1,i), x(1,j), 'og'); hold on;
      for k = 1:size(x,1)
	plot(x(k,i), x(k,j), '.', 'MarkerSize', 10+10*log10(k));
      end
      plot(x(end,i), x(end,j), 'or');
      grid on;
      title('CAVE ASPECT RATIO');
    end
  end
end

% ---------------------------------------------------------------  
% ---------------------------------------------------------------  
function res = isoctave
% any hack to find out whether we are under octave
  s = version;
  res = 0;
  if exist('fflush', 'builtin')
    res = 1;
  end

