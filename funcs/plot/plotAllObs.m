%%% =======================================================================
%%% = plotAllObs.m
%%% = Alex Turner
%%% = 06/02/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Plot the raw observations and the box-model output.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): obs  -- Structure containing the observations.
%%% =  ( 2): type -- String indicating the type of observations.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = plotAllObs( St, avg_obs, obs, tAvg, type, baseName, deseas, plot_all )

%%% Get file extension
fExten = strsplit(baseName,'.');
fExten = fExten(end);
if strcmp(fExten,'tif')
    printOpts = {'-dtiff','-r300'};
elseif strcmp(fExten,'eps')
    printOpts = {'-depsc2'};
elseif strcmp(fExten,'pdf')
    printOpts = {'-dpdf'};
elseif strcmp(fExten,'png')
    printOpts = {'-dpng'};
end

%%% Get info
sNames = fieldnames(obs.obs);

%%% Define parameters for the block averaging
fDays = 365.25; % Number of days in the block average
if strcmp(tAvg,'month') || strcmp(tAvg,'MONTH') || strcmp(tAvg,'monthly')
    fDays = fDays / 12;
end

%%% What type of simulation are we doing?
if strcmp(type,'ch4')
    yLims  = [ 1500 : 100 : 2000]';
    yTitle = 'CH_4 (ppb)';
    inputN = 'ch4';
elseif strcmp(type,'d13C')
    yLims  = [-48.0 : .20 : -46.8]';
    yTitle = sprintf('\\delta^{13}CH_{4} (%s)',char(8240));
    inputN = 'ch4c13';
elseif strcmp(type,'mcf')
    yLims  = [ 0 : 50 : 200]';
    yTitle = 'CH_3CCl_3 (ppt)';
    inputN = 'mcf';
elseif strcmp(type,'c2h6')
    yLims  = [ 0 : 1000 : 9000]';
    yTitle = 'C_2H_6 (ppt)';
    inputN = 'c2h6';
end
yrs   = datevec(St);
xLims = [datenum(yrs(1,1),1,1),datenum(yrs(end,1),1,1)]';

%%% Define colors
nhCol    = [204, 179, 102]./256;
shCol    = [ 58, 106, 176]./256;
gloCol   = [ 29, 176,  80]./256;
obsCol   = [  0,   0,   0]./256;

%%% Set the plot options
lOpts = {'LineWidth', 2};
pOpts = {'LineWidth',2,'FontName','Helvetica','FontWeight','Bold',...
           'FontSize',16,'YGrid','on','XMinorTick','on','YMinorTick','on'};
tOpts = {'FontSize',20};
cmap  = parula(length(sNames));

%%% Plot
% All together
h = figure();
set(gca,pOpts{:},'YTick',yLims)
box on
ylabel(yTitle,tOpts{:});
xlabel('Year',tOpts{:});
xlim(xLims)
ylim([yLims(1),yLims(end)])
hold on
for i = 1:length(sNames)
    tDat = obs.tim.(sNames{i});
    yDat = obs.obs.(sNames{i});
    if deseas
        yDat = DeseasonalizeData(tDat,yDat,fDays);
    end
    plot(tDat,yDat,'.',lOpts{:},'Color',cmap(i,:));
end
datetick('x','yyyy','keeplimits')
print(h,printOpts{:},sprintf(baseName,type,'IND'))
% All together (colored by NH/SH)
h = figure();
set(gca,pOpts{:},'YTick',yLims)
box on
ylabel(yTitle,tOpts{:});
xlabel('Year',tOpts{:});
xlim(xLims)
ylim([yLims(1),yLims(end)])
hold on
for i = 1:length(sNames)
    if obs.lat.(sNames{i}) > 0
        col = nhCol;
    else
        col = shCol;
    end
    tDat = obs.tim.(sNames{i});
    yDat = obs.obs.(sNames{i});
    if deseas
        yDat = DeseasonalizeData(tDat,yDat,fDays);
    end
    plot(tDat,yDat,'-',lOpts{:},'Color',col);
end
datetick('x','yyyy','keeplimits')
print(h,printOpts{:},sprintf(baseName,type,'NHSH'))
% All together (colored by NH/SH)
h = figure();
set(gca,pOpts{:},'YTick',yLims)
box on
ylabel(yTitle,tOpts{:});
xlabel('Year',tOpts{:});
xlim(xLims)
ylim([yLims(1),yLims(end)])
hold on
for i = 1:length(sNames)
    if obs.lat.(sNames{i}) > 0
        col = nhCol;
    else
        col = shCol;
    end
    tDat = obs.tim.(sNames{i});
    yDat = obs.obs.(sNames{i});
    if deseas
        yDat = DeseasonalizeData(tDat,yDat,fDays);
    end
    plot(tDat,yDat,'-',lOpts{:},'Color',col);
end
if ~strcmp(type,'c2h6')% We don't actually have the hemispheric averages for C2H6
    for i = 1:length(St)
        yV = avg_obs.(sprintf('nh_%s',inputN))(i); yE = avg_obs.(sprintf('nh_%s_err',inputN))(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],'k-','Color',obsCol)
        end
        yV = avg_obs.(sprintf('sh_%s',inputN))(i); yE = avg_obs.(sprintf('sh_%s_err',inputN))(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],'k-','Color',obsCol)
        end
    end
    plot(St,avg_obs.(sprintf('nh_%s',inputN)),'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','none')
    plot(St,avg_obs.(sprintf('sh_%s',inputN)),'v','Color',obsCol,'MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor','none')
end
datetick('x','yyyy','keeplimits')
print(h,printOpts{:},sprintf(baseName,type,'MEAN'))
% Individually
if plot_all
    for i = 1:length(sNames)
        figure();
        set(gca,pOpts{:},'YTick',yLims)
        box on
        ylabel(yTitle,tOpts{:});
        xlabel('Year',tOpts{:});
        h_tit = title(sprintf('Site: %s',sNames{i}),tOpts{:});
        set(h_tit,'interpreter','none')
        xlim(xLims)
        ylim([yLims(1),yLims(end)])
        hold on
        tDat = obs.tim.(sNames{i});
        yDat = obs.obs.(sNames{i});
        if deseas
            yDatO = yDat;
            yDat  = DeseasonalizeData(tDat,yDat,fDays);
            plot(tDat,yDatO,'.',lOpts{:},'Color',cmap(i,:),'MarkerSize',8);
            plot(tDat,yDat,'-',lOpts{:},'Color',cmap(i,:));
        else
            plot(tDat,yDat,'-',lOpts{:},'Color',cmap(i,:));
        end
        datetick('x','yyyy','keeplimits')
    end
end


end


%%% =======================================================================
%%% = END
%%% =======================================================================
