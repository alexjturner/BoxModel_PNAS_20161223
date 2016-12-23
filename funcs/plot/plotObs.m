%%% =======================================================================
%%% = plotObs.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Plot the observations and the box-model output.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): model    -- Structure containing the model results.
%%% =  ( 3): data     -- Structure containing the observations.
%%% =  ( 4): baseName -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = plotObs( St, model, data, baseName )

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

%%% Add errorbars?
add_error = true;

%%% Set the axes limits 
yrs           = datevec(St);
xLims         = [datenum(yrs(1,1),1,1),datenum(yrs(end,1),1,1)]';
yLims_ch4     = [ 1500 : 100 :  1900]';
yLims_ch4c13  = [-47.8 : .20 : -46.8]';
yLims_mcf     = [    0 :  50 :   200]';
yLims_Rch4    = [  -10 :   2 :    10]';
yLims_Rch4c13 = [ -0.1 : .05 :   0.1]';
yLims_Rmcf    = [   -3 :   1 :     3]';
yLims_Lmcf    = [2,5,10,20,50,100,200];
yLims_Ich4    = [   80 :   5 :   115]';
yLims_Ich4c13 = [ -0.8 : 0.2 :   0.3]';
% Get the labels
yLims_ch4_lab     = cell(size(yLims_ch4));
yLims_ch4c13_lab  = cell(size(yLims_ch4c13));
yLims_mcf_lab     = cell(size(yLims_mcf));
yLims_Rch4_lab    = cell(size(yLims_Rch4));
yLims_Rch4c13_lab = cell(size(yLims_Rch4c13));
yLims_Rmcf_lab    = cell(size(yLims_Rmcf));
yLims_Lmcf_lab    = cell(size(yLims_Lmcf));
yLims_Ich4_lab    = cell(size(yLims_Rch4));
yLims_Ich4c13_lab = cell(size(yLims_Rch4c13));
for i = 1:length(yLims_ch4);     yLims_ch4_lab{i}     = sprintf('%4.0f',yLims_ch4(i));     end
for i = 1:length(yLims_ch4c13);  yLims_ch4c13_lab{i}  = sprintf('%0.1f',yLims_ch4c13(i));  end
for i = 1:length(yLims_mcf);     yLims_mcf_lab{i}     = sprintf('%1.0f',yLims_mcf(i));     end
for i = 1:length(yLims_Rch4);    yLims_Rch4_lab{i}    = sprintf('%4.0f',yLims_Rch4(i));    end
for i = 1:length(yLims_Rch4c13); yLims_Rch4c13_lab{i} = sprintf('%0.1f',yLims_Rch4c13(i)); end
for i = 1:length(yLims_Rmcf);    yLims_Rmcf_lab{i}    = sprintf('%1.0f',yLims_Rmcf(i));    end
for i = 1:length(yLims_Lmcf);    yLims_Lmcf_lab{i}    = sprintf('%1.0f',yLims_Lmcf(i));    end
for i = 1:length(yLims_Ich4);    yLims_Ich4_lab{i}    = sprintf('%4.0f',yLims_Ich4(i));    end
for i = 1:length(yLims_Ich4c13); yLims_Ich4c13_lab{i} = sprintf('%0.1f',yLims_Ich4c13(i)); end

%%% Make the titles
title_ch4    = 'CH_4 (ppb)';
title_ch4c13 = sprintf('\\delta^{13}CH_{4} (%s)',char(8240));
title_mcf    = 'CH_3CCl_3 (ppt)';

%%% Set the plot options
nhCol    = [204, 179, 102]./256;
shCol    = [ 58, 106, 176]./256;
gloCol   = [ 29, 176,  80]./256;
obsCol   = [  0,   0,   0]./256;
nhOpts  = {'-','Color', nhCol, 'LineWidth', 2};
shOpts  = {'-','Color', shCol, 'LineWidth', 2};
eOpts   = {'-','Color',obsCol, 'LineWidth', 1};
onhOpts = {'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none'};
oshOpts = {'v','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none'};
rnhOpts = {'^','Color',nhCol,'MarkerSize',6,'MarkerFaceColor',nhCol,'MarkerEdgeColor','none'};
rshOpts = {'v','Color',shCol,'MarkerSize',6,'MarkerFaceColor',shCol,'MarkerEdgeColor','none'};
ihdOpts = {'-','Color',obsCol, 'LineWidth', 2};
pOpts   = {'LineWidth',2,'FontName','Helvetica','FontWeight','Bold',...
           'FontSize',16,'YGrid','on','XMinorTick','on','YMinorTick','on'};
tOpts   = {'FontSize',20};
lOpts   = {'HorizontalAlignment','Right','FontSize',18,'FontName','Helvetica','FontWeight','Bold'};
xloc    = .975*(    xLims(end) -     xLims(1)) +     xLims(1);
yloc    = .135*(yLims_ch4(end) - yLims_ch4(1)) + yLims_ch4(1);
xlocO   = .025*(    xLims(end) -     xLims(1)) +     xLims(1);
ylocO   = .875*(yLims_ch4(end) - yLims_ch4(1)) + yLims_ch4(1);
ylocR   = .125*(yLims_Rmcf(end) - yLims_Rmcf(1)) + yLims_Rmcf(1);
spac    = .250*(yLims_ch4(end) - yLims_ch4(1));
spacR   = .250*(yLims_Rmcf(end) - yLims_Rmcf(1));


%%% Plot the observations
h = figure();
% CH4
ax(1) = subplot(3,1,1);p = get(ax(1),'pos');
set(ax(1),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_ch4)
box on
set(ax(1),'yaxislocation','left','YTick',yLims_ch4,'YTickLabel',yLims_ch4_lab)
ylabel(ax(1),title_ch4,tOpts{:})
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4(i); yE = data.nh_ch4_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_ch4(i); yE = data.sh_ch4_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
plot(St, data.nh_ch4,onhOpts{:})
plot(St, data.sh_ch4,oshOpts{:})
plot(St,model.nh_ch4, nhOpts{:})
plot(St,model.sh_ch4, shOpts{:})
text(xlocO,ylocO,            'Observations','Color',obsCol,lOpts{:},'HorizontalAlignment','Left')
text(xloc,yloc+1*spac,'Northern Hemisphere','Color', nhCol,lOpts{:})
text(xloc,yloc+0*spac,'Southern Hemisphere','Color', shCol,lOpts{:})
xlim(xLims)
ylim([yLims_ch4(1),yLims_ch4(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% delta13C
ax(2) = subplot(3,1,2);p = get(ax(2),'pos');
set(ax(2),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_ch4c13)
box on
set(ax(2),'yaxislocation','right','YTick',yLims_ch4c13,'YTickLabel',yLims_ch4c13_lab)
ylabel(ax(2),title_ch4c13,tOpts{:})
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4c13(i); yE = data.nh_ch4c13_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_ch4c13(i); yE = data.sh_ch4c13_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
plot(St, data.nh_ch4c13,onhOpts{:})
plot(St, data.sh_ch4c13,oshOpts{:})
plot(St,model.nh_ch4c13, nhOpts{:})
plot(St,model.sh_ch4c13, shOpts{:})
xlim(xLims)
ylim([yLims_ch4c13(1),yLims_ch4c13(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% MCF
ax(3) = subplot(3,1,3);p = get(ax(3),'pos');
set(ax(3),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
semilogy([NaN,NaN],[NaN,NaN]);
set(gca,pOpts{:},'YTick',yLims_Lmcf,'YTickLabel',yLims_Lmcf_lab)
hold on;
box on
ylabel(title_mcf,tOpts{:})
if add_error
    for i = 1:length(St)
        yV = data.nh_mcf(i); yE = data.nh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            semilogy([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_mcf(i); yE = data.sh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            semilogy([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
semilogy(St, data.nh_mcf,onhOpts{:})
semilogy(St, data.sh_mcf,oshOpts{:})
semilogy(St,model.nh_mcf, nhOpts{:})
semilogy(St,model.sh_mcf, shOpts{:})
xlim(xLims)
ylim([yLims_Lmcf(1),yLims_Lmcf(end)])
datetick('x','yyyy','keeplimits')
% Save the plot
print(h,printOpts{:},sprintf(baseName,'allSpecies'))

%%% Plot the residuals
h = figure();
% CH4
ax(1) = subplot(3,1,1);p = get(ax(1),'pos');
set(ax(1),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_Rch4)
box on
set(ax(1),'yaxislocation','left','YTick',yLims_Rch4,'YTickLabel',yLims_Rch4_lab)
ylabel(ax(1),title_ch4,tOpts{:})
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4(i) - model.nh_ch4(i); yE = data.nh_ch4_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)]-40,[yV-yE,yV+yE],eOpts{:},'Color',nhCol)
        end
        yV = data.sh_ch4(i) - model.sh_ch4(i); yE = data.sh_ch4_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)]+40,[yV-yE,yV+yE],eOpts{:},'Color',shCol)
        end
    end
end
plot(St-40,data.nh_ch4 - model.nh_ch4,rnhOpts{:})
plot(St+40,data.sh_ch4 - model.sh_ch4,rshOpts{:})
xlim(xLims)
ylim([yLims_Rch4(1),yLims_Rch4(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% delta13C
ax(2) = subplot(3,1,2);p = get(ax(2),'pos');
set(ax(2),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_Rch4c13)
box on
set(ax(2),'yaxislocation','right','YTick',yLims_Rch4c13,'YTickLabel',yLims_Rch4c13_lab)
ylabel(ax(2),title_ch4c13,tOpts{:})
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4c13(i) - model.nh_ch4c13(i); yE = data.nh_ch4c13_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)]-40,[yV-yE,yV+yE],eOpts{:},'Color',nhCol)
        end
        yV = data.sh_ch4c13(i) - model.sh_ch4c13(i); yE = data.sh_ch4c13_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)]+40,[yV-yE,yV+yE],eOpts{:},'Color',shCol)
        end
    end
end
plot(St,data.nh_ch4c13 - model.nh_ch4c13,rnhOpts{:})
plot(St,data.sh_ch4c13 - model.sh_ch4c13,rshOpts{:})
xlim(xLims)
ylim([yLims_Rch4c13(1),yLims_Rch4c13(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% MCF
ax(3) = subplot(3,1,3);p = get(ax(3),'pos');
set(ax(3),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_Rmcf)
box on
set(ax(3),'yaxislocation','left','YTick',yLims_Rmcf,'YTickLabel',yLims_Rmcf_lab)
ylabel(ax(3),title_mcf,tOpts{:})
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_mcf(i) - model.nh_mcf(i); yE = data.nh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)]-40,[yV-yE,yV+yE],eOpts{:},'Color',nhCol)
        end
        yV = data.sh_mcf(i) - model.sh_mcf(i); yE = data.sh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)]+40,[yV-yE,yV+yE],eOpts{:},'Color',shCol)
        end
    end
end
plot(St,data.nh_mcf - model.nh_mcf,rnhOpts{:})
plot(St,data.sh_mcf - model.sh_mcf,rshOpts{:})
xlim(xLims)
ylim([yLims_Rmcf(1),yLims_Rmcf(end)])
datetick('x','yyyy','keeplimits')
% Save the plot
print(h,printOpts{:},sprintf(baseName,'Residuals'))

%%% Plot Methylchloroform
h = figure();
semilogy([NaN,NaN],[NaN,NaN]);
hold on;
box on
set(gca,pOpts{:},'YTick',yLims_Lmcf,'YTickLabel',yLims_Lmcf_lab)
ylabel(title_mcf,tOpts{:})
if add_error
    for i = 1:length(St)
        yV = data.nh_mcf(i); yE = data.nh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            semilogy([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_mcf(i); yE = data.sh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            semilogy([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
semilogy(St, data.nh_mcf,onhOpts{:})
semilogy(St, data.sh_mcf,oshOpts{:})
semilogy(St,model.nh_mcf, nhOpts{:})
semilogy(St,model.sh_mcf, shOpts{:})
xlim(xLims)
ylim([yLims_Lmcf(1),yLims_Lmcf(end)])
datetick('x','yyyy','keeplimits')
% Save the plot
print(h,printOpts{:},sprintf(baseName,'MCFobs'))


%%% Plot the interhemispheric differences
h = figure();
% CH4
ax(1) = subplot(2,1,1);p = get(ax(1),'pos');
set(ax(1),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_Ich4)
box on
set(ax(1),'yaxislocation','left','YTick',yLims_Ich4,'YTickLabel',yLims_Ich4_lab)
ylabel(ax(1),sprintf('IHD of %s',title_ch4),tOpts{:})
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4(i) - data.sh_ch4(i);
        yE = sqrt( (data.nh_ch4_err(i)) + (data.sh_ch4_err(i)) );
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
plot(St, data.nh_ch4 -  data.sh_ch4,onhOpts{:})
plot(St, data.nh_ch4 -  data.sh_ch4,oshOpts{:})
plot(St,model.nh_ch4 - model.sh_ch4,ihdOpts{:})
xlim(xLims)
ylim([yLims_Ich4(1),yLims_Ich4(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% delta13C
ax(2) = subplot(2,1,2);p = get(ax(2),'pos');
set(ax(2),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_Ich4c13)
box on
set(ax(2),'yaxislocation','right','YTick',yLims_Ich4c13,'YTickLabel',yLims_Ich4c13_lab)
ylabel(ax(2),sprintf('IHD of %s',title_ch4c13),tOpts{:})
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4c13(i) - data.sh_ch4c13(i);
        yE = sqrt( (data.nh_ch4c13_err(i)) + (data.sh_ch4c13_err(i)) );
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
plot(St, data.nh_ch4c13 -  data.sh_ch4c13,onhOpts{:})
plot(St, data.nh_ch4c13 -  data.sh_ch4c13,oshOpts{:})
plot(St,model.nh_ch4c13 - model.sh_ch4c13,ihdOpts{:})
xlim(xLims)
ylim([yLims_Ich4c13(1),yLims_Ich4c13(end)])
datetick('x','yyyy','keeplimits')
% Save the plot
print(h,printOpts{:},sprintf(baseName,'IHD'))



end


%%% =======================================================================
%%% = END
%%% =======================================================================
