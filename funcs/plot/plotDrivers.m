%%% =======================================================================
%%% = plotDrivers.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Plot the drivers for the box-model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): ems_post -- Matrix with the posterior emissions.
%%% =  ( 3): ems_pri  -- Matrix with the prior emissions.
%%% =  ( 4): baseName -- Prefix for the plots.
%%% =  ( 5): dataDir  -- Directory containing the data.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = plotDrivers( St, ems_post, ems_pri, baseName, dataDir )

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

%%% Load the data from McNorton et al.
fname   = sprintf('%sobs/oh/%s',dataDir,'OH_anons.csv');
oh_dat  = csvread(fname,1,0);
oh_dat(oh_dat == 0) = NaN;
oh_time = datenum(oh_dat(:,1),ones(size(oh_dat,1),1),ones(size(oh_dat,1),1));
% Get the OH offset
ind       = oh_time(1) <= St & St <= oh_time(end);
oh_offset = mean((mean([ems_post(ind,7),ems_post(ind,8)],2)-1)*100);

%%% Turn 
yrs          = datevec(St);
xLims        = [datenum(yrs(1,1),1,1),datenum(yrs(end,1),1,1)]';
yLims_ch4    = [ 0    : 200 :   600]';
yLims_ch4c13 = [-55.0 : 0.5 : -48.0]';
yLims_d13C   = [-53.5 : 1.0 : -48.5]';
yLims_mcf    = [    0 : 200 :   800]';
yLims_oh     = [  -15 :   5 :    15]';
yLims_Dch4   = [  -20 :  10 :    40]';
% Get the labels
yLims_Dch4_lab = cell(size(yLims_Dch4));
yLims_d13C_lab = cell(size(yLims_d13C));
yLims_oh_lab   = cell(size(yLims_mcf));
for i = 1:length(yLims_Dch4); yLims_Dch4_lab{i} = sprintf('%4.0f',yLims_Dch4(i)); end
for i = 1:length(yLims_d13C); yLims_d13C_lab{i} = sprintf('%0.1f',yLims_d13C(i)); end
for i = 1:length(yLims_oh);   yLims_oh_lab{i}   = sprintf('%1.0f',yLims_oh(i));   end

title_ch4    = 'CH_4 Emissions (Tg/yr)';
title_ch4c13 = sprintf('\\delta^{13}CH_{4} Composition (%s)',char(8240));
title_mcf    = 'CH_3CCl_3 Emissions (Gg/yr)';
title_oh     = '[OH] anomaly (%)';

%%% Set the plot options
nhCol    = [204, 179, 102]./256;
shCol    = [ 58, 106, 176]./256;
gloCol   = [ 29, 176,  80]./256;
obsCol   = [  0,   0,   0]./256;
nhOptsA  = {'-','Color', nhCol, 'LineWidth', 2};
shOptsA  = {'-','Color', shCol, 'LineWidth', 2};
gloOptsA = {'-','Color',gloCol, 'LineWidth', 2};
nhOptsB  = {'--','Color', nhCol, 'LineWidth', 3};
shOptsB  = {'--','Color', shCol, 'LineWidth', 3};
gloOptsB = {'--','Color',gloCol, 'LineWidth', 3};
pOpts    = {'LineWidth',2,'FontName','Helvetica','FontWeight','Bold',...
            'FontSize',16,'YGrid','on','XMinorTick','on','YMinorTick','on'};
tOpts    = {'FontSize',20};
lOpts    = {'HorizontalAlignment','Left','FontSize',18,'FontName','Helvetica','FontWeight','Bold'};


% Get the x and y locations for the labels
xloc        = .015*(xLims(end) - xLims(1)) + xLims(1);
yLims       = yLims_ch4;
yloc_ch4    = .04*(yLims(end) - yLims(1)) + yLims(1);
spac_ch4    = .075*(yLims(end) - yLims(1));
yLims       = yLims_ch4c13;
yloc_ch4c13 = .03*(yLims(end) - yLims(1)) + yLims(1);
spac_ch4c13 = .075*(yLims(end) - yLims(1));
yLims       = yLims_mcf;
yloc_mcf    = .04*(yLims(end) - yLims(1)) + yLims(1);
spac_mcf    = .075*(yLims(end) - yLims(1));
yLims       = yLims_oh;
yloc_oh     = .04*(yLims(end) - yLims(1)) + yLims(1);
spac_oh     = .075*(yLims(end) - yLims(1));
yLims       = yLims_oh;
yloc_ALL    = .04*(yLims(end) - yLims(1)) + yLims(1);
spac_ALL    = .075*(yLims(end) - yLims(1));

%%% Plot the 3 main drivers (CH4, isotopes, & OH)
h = figure();
% CH4
ax(1) = subplot(3,1,1);p = get(ax(1),'pos');
set(ax(1),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_Dch4)
box on
set(ax(1),'yaxislocation','left','YTick',yLims_Dch4,'YTickLabel',yLims_Dch4_lab)
ylabel(ax(1),sprintf('\\Delta %s',title_ch4),tOpts{:})
hold on
plot(St,ems_post(:,1) - ems_pri(:,1),nhOptsA{:})
plot(St,ems_post(:,4) - ems_pri(:,4),shOptsA{:})
plot(St,(ems_post(:,1)+ems_post(:,4)) - (ems_pri(:,1)+ems_pri(:,4)),gloOptsA{:})
xlim(xLims)
ylim([yLims_Dch4(1),yLims_Dch4(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% delta13C
ax(2) = subplot(3,1,2);p = get(ax(2),'pos');
set(ax(2),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_d13C)
box on
set(ax(2),'yaxislocation','right','YTick',yLims_d13C,'YTickLabel',yLims_d13C_lab)
ylabel(ax(2),title_ch4c13,tOpts{:})
hold on
plot(St,ems_post(:,2),nhOptsA{:})
plot(St,ems_post(:,5),shOptsA{:})
xlim(xLims)
ylim([yLims_d13C(1),yLims_d13C(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% OH
ax(3) = subplot(3,1,3);p = get(ax(3),'pos');
set(ax(3),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_oh)
box on
set(ax(3),'yaxislocation','left','YTick',yLims_oh,'YTickLabel',yLims_oh_lab)
ylabel(ax(3),title_oh,tOpts{:})
hold on
if any(~isnan(ems_pri));
    for i = 2:size(oh_dat,2)
        plot(oh_time,oh_offset + oh_dat(:,i),'k-','LineWidth',1,'Color',obsCol)
    end
end
plot(St,(ems_post(:,7)-1)*100,nhOptsA{:})
plot(St,(ems_post(:,8)-1)*100,shOptsA{:})
plot(St,(mean([ems_post(:,7),ems_post(:,8)],2)-1)*100,gloOptsA{:})
text(xloc,yloc_ALL-0*spac_ALL,'Northern Hemisphere','Color', nhCol,lOpts{:})
text(xloc,yloc_ALL-1*spac_ALL,'Southern Hemisphere','Color', shCol,lOpts{:})
text(xloc,yloc_ALL-2*spac_ALL,             'Global','Color',gloCol,lOpts{:})
xlim(xLims)
ylim([yLims_oh(1),yLims_oh(end)])
datetick('x','yyyy','keeplimits')
% Save the plot
print(h,printOpts{:},sprintf(baseName,'mainDrivers'))

% CH4
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',yLims_ch4)
ylabel(title_ch4,tOpts{:})
plot(St,ems_post(:,1),nhOptsA{:})
plot(St,ems_post(:,4),shOptsA{:})
plot(St,ems_post(:,1)+ems_post(:,4),gloOptsA{:})
plot(St,ems_pri(:,1),nhOptsB{:})
plot(St,ems_pri(:,4),shOptsB{:})
plot(St,ems_pri(:,1)+ems_pri(:,4),gloOptsB{:})
text(xloc,yloc_ch4+2*spac_ch4,             'Global','Color',gloCol,lOpts{:})
text(xloc,yloc_ch4+1*spac_ch4,'Northern Hemisphere','Color', nhCol,lOpts{:})
text(xloc,yloc_ch4+0*spac_ch4,'Southern Hemisphere','Color', shCol,lOpts{:})
xlim(xLims)
ylim([yLims_ch4(1),yLims_ch4(end)])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'CH4emissions'))

% delta13C
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',yLims_ch4c13)
ylabel(title_ch4c13,tOpts{:})
plot(St,ems_post(:,2),nhOptsA{:})
plot(St,ems_post(:,5),shOptsA{:})
plot(St,ems_pri(:,2),nhOptsB{:})
plot(St,ems_pri(:,5),shOptsB{:})
text(xloc,yloc_ch4c13+1*spac_ch4c13,'Northern Hemisphere','Color', nhCol,lOpts{:})
text(xloc,yloc_ch4c13+0*spac_ch4c13,'Southern Hemisphere','Color', shCol,lOpts{:})
xlim(xLims)
ylim([yLims_ch4c13(1),yLims_ch4c13(end)])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'CH4C13composition'))

% MCF
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',yLims_mcf)
ylabel(title_mcf,tOpts{:})
plot(St,ems_post(:,3),nhOptsA{:})
plot(St,ems_post(:,6),shOptsA{:})
plot(St,ems_pri(:,3),nhOptsB{:})
plot(St,ems_pri(:,6),shOptsB{:})
text(xloc,yloc_mcf+1*spac_mcf,'Northern Hemisphere','Color', nhCol,lOpts{:})
text(xloc,yloc_mcf+0*spac_mcf,'Southern Hemisphere','Color', shCol,lOpts{:})
xlim(xLims)
ylim([yLims_mcf(1),yLims_mcf(end)])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'MCFemissions'))

% MCF (log)
h = figure();
box on
hold on
set(gca,pOpts{:},'YScale','log')
ylabel(title_mcf,tOpts{:})
mcf_A = ems_post(:,3);mcf_A(mcf_A<=0) = NaN;
mcf_B = ems_pri(:,3);mcf_B(mcf_B<=0) = NaN;
plot(St,mcf_A,nhOptsA{:})
plot(St,mcf_B,nhOptsB{:})
xlim(xLims)
ylim([1e-1,1e3])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'MCFemissions_Log'))

% OH
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',yLims_oh)
ylabel(title_oh,tOpts{:})
if any(~isnan(ems_pri));
    for i = 2:size(oh_dat,2)
        plot(oh_time,oh_offset + oh_dat(:,i),'k-','LineWidth',1,'Color',obsCol)
    end
end
plot(St,(ems_post(:,7)-1)*100,nhOptsA{:})
plot(St,(ems_pri(:,7)-1)*100,nhOptsB{:})
plot(St,(ems_post(:,8)-1)*100,shOptsA{:})
plot(St,(ems_pri(:,8)-1)*100,shOptsB{:})
plot(St,(mean([ems_post(:,7),ems_post(:,8)],2)-1)*100,gloOptsA{:})
plot(St,(mean([ems_pri(:,7),ems_pri(:,8)],2)-1)*100,gloOptsB{:})
text(xloc,yloc_oh+2*spac_oh,             'Global','Color',gloCol,lOpts{:})
text(xloc,yloc_oh+1*spac_oh,'Northern Hemisphere','Color', nhCol,lOpts{:})
text(xloc,yloc_oh+0*spac_oh,'Southern Hemisphere','Color', shCol,lOpts{:})
xlim(xLims)
ylim([yLims_oh(1),yLims_oh(end)])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'OHanomaly'))

% OH ratio
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',[.85:.05:1.15])
ylabel('NH/SH OH Ratio',tOpts{:})
plot(St,ems_post(:,7)./ems_post(:,8),gloOptsA{:})
plot(St,ems_pri(:,7)./ems_pri(:,8),gloOptsB{:})
% Patra
plot(xLims,.97*[1,1],'k-','LineWidth',1,'Color',obsCol)
xlim(xLims)
ylim([0.85,1.15])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'OHratio'))

% Methane lifetime
tau_post = getLifetime(St,ems_post);
tau_pri  = getLifetime(St,ems_pri);
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',[7.5:0.5:10.5])
ylabel('CH_4 Lifetime (yr)',tOpts{:})
plot(St,tau_post.nh,nhOptsA{:})
plot(St,tau_pri.nh,nhOptsB{:})
plot(St,tau_post.sh,shOptsA{:})
plot(St,tau_pri.sh,shOptsB{:})
plot(St,tau_post.glo,gloOptsA{:})
plot(St,tau_pri.glo,gloOptsB{:})
xlim(xLims)
ylim([7.5,10.5])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'CH4lifetime'))

% MCF lifetime
tau_post = getMCFLifetime(St,ems_post);
tau_pri  = getMCFLifetime(St,ems_pri);
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',[4.0:0.5:7.0])
ylabel('MCF Lifetime (yr)',tOpts{:})
plot(St,tau_post.nh,nhOptsA{:})
plot(St,tau_pri.nh,nhOptsB{:})
plot(St,tau_post.sh,shOptsA{:})
plot(St,tau_pri.sh,shOptsB{:})
plot(St,tau_post.glo,gloOptsA{:})
plot(St,tau_pri.glo,gloOptsB{:})
xlim(xLims)
ylim([4.0,7.0])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'MCFlifetime'))

% Delta CH4 emissions NH/SH
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',yLims_Dch4)
ylabel(sprintf('\\Delta %s',title_ch4),tOpts{:})
plot(St,ems_post(:,1) - ems_pri(:,1),nhOptsA{:})
plot(St,ems_post(:,4) - ems_pri(:,4),shOptsA{:})
plot(St,(ems_post(:,1)+ems_post(:,4)) - (ems_pri(:,1)+ems_pri(:,4)),gloOptsA{:})
xlim(xLims)
ylim([yLims_Dch4(1),yLims_Dch4(end)])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'CH4emissions_Delta'))

% CH4 emissions total
h = figure();
box on
hold on
set(gca,pOpts{:},'YTick',[500:10:620])
ylabel(title_ch4,tOpts{:})
plot(St,ems_post(:,1)+ems_post(:,4),gloOptsA{:})
xlim(xLims)
ylim([500,600])
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'CH4emissionsGlobal'))

%%% Plot the 3 OH analysis
h = figure();
% OH
ax(1) = subplot(3,1,1);p = get(ax(1),'pos');
set(ax(1),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_oh)
box on
set(ax(1),'yaxislocation','left','YTick',yLims_oh,'YTickLabel',yLims_oh_lab)
ylabel(ax(1),title_oh,tOpts{:})
hold on
if any(~isnan(ems_pri));
    for i = 2:size(oh_dat,2)
        plot(oh_time,oh_offset + oh_dat(:,i),'k-','LineWidth',1,'Color',[0,0,0])
    end
end
plot(St,(ems_post(:,7)-1)*100,nhOptsA{:})
plot(St,(ems_post(:,8)-1)*100,shOptsA{:})
plot(St,(mean([ems_post(:,7),ems_post(:,8)],2)-1)*100,gloOptsA{:})
xlim(xLims)
ylim([yLims_oh(1),yLims_oh(end)])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% OH Ratio
ax(2) = subplot(3,1,2);p = get(ax(2),'pos');
set(ax(2),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_d13C)
box on
set(ax(2),'yaxislocation','right','YTick',[.9:.05:1.1])
ylabel(ax(2),'NH/SH OH Ratio',tOpts{:})
hold on
set(gca,pOpts{:},'YTick',[.9:.05:1.1])
plot(St,ems_post(:,7)./ems_post(:,8),gloOptsA{:})
plot(xLims,.97*[1,1],'k-','LineWidth',1,'Color',obsCol)% Patra
xlim(xLims)
ylim([0.9,1.1])
datetick('x','yyyy','keeplimits')
set(gca,'XTickLabel',{})
% Methane lifetime
tau_post = getLifetime(St,ems_post);
ax(3) = subplot(3,1,3);p = get(ax(3),'pos');
set(ax(3),'pos',[p(1),p(2)-0.04625,p(3)-.04,p(4)+0.06625])
set(gca,pOpts{:},'YTick',yLims_oh)
box on
set(ax(3),'yaxislocation','left','YTick',[8:0.5:10])
ylabel(ax(3),'CH_4 Lifetime (yr)',tOpts{:})
hold on
plot(St,tau_post.nh,nhOptsA{:})
plot(St,tau_post.sh,shOptsA{:})
plot(St,tau_post.glo,gloOptsA{:})
plot(xLims,9.1*[1,1],'k-','LineWidth',1,'Color',obsCol)% Prather
xlim(xLims)
ylim([8,10])
datetick('x','yyyy','keeplimits')
% Save the plot
print(h,printOpts{:},sprintf(baseName,'OHanalysis'))


end


%%% =======================================================================
%%% = END
%%% =======================================================================
