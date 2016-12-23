%%% =======================================================================
%%% = plotEdObs.m
%%% = Alex Turner
%%% = 06/02/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Plot the observations using Ed's averages.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): obs  -- Structure containing the observations.
%%% =  ( 2): type -- String indicating the type of observations.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = plotEdObs( St, ajt_obs, ed_obs, baseName )

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


figure();plot(St,ajt_obs.nh_ch4,'b-',St,ed_obs.nh_ch4,'b--',St,ajt_obs.sh_ch4,'r-',St,ed_obs.sh_ch4,'r--');
figure();plot(St,ajt_obs.nh_ch4c13,'b-',St,ed_obs.nh_ch4c13,'b--',St,ajt_obs.sh_ch4c13,'r-',St,ed_obs.sh_ch4c13,'r--');
figure();plot(St,ajt_obs.nh_mcf,'b-',St,ed_obs.nh_mcf,'b--',St,ajt_obs.sh_mcf,'r-',St,ed_obs.sh_mcf,'r--');

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
obsCol   = [  0,   0,   0]./256;
EnhOpts  = {'--','Color', nhCol, 'LineWidth', 3};
EshOpts  = {'--','Color', shCol, 'LineWidth', 3};
AnhOpts  = {'-', 'Color', nhCol, 'LineWidth', 2};
AshOpts  = {'-', 'Color', shCol, 'LineWidth', 2};
pOpts    = {'LineWidth',2,'FontName','Helvetica','FontWeight','Bold',...
           'FontSize',16,'YGrid','on','XMinorTick','on','YMinorTick','on'};
tOpts   = {'FontSize',20};
lOpts   = {'HorizontalAlignment','Right','FontSize',18,'FontName','Helvetica','FontWeight','Bold'};
xloc    = .975*(    xLims(end) -     xLims(1)) +     xLims(1);
yloc    = .135*(yLims_ch4(end) - yLims_ch4(1)) + yLims_ch4(1);
xlocO   = .025*(    xLims(end) -     xLims(1)) +     xLims(1);
ylocO   = .875*(yLims_ch4(end) - yLims_ch4(1)) + yLims_ch4(1);
spac    = .250*(yLims_ch4(end) - yLims_ch4(1));


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
plot(St, ed_obs.nh_ch4,EnhOpts{:})
plot(St, ed_obs.sh_ch4,EshOpts{:})
plot(St,ajt_obs.nh_ch4,AnhOpts{:})
plot(St,ajt_obs.sh_ch4,AshOpts{:})
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
plot(St, ed_obs.nh_ch4c13,EnhOpts{:})
plot(St, ed_obs.sh_ch4c13,EshOpts{:})
plot(St,ajt_obs.nh_ch4c13,AnhOpts{:})
plot(St,ajt_obs.sh_ch4c13,AshOpts{:})
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
semilogy(St, ed_obs.nh_mcf,EnhOpts{:})
semilogy(St, ed_obs.sh_mcf,EshOpts{:})
semilogy(St,ajt_obs.nh_mcf,AnhOpts{:})
semilogy(St,ajt_obs.sh_mcf,AshOpts{:})
xlim(xLims)
ylim([yLims_Lmcf(1),yLims_Lmcf(end)])
datetick('x','yyyy','keeplimits')
% Save the plot
print(h,printOpts{:},sprintf(baseName,'allSpecies'))


end


%%% =======================================================================
%%% = END
%%% =======================================================================
