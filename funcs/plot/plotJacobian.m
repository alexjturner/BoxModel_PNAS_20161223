%%% =======================================================================
%%% = plotJacobian.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Plot different rows of the Jacobian for the box-model.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St           -- Our time vector.
%%% =  ( 2): jacobian_ems -- Jacobian matrix for the sources.
%%% =  ( 2): tRes         -- String containing the temporal resolution.
%%% =  ( 4): baseName     -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = plotJacobian( St, jacobian_ems, tRes, baseName )

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

%%% Are we perturbing a month or a year?
nT    = length(St);
yrFac = 12;
pTime = '(1 MONTH)';
if strcmp(tRes,'year') || strcmp(tRes,'YEAR') || strcmp(tRes,'yearly')
    yrFac = 1;
    pTime = '(1 YEAR)';
end

%%% Titles and limits
yrs       = datevec(St);
xLims     = [datenum(yrs(1,1),1,1),datenum(yrs(end,1),1,1)]';
title_ch4 = sprintf('S_x = NH EMISSION CHANGE AT 1988 %s',pTime);
title_oh  = sprintf('S_x = NH OH CHANGE AT 1988 %s',pTime);
yLab_ch4a = 'dCH_4/dS_x (ppb/Tg)';
yLab_ch4b = sprintf('d\\delta^{13}CH_4/dS_x (%s/Tg)',char(8240));
yLab_oha  = 'dCH_4/dS_x (ppb/%OH)';
yLab_ohb  = sprintf('d\\delta^{13}CH_4/dS_x (%s/%%OH)',char(8240));

%%% Set the plot options
nhCol    = [204, 179, 102]./256;
shCol    = [ 58, 106, 176]./256;
gloCol   = [ 29, 176,  80]./256;
obsCol   = [  0,   0,   0]./256;
pOpts    = {'LineWidth',2,'FontName','Helvetica','FontWeight','Bold',...
            'FontSize',16,'YGrid','on','XMinorTick','on','YMinorTick','on'};
tOpts    = {'FontSize',20};
lOpts = {'LineWidth',2,'Color'};


%%% NH CH4 emission change on CH4
h = figure();
box on
hold on
set(gca,pOpts{:})
ylabel(yLab_ch4a,tOpts{:})
title(title_ch4,tOpts{:})
plot(St,jacobian_ems(0*nT+0*nT+[1:nT],9*yrFac,1),lOpts{:},nhCol)
plot(St,jacobian_ems(0*nT+3*nT+[1:nT],9*yrFac,1),lOpts{:},shCol)
xlim(xLims)
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'CH4ems_CH4'))

%%% NH CH4 emission change on CH4C13
h = figure();
box on
hold on
set(gca,pOpts{:})
ylabel(yLab_ch4b,tOpts{:})
title(title_ch4,tOpts{:})
plot(St,jacobian_ems(1*nT+0*nT+[1:nT],9*yrFac,1),lOpts{:},nhCol)
plot(St,jacobian_ems(1*nT+3*nT+[1:nT],9*yrFac,1),lOpts{:},shCol)
xlim(xLims)
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'CH4ems_CH4C13'))

%%% NH OH change on CH4
h = figure();
box on
hold on
set(gca,pOpts{:})
ylabel(yLab_oha,tOpts{:})
title(title_oh,tOpts{:})
plot(St,jacobian_ems(0*nT+0*nT+[1:nT],9*yrFac,7),lOpts{:},nhCol)
plot(St,jacobian_ems(0*nT+3*nT+[1:nT],9*yrFac,7),lOpts{:},shCol)
xlim(xLims)
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'OH_CH4'))

%%% NH OH change on CH4C13
h = figure();
box on
hold on
set(gca,pOpts{:})
ylabel(yLab_ohb,tOpts{:})
title(title_oh,tOpts{:})
plot(St,jacobian_ems(1*nT+0*nT+[1:nT],9*yrFac,7),lOpts{:},nhCol)
plot(St,jacobian_ems(1*nT+3*nT+[1:nT],9*yrFac,7),lOpts{:},shCol)
xlim(xLims)
datetick('x','yyyy','keeplimits')
print(h,printOpts{:}, sprintf(baseName,'OH_CH4C13'))


end


%%% =======================================================================
%%% = END
%%% =======================================================================
