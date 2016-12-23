%%% =======================================================================
%%% = writeData.m
%%% = Alex Turner
%%% = 06/06/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Saves the output as a csv file.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): obs      -- Structure containing the observations.
%%% =  ( 3): ems      -- Emission sources (and OH) for the box model.
%%% =  ( 4): mod      -- Structure containing the model results.
%%% =  ( 5): IC       -- Initial conditions.
%%% =  ( 6): baseName -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = writeData( St, obs, mod, ems, IC, baseName )

%%% Assemble all the data into a single matrix
% Initialized
nObs = 6;           % NH/SH CH4, NH/SH d13C, NH/SH MCF 
nEms = 8;           % NH CH4, NH d13C, NH MCF, SH CH4, SH d13C, SH MCF, NH OH, SH OH
nSt  = length(St);
dat  = nan(nSt,nObs+nObs+nObs+nEms+1);
% Fill the years
dat(:,1)  = St;
% Fill obs
dat(:,2)  = obs.nh_ch4;
dat(:,3)  = obs.nh_ch4c13;
dat(:,4)  = obs.nh_mcf;
dat(:,5)  = obs.sh_ch4;
dat(:,6)  = obs.sh_ch4c13;
dat(:,7)  = obs.sh_mcf;
% Fill obs error
dat(:,8)  = obs.nh_ch4_err;
dat(:,9)  = obs.nh_ch4c13_err;
dat(:,10)  = obs.nh_mcf_err;
dat(:,11) = obs.sh_ch4_err;
dat(:,12) = obs.sh_ch4c13_err;
dat(:,13) = obs.sh_mcf_err;
% Fill model output
dat(:,14) = mod.nh_ch4;
dat(:,15) = mod.nh_ch4c13;
dat(:,16) = mod.nh_mcf;
dat(:,17) = mod.sh_ch4;
dat(:,18) = mod.sh_ch4c13;
dat(:,19) = mod.sh_mcf;
% Fill emissions
dat(:,20) = ems(:,1);           % NH CH4
dat(:,21) = ems(:,2);           % NH CH4C13
dat(:,22) = ems(:,3);           % NH MCF
dat(:,23) = ems(:,4);           % SH CH4
dat(:,24) = ems(:,5);           % SH CH4C13
dat(:,25) = ems(:,6);           % SH MCF
dat(:,26) = (ems(:,7)-1)*100;   % NH OH
dat(:,27) = (ems(:,8)-1)*100;   % SH OH
% Make the strings
dat_Head = {'JULIAN DATE',...
            'NH CH4 (obs)','NH d13C (obs)','NH MCF (obs)','SH CH4 (obs)','SH d13C (obs)','SH MCF (obs)',...
            'NH CH4 (err)','NH d13C (err)','NH MCF (err)','SH CH4 (err)','SH d13C (err)','SH MCF (err)',...
            'NH CH4 (mod)','NH d13C (mod)','NH MCF (mod)','SH CH4 (mod)','SH d13C (mod)','SH MCF (mod)',...
            'NH CH4 (ems)','NH d13C (ems)','NH MCF (ems)','SH CH4 (ems)','SH d13C (ems)','SH MCF (ems)','NH OH (ems)','SH OH (ems)'};
dat_Unit = {'days since Jan-1-0000',...
                     'ppb',       'permil',         'ppt',         'ppb',       'permil',         'ppt',...
                     'ppb',       'permil',         'ppt',         'ppb',       'permil',         'ppt',...
                     'ppb',       'permil',         'ppt',         'ppb',       'permil',         'ppt',...
                   'Tg/yr',       'permil',       'Gg/yr',       'Tg/yr',       'permil',       'Gg/yr',          '%',          '%'};
fspecA = '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\n';
fspecB = '%9i,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,\n';

%%% Write the data
fid = fopen(sprintf(baseName,'Data'),'w');
fprintf(fid,fspecA,dat_Head{:}); % Header
fprintf(fid,fspecA,dat_Unit{:}); % Units
% Data
for i = 1:nSt
    fprintf(fid,fspecB,dat(i,:));
end
fclose(fid);

%%% Write the ICs
fid = fopen(sprintf(baseName,'ICs'),'w');
fprintf(fid,'NH 12CH4,NH 13CH4,NH MCF,SH 12CH4,SH 13CH4,SH MCF\n'); % Header
fprintf(fid,'ppb,permil,ppt,ppb,permil,ppt\n');                     % Units
fprintf(fid,'%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,%9.3f,\n',IC);           % Data
fclose(fid);

end


%%% =======================================================================
%%% = END
%%% =======================================================================
