%%% =======================================================================
%%% = getTime.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Evaluates a distribution over a given domain.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): sYear -- Start year for the simulation.
%%% =  ( 2): eYear -- End year for the simulation.
%%% =  ( 3): tRes  -- String containing the temporal resolution.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): St -- Our time vector.
%%% =======================================================================

function [ St ] = getTime( sYear, eYear, tRes )

if strcmp(tRes,'year') || strcmp(tRes,'YEAR') || strcmp(tRes,'yearly')
    St = [sYear : eYear]';
    St = datenum(St,ones(size(St)),ones(size(St)));
else % Otherwise use monthly
    yr  = [sYear : eYear];
    mon = [1:12];
    St = nan(length(yr)*length(mon),1);
    k = 1;
    for i = 1:length(yr)
        for j = 1:length(mon)
            St(k) = datenum(yr(i),mon(j),1);
            k    = k + 1;
        end
    end
end

end


%%% =======================================================================
%%% = END
%%% =======================================================================
