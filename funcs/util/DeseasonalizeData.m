%%% =======================================================================
%%% = DeseasonalizeData.m
%%% = Alex Turner
%%% = 06/02/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Deseasonalizes a timeseries with a stable seasonal filter.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): t = Julian date for the data.
%%% =  ( 2): y = y-data to deseasonalize.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): yO = Deseasonalized data.
%%% =======================================================================

function [ yO ] = DeseasonalizeData(t,y,fDays)

%%% Step 0: Get the seasonal indices and specify filter time period
Sind  = datevec(t);
Sind  = Sind(:,2);

%%% Do two passes (second estimate will improve estimate of trend)
yO = y;
for k = 1:2
    % Step 1: Obtain a first estimate of the trend component
    [tT, yT] = BlockAverage(t,yO,ones(size(t)),fDays);
    % Step 2: Detrend the original series
    yDetrend = yO - interp1(tT,yT,t);
    % Step 3: Apply a stable seasonal filter
    sst = nan(size(t));
    for i = 1:length(t)
        ind    = Sind == Sind(i);
        sst(i) = nanmean(yDetrend(ind));
    end
    Sbar      = nanmean(sst);
    SeasCycle = sst - Sbar;
    % Step 4: Deseasonalize the time series
    yO = yO - SeasCycle;
end

%%% Just get the trend component without the tails
[tT, yT] = BlockAverage(t,yO,ones(size(t)),fDays,'NoTail');
yO       = interp1(tT,yT,t);

end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================