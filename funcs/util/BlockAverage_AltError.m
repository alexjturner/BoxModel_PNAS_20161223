%%% =======================================================================
%%% = BlockAverage.m
%%% = Alex Turner
%%% = 08/14/2013
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Gives a weighted moving average with an alternate error calc.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): x  = Date in datenum format.
%%% =  ( 2): y  = y-data to filter.
%%% =  ( 3): pI = Pressure.
%%% =  ( 4): nD = Number of days for the filter
%%% =  ( 5): Flag for removing tails (any input = remove tails).
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): xO = Date in datenum format.
%%% =  ( 2): yO = Filtered data.
%%% =  ( 3): E  = Error.
%%% =======================================================================

function [ xO, yO, E ] = BlockAverage_AltError(x,y,p,nD,varargin)
% Initialize matricies
xO = [x(1) : x(end)]'; 
yO = zeros(size(xO));
E  = zeros(size(xO));
% Sort the data
[~,I] = sort(x); x = x(I); y = y(I); p = p(I);
% Loop over the x's
for i = 1:length(xO)
    % Get the indicies in our time period
    ind = ( (xO(i) - nD/2 < x) & (x < xO(i) + nD/2) );
    % Make sure we've got obs
    if sum(ind) < 1
        yO(i) = -9999;
        E(i)  = -9999;
    else
        % Get the weighted average
        tmpy = y(ind);
        tmpp = p(ind);
        w    = tmpp ./ sum(tmpp);
        for j = 1:sum(ind)
            yO(i) = yO(i) + w(j).*tmpy(j);
        end
        E(i) = sqrt( sum( w .* (tmpy - yO(i)) .^ 2 ) );
    end
end
% Remove NaNs
ind = ~isnan(yO); xO = xO(ind); yO = yO(ind); E = E(ind);
% Replace -9999's
yO(yO == -9999) = NaN;
E(E == -9999)   = NaN;
% Set the beginning and end of the filtered data to NaN
if (nargin == 5)
    ind      = ( (xO(1) + nD/2 < xO) & (xO < xO(end) - nD/2) );
    yO(~ind) = NaN;
end
end


%%% =======================================================================
%%% =                             E N D                                   =
%%% =======================================================================