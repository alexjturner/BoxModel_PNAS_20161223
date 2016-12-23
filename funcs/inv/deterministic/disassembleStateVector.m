%%% =======================================================================
%%% = getTime.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Break the state vector into a matrix of emission sources and a
%%% =        vector of ICs.
%%% =  ( 2): "assembleStateVector" will do the opposite.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): x  -- State vector.
%%% =  ( 2): nT -- Number of time steps in our time vector.
%%% =  ( 3): nE -- Number of different emission sources.
%%% =  ( 4): nI -- Number of ICs.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): ems -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC  -- Vector with the initial conditions.
%%% =======================================================================

function [ ems, IC ] = disassembleStateVector( x, nT, nE, nI )
% Construct the sources and IC vectors
ems = nan(nT,nE);
IC  = nan(1,nI);
ii  = 1;
for i = 1:nE
    for j = 1:nT
        ems(j,i) = x(ii);
        ii       = ii + 1;
    end
end
IC(:) = x(ii:end);
end


%%% =======================================================================
%%% = END
%%% =======================================================================
