%%% =======================================================================
%%% = assembleStateVector.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Transform the emission sources and ICs from matrice/vectors to
%%% =        a state vector for the inversion.
%%% =  ( 2): "disassembleStateVector" will do the opposite.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): ems -- Matrix with the emission sources (and OH).
%%% =  ( 2): IC  -- Vector with the initial conditions.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =  ( 1): x -- State vector.
%%% =======================================================================

function [ x ] = assembleStateVector( ems, IC )
% Get sizes
nT = size(ems,1);
nE = size(ems,2);
nI = length(IC);
% Construct the state vector
x_ems = nan(nT*nE,1);
x_IC  = nan(nI,1);
% Sources
ii = 1;
for i = 1:nE
    for j = 1:nT
        x_ems(ii) = ems(j,i);
        ii        = ii + 1;
    end
end
% ICs
x_IC(:) = IC(:);
% Full state vector
x = [x_ems;x_IC];
end


%%% =======================================================================
%%% = END
%%% =======================================================================
