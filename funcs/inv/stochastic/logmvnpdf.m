%%% =======================================================================
%%% = logmvnpdf.m
%%% = Alex Turner
%%% = Downloaded on: 03/30/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Computes a log-likelihood array for observations x  where 
%%% =        x_n ~ N(mu,Sigma).
%%% =  ( 2): x is NxD, mu is 1xD, Sigma is DxD.
%%% =  ( 3): From: "http://www.mathworks.com/matlabcentral/fileexchange/...
%%% =               34064-log-multivariate-normal-distribution-function/...
%%% =               content/logmvnpdf.m".
%%% =  ( 4): Author:  Benjamin Dichter.
%%% =======================================================================

function [logp] = logmvnpdf(x,mu,Sigma)

[N,D] = size(x);
const = -0.5 * D * log(2*pi);

xc = bsxfun(@minus,x,mu)';

term1 = -0.5 * sum((xc / Sigma) .* xc, 2);  % N x 1
term2 = const - 0.5 * logdet(Sigma);        % scalar
logp = term1' + term2;

end

function y = logdet(A)

U = chol(A);
y = 2*sum(log(diag(U)));

end


%%% =======================================================================
%%% =                              E N D                                  =
%%% =======================================================================