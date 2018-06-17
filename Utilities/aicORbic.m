function [ order ] = aicORbic(L,k,n,criterion)
% Get the AIC and BIC for a log-likelihood relation.
%
% Input : L    -> Log-likelihood
%         k    -> Number of free parameters
%         n    -> Number of observations
%         criterion -> Choose between 'aic' and 'bic' (default)
%
% Output  aic  -> Aikaike information criterion
%         bic  -> Bayesian information criterion
%
% Author : Maxime Tremblay, Patrick Desrosiers
% Last update : September 2015

    if nargin < 4
        criterion = 'bic';
    end
    
    if strcmp(criterion,'bic')
        order = -2*L + k*log(n);
        
    elseif strcmp(criterion,'aic')
        % Cannot compute AIC if number of free parameters is greater or 
        % equal to number of observations
        if k >= n
            order = nan;
        else
            order = -2*L + 2*k*n/(n-k-1);
        end
    else
        order = nan;
        fprinf('Invalid criterion');
    end
    

end

