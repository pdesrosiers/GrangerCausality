function [ connectMat,modelOrder ] = granger(X, maxModelOrder, ...
                                                  infoCriterion, normalize)
% This function computes an auto-regressive model for the time series with 
% the least-square method. The best order of the auto-regression is 
% selected by minimising an information criterion.  
%
% Input:  X                -> matrix of size m*n (m variables,n samples)
%                             each row is a different time series
%                             each column is a time step
%         maxModelOrder    -> A positive integer for the maximum model
%                             order of the MVAR (default 8)
%         infoCriterion    -> Information criterion used to 
%                             select model order. Either 'aic'(default)
%                             or 'bic' 
%         normalize        -> Boolean.  If yes max connectivity is 1.
% 
% Output: connectMat        -> m*m connectivity matrix
%                             values based on Granger causaliy
%         modelOrder        ->int, order choose for the linear regression
% Example:
%         X=rand(3,100);
%         X(2,3:end)=2*X(1,2:end-1)-X(1,1:end-2);
%         X(3,3:end)=-3*X(1,1:end-2);
%         X = X + 0.2*randn(3,100);
%         M = granger(X,10)
%
% References :
%   The MVGC multivariate Granger causality toolbox: a new approach to Granger-causal inference.
%   L. Barnett, A. K. Seth (2014)
%   
%   Comparative performance evaluation of data-driven causality measures applied to brain networks.
%   A. Fasoula, Y. Attal, D. Schwartz (2013)
%
%
% Authors : Maxime Tremblay, Patrick Desrosiers
% Last update : September 2015
    
    [m,n] = size(X);

    % Set default values
    if nargin < 4, normalize = false; end
    if nargin < 3, infoCriterion = 'aic'; end
    if nargin < 2, maxModelOrder = 8; end
    
    % Verify error
    
    if maxModelOrder < 1
        error('maxModelOrder must be a positive integer');
    end
    
    if and(~strcmp(infoCriterion,'aic'),~strcmp(infoCriterion,'bic'))
        error('Invalid information criterion');
    end
    
    if maxModelOrder >= n
        maxModelOrder = n-1;
    end
    
    
    %Find the best model order the information criterion
    
    modelOrder = nan(1,maxModelOrder);
    
    for l=1:maxModelOrder
        % Complete a least-square AR of order l and get the error
        [~,E] = linearAutoRegression(X,l);
        % Compute the estimate residual covariance matrix
        sigma = cov(E');
        % Compute the log-likelihood of estimate value for order l
        loglike = (l-n) * log(det(sigma)) / 2;
        % Compute AIC or BIC for order l
        modelOrder(l) = aicORbic(loglike,m*m*l,n,infoCriterion); 
    end
    
    
    [~,modelOrder] = min(modelOrder); % Best model order
    
    textAboutOrder = [newline, ...
        'Best regression order found: ',num2str(modelOrder)];
    disp(textAboutOrder);
 
    % Connectivity Matrix
    connectMat = zeros(m,m);
    
    for i=1:m
        for j=1:m
            if i~=j
                connectMat(i,j) = granger2D(X(i,:),X(j,:),modelOrder);
            end
        end
    end
    
    connectMat(connectMat<0) = 0;
    
    if normalize
        maxVal = max(max(connectMat));
        connectMat = connectMat./maxVal;
    end
    
end