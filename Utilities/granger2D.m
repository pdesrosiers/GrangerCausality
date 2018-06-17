function [ coefG ] = granger2D(X, Y, modelOrder)
% This function computes an auto-regressive model for two time series with 
% the least-square method, and then infers the Granger's causality X->Y
%
% Input:  X                -> matrix of size 1 x n 
%                             first time series with n time steps
%         Y                -> matrix of size 1 x n 
%                             a time series with n time steps
%         modelOrder       -> A positive integer 
%                             number of time steps used in regression

% 
% Output: coefG            -> real number
%                             Granger's causality coefficient for X->Y
% References :
%   The MVGC multivariate Granger causality toolbox: a new approach to Granger-causal inference.
%   L. Barnett, A. K. Seth (2014)
%   
%   Comparative performance evaluation of data-driven causality measures applied to brain networks.
%   A. Fasoula, Y. Attal, D. Schwartz (2013)
%
%
% Authors : Maxime Tremblay, Patrick Desrosiers
% Last update : Seotember 2015
    
    [~,nx] = size(X);
    [~,ny] = size(Y);
    n = min(nx,ny);

    % Auto-regression of Y with itself
    
    [~,E] = linearAutoRegression(Y, modelOrder);
    
    
    % X is now divided into modelOrder new variables 
    
    newX = zeros(modelOrder,n-modelOrder);
    
    for i = 1:modelOrder
        newX(i,:) = X(:,i:n-modelOrder+i-1);
    end
    
    
    % Now linear regression of E with new data from X
    
    [~, newE] = linearRegression(newX,E);
    
    % Granger's causality coefficient
    
    coefG = log( var(E)/var(newE) );
    
    