function [ network ] = inferNetworkWithGranger( timeSeries, parameters)
% Infers a network from times series by using Granger causality.  Relies on
% file granger.m.
%
% INPUT
%   timeSeries  matrix of dim (p,n)
%               p is the number of variables (number of nodes or cells)
%               n is the number of samples (number of time steps)
%   parameters  struct whose fields are the values of the parameters
%               field names are: 
%                   alpha:  the significance level;
%                   maxLag: non neg integer, number of time steps used 
%                           in the past to compute the cross-correlations
% OUPUT
%   network     Object of class graph or digraph
%
% EXAMPLE
%   X = sin(12*rand(4,1000));
%   X(2,2:end) = 2*X(1,1:end-1)+0.1*randn(1,999); 
%   X(3,3:end) = -X(2,2:end-1) + X(1,1:end-2);
%   G = inferNetworkWithGranger(X);
%   plot(G);
%   disp(G.Edges)
%
% Patrick Desrosiers, DCClab.  Created 2016-09.  Modified 2017-07

%% Check input
    if ~ismatrix(timeSeries)
        error('The 1st argument must be a matrix')
    elseif nargin<2
        parameters = struct;
        parameters.alpha = 0.01;
        parameters.maxLag = 3;
    end
%% Set global variables
    alpha = parameters.alpha;
    maxLag = parameters.maxLag;
    [~, nbTimeSteps] = size(timeSeries);
    X = double(timeSeries);
    
%% Compute Granger coeffs
    [connectMatGranger,modelOrder] = granger(X, maxLag);

%% Keep only significant correlation
    pValues = pValueGranger(connectMatGranger,modelOrder,nbTimeSteps);
    connectMatGranger = connectMatGranger.*(pValues<alpha);

%% Check for significant Granger coeffs
    maxGrangerCoeff = max(max(connectMatGranger));
    
    if maxGrangerCoeff>0
        % Granger analysis ok
        connectMatGranger = connectMatGranger/maxGrangerCoeff;
    else
        % Granger analysis has failed.  Parameters are changed
        minPValue = min(min(pValues));
        warning(['No significant Granger coefficient detected.  '....
                'The significance level should be at least ', ...
                num2str(minPValue),...
                '.  We now assume that the number of time steps is 1e4']);
        pValues = pValueGranger(connectMatGranger,modelOrder,1e4);
        minPValue = min(min(pValues));
        disp(['p-value min: ', num2str(minPValue)]);
        connectMatGranger = connectMatGranger.*(pValues<alpha);
    end
    
%% Remove the diagonal elelements    
    connectMatGranger = connectMatGranger-diag(diag(connectMatGranger));

%% Return a digraph object
    network = digraph(double(connectMatGranger));
end

