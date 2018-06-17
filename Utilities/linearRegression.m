function [B,E] = linearRegression(X,Y)
% We want to linearly modelize k random variables (Y_1,..,Y_k)
% with m other variables (X_1,...X_m) by using n instances of these
% variables.
%
% We solve for B:    Y = B * X, B is of size (k x m)
%     
% Input :  
%          X    -> data matrix of size (m x n)
%                   m variables, n sample points for each variable
%          Y    -> data matrix of size (k x n)
%                   k variables, n sample points for each variable                   
%
% Output : 
%          B    -> matrix of size (k x m) where the element 
%                  i,j in matrix k is the regression coefficient describing
%                  the influence of X_j on Y_i
%                
%          E    -> matrix of size (k x n) 
%                  error defined by E = Y - BX
%
%  Example:
%           X=rand(3,1000);
%           B=[1,2,3;-4,-2,0];
%           Y = B*X+0.2*rand(2,1000);
%           [BB,E] = linearRegression(X,Y);
%           [B,BB]
%
% Author : Maxime Tremblay, Patrick Desrosiers
% Last update : September 2015

    if nargin < 2 
        error('At least two arguments. See help for info');
    end
        
    if ~ismatrix(X)||~ismatrix(X)
        error('inputs arguments must be matrices');
    end
       
    [~,nX] = size(X);
    [~,nY] = size(Y);
    n = min(nX,nY);
    
    newX = X(:,1:n);
    newY = Y(:,1:n);
    

    % Solve Y = B*K (least square solution) with the pseudo-inverse of X 
    B = newY/newX;
    
    %Compute the error E(t) 
    E = (newY - B*newX);
        
end

