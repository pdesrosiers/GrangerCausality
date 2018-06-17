function [B,E] = linearAutoRegression(X, l)
% Solve for B :    X(t) = B * K(t)
%                  K(t) = [X(t-1)  X(t-2)  X(t-3) ... X(t-l)]^T   
%     
% Input :  
%          X    -> matrix in which each row is random variable (a time 
%                  series) and each column sample (a time step)
%
%          l    -> Maximum lag
%
%
% Output : 
%          B    -> array of l matrices of size (m x m) where the element 
%                  i,j in matrix k is the regression coefficient from  
%                  j to i at with lag k
%                
%          E    -> error defined by E(t) = ( X(t) - B*K(t) ).
%                  E is (m x (n-l)) and the element i,j is the error of node
%                  i at time j.
%
%  Example:
%           X=rand(3,100);
%           X(2,3:end)=2*X(1,2:end-1)-X(1,1:end-2);
%           X(3,3:end)=-3*X(1,1:end-2);
%           [B,E] = linearAutoRegression(X,2);
%
% Author : Maxime Tremblay, Patrick Desrosiers
% Last update : September 2015

    if nargin < 2 
        error('At least two input parameters are needed. See help for info');
    end   
    if ~ismatrix(X)
        error('X must be at most bidimensional arrays');
    end
       
    [m,n] = size(X);
    
    % Compute Y = X for t = l+1 : n
    Y = X(:,l+1:n);
    
    % Get the index where to take value in X to obtain K
    index = repmat(l+1:n,l,1) - repmat((1:l)',1,n-l);
    
    % Get the value in X or spk from index
    K = reshape(permute(reshape(X(:,index)',l,n-l,m),[1 3 2]),m*l,n-l);

    % We can now solve Y = B*K (least square solution)
    B = Y/K;
    
    %Compute the error E(t) 
    E = (Y - B*K);
    
    % Reshape B
    B = permute(reshape(B',l,m,m) ,[3 2 1]);
        
end

