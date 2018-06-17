function [ pval ] = pValueGranger(g,tau,n)
% Compute the p-values for all granger coefficients in an array
%
%     input: g  ==> array whose elements are the Granger coefficients
%            tau==> positive int, number of time steps used in regression
%            n  ==> int, the number of data used when computing g
% 
%     output: pval ==> matrix with floats, each element is a p-value  

    %Parameters
    
    m=size(g,1);
    nx=1;
    ny=1;
    %nz=m-nx-ny;% We consider that every other nodes may have an effect on 
                % the Granger coef even if we do the pairewise method
    nz=0;
    d1=tau*ny;
    d2=n-tau*(1+ny+nz);
    if d2<0
        disp('PROB');
        warning(['The number of time steps is insufficient.', ...
            'The p-values are calculated as if the sample size were 1e5']);
        d2=1e5-tau*(1+ny+nz);
    end
    
    % Evaluate the pValue 
    xx = exp(g)-ones(size(g));
    pval = 1-fcdf((d2/d1)*xx,d1,d2);
end