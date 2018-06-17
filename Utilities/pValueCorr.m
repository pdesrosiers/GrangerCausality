function [ pval ] = pValueCorr(r,n)
% Compute the p-values for all correlation coefficients and n the number of data use to compute r.
% 
%     input: r ==> float, the correlation coefficient
%            n==> int, the number of data used when computing r
% 
%     output: pval ==> m*m matrix with floats, each element is a p-value  
%
% We suppose that the data used in computing r are 
%     distributed like a binomial law. we use:
% 
%         p-value= 2*P( T_{n-2} >= |t| )
% 
%     For fast computation we use : 
% 
%         t^2 = r*r * ((n-2) / ((1.0 - r) * (1.0 + r)))
%         p-value = beta_inc(0.5*(n-2), 0.5, (n-2) / ((n-2) + t^2)) , where
% beta_inc is the incomplete beta function ==> \Beta(x;\,a,b) = \int_0^x t^{a-1}*(1-t)^{b-1} dt     
% 
% 
    % Degrees of freedom
    df = (n-2)*ones(size(r));
    % t-squared variable
    tSquared = (r.*r).*(df./(ones(size(r))-r.*r)) ;
    % new variable
    new = rdivide(df,df + tSquared);
    new(new<0)=0;
    new(new>1)=1;
    % pvalue
    pval = betainc(new, 0.5*df, 0.5*ones(size(r)));
end