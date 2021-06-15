function [binCoeff] = calcBinCoeff(N,s,M,i)
% based on the power Matrix M and the requested Moment s-1, the binomial
% coefficient for the i-th combination of powers in M is calculated
binCoeff = nchoosek(s-1,M(i,1));
n=s-1;
for j=2:N
    n = n-M(i,j-1);
    binCoeff = binCoeff*nchoosek(n,M(i,j));
end
end