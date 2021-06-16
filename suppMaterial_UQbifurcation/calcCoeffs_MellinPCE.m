function [c_MellinPCE] = calcCoeffs_MellinPCE(c,N,s)
% Calculates the coefficients of the Mellin transform of the PCE expansion of a random variable
% X with coefficients coeffs c_n via U~U(0,1) and shifted Legendre polynomials,
% i.e. X = \sum_{n=0}^{N} c_nP_n(U) in terms of coefficients in front of
% monomials

%% flexible implementation for arbitrary value of N_trunc
[numRows,M] = build_PowerMatrix(N,s);
sumM = sum(M,2);
c_MellinPCE = zeros(1,N*(s-1)+1);
for i=1:numRows
    binCoeff = calcBinCoeff(N,s,M,i);
    cProd = c(1)^(s-1-sumM(i));
    for n=1:N
        cProd = cProd*c(n+1)^(M(i,end-n+1));
    end
    u_exp = sum((N:-1:1).*M(i,:));
    c_MellinPCE(u_exp+1) = c_MellinPCE(u_exp+1) + binCoeff*cProd;
end

end
