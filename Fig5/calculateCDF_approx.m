function [cdf_approx] = calculateCDF_approx(moms_tN,supp,N_trunc,x)
%% Calculation of the cdf approximation according to "Hausdorff moment problem: Reconstruction of distributions", Mnatsakanov(2008), formula (2)
doubleSum = zeros(1,N_trunc+1);
for k=0:N_trunc
    outerAux = 0;
for m=0:k
    aux = 0;
    for j=m:N_trunc
        aux = aux + nchoosek(N_trunc,j)*nchoosek(j,m)*(-1)^(j-m)/supp(2)^j*moms_tN(j+1);
    end
    outerAux = outerAux + aux;
end
    doubleSum(k+1) = outerAux;
end
% Determine upper summation bound from given x
upSumInd = floor(N_trunc*x/supp(2));
cdf_approx = doubleSum(upSumInd+1);


% cdf_approx = 0;
% for k=0:upSumInd
%     aux = 0;
%     for j=k:N_trunc
%         aux = aux + nchoosek(N_trunc,j)*nchoosek(j,k)*(-1)^(j-k)/supp(2)^j*moms_tN(j+1);
%     end
%     cdf_approx = cdf_approx + aux;
% end
end

