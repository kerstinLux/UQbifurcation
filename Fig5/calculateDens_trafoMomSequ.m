function [dens_approx] = calculateDens_trafoMomSequ(moms_tN,supp,N_trunc,x)
%% Calculation of the pdf approximation according to "Hausdorff moment problem: Reconstruction of PDFs", Mnatsakanov(2008b), formula (6)
sumInDens = zeros(1,N_trunc+1);
for k=0:N_trunc
    aux = 0;
    for m=0:N_trunc-k
        aux = aux + (-1/supp(2))^m*moms_tN(m+k+1)/(factorial(m)*factorial(N_trunc-k-m));
    end
    sumInDens(k+1) = aux;
end
% Determine upper summation bound for the given x
upSumInd = floor(N_trunc*x/supp(2));

dens_approx = gamma(N_trunc+2)*gamma(upSumInd+1).^(-1).*(supp(2).^(upSumInd+1)).^(-1).*sumInDens(upSumInd+1);
% % Determine upper summation bound from given x
% upSumInd = floor(N_trunc*x/supp(2));
% dens_approx = 0;
% for m=0:N_trunc-upSumInd
%     dens_approx = dens_approx + (-1/supp(2))^m*moms_tN(m+upSumInd+1)/(factorial(m)*factorial(N_trunc-upSumInd-m));
% end
% dens_approx = gamma(N_trunc+2)/gamma(upSumInd+1)*1/(supp(2)^(upSumInd+1))*dens_approx;


