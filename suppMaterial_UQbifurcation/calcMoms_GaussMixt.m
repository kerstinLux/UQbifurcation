function [moms_GaussMixt] = calcMoms_GaussMixt(pars,n_Gauss,num_moms)
% Calculates the k-th moment of a Gaussian mixture distribution with n_Gauss components specified
% by pars=[weights mean sigma], where weights is a 1xn_Gauss-1 vector, and
% meand and sigma are 1xn_Gauss vectors
mu = pars(n_Gauss:2*n_Gauss-1);
sigma = pars(2*n_Gauss:end);
compMoms_GMM = ones(num_moms,n_Gauss);
compMoms_GMM(2,:) = mu;
compMoms_GMM(3,:) = sigma.^2 + mu.^2;
for n=1:n_Gauss
    for i=3:num_moms
        compMoms_GMM(i+1,n) = calcMomGaussian(compMoms_GMM(2,n),sigma(n),i);
    end
end

weights = [pars(1:n_Gauss-1) 1-sum(pars(1:n_Gauss-1))];
mean_GMM = weights * compMoms_GMM(2,:)';
mom2_GMM = weights * compMoms_GMM(3,:)';
moms_GaussMixt = ones(num_moms+1,1);
moms_GaussMixt(2) = mean_GMM;
moms_GaussMixt(3) = mom2_GMM;
for i=3:num_moms
    moms_GaussMixt(i+1) = weights * compMoms_GMM(i+1,:)';
end
% avoid zero-th moment in optimization
moms_GaussMixt = moms_GaussMixt(2:end);
end

