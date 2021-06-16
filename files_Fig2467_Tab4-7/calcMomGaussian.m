function [mom_k] = calcMomGaussian(mu,sigma,k)
% Calculates the k-th moment E[x^k] of a random variable X~N(mu,sigma^2)
if mod(k,2)==0 % k is even
    mom_k = sigma^k*2^(k/2)*gamma((k+1)/2)/sqrt(pi)*hypergeom(-k/2,0.5,-mu^2/(2*sigma^2));
else
    mom_k = mu*sigma^(k-1)*2^((k+1)/2)*gamma(k/2+1)/sqrt(pi)*hypergeom((1-k)/2,1.5,-mu^2/(2*sigma^2));
end

end

