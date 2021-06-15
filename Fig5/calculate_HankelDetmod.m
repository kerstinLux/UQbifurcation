function [d] = calculate_HankelDetmod(m_Beta,i,k)
% Calculates the modified Hankel determinant as in "DistApproxandModellingViaOrthogonalPolynomialSequences", Provost (2015)
M = zeros(i,i);
for r=0:i-1
    M(r+1,:) = [m_Beta(r+1:k-1+r+1)' m_Beta(k+1+r+1:i+r+1)'];
end
d = det(M);
end

