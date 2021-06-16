function [d] = calculate_HankelDet(m_Beta,i)
% Calculates the Hankel determinant as in "DistApproxandModellingViaOrthogonalPolynomialSequences", Provost et al. (2015)
M = zeros(i,i);
for r=0:i-1
    M(r+1,:) = m_Beta(r+1:r+i);
end
d = det(M);
end