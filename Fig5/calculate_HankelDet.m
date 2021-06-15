function [d] = calculate_HankelDet(m_Beta,i)
% Calculates the modified Hankel determinant as in "DistApproxandModellingViaOrthogonalPolynomialSequences", Provost (2015)
M = zeros(i,i);
for r=0:i-1
    M(r+1,:) = m_Beta(r+1:r+i);
end
d = det(M);
end