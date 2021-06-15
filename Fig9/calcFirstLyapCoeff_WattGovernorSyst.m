function [lyap_coeff] = calcFirstLyapCoeff_WattGovernorSyst(r)
% Calculates the first Lyapunov coefficient of the Watt governor system according
% to "Bifurcation anallysis of the Watt governor system", Sotomayor et al. (2007), eq. (60)
% lyap_coeff = -0.5*((r(1,:).*r(2,:).^1.5.*(1-r(2,:).^2).*(3+(r(1,:).^2-5).*r(2,:).^2+r(1,:).^4.*r(2,:).^6))./((1-r(2,:).^2+r(1,:).^2.*r(2,:).^4).*(1-r(2,:).^2+4*r(1,:).^2.*r(2,:).^4)));
lyap_coeff = -(3+(r(1,:).^2-5).*r(2,:).^2+r(1,:).^4.*r(2,:).^6);
if abs(3+(r(1,:).^2-5).*r(2,:).^2+r(1,:).^4.*r(2,:).^6) < eps
    disp('Lyapunov coefficient equal to zero');
end
end

