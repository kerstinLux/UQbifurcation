function [PDF_approx] = approxPDFviaMoments(PDF_method,para,N_trunc,moms)
% Approximates according to the chosen PDF_method the PDF of the chosen
% distribution (dist_Choice) with corresponding parameters para based on
% the sequence of moments moms of length N_trunc+1
switch PDF_method
    case 'Legendre'
        %% Calculate moment-based density approximation according to Provost(2005) (same normalization Legendre polys as Matlab)
        syms x
        aux = 0;
        for k=0:N_trunc
            c = flip(sym2poly(legendreP(k,(para.supp(1)+para.supp(2)-2*x)/(para.supp(1)-para.supp(2)))));
            aux = aux + ((2*k+1)/(para.supp(2)-para.supp(1))*c*moms(1:length(c)))*legendreP(k,(para.supp(1)+para.supp(2)-2*x)/(para.supp(1)-para.supp(2)));
        end
        PDF_approx = matlabFunction(aux);
    case 'monic'
        %% Calculate monic orthogonal polynomial density approximation as in "DistApproxandModellingViaOrthogonalPolynomialSequences", Provost et al. (2015)
        syms x
        a = para.supp(1);
        b = para.supp(2);
        %% Calculate alpha,beta from B~Beta(alpha,beta) s.th. c_w*E[B]=E[tN], c_w*E[B^2]=E[tN^2], where c_w = 1/integral(w,a,b) is a normalizing constant with
        alpha = (moms(2)-a)*(moms(3)+a*b - (a+b)*moms(2))/((b-a)*(moms(2)^2-moms(3)));
        eta = (b-moms(2))*(moms(3)+a*b - (a+b)*moms(2))/((b-a)*(moms(2)^2-moms(3)));
        w = @(x) (x-a).^(alpha-1).*(b-x).^(eta-1);
        c_w = 1/integral(w,a,b); % normalizing constant
        %% Calculate monic orthogonal polynomial coefficients
        % Calculate moments of correponding Beta distribution
        paraMonic = para;
        paraMonic.alpha = alpha;
        paraMonic.beta = eta;
        moms_Beta = calcMoms('Beta', paraMonic, 2*N_trunc-1);
        m_Beta = 1/c_w*moms_Beta;
        % Calculate Hankel determinant
        D_Hankel = zeros(N_trunc+1,1);
        D_Hankel(1,1) = 1;
        for i=1:N_trunc
            D_Hankel(i+1,1) = calculate_HankelDet(m_Beta,i);
        end
        % Calculate modified Hankel determinant
        D_Hankelmod = zeros(N_trunc+1,N_trunc);
        D_Hankelmod(1,1) = 0;
        D_Hankelmod(2,1) = 1;
        for i=2:N_trunc
            for k=0:i-1
                D_Hankelmod(i+1,k+1) = calculate_HankelDetmod(m_Beta,i,k);
            end
        end
        %% Calculate corresponding monic orthogonal polynomials
        d_monic = zeros(N_trunc+1,N_trunc+1);
        for i=3:N_trunc
            for k=0:i-1
                d_monic(i+1,k+1) = (-1)^(i-k)*D_Hankelmod(i+1,k+1)/D_Hankel(i+1,1);
            end
            d_monic(i+1,i+1) = 1;
        end

        %% Calculate coefficients lambda_i in density approximation (14) in Provost et al. (2015)
        lambda = zeros(1,N_trunc+1);
        aux = 0;
        for i=3:N_trunc
            integrand = @(x) w(x).*(polyval(flip(d_monic(i+1,1:i+1)),x)).^2;
            intVal = integral(integrand,a,b);
            for k=0:i
                aux = aux + d_monic(i+1,k+1)*moms(k+1);
            end
            lambda(i+1) = aux/intVal;
            aux = 0;
        end
        %% Calculate density approximation according to formula (14) in Provost et al. (2015)
        aux = @(x) 0;
        for i=3:N_trunc
            aux = @(x) aux(x) + lambda(i+1)*w(x).*polyval(flip(d_monic(i+1,1:i+1)),x);
        end
        PDF_approx = @(x) c_w*w(x).*(x>a & x<b) + aux(x).*(x>a & x<b);
        
    case 'transformedMoments'
        %% Calculate moment-based PDF approximation according to "Hausdorff moment problem: Reconstruction of PDFs", Mnatsakanov(2008b), formula (6)
        PDF_approx = @(x) calculateDens_trafoMomSequ(moms,para.supp,N_trunc,x);
end
end