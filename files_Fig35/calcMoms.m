function [moms] = calcMoms(dist_Choice, para, N_trunc)
% Calculate the first N_trunc+1 moments of the chosen distribution
% (dist_Choice) with corresponding parameters para
switch dist_Choice
    case 'Uniform'
        %% test random variable X~U(para.supp(1),para.supp(2))
        moms = zeros(N_trunc+1,1);
        moms(1) = 1; % Set moment E[X^0]=1
        for i=1:N_trunc
            moms(i+1) = (para.supp(2)^(i+1)-para.supp(1)^(i+1))/(i+1); % Mellin transform of U~\UU(a,b)
        end
    case 'truncatedNormal'
        %% test random variable X~N_{para.supp(1)}^{para.supp(2)}(para.meanN,para.stdN^2)
        alpha = (para.supp(1)-para.meanN)/para.stdN;
        beta = (para.supp(2)-para.meanN)/para.stdN;
        moms = zeros(N_trunc,1);
        EW_tN = para.meanN+(normpdf((para.supp(1)-para.meanN)/para.stdN)-normpdf((para.supp(2)-para.meanN)/para.stdN))/(normcdf((para.supp(2)-para.meanN)/para.stdN)-normcdf((para.supp(1)-para.meanN)/para.stdN))*para.stdN;
        Var_tN = para.stdN^2*(1+ (alpha*normpdf(alpha)-beta*normpdf(beta))/(normcdf(beta)-normcdf(alpha)) -((normpdf(alpha)-normpdf(beta))/(normcdf(beta)-normcdf(alpha)))^2);
        Mom2_tN = Var_tN + EW_tN^2;
        
        moms(1) = 1; % Set moment E[X^0]=1
        moms(2) = EW_tN;
        moms(3) = Mom2_tN;
        % Compute moments of truncated Gaussian according to https://people.smp.uq.edu.au/YoniNazarathy/teaching_projects/studentWork/EricOrjebin_TruncatedNormalMoments.pdf
        for i=3:N_trunc
            moms(i+1) = (i-1)*para.stdN^2*moms(i+1-2) + para.meanN*moms(i+1-1) - para.stdN*(para.supp(2)^(i-1)*normpdf(beta)-para.supp(1)^(i-1)*normpdf(alpha))/(normcdf(beta)-normcdf(alpha));
        end
    case 'Beta'
        %% test random variable X~Beta(alpha,beta)
        moms = ones(N_trunc+1,1);
        moms(1) = 1; % Set moment E[X^0]=1
        for i=1:N_trunc
            aux = 0;
            for j=0:i
                aux = aux + nchoosek(i,j)*(para.supp(1)^(i-j)*(para.supp(2)-para.supp(1))^j*gamma(para.alpha+para.beta)*gamma(para.alpha+j))/(gamma(para.alpha+para.beta+j)*gamma(para.alpha));
            end
            moms(i+1) = aux;
        end
end

