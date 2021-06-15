%% Inverse Mellin transform
clear; clc;
%% Choose distribution to approximate
dist_Choice = 'Beta';
%% Define parameters according to chosen distribution
para.supp = [0 1];
switch dist_Choice
    case 'Uniform'
    case 'truncatedNormal'
        para.meanN = 0.65;
        para.stdN = 0.2;
    case 'Beta'
        % Choose first two moments to match up to normalizing constant
        para.mean = 0.6;
        para.mom2 = 0.4;
        % Calculate alpha,beta from B~gBeta(alpha,beta) s.th. c_w*E[B]=para.mean,
        % c_w*E[B^2]=para.mom2, where c_w is a normalizing constant w = @(x) x.^(alpha-1).*(1-x).^(beta-1);
        para.alpha = (para.mean-para.supp(1))*(para.mom2+para.supp(1)*para.supp(2) - (para.supp(1)+para.supp(2))*para.mean)/((para.supp(2)-para.supp(1))*(para.mean^2-para.mom2));
        para.beta = (para.supp(2)-para.mean)*(para.mom2+para.supp(1)*para.supp(2) - (para.supp(1)+para.supp(2))*para.mean)/((para.supp(2)-para.supp(1))*(para.mean^2-para.mom2));
end

%% Calculate moments of chosen distribution up to order N_trunc
N_trunc = 4;
moms = calcMoms(dist_Choice, para, N_trunc);

%% Choose method to approximate probability density function
PDF_Legendre = approxPDFviaMoments('Legendre',para,N_trunc,moms);
PDF_monic = approxPDFviaMoments('monic',para,N_trunc,moms);
PDF_trafoMom = approxPDFviaMoments('transformedMoments',para,N_trunc,moms);


%% Plot density and density approximations
fig1 = figure(1);
hold on
plot(para.supp(1):0.01:para.supp(2),PDF_monic(para.supp(1):0.01:para.supp(2)),'m*','MarkerSize',6)
plot(para.supp(1):0.01:para.supp(2),PDF_trafoMom(para.supp(1):0.01:para.supp(2)),'k--','LineWidth',2)
switch dist_Choice
    case 'Uniform'
        plot(para.supp(1):0.01:para.supp(2),ones(size(para.supp(1):0.01:para.supp(2)))*PDF_Legendre(),'ro','MarkerSize',4,'LineWidth',1);
        plot(para.supp(1):10^(-3):para.supp(2),pdf('Uniform',para.supp(1):10^(-3):para.supp(2),para.supp(1),para.supp(2)),'b','LineWidth',1)
        leg = legend('monic orthogonal approx truncated','transformed moments approx','Legendre approx truncated',strcat(dist_Choice,' density'));
        leg.Location = 'Southeast';
        leg.Interpreter = 'latex';
        leg.FontSize= 12;
        xlabel('Support of distribution');
        ylim([0 1.2])
%         % Save figure(1)
%         savefig(fig1,strcat('PDF_approx_',dist_Choice,'supp',num2str(para.supp(1)),num2str(para.supp(2)),'N',num2str(N_trunc),'.fig'));
%         saveas(fig1,strcat('PDF_approx_',dist_Choice,'supp',num2str(para.supp(1)),num2str(para.supp(2)),'N',num2str(N_trunc),'.eps'),'epsc');
    case 'truncatedNormal'
        plot(para.supp(1):0.01:para.supp(2),PDF_Legendre(para.supp(1):0.01:para.supp(2)),'ro','MarkerSize',4,'LineWidth',1);
        pd = makedist('Normal','mu',para.meanN,'sigma',para.stdN);
        pd_tN = truncate(pd,para.supp(1),para.supp(2));
        plot(para.supp(1):10^(-3):para.supp(2),pdf(pd_tN,para.supp(1):10^(-3):para.supp(2)),'b','LineWidth',2)
        leg = legend('monic orthogonal approx truncated','transformed moments approx','Legendre approx truncated',strcat(dist_Choice,' density'));
        leg.Location = 'Southeast';
        xlabel('Support of distribution');
        leg.Interpreter = 'latex';
        leg.FontSize= 12;
        ylim([-0.5 4])
%         % Save figure(1)
%         savefig(fig1,strcat('PDF_approx_',dist_Choice,'supp',num2str(para.supp(1)),num2str(para.supp(2)),'mu',strrep(num2str(para.meanN),'.','K'),'std',strrep(num2str(para.stdN),'.','K'),'N',num2str(N_trunc),'.fig'));
%         saveas(fig1,strcat('PDF_approx_',dist_Choice,'supp',num2str(para.supp(1)),num2str(para.supp(2)),'mu',strrep(num2str(para.meanN),'.','K'),'std',strrep(num2str(para.stdN),'.','K'),'N',num2str(N_trunc),'.eps'),'epsc');
    case 'Beta'
        plot(para.supp(1):0.01:para.supp(2),PDF_Legendre(para.supp(1):0.01:para.supp(2)),'ro','MarkerSize',4,'LineWidth',1);
        plot(para.supp(1):10^(-3):para.supp(2),pdf('Beta',para.supp(1):10^(-3):para.supp(2),para.alpha,para.beta),'b','LineWidth',2)
        leg = legend('monic orthogonal approximation','transformed moments approximation','Legendre approximation',strcat(dist_Choice,' density'));
        leg.Location = 'Southeast';
        xlabel('Support of distribution');
        leg.Interpreter = 'latex';
        leg.FontSize= 12;
        ylim([-0.5 2.5])
%         % Save figure(1)
%         savefig(fig1,strcat('PDF_approx_',dist_Choice,'supp',num2str(para.supp(1)),num2str(para.supp(2)),'alpha',strrep(num2str(para.alpha),'.','K'),'beta',strrep(num2str(para.beta),'.','K'),'N',num2str(N_trunc),'.fig'));
%         saveas(fig1,strcat('PDF_approx_',dist_Choice,'supp',num2str(para.supp(1)),num2str(para.supp(2)),'alpha',strrep(num2str(para.alpha),'.','K'),'beta',strrep(num2str(para.beta),'.','K'),'N',num2str(N_trunc),'.eps'),'epsc');
end