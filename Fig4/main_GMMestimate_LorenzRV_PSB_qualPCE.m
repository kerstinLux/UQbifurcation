%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the framework ( if not already started)
clear; clc;
uqlab
%% Define the model
% Define the model options
modelopts.mFile = 'LorenzCMaux_function';
% Create and add the model to UQLab
myLorenzCM_Model = uq_createModel(modelopts);

%% Specify parameters
% Choose first two moments to match up to normalizing constant
%         para.mean_alpha = 0.4;
%         para.mom2_alpha = 0.2;
%         % Calculate alpha,beta from B~gBeta(alpha,beta) s.th. c_w*E[B]=para.mean,
%         para.alpha1 = para.mean_alpha^2*(1-para.mean_alpha)/(para.mom2_alpha-para.mean_alpha^2)-para.mean_alpha;
%         para.beta1 = (para.mean_alpha*(1-para.mean_alpha)/(para.mom2_alpha-para.mean_alpha^2)-1)*(1-para.mean_alpha);
para.alpha1 = 2;
para.beta1 = 2;
%% Define numInputRV independent random variables
numInputRV = 1;
% supp_RV = [4 6];
% supp_RV = [-0.5 0.5];
supp_RV = [0 1];
for i = 1 : 1
%     IOpts.Marginals(i).Type = 'Uniform' ;
IOpts.Marginals(i).Type = 'Beta' ;
% IOpts.Marginals(i).Type = 'Gaussian' ;
    IOpts.Marginals(i).Parameters = supp_RV ;
IOpts.Marginals(i).Parameters = [para.alpha1 para.beta1 supp_RV] ;
end
myInput = uq_createInput(IOpts);
%% Setup of PCE
MetaOpts.Type = 'uq_metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.PolyTypes = {'Legendre'};

% Specification of PCE degree, Gaussian quadrature based projection
MetaOpts.Degree = 2;
N = MetaOpts.Degree;
MetaOpts.Method = 'quadrature';
% Creation of the metamodel :
myPCE_Quadrature = uq_createModel(MetaOpts);
%% Report basic information
uq_print(myPCE_Quadrature);
% moments
Mom_PCE = myPCE_Quadrature.PCE.Moments;
PCEcoeffs = myPCE_Quadrature.PCE.Coefficients;
%% Graphical representation of spectrum of non-zero coefficients
% uq_display(myPCE_Quadrature);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_moms = 5;

%%%%%%%%%%%%% MC simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plain MC countercheck
M=10^6;
rng default
% % Uniform Input
% a = IOpts.Marginals(i).Parameters(1); b = IOpts.Marginals(i).Parameters(2);
% U = a+(b-a)*rand(M,1);
% Beta Input
% U = betarnd(para.alpha1, para.beta1,[M 1]);
U_test = uq_getSample(M);
U = U_test;
% % Gaussian input
% mu = IOpts.Marginals(ii).Parameters(1); sigma = IOpts.Marginals(ii).Parameters(2);
% U = mu + sigma*randn(M,1);
% QoI = LorenzCMaux_function(U_test);
QoI = LorenzCMaux_function(U);
EW = mean(QoI);
Var = std(QoI).^2;
%% test whether degree of PCE is reasonable
U_test = uq_getSample(M);
QoI_test = uq_evalModel(myPCE_Quadrature,U_test);
%% Plot comparison MC-based histogram of QoI vs. PCE-based histogram of QoI
fig1 = figure(1);
h_MC = histogram(QoI);
hold on;
h_PCE = histogram(QoI_test);
xlabel('QoI realization');
ylabel('Frequency');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend([h_MC h_PCE],{'MC-based samples QoI','PCE-based samples QoI'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northwest';
set(gca,'FontSize',12);
% savefig(fig1,'hist_X_div_1plusX_Beta22_LegBasis_N2_M10P6.fig');
% saveas(fig1,'hist_X_div_1plusX_Beta22_LegBasis_N2_M10P6.eps','epsc');