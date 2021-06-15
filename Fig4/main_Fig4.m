%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
%% Download and installation of the framework UQLab from https://www.uqlab.com/ needed for this file to run
% UQLab: A Framework for Uncertainty Quantification in MATLAB,
% Stefano Marelli and Bruno Sudret,
% In The 2nd International Conference on Vulnerability and Risk Analysis and Management (ICVRAM 2014),
% University of Liverpool, United Kingdom, July 13-16, 2014, pp. 2554–2563.
% DOI: 10.1061/9780784413609.257
uqlab
% structure of PCE calculation taken from:
% @TechReport{Sudret.2019,
% 	author = {Marelli, S. and Sudret, B.},
% 	title = {{UQLab user manual -- Polynomial chaos expansions}},
% 	institution = {Chair of Risk, Safety and Uncertainty Quantification, ETH Zurich,
% 	Switzerland},
% 	year = {2019},
% 	note = {Report \# UQLab-V1.3-104}
% }
%% Define the model
% Define the model options
modelopts.mFile = 'LorenzCMaux_function';
% Create and add the model to UQLab
myLorenzCM_Model = uq_createModel(modelopts);

%% Specify parameters
para.alpha1 = 2;
para.beta1 = 2;
%% Define numInputRV independent random variables
numInputRV = 1;
supp_RV = [0 1];
for i = 1 : 1
    IOpts.Marginals(i).Type = 'Beta' ;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_moms = 5;

%%%%%%%%%%%%% MC simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plain MC countercheck
M=10^6;
rng default
%% Beta Input
U_test = uq_getSample(M);
U = U_test;
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