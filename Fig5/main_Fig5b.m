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
%% Define independent random variables
for ii = 1 : 1
    IOpts.Marginals(ii).Type = 'Uniform' ;
    IOpts.Marginals(ii).Parameters = [4, 6] ;
end
myInput = uq_createInput(IOpts);

%% Setup of PCE
MetaOpts.Type = 'uq_metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.PolyTypes = {'Legendre'};

% Specification of 14th degree, Gaussian quadrature based projection
MetaOpts.Degree = 2;
N=MetaOpts.Degree;
MetaOpts.Method = 'quadrature';
% Creation of the metamodel :
myPCE_Quadrature = uq_createModel(MetaOpts);
%% Report basic information
uq_print(myPCE_Quadrature);
Mom_PCE = myPCE_Quadrature.PCE.Moments;
c = myPCE_Quadrature.PCE.Coefficients;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_trunc = 7;
%% Plain MC approach
M=10^6;
rng default
% Uniform Input
a = IOpts.Marginals(ii).Parameters(1); b = IOpts.Marginals(ii).Parameters(2);
U = a+(b-a)*rand(M,1);
QoI = LorenzCMaux_function(U);
EW = mean(QoI);
Var = std(QoI).^2;

%% Plain MC approach for zeta/(Beta*(1+zeta))
%% Simulate gamma(8,1)
% rng default
gamma_val = 8;
Beta = gamrnd(gamma_val,1,[M 1]);
LorenzRV = Beta.^(-1).*QoI;
EW_LorenzRV = mean(LorenzRV);
Var_LorenzRV = std(LorenzRV)^2;
moms_MC = zeros(N_trunc+1,1);
moms_MC(1) = 1;
for k=1:N_trunc
    moms_MC(k+1) = 1/M*sum(LorenzRV.^k);
end
histfit(LorenzRV);

%% Generate approximating PCE polynomial in u
syms u;
expPoly = 0;
%% for u \in [0,1]
for n=1:N+1
    expPoly = expPoly + c(n)*1/(sqrt(1/(2*(n-1)+1)))*legendreP(n-1,2*u-1);
end
cumCoeffs = flip(sym2poly(expPoly));

PCE = 0;
for n=1:N+1
    PCE = PCE + c(n)*1/(sqrt(1/(2*(n-1)+1)))*legendreP(n-1,u);
end

moms_MellinPCE = zeros(N_trunc+1,1);
moms_MellinPCE(1) = 1;
for n=1:N_trunc
    c_MellinPCE = calcCoeffs_MellinPCE(cumCoeffs,N,n+1);
    % for uniform distribution U~U(0,1) and shifted Legendre polynomials
    counter = 1:1:N*n+1;
    moms_MellinPCE(n+1) = sum(c_MellinPCE./counter);
end

%% Theoretical moments of Lorenz random variable via Mellin transform
moms_LorenzRV_Mellin = zeros(N_trunc+1,1);
moms_LorenzRV_Mellin(1) = 1;
for n=1:N_trunc
    moms_LorenzRV_Mellin(n+1) = gamma(gamma_val-(n+1)+1)/gamma(gamma_val)*moms_MellinPCE(n+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate moment-based density approximation according to "Moment-based Density Approximants", Provost(2005)
syms x
aux = 0;
supp = [0,0.5];
for k=0:N_trunc
    legCoeffs = flip(sym2poly(legendreP(k,(supp(1)+supp(2)-2*x)/(supp(1)-supp(2)))));
    aux = aux + ((2*k+1)/(supp(2)-supp(1))*legCoeffs*moms_LorenzRV_Mellin(1:length(legCoeffs)))*legendreP(k,(supp(1)+supp(2)-2*x)/(supp(1)-supp(2)));
end
dens_Mellin = aux;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots
fig2 = figure(2)
h=histogram(LorenzRV);
counts = h.Values;
binCenters = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;
binWidth = (h.BinEdges(2)-h.BinEdges(1));
barLorenzRV = bar(binCenters,counts/(M*binWidth),'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
hold on;
l1 = fplot(dens_Mellin,[0 max(binCenters)],'r--','LineWidth',2);
xlabel('Value of bifurcation coefficient');
ylabel('PDF approximation');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend([barLorenzRV l1],{'Normalized histogram','Legendre polynomial based'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
% savefig(fig2,'histVsPDFapprox_LorenzRV_Uniform46_Gamma81_PolyApprox_chosenSupp0_0K5_N2_numMoms7.fig');
% saveas(fig2,'histVsPDFapprox_LorenzRV_Uniform46_Gamma81_PolyApprox_chosenSupp0_0K5_N2_numMoms7.eps','epsc');