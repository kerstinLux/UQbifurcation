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

%% Plain MC approach for sigma/(theta*(1+sigma))
%% Simulate Gamma(gamma_val,1)
rng default
gamma_val = 8;
Theta = gamrnd(gamma_val,1,[M 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LorenzRV = Theta.^(-1).*QoI;
EW_LorenzRV = mean(LorenzRV);
Var_LorenzRV = std(LorenzRV)^2;
moms_MC = zeros(num_moms+1,1);
moms_MC(1) = 1;
for k=1:num_moms
    moms_MC(k+1) = 1/M*sum(LorenzRV.^k);
end
histogram(LorenzRV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compose PCE polynomial PCE from PCEcoeffs and chosen orthogonal polynomial MetaOpts.PolyTypes
syms u;
PCE = 0;
for n=1:N+1
    PCE = PCE + PCEcoeffs(n)*1/(sqrt(1/(2*(n-1)+1)))*legendreP(n-1,u);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate approximating PCE polynomial in u
expPoly = 0;
%% for u \in [0,1]
for n=1:N+1
    expPoly = expPoly + PCEcoeffs(n)*1/(sqrt(1/(2*(n-1)+1)))*legendreP(n-1,2*u-1);
end
cumCoeffs = flip(sym2poly(expPoly));

%% Calculate num_moms moments of PCE approximated rv via Mellin transform
moms_MellinPCE = zeros(num_moms+1,1);
moms_MellinPCE(1) = 1;
for n=1:num_moms
    c_MellinPCE = calcCoeffs_MellinPCE(cumCoeffs,N,n+1);
    % for uniform distribution U~U(0,1) and shifted Legendre polynomials
    counter = 1:1:N*n+1;
    moms_MellinPCE(n+1) = sum(c_MellinPCE./counter);
end
%% Theoretical moments of Lorenz random variable via Mellin transform
moms_LorenzRV_Mellin = zeros(num_moms+1,1);
moms_LorenzRV_Mellin(1) = 1;
for n=1:num_moms
    moms_LorenzRV_Mellin(n+1) = gamma(gamma_val-(n+1)+1)/gamma(gamma_val)*moms_MellinPCE(n+1);
end

%%%%%%%%%%%%%%% ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use Gaussian mixture model
numComponents = 2;
GMM = fitgmdist(LorenzRV,numComponents);

%% Calculate num_moms moments of GMM Matlab fit
compMoms_GMM = ones(num_moms,numComponents);
compMoms_GMM(2,:) = GMM.mu';
compMoms_GMM(3,:) = reshape(GMM.Sigma,1,numComponents) + GMM.mu.^2';
compStd_GMM = sqrt(reshape(GMM.Sigma,1,numComponents));
for n=1:numComponents
    for i=3:num_moms
        compMoms_GMM(i+1,n) = calcMomGaussian(compMoms_GMM(2,n),compStd_GMM(n),i);
    end
end

pars_GMM = [GMM.ComponentProportion(1) GMM.mu' compStd_GMM];

mean_GMM = GMM.ComponentProportion*compMoms_GMM(2,:)';
mom2_GMM = GMM.ComponentProportion*compMoms_GMM(3,:)';
moms_GMM = ones(num_moms+1,1);
moms_GMM(2) = mean_GMM;
moms_GMM(3) = mom2_GMM;
for i=3:num_moms
    moms_GMM(i+1) = GMM.ComponentProportion*compMoms_GMM(i+1,:)';
end
% Sample from GMM
samplesGMM = random(GMM,M);
meanGMM_sampled = mean(samplesGMM);
stdGMM_sampled = std(samplesGMM);
mom2GMM_sampled = std(samplesGMM)^2 + mean(samplesGMM)^2;
if abs(mom2GMM_sampled-moms_GMM(3))>10^(-3)
    disp('Attention: sample second moment GMM not close to calculated second moment.')
end
plot((0:0.01:2)',pdf(GMM,(0:0.01:2)'),'k');

%% Estimate Gaussian mixture model from Mellin-based moments
%% Choose number of Gaussians
n_Gauss = numComponents;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define moment conditions
moms_GaussMixt = @(pars) calcMoms_GaussMixt(pars,n_Gauss,num_moms);
mom_conds = @(pars) moms_LorenzRV_Mellin(2:end) - moms_GaussMixt(pars);

%% Define quadratic form to minimize
% W = eye(num_moms); % positive definite weighting matrix if zero moment appears
W = diag(abs(moms_LorenzRV_Mellin(2:end)).^(-1)); % positive definite weighting matrix
qForm = @(pars) mom_conds(pars)'*W*mom_conds(pars);
init_weights = zeros(1,n_Gauss-1);
upperBound = 1;
rng(1)
for i=1:n_Gauss-1
    init_weights(i) = upperBound*rand(1);
    upperBound = sum(init_weights);
end
init_means = [0.2 0.8];
init_std = [0.01 0.1];
pars_init = [init_weights init_means init_std];
options = optimset('Display','iter','MaxFunEvals',10^4,'MaxIter',5*10^3,'TolX',1e-12);
[pars,qForm_val] = fmincon(qForm,pars_init,[],[],[],[],[zeros(1,n_Gauss-1) -Inf*ones(1,n_Gauss) zeros(1,n_Gauss)],[ones(1,n_Gauss-1) Inf*ones(1,2*n_Gauss)],[],options);
pars_new = pars;

% Define density of Gaussian mixture model
dens_Gmix = @(x) 0;
for i=1:n_Gauss-1
    dens_Gmix = @(x) dens_Gmix(x) + pars_new(i)*normpdf(x,pars_new(n_Gauss+i-1),pars_new(2*n_Gauss+i-1));
end
dens_Gmix = @(x) dens_Gmix(x) + (1-sum(pars_new(1:n_Gauss-1)))*normpdf(x,pars_new(n_Gauss+n_Gauss-1),pars_new(2*n_Gauss+n_Gauss-1));
% Define CDF of Gaussian mixture model
cdf_Gmix = @(x) 0;
for i=1:n_Gauss-1
    cdf_Gmix = @(x) cdf_Gmix(x) + pars_new(i)*normcdf(x,pars_new(n_Gauss+i-1),pars_new(2*n_Gauss+i-1));
end
cdf_Gmix = @(x) cdf_Gmix(x) + (1-sum(pars_new(1:n_Gauss-1)))*normcdf(x,pars_new(n_Gauss+n_Gauss-1),pars_new(2*n_Gauss+n_Gauss-1));

%% Calculate MC CDF values
rangeLorenzRV = -0.1:0.001:0.4;
rangeLorenzRV_rough = rangeLorenzRV(1:10:end);
cdf_MC = zeros(size(rangeLorenzRV_rough));
for i=1:length(rangeLorenzRV_rough)
    cdf_MC(i) = sum(LorenzRV<=rangeLorenzRV_rough(i))/M;
end

mean_Gmix = 0;
for i=1:n_Gauss-1
    mean_Gmix = mean_Gmix + pars_new(i)*pars_new(n_Gauss+i-1);
end
mean_Gmix = mean_Gmix + (1-sum(pars_new(1:n_Gauss-1)))*pars_new(n_Gauss+n_Gauss-1);

moms_ownGMM = calcMoms_GaussMixt(pars,n_Gauss,num_moms);

%% plot results
fig2 = figure(2);
h=histogram(LorenzRV);
counts = h.Values;
binCenters = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;
binWidth = (h.BinEdges(2)-h.BinEdges(1));
barLorenzRV = bar(binCenters,counts/(M*binWidth),'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
hold on;
l1 = plot(rangeLorenzRV, dens_Gmix(rangeLorenzRV),'c','LineWidth',2);
l2 = plot((rangeLorenzRV)',pdf(GMM,(rangeLorenzRV)'),'m--','LineWidth',2);
xlabel('Value of bifurcation coefficient');
ylabel('PDF approximation');
ax.FontSize = 12;
ax.Interpreter = 'latex';
xlim([-0.1 0.4])
leg = legend([barLorenzRV l1 l2],{'Normalized histogram','Moment-based Gmix','EM Gmix'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
% savefig(fig2,'histVsPDFapprox_LorenzRV_Beta22_Gamma81_LegBasis_N2_numMoms5_numComp2_maxIter5mal1e+3_StepTol1e-12_stopMaxfeval1e+4.fig');
% saveas(fig2,'histVsPDFapprox_LorenzRV_Beta22_Gamma81_LegBasis_N2_numMoms5_numComp2_maxIter5mal1e+3_StepTol1e-12_stopMaxfeval1e+4.eps','epsc');

fig3 = figure(3);
l1 = plot(rangeLorenzRV, cdf_Gmix(rangeLorenzRV),'c','LineWidth',2);
hold on;
l2 = plot((rangeLorenzRV)',cdf(GMM,(rangeLorenzRV)'),'m--','LineWidth',2);
l3 = plot(rangeLorenzRV_rough,cdf_MC,'b*','MarkerSize',8);
xlabel('Value of bifurcation coefficient');
ylabel('CDF approximation');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend([l1 l2 l3],{'Moment-based Gmix','EM Gmix','MC-based'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Southeast';
set(gca,'FontSize',12);
% savefig(fig3,'CDFapprox_LorenzRV_Beta22_Gamma81_LegBasis_N2_numMoms5_numComp2_maxIter5mal1e+3_StepTol1e-12_stopMaxfeval1e+4.fig');
% saveas(fig3,'CDFapprox_LorenzRV_Beta22_Gamma81_LegBasis_N2_numMoms5_numComp2_maxIter5mal1e+3_StepTol1e-12_stopMaxfeval1e+4.eps','epsc');

%% For calculation of probability of subcritical Hopf bifurcation
y = [0 0.05 0.1 0.35];
prob_less_y_Gmix = zeros(1,length(y)); prob_less_y_GmixMat = zeros(1,length(y));
prob_less_y_MC = zeros(1,length(y));
for i=1:length(y)
    prob_less_y_Gmix(i) = cdf_Gmix(y(i));
    prob_less_y_GmixMat(i) = cdf(GMM,y(i));
    prob_less_y_MC(i) = sum(LorenzRV<=y(i))/M;
end

%% Table of probabilities
MethodName = {'y';'Moment-based Gmix';'EM Gmix';'MC'};
prob_bifCoeffLess = [y; prob_less_y_Gmix; prob_less_y_GmixMat; prob_less_y_MC];

myTable_prob = table(prob_bifCoeffLess,...
    'RowNames',MethodName)

%% Table of moments
mom_method = {'Moment-based Gmix';'Mellin-based';'EM Gmix';'MC'};
calcMoms = [moms_ownGMM'; moms_LorenzRV_Mellin(2:end)'; moms_GMM(2:end)'; moms_MC(2:end)'];

myTable_moms = table(calcMoms,...
    'RowNames',mom_method)

%% Table of optimized parameters
parsComp = {'Moment-based Gmix';'EM Gmix'};
pars_opt = [pars; pars_GMM];

myTable_pars = table(pars_opt,...
    'RowNames',parsComp)