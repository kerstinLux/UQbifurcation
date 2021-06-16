clear; clc;
%%%%%%%%%%%%% MC simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plain MC countercheck
M=10^6;
num_moms = 5;
rng default
% Uniform Input
ul = 0.95; ur = 1.2;
b = ul+(ur-ul)*rand(M,1);
a=0.3; d=0.7;
cumCoeffs = [-a -2/3*d 2/3];
testfun = @(b) cumCoeffs(3)*b.^2 + cumCoeffs(2)*b +cumCoeffs(1);
QoI = cumCoeffs(3)*b.^2 + cumCoeffs(2)*b +cumCoeffs(1); % Hindmarsh-Rose Lyapunov coeff. part decisive for sign
EW_QoI = mean(QoI);
Var_QoI = std(QoI).^2;
moms_MC = zeros(num_moms+1,1);
moms_MC(1) = 1;
for k=1:num_moms
    moms_MC(k+1) = 1/M*sum(QoI.^k);
end
histogram(QoI); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(cumCoeffs)-1;
%% Theoretical Moments of Lorenz random variable via Mellin Transform
%% Calculate num_moms moments of rv b via Mellin transform
moms_MellinPCE = zeros(num_moms+1,1);
moms_MellinPCE(1) = 1;
for n=1:num_moms
    c_MellinPCE = calcCoeffs_MellinPCE(cumCoeffs,N,n+1);
    % for uniform distribution U~U(0,1) and shifted Legendre polynomials
    counter = 1:1:N*n+1;
    moms_MellinPCE(n+1) = sum(c_MellinPCE./counter.*(ur.^counter-ul.^counter));
end

%%%%%%%%%%%%%%% ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian mixture model
numComponents = 3;
options = statset('MaxIter',1000);
GMM = fitgmdist(QoI,numComponents,'Options',options);
% options = statset('MaxIter',1000);
% GMM = fitgmdist(LorenzRV,numComponents,'Options',options);

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

pars_GMM = [GMM.ComponentProportion(1:end-1) GMM.mu' compStd_GMM];

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
mom_conds = @(pars) moms_MellinPCE(2:end) - moms_GaussMixt(pars);

%% Define quadratic form to minimize
% W = eye(num_moms); % positive definite weighting matrix
W = diag(abs(moms_MellinPCE(2:end)).^(-1)); % positive definite weighting matrix
qForm = @(pars) mom_conds(pars)'*W*mom_conds(pars);
init_weights = zeros(1,n_Gauss-1);
upperBound = 1;
rng(1)
for i=1:n_Gauss-1
    init_weights(i) = upperBound*rand(1);
    upperBound = sum(init_weights);
end
init_means = [0 0 0];
init_std = [1 1 1];
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

% Calculate MC CDF values
rangeLorenzRV = -0.6:0.001:1;
rangeLorenzRV_rough = rangeLorenzRV(1:10:end);
cdf_MC = zeros(size(rangeLorenzRV_rough));
for i=1:length(rangeLorenzRV_rough)
    cdf_MC(i) = sum(QoI<=rangeLorenzRV_rough(i))/M;
end

mean_Gmix = 0;
for i=1:n_Gauss-1
    mean_Gmix = mean_Gmix + pars_new(i)*pars_new(n_Gauss+i-1);
end
mean_Gmix = mean_Gmix + (1-sum(pars_new(1:n_Gauss-1)))*pars_new(n_Gauss+n_Gauss-1);

moms_ownGMM = calcMoms_GaussMixt(pars,n_Gauss,num_moms);

%% plot results
fig2 = figure(2);
h=histogram(QoI);
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
leg = legend([barLorenzRV l1 l2],{'Normalized histogram','Moment-based Gmix','EM Gmix'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
% savefig(fig2,'PDFapprox_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.fig');
% saveas(fig2,'PDFapprox_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.eps','epsc');

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
% savefig(fig3,'CDFapprox_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.fig');
% saveas(fig3,'CDFapprox_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.eps','epsc');

%% Calculate probability of supercritical Hopf bifurcation
y = [0 0.05 0.1 0.35];
prob_less_y_Gmix = zeros(1,length(y)); prob_less_y_GmixMat = zeros(1,length(y));
prob_less_y_MC = zeros(1,length(y));
for i=1:length(y)
    prob_less_y_Gmix(i) = cdf_Gmix(y(i)); % probability of supercritical Hopf bifurcation
    prob_less_y_GmixMat(i) = cdf(GMM,y(i));
    prob_less_y_MC(i) = sum(QoI<=y(i))/M;
end

%% Table of probabilities
MethodName = {'y';'Moment-based Gmix';'EM Gmix';'MC'};
prob_bifCoeffLess = [y; prob_less_y_Gmix; prob_less_y_GmixMat; prob_less_y_MC];

myTable_prob = table(prob_bifCoeffLess,...
    'RowNames',MethodName)
% writetable(myTable_prob,'prob_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.xlsx','WriteRowNames',true);
% writetable(myTable_prob,'prob_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.xlsx','WriteRowNames',true);

%% Table of moments
mom_method = {'Moment-based Gmix';'Mellin-based';'EM Gmix';'MC'};
calcMoms = [moms_ownGMM'; moms_MellinPCE(2:end)'; moms_GMM(2:end)'; moms_MC(2:end)'];

myTable_moms = table(calcMoms,...
    'RowNames',mom_method)
% writetable(myTable_moms,'moms_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.xlsx','WriteRowNames',true);
% writetable(myTable_moms,'moms_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.xlsx','WriteRowNames',true);

%% Table of optimized parameters
parsComp = {'Moment-based Gmix';'EM Gmix'};
pars_opt = [pars; pars_GMM];

myTable_pars = table(pars_opt,...
    'RowNames',parsComp)
% writetable(myTable_pars,'GMMpars_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.xlsx','WriteRowNames',true)
% writetable(myTable_pars,'GMMpars_HindmarshRose_Uni0K95_1K2_a_0K3_d0K7_numMoms5_numComp3_Maxfeval1e+4_maxIter5mal1e+3_stopStepTol1e-12.xlsx','WriteRowNames',true)