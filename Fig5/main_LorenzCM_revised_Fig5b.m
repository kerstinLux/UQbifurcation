%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the framework ( if not already started)
uqlab
clear; clc;
%% Define the model
% Define the model options
modelopts.mFile = 'LorenzCMaux_function';
% Create and add the model to UQLab
myLorenzCM_Model = uq_createModel(modelopts);
% X = (-5:0.01:-1.01)';
% Y = uq_evalModel(X);
%% Define three independent random variables
for ii = 1 : 1
IOpts.Marginals(ii).Type = 'Uniform' ;
% IOpts.Marginals(ii).Type = 'Gaussian' ;
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
%% Graphical representation of spectrum of non-zero coefficients
% uq_display(myPCE_Quadrature);
Mom_PCE = myPCE_Quadrature.PCE.Moments;
c = myPCE_Quadrature.PCE.Coefficients

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N_trunc = MetaOpts.Degree;
N_trunc = 7;
%% Plain MC approach
M=10^6;
rng default
% Uniform Input
a = IOpts.Marginals(ii).Parameters(1); b = IOpts.Marginals(ii).Parameters(2);
U = a+(b-a)*rand(M,1);
% % Gaussian input
% mu = IOpts.Marginals(ii).Parameters(1); sigma = IOpts.Marginals(ii).Parameters(2);
% U = mu + sigma*randn(M,1);
QoI = LorenzCMaux_function(U);
EW = mean(QoI);
Var = std(QoI).^2;

%% Plain MC approach for sigma/(beta*(1+sigma))
%% Simulate beta(8,1)
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
% %% for u \in [-1,1]
% for n=1:N+1
%     expPoly = expPoly + c(n)*legendreP(n-1,u);
% end
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

%% Theoretical Moments of Lorenz randoom variable via Mellin Transform
moms_LorenzRV_Mellin = zeros(N_trunc+1,1);
moms_LorenzRV_Mellin(1) = 1;
for n=1:N_trunc
    moms_LorenzRV_Mellin(n+1) = gamma(gamma_val-(n+1)+1)/gamma(gamma_val)*moms_MellinPCE(n+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate moment-based density approximation according to "Moment-based Density Approximants", Provost(2005)
syms x
aux = 0;
supp = [0,0.5]; % Change support here for Figure 5(b)
for k=0:N_trunc
    legCoeffs = flip(sym2poly(legendreP(k,(supp(1)+supp(2)-2*x)/(supp(1)-supp(2)))));
    aux = aux + ((2*k+1)/(supp(2)-supp(1))*legCoeffs*moms_LorenzRV_Mellin(1:length(legCoeffs)))*legendreP(k,(supp(1)+supp(2)-2*x)/(supp(1)-supp(2)));
end
dens_Mellin = aux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate monic orthogonal polynomial density approximation as in "DistApproxandModellingViaOrthogonalPolynomialSequences", Provost (2015)
syms x
%% Calculate alpha,beta from B~gBeta(alpha,beta) s.th. c_w*E[B]=E[LorenzRV], c_w*E[B^2]=E[LorenzRV^2], where c_w is a normalizing constant
alpha = (moms_LorenzRV_Mellin(2)-supp(1))*(moms_LorenzRV_Mellin(3)+supp(1)*supp(2) - (supp(1)+supp(2))*moms_LorenzRV_Mellin(2))/((supp(2)-supp(1))*(moms_LorenzRV_Mellin(2)^2-moms_LorenzRV_Mellin(3)));
beta = (supp(2)-moms_LorenzRV_Mellin(2))*(moms_LorenzRV_Mellin(3)+supp(1)*supp(2) - (supp(1)+supp(2))*moms_LorenzRV_Mellin(2))/((supp(2)-supp(1))*(moms_LorenzRV_Mellin(2)^2-moms_LorenzRV_Mellin(3)));
w = @(x) x.^(alpha-1).*(1-x).^(beta-1);
c_w = 1/integral(w,supp(1),supp(2)); % normalizing constant

%% Calculate monic orthogonal polynomial coefficients
moms_Beta = ones(2*N_trunc-1,1);
for i=1:2*N_trunc-1
    aux = 0;
    for j=0:i
    aux = aux + nchoosek(i,j)*(supp(1)^(i-j)*(supp(2)-supp(1))^j*gamma(alpha+beta)*gamma(alpha+j))/(gamma(alpha+beta+j)*gamma(alpha));
    end
    moms_Beta(i+1) = aux;
end
m_Beta = 1/c_w*moms_Beta;
% Calculate Hankel determinant
D_Hankel = zeros(N_trunc+1,1);
D_Hankel(1,1) = 1;
for i=1:N_trunc
    D_Hankel(i+1,1) = calculate_HankelDet(m_Beta,i);
end

% Calculate modified Hankel determinant
D_Hankelmod = zeros(N_trunc+1,1);
D_Hankelmod(1,1) = 0;
D_Hankelmod(2,1) = 1;
for i=2:N_trunc
    for k=0:i-1
        D_Hankelmod(i+1,k+1) = calculate_HankelDetmod(m_Beta,i,k);
    end
end
%% Calculate corresponding monic orthogonal polynomials
d_monic = zeros(N_trunc-2,N_trunc+1);
for i=1:N_trunc-2
for k=0:i+2-1
    d_monic(i,k+1) = (-1)^(i+2-k)*D_Hankelmod(i+2+1,k+1)/D_Hankel(i+2+1,1);
end
d_monic(i,k+2) = 1;
end

%% ATTENTION: Coefficients of pi_i(x) do not coincide with example 3.1 from Provost(2015) --> Debug code for Hankel matrices and their determinants
% d_monic = [-0.1337 0.8861 -1.7096 1 0; 0.0529 -0.5313 1.722 -2.2302 1];

%% Calculate coefficients lambda_i in density approximation (14) in Provost (2015)
lambda = zeros(1,N_trunc+1);
aux = 0;
for i=1:N_trunc-2
    integrand = @(x) w(x).*(polyval(flip(d_monic(i,:)),x)).^2;
    intVal = integral(integrand,supp(1),supp(2));
    for k=0:i+2
        aux = aux + d_monic(i,k+1)*moms_LorenzRV_Mellin(k+1)/intVal;
    end
    lambda(i) = aux;
end

%% Calculate density approximation according to formula (14) in "2015DistApproxandModellingViaOrthogonalPolynomialSequences"
aux = @(x) 0;
for i=1:N_trunc-2
    aux = @(x) aux(x) + lambda(i)*polyval(flip(d_monic(i,:)),x);
end
dens_monicMellin = @(x) c_w*w(x) + w(x).*aux(x);


%% Calculate moment-based density approximation according to "Hausdorff moment problem: Reconstruction of PDFs", Mnatsakanov(2008b), formula (6)
% dens_trafoMomSequ = @(x) calculateDens_trafoMomSequ(moms_LorenzRV_Mellin,supp,N_trunc,x);

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
% fplot(dens_Mellin,[0 0.4],'r--');
% hold on
% plot(0:0.01:0.4,dens_monicMellin(0:0.01:0.4),'m*');
% l2 = plot(0:0.05:max(binCenters),dens_trafoMomSequ(0:0.05:max(binCenters)),'k--*','LineWidth',1);
xlabel('Value of bifurcation coefficient');
ylabel('PDF approximation');
ax.FontSize = 12;
ax.Interpreter = 'latex';
% leg = legend([barLorenzRV l1 l2],{'Normalized histogram','Legendre polynomial based','Transformed Moments based'});
leg = legend([barLorenzRV l1],{'Normalized histogram','Legendre polynomial based'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
% savefig(fig2,'histVsPDFapprox_LorenzRV_Uniform46_Gamma81_PolyApprox_chosenSupp0_0K5_N2_numMoms7.fig');
% saveas(fig2,'histVsPDFapprox_LorenzRV_Uniform46_Gamma81_PolyApprox_chosenSupp0_0K5_N2_numMoms7.eps','epsc');


%% Approximated cumulative distribution function
% Calculate moment-based cdf approximation according to "Hausdorff moment problem: Reconstruction of distributions", Mnatsakanov(2008), formula (2)
aux = 0;
cdf_approx = @(x) calculateCDF_approx(moms_LorenzRV_Mellin,supp,N_trunc,x);
y = 0.05;
% y = 0.1;
if y<supp(1) | y>supp(2)
    disp('Choose y value in within support of distribution.');
end
prob_less_y_Mellin = double(int(dens_Mellin,x,supp(1),y));
prob_less_y_monicMellin = integral(dens_monicMellin,supp(1),y);
prob_less_y_trafoMomSequ = cdf_approx(y);
prob_less_y_MC = sum(LorenzRV<=y)/M;
