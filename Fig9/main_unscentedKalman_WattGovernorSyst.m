%% Input distribution 1: X1~U(a,b)
a = 0;
b = 0.2;
meanX1 = (a+b)/2;
varX1 = 1/12*(b-a)^2;
mom3X1 = ((b-meanX1)^(3+1)-(a-meanX1)^(3+1))/(3+1);
mom4X1 = ((b-meanX1)^(4+1)-(a-meanX1)^(4+1))/(4+1);
% %% Input distribution X~N(mu,sigma)
% mu = 2;
% sigma = 0.5;
% meanX = mu;
% varX = sigma^2;
% mom3X = 0;
% mom4X = 3*sigma^4;

%% Input distribution 2: X2~U(aa,bb)
aa = 0.7;
bb = 0.95;
meanX2 = (aa+bb)/2;
varX2 = 1/12*(bb-aa)^2;
mom3X2 = ((bb-meanX2)^(3+1)-(aa-meanX2)^(3+1))/(3+1);
mom4X2 = ((bb-meanX2)^(4+1)-(aa-meanX2)^(4+1))/(4+1);

%% Mean and covariance matrix assuming that X1 and X2 are independent
meanX = [meanX1; meanX2];
CX = diag([varX1 varX2]);
rootCX = chol(CX);

%% Sample from input distribution and transform samples
M = 10^5;
rng('default')
%% Simulate X1 and calculate empirical moments
X1 = a + (b-a)*rand(1,M);
% X = mu + sigma*randn(M,1);
meanX1_MC = mean(X1);
varX1_MC = var(X1);
mom3X1_MC = 1/M*sum((X1-meanX1).^3);
mom4X1_MC = 1/M*sum((X1-meanX1).^4);

%% Simulate X2 and calculate empirical moments
X2 = aa + (bb-aa)*rand(1,M);
% X = mu + sigma*randn(M,1);
meanX2_MC = mean(X2);
varX2_MC = var(X2);
mom3X2_MC = 1/M*sum((X2-meanX2).^3);
mom4X2_MC = 1/M*sum((X2-meanX2).^4);

% %% calculate output statistics for X1
% Y1 = g(X1);
% meanY_MC = mean(Y1);
% varY_MC = var(Y1);

%% calculate output statistics for X1 and X2
Y = calcFirstLyapCoeff_WattGovernorSyst([X1;X2]);
prob_subHopf_MC = sum(Y>0)/M;
meanY_MC = mean(Y);
varY_MC = var(Y);
mom3_MC = 1/M*sum(Y.^3);
mom4_MC = 1/M*sum(Y.^4);
mom5_MC = 1/M*sum(Y.^5);
mom6_MC = 1/M*sum(Y.^6);
mom7_MC = 1/M*sum(Y.^7);
mom8_MC = 1/M*sum(Y.^8);
moms_MC = [1; meanY_MC; varY_MC+meanY_MC^2; mom3_MC; mom4_MC; mom5_MC; mom6_MC; mom7_MC; mom8_MC];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unscented Kalman filter
precision = 5;
switch precision
    case 3
        %% Compute sigma points via eq. (12) from "Unscented Filtering and Nonlinear Estimation", Julier, Uhlmann (2004)
        % Use Cholesky decomposition for matrix square root in calculation of sigma
        % points
        dim = 2;
        numPoints = 2*dim+1;
        % w0 = 1/(2*dim+1);
        w0 = 1- dim/3;
        x = [meanX meanX+sqrt(dim/(1-w0))*rootCX(1,:)' meanX+sqrt(dim/(1-w0))*rootCX(2,:)' meanX-sqrt(dim/(1-w0))*rootCX(1,:)' meanX-sqrt(dim/(1-w0))*rootCX(2,:)'];
        %% Assign a weight to each sigma point
        kappa = 2; % for Gaussian dist matches first 4 moms if dim+kappa=3
        % w0 = kappa/(kappa+n);
        % w1 = 1/(2*(n+kappa));
        % w2 = 1/(2*(n+kappa));
        w1 = (1-w0)/(2*dim);
        w2 = (1-w0)/(2*dim);
        w3 = (1-w0)/(2*dim);
        w4 = (1-w0)/(2*dim);
        w = [w0; w1; w2; w3; w4];
        % w = 1/length(x)*ones(size(x));
        if sum(w)-1>eps
            disp('Weights must sum up to 1.')
        end
    case 5
        %% Compute sigma points via fifth order set of [37] from "A Systematization of the Unscented
        %% Kalman Filter Theory", Menegaz et al. (2015)
        dim = 2;
        numPoints = 2*dim^2+1;
        x = zeros(2,numPoints);
        w = zeros(1,numPoints);
        % set up weights of sigma points
        w0 = (dim^2-7*dim)/18 + 1;
        w(1) = w0;
        w(2:2*dim+1) = (4-dim)/18;
        w(2*dim+2:2*dim^2+1) = 1/36;
        if sum(w)-1>eps
            disp('Weights must sum up to 1.')
        end
        % set up sigma points
        x(:,1) = meanX; % first point is mean
        % type 1 and type 2 points
        genPoint = sqrt(3);
        xi = zeros(2,numPoints);
        % type 1 point
        xi(:,2) = [genPoint;0];
        xi(:,3) = [0;genPoint];
        xi(:,4) = [-genPoint;0];
        xi(:,5) = [0;-genPoint];
        % type 2 point
        xi(:,6) = [genPoint;genPoint];
        xi(:,7) = [genPoint;-genPoint];
        xi(:,8) = [-genPoint;genPoint];
        xi(:,9) = [-genPoint;-genPoint];
        for i=2:numPoints
            x(:,i) = meanX + rootCX*xi(:,i);
        end
    otherwise
        disp('Precision not specified.');
end

%% Can sigma points capture mean and var of X?
meanX1_UT = sum(w.*x(1,:));
meanX2_UT = sum(w.*x(2,:));
meanX_UT = [meanX1_UT; meanX2_UT];
cMeanX1 = x(1,:)-meanX1;
cMeanX2 = x(2,:)-meanX2;
varX1X2_UT = sum(w.*(cMeanX1)*(cMeanX2)');
varX1_UT = sum(w.*(cMeanX1)*(cMeanX1)');
varX2_UT = sum(w.*(cMeanX2)*(cMeanX2)');
% varX_UT = sum(w.*(x-meanX1).^2);
CX_UT = [varX1_UT varX1X2_UT; varX1X2_UT varX2_UT];
% cMeanX1 = X1-meanX1;
% cMeanX2 = X2-meanX2;
% varX1X2_UT = sum(w.*(cMeanX1)*(cMeanX2)');
% varX1_UT = sum(w.*(cMeanX1)*(cMeanX1)');
% varX2_UT = sum(w.*(cMeanX2)*(cMeanX2)');
% % varX_UT = sum(w.*(x-meanX1).^2);
% CX_UT = [varX1_UT varX1X2_UT; varX1X2_UT varX2_UT];

%% Transform sigma points through nonlinear transformation
y = calcFirstLyapCoeff_WattGovernorSyst(x);
%% Compute distribution from weighted points
prob_subHopf_UT = sum(y>0)/numPoints;
% prob_subHopf_UT = sum(y>0)/(2*dim+1);
meanY_UT = sum(w.*y(1,:));
varY_UT = sum(w.*(y-meanY_UT).^2);

% %% Approximate density function via characteristic function according to DURRETT probability and measures, p. 118
% syms t
% charfunY = 1 +1i*t*meanY_UT - t^2*(varY_UT+meanY_UT^2)/2; % approximation is o(t^2)
% % From characteristic function phi(t) to probability density rhoY if
% % int|phi(t)| dt < \infty, then probability measure mu has bounded
% % continuous density rhoY  = 1/(2*pi)int e^(-ity)phi(t) dt
% syms y
% integrand = 1/(2*pi)*exp(-1i*t*y)*charfunY;
% LargeNum = 10^5;
% intVal = int(integrand,t,-LargeNum,LargeNum);
% 
fig1 = figure(1)
h = histfit(Y,[],'kernel');
ax=gca;
yl = ylim;
hold on;
plot(y,0,'mo','MarkerSize',8,'MarkerFaceColor','magenta','MarkerEdgeColor','magenta');
plot([0 0],yl,'k', 'LineWidth',2)
xlabel('$s(\beta,\alpha)$','interpreter','latex');
% savefig(fig1,'LyapunovCoeff_WattGovernor_alphaU0_0K2_beta_U0K7_0K95_prec5.fig');
% saveas(fig1,'LyapunovCoeff_WattGovernor_alphaU0_0K2_beta_U0K7_0K95_prec5.eps','epsc');

