%% PDF of product of uniform RV with negative part of support and gamma RV via generalized Mellin transform
clear; clc;
%% Sample from product distribution
M=10^6;
rng default
X=-1+4*rand(1,M);
Y=gamrnd(3,1,1,M);
QoI = X.*Y;
%% According to procedure in , "Some applications of the Mellin transform in Statistics", Epstein (1948)
h1 = @(y) 1/(4*gamma(3))*(1/3*y.*exp(-1/3*y)+exp(-1/3*y)).*(y>=0);
h3 = @(y) 1/(4*gamma(3))*(-y.*exp(y)+exp(y)).*(y<0);
%% plot results
fig1 = figure(1);
h=histogram(QoI);
counts = h.Values;
binCenters = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;
binWidth = (h.BinEdges(2)-h.BinEdges(1));
barProdRV = bar(binCenters,counts/(M*binWidth),'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
hold on;
xrange = -10:0.01:30;
l1 = plot(xrange,h1(xrange)+h3(xrange),'r','LineWidth',2);
plot([0 0],[0 0.14],'k','LineWidth',2);
xlabel('Value of bifurcation coefficient');
ylabel('PDF of bifurcation coefficient');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend([barProdRV l1],{'Normalized histogram','Mellin-based PDF'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
% savefig(fig1,'PDF_ana_Prod_UniformMinus1_3_Gamma3_1.fig');
% saveas(fig1,'PDF_ana_Prod_UniformMinus1_3_Gamma3_1.eps','epsc');


%% Calculate CDF of product of uniform RV with negative part of support and gamma RV
H = @(y) (1/(4*gamma(3))*((2-y).*exp(y))).*(y<0) + (1/(4*gamma(3))*(-y.*exp(-1/3*y)-6*exp(-1/3*y)+8)).*(y>=0);
xrange_rough = xrange(1:100:end);
CDF_MC = zeros(size(xrange_rough));
for i=1:length(xrange_rough)
    CDF_MC(i) = sum(QoI<=xrange_rough(i))/M;
end
fig2 = figure(2);
l1 = plot(xrange,H(xrange),'r','LineWidth',2);
hold on
l2 = plot(xrange_rough,CDF_MC,'b*','MarkerSize',8);
plot([0 0],[0 1],'k','LineWidth',2);
xlabel('Value of bifurcation coefficient');
ylabel('CDF of bifurcation coefficient');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend([l1 l2],{'Mellin-based CDF','MC-based CDF approx'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Southeast';
set(gca,'FontSize',12);
% savefig(fig2,'CDF_ana_Prod_UniformMinus1_3_Gamma3_1.fig');
% saveas(fig2,'CDF_ana_Prod_UniformMinus1_3_Gamma3_1.eps','epsc');

%% Calculate probability of subcritical bifurcation P(-r_1r_2>0)=P(r_1r_2<0)
probSubBif = H(0);
probSubBif_MC = sum(QoI<=0)/M;
probSubBifMessage = ['The probability of observing a subcritical bifurcation is',num2str(probSubBif), 'and the MC probability is',num2str(probSubBif_MC)];
disp(probSubBifMessage);