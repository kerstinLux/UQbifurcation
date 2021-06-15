%% Create input density plot for LorenzRV with r_1 Beta and r_2 Gamma distributed
clear; clc;
% generalized support
a=-0.5;
b=0.5;
xrangeBeta = a:0.01:b;
% for Beta distribution of r_1
alpha = 2; beta = 5;
% for Gamma distribution of r_2
gamma = 8;
shape = 1; % shape in gampdf and not rate (shape = 1/rate)
xrangeGamma = 0:0.01:25;
% %% Plot densities of r_1 and r_2
% % % Open Figure hist_genBeta25.fig if countercheck with histogram von
% % % "uq_getSample(M)" desired
% % h=get(gca, 'children');
% % counts = h.Values;
% % binCenters = (h.BinEdges(2:end)+h.BinEdges(1:end-1))/2;
% % binWidth = (h.BinEdges(2)-h.BinEdges(1));
% % barGenBeta25 = bar(binCenters,counts/(sum(counts)*binWidth),'FaceColor',[0 0 1],'EdgeColor',[0 0 1]);
% hold on;
% l1 = plot(xrangeBeta,betapdf((xrangeBeta-a)/(b-a),2,5),'b','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b');
% l2 = plot(xrangeGamma, gampdf(xrangeGamma,8,1),'b--','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b');
% xlabel('Input value');
% ylabel('Input PDF');
% ax.FontSize = 12;
% ax.Interpreter = 'latex';
% % % if counter check used
% % leg = legend([barGenBeta25 l1 l2],{'Normalized histogram','Input density of $r_1$','Input density of $r_2$'});
% leg = legend([l1 l2],{'Input density of $r_1$','Input density of $r_2$'});
% leg.Interpreter = 'latex';
% leg.FontSize= 12;
% leg.Location = 'Northeast';
% set(gca,'FontSize',12);
%% Plot input density of r_2
fig1 = figure(1)
l1 = plot(xrangeBeta,betapdf((xrangeBeta-a)/(b-a),2,5),'b','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b');
xlabel('Input value');
ylabel('Input PDF');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend([l1],{'Input density of $r_1$'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
savefig(fig1,'inputPDF_genBeta25_suppMinus0K5_0K5.fig');
saveas(fig1,'inputPDF_genBeta25_suppMinus0K5_0K5.eps','epsc');
%% Plot input density of r_2
fig2 = figure(2)
l2 = plot(xrangeGamma, gampdf(xrangeGamma,8,1),'b','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b');
xlabel('Input value');
ylabel('Input PDF');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend(l2,{'Input density of $r_2$'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
savefig(fig2,'inputPDF_Gamma81.fig');
saveas(fig2,'inputPDF_Gamma81.eps','epsc');