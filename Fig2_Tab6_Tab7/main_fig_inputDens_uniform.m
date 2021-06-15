%% Create input density plot for LorenzRV with r_1 Beta and r_2 Gamma distributed
clear; clc;
% generalized support
a=4;
b=6;
xrangeUniform = a-0.5:0.01:b+0.5;
%% Plot input density of r_1
fig1 = figure(1)
l1 = plot(xrangeUniform, unifpdf(xrangeUniform,a,b),'b','LineWidth',2,'MarkerSize',2,'MarkerEdgeColor','b','MarkerFaceColor','b');
xlabel('Input value');
ylabel('Input PDF');
ax.FontSize = 12;
ax.Interpreter = 'latex';
leg = legend([l1],{'Input density of $r_1$'});
leg.Interpreter = 'latex';
leg.FontSize= 12;
leg.Location = 'Northeast';
set(gca,'FontSize',12);
ylim([0 0.7])
savefig(fig1,'inputPDF_uniform46.fig');
saveas(fig1,'inputPDF_uniform46.eps','epsc');