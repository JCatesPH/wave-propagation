%% Plots data from the new function and struct
%% Read in the parameters
pathhead = "PMLtest/UPML_om01th01_";
tab = readtable(strcat(pathhead, "params.txt"));

lamb0 = 6e8 .* pi ./ tab.omeg;
lambx = lamb0 .* sin(tab.thetai);
lamby = lamb0 .* cos(tab.thetai);
sfr = tab.Lx + tab.bufsize;
p = floor(0.5 * (tab.Nx-sfr) + sfr);

%% Load in data
xfile = strcat(pathhead, "X.csv");
yfile = strcat(pathhead, "Y.csv");
Efile = strcat(pathhead, "Ez.csv");
x = 1e6*readmatrix(xfile);
y = 1e6*readmatrix(yfile);
Ez = readmatrix(Efile);
Ez = reshape(Ez, length(x), length(y));

%% Plot filled contour and lines
f = figure;
[M,c] = contourf(x, y, real(Ez'), 64);
set(c,'LineColor','none')
xlabel('x [{\mu}m]');
ylabel('y [{\mu}m]');
colorbar("EastOutside");
% ylim([y(1) y(end)]);
% xlim([x(1) x(end)]);
hold on;
%lims = max(Ez(i,:));
%lims = 1;
% Add lines for the SF region and region boundary
plot([x(sfr),x(sfr)], [y(1),y(end)], '--r', "linewidth", 2.5)
plot([x(p),x(p)], [y(1),y(end)], '-k', "linewidth", 2.5)
plot([x(tab.Rpind),x(tab.Rpind)], [y(1),y(end)], '--k', "linewidth", 2.5)
plot([x(tab.Rpind),x(tab.Rpind)], [y(1),y(end)], '--k', "linewidth", 2.5)

quiver(x(sfr+20), y(20), lamby*1e6, lambx*1e6, '-r', "linewidth", 2.5)

freq_text = sprintf("%.2e", tab.omeg);
title(strcat('E_z(x,y) for \omega = ', freq_text));
%text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
%text(z(sfr+10), -0.9, '\epsilon = \epsilon_0');
%text(z(p+50), -0.9, '\epsilon = 5\epsilon_0')
%text(z(q+10), -0.9, '\epsilon = \epsilon_0')

%%
figfile = strcat(pathhead, "Ez.png");
saveas(f, figfile)
%system ("wslview figs/test.png")
