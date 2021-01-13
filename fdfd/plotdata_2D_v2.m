R = csvread("data/reflectivity.csv");
omeg = R(1,1);

graphics_toolkit gnuplot;
f = figure('visible','off');

xfile = sprintf("data/X.csv");
yfile = sprintf("data/Y.csv");
efile = sprintf("data/Ez.csv");
x = 1e6*csvread(xfile);
y = 1e6*csvread(yfile);
Ez = csvread(efile);
Ez = reshape(Ez, length(x), length(y));
clf;
imagesc(x, y, real(Ez'));
set(gca, 'YDir', 'normal');
xlabel('x [{\mu}m]');
ylabel('y [{\mu}m]');
colorbar("EastOutside");
ylim([y(1) y(end)]);
xlim([x(1) x(end)]);
hold on;
%lims = max(Ez(i,:));
%lims = 1;
plot([x(1),x(end)], [y(SFy),y(SFy)], '-k', "linewidth", 2.5)
plot([x(1),x(end)], [y(p),y(p)], '-r', "linewidth", 2.5)
%plot([x(q),x(q)], [y(1),y(end)], '-k', "linewidth", 2.5)
%plot([z(p),z(p)], [-lims,lims], '-k', "linewidth", 2.5)
%plot([z(q),z(q)], [-lims,lims], '-k', "linewidth", 2.5)
%plot(z,Q(1:floor(end/2),1:floor(end/2)) * real(fsrc(z))', "-r;Re[f(z)];");
%legend("hide");
%plot([x(sfr+20), x(sfr+20)+2*pi/ki_x*1e6], [y(20), y(20)+2*pi/ki_x*1e6], '-r', "linewidth", 2.5)
quiver(x(2), y(SFy+20), lambx*1e6, lamby*1e6, '-r', "linewidth", 2.5)

freq_text = sprintf("%.2e", omeg);
title(['Re[E_z(x,y)] for \omega = ' freq_text]);
%text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
%text(z(sfr+10), -0.9, '\epsilon = \epsilon_0');
%text(z(p+50), -0.9, '\epsilon = 5\epsilon_0')
%text(z(q+10), -0.9, '\epsilon = \epsilon_0')

file_text=sprintf("figs/Ez.png");
saveas (gca, file_text);
%system ("wslview figs/test.png")

sprintf("\n Plotting Finished \n")