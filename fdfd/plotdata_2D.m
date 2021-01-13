R = csvread("data/reflectivity.csv");
omegs = R(:,1);
R = R(:,2);

graphics_toolkit gnuplot;
f = figure('visible','off');

%Rf = (((1-sqrt(5))/(1+sqrt(5))))^2 * ones(length(omegs),1);
clf;
plot(omegs, R, "-b;R;");
xlabel('\omega');
ylabel('R');
print("figs/reflection.png");

for n = 1:length(omegs)
  omeg = omegs(n);
  xfile = sprintf("data/X_%02d.csv", n);
  yfile = sprintf("data/Y_%02d.csv", n);
  efile = sprintf("data/Ez_%02d.csv", n);
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
  hold on;
  %lims = max(Ez(i,:));
  %lims = 1;
  plot([x(sfr),x(sfr)], [y(1),y(end)], '-k', "linewidth", 2.5)
  plot([x(p),x(p)], [y(1),y(end)], '-r', "linewidth", 2.5)
  %plot([x(q),x(q)], [y(1),y(end)], '-k', "linewidth", 2.5)
  %plot([z(p),z(p)], [-lims,lims], '-k', "linewidth", 2.5)
  %plot([z(q),z(q)], [-lims,lims], '-k', "linewidth", 2.5)
  %plot(z,Q(1:floor(end/2),1:floor(end/2)) * real(fsrc(z))', "-r;Re[f(z)];");
  %legend("hide");
  
  freq_text = sprintf("%.2e", omeg)
  title(['Re[E_z(x,y)] for \omega = ' freq_text]);
  %text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
  %text(z(sfr+10), -0.9, '\epsilon = \epsilon_0');
  %text(z(p+50), -0.9, '\epsilon = 5\epsilon_0')
  %text(z(q+10), -0.9, '\epsilon = \epsilon_0')
  
  file_text=sprintf("figs/Ez_%02d.png",n);
  saveas (gca, file_text);
  %system ("wslview figs/test.png")
endfor