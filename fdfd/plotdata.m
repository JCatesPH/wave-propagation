Gamm = csvread("data/Gamm.csv");
Sz = csvread("data/sz.csv");
omegs = Gamm(:,1);
Gamma = Gamm(:,2);
z = Sz(:,1);
sz = Sz(:,2);

graphics_toolkit gnuplot;
f = figure('visible','off');

clf;
plot(omegs, Gamma);
xlabel('\omega');
ylabel('|\Gamma|');
print("figs/reflection.png");

clf;
plotyy(z, real(sz), z, imag(sz));
xlim([z(1) z(end)]);
xlabel('z [m]');
ylabel('s_z');
print("figs/sig.png");

sfr = 5000;
p = 12500;
q = 18750;

for i = 1:length(omegs)
  omeg = omegs(i);
  outfile = sprintf("data/ex_%02d.csv", i);
  Ex = csvread(outfile);
  Ex = Ex(:,2);
  clf;
  plot(z, real(Ex), '-b;Re[E(\omega,z)];');
  xlabel('z [m]');
  %xlim([z(1) z(end)]);
  hold on;
  %lims = max(Ez(i,:));
  lims = 1;
  plot([z(sfr),z(sfr)], [-lims,lims], '--k', "linewidth", 2.5)
  plot([z(p),z(p)], [-lims,lims], '-k', "linewidth", 2.5)
  plot([z(q),z(q)], [-lims,lims], '-k', "linewidth", 2.5)
  %plot(z,Q(1:floor(end/2),1:floor(end/2)) * real(fsrc(z))', "-r;Re[f(z)];");
  legend("hide");
  
  freqtext = sprintf("%.2e", omeg);
  text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
  text(z(sfr+10), -0.9, '\epsilon = \epsilon_0');
  text(z(p+10), -0.9, '\epsilon = 2\epsilon_0')
  text(z(q+10), -0.9, '\epsilon = \epsilon_0')
  
  file_text=sprintf("figs/output%d.png",i);
  saveas (gca, file_text);
  %system ("wslview figs/test.png")
endfor