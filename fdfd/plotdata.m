Gamm = csvread("data/Gamm.csv");
Sz = csvread("data/sz.csv");

graphics_toolkit gnuplot;
f = figure('visible','off');

clf;
plot(Gamm(:,1), Gamm(:,2));
xlabel('\omega');
ylabel('|\Gamma|');
print("figs/reflection.png");

clf;
plot(Sz(:,1), Sz(:,2));
xlim([Sz(1,1) Sz(end,1)]);
xlabel('z [m]');
ylabel('s_z');
print("figs/sig.png");