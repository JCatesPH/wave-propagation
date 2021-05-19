%% Plots output from *.c FDTD simulations
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%%
X = readmatrix("../fdtd_c/data/X.csv");
Y = readmatrix("../fdtd_c/data/Y.csv");
Ez = readmatrix("../fdtd_c/data/Ez0000.csv");

%%
f = figure(1);
clf;
nlevels = 30;
colormap parula;
[C,h] = contourf(X, Y, Ez);
set(h,'LineColor','none');
set(h,'ZDataSource','Ez');
xlabel('$x$ [m]');
ylabel('$y$ [m]');
caxis([-0.4 1.0]);
colorbar();

for n = 1:10
    T = (n-1)*10;
    matfile = sprintf("../fdtd_c/data/Ez%04d.csv", T);
    Ez = readmatrix(matfile);
    
    titletxt = sprintf("$E_z(x,y)$ at $T=%d$", T);
    title(titletxt);
    
    refreshdata
    
    pngfile = sprintf("figs/output%d.png", n);
    saveas(f, pngfile);
end

