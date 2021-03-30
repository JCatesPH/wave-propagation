%% Tests the FDFD with a box in TF region.
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
%% Mesh specifications
% Number of cells to split vacuum wavelengths into
param2D.dx = 5e-8; 
param2D.dy = 5e-8;
% Number of cells for PML
param2D.Lx = 20; 
param2D.Ly = 20; 
% Enforce distance of PML from TF region
param2D.bx = 80; 
param2D.by = 80; 
% Total number of cells for x and y axis
param2D.Nx = 800 + 2*(param2D.Lx + param2D.bx);
param2D.Ny = 600 + 2*(param2D.Ly + param2D.by);
Nc = param2D.Nx*param2D.Ny

%% PML parameters
% Set the polynomial order. (1 <= p < 5)
param2D.sx_m = 3;
param2D.sy_m = 3;
% Desired reflectivity coefficient
param2D.sx_R = 1e-9;
param2D.sy_R = 1e-9;
% Set s_max. (0 <= s_max <= 5)
param2D.sx_smax = 3;
param2D.sy_smax = 3;

%% Space parameters
% Set the boundaries 
% -- Currently only for region split into two halves (ADD FUNCTIONALITY LATER)
% param2D.xbounds = [];
% Set the relative permittivities and permeabilities in each region
%   Only for second region in current implementation.
epsr1 = 2.25;
param2D.epsr = epsr1;
param2D.mur = 1

%% Domain parameters
ff = 300;
bf = 500;
ls = param2D.Ly + param2D.by;
rs = param2D.Ny - (param2D.Ly + param2D.by);

%%
eps_zz = ones(param2D.Nx, param2D.Ny);
eps_zz(ff:bf-1,ls:rs) = epsr1 * eps_zz(ff:bf-1,ls:rs);
eps_zz(bf,ls:rs-1) = epsr1 * eps_zz(bf,ls:rs-1);

Q = zeros(param2D.Nx, param2D.Ny);
Q(ff:bf-1,ls:rs) = Q(ff:bf-1,ls:rs) + 1;
Q(bf,ls:rs-1) = Q(bf,ls:rs-1) + 1;
q = Q(:);

%% Set the omeg sampling and and angle of incidence
omeg = 1000e12;
thetai = 0.0;

%% Test
chi3 = 1e-16;
dP = spalloc(Nc, Nc, sum(q));
PNL = @(x) 3*chi3/4*(q.*abs(x)).^2 .*x;
dPNL = @(x) spdiags(9*chi3/4*(q.*abs(x)).^2, 0, dP);
pathhead = sprintf("nltest/kerr");
param2D.dir = pathhead;
[x, y, Ez] = fdfd2D_nonlinear(omeg, thetai, eps_zz, PNL, dPNL, param2D);


%% Plots data
lamb0 = 6e8 .* pi ./ omeg;
lambx = lamb0 .* sin(thetai);
lamby = lamb0 .* cos(thetai);

f = figure(1);
clf;
[M,c] = contourf(x, y, real(Ez'), 64);
set(c,'LineColor','none')
xlabel('$x$ [m]');
ylabel('$y$ [m]');
colorbar("EastOutside");
% ylim([y(1) y(end)]);
% xlim([x(1) x(end)]);
hold on;
%lims = max(Ez(i,:));
%lims = 1;

% Add PML boundaries
plot([x(param2D.Lx),x(param2D.Lx)], [y(1),y(end)], '--r', "linewidth", 1.5)
plot([x(end-param2D.Lx),x(end-param2D.Lx)], [y(1),y(end)], '--r', "linewidth", 1.5)
plot([x(1),x(end)], [y(param2D.Ly),y(param2D.Ly)], '--r', "linewidth", 1.5)
plot([x(1),x(end)], [y(end-param2D.Ly),y(end-param2D.Ly)], '--r', "linewidth", 1.5)

% Add outline of Total-Field Region
plot([x(param2D.Lx+param2D.bx),x(param2D.Lx+param2D.bx)], ...
    [y(param2D.Ly+param2D.by),y(end-param2D.Ly-param2D.by)], '-k', "linewidth", 1.5)
plot([x(end-param2D.Lx-param2D.bx),x(end-param2D.Lx-param2D.bx)], ...
    [y(param2D.Ly+param2D.by),y(end-param2D.Ly-param2D.by)], '-k', "linewidth", 1.5)
plot([x(param2D.Lx+param2D.bx),x(end-param2D.Lx-param2D.bx)], ...
    [y(param2D.Ly+param2D.by),y(param2D.Ly+param2D.by)], '-k', "linewidth", 1.5)
plot([x(param2D.Lx+param2D.bx),x(end-param2D.Lx-param2D.bx)], ...
    [y(end-param2D.Ly-param2D.by),y(end-param2D.Ly-param2D.by)], '-k', "linewidth", 1.5)

% Add outline of 'box'
plot([x(ff),x(ff)], [y(ls),y(rs)], '--k', "linewidth", 1.5)
plot([x(bf),x(bf)], [y(ls),y(rs)], '--k', "linewidth", 1.5)
plot([x(bf),x(ff)], [y(ls),y(ls)], '--k', "linewidth", 1.5)
plot([x(bf),x(ff)], [y(rs),y(rs)], '--k', "linewidth", 1.5)

%quiver(x(sfr+20), y(20), lamby*1e6, lambx*1e6, '-r', "linewidth", 2.5)

freq_text = sprintf(" %.2e THz", omeg*1e12);
title(strcat('$E_z(x,y)$ for $\omega$ = ', freq_text));
%text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
%text(z(sfr+10), -0.9, '\epsilon = \epsilon_0');
%text(z(p+50), -0.9, '\epsilon = 5\epsilon_0')
%text(z(q+10), -0.9, '\epsilon = \epsilon_0')

%%
figfile = strcat(pathhead, "Ez.png");
saveas(f, figfile)
%system ("wslview figs/test.png")

%%
figure(2);
plot(x, real(Ez(:,ceil(end/2))));
