%% Tests the FDFD with a box in TF region.
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Select frequency
%lamb0 = 6e8 .* pi ./ omeg;
lamb0 = 2e-6
omeg = 6e8 * pi ./ lamb0
thetai = 0.0;

%% Mesh specifications
% Number of cells to split vacuum wavelengths into
param2D.dx = 5e-8; 
param2D.dy = 5e-8;
% Number of cells for PML
param2D.Lx = 40; 
% Enforce distance of PML from TF region
param2D.bx = 120; 
% Total number of cells for x and y axis
param2D.Nx = 800 + 2*param2D.Lx + param2D.bx;
param2D.Ny = 201;
Nc = param2D.Nx*param2D.Ny

param2D.Floquet = 1;

%% PML parameters
% Set the polynomial order. (1 <= p < 5)
param2D.sx_m = 3;
% Desired reflectivity coefficient
param2D.sx_R = 1e-10;
% Set s_max. (0 <= s_max <= 5)
param2D.sx_smax = 4;

%% Space parameters
% Set the boundaries 
% -- Currently only for region split into two halves (ADD FUNCTIONALITY LATER)
% param2D.xbounds = [];
% Set the relative permittivities and permeabilities in each region
%   Only for second region in current implementation.
epsr1 = 1.0;
param2D.epsr = epsr1;
param2D.mur = 1

%% Domain parameters
% Specify region that is not vacuum.
%   ff : front face
%   bf : back face
ff = 400;
bf = 600;

%%
eps_zz = ones(param2D.Nx, param2D.Ny);
eps_zz(ff:bf-1,:) = epsr1 * eps_zz(ff:bf-1,:);

Q = zeros(param2D.Nx, param2D.Ny);
Q(ff:bf-1,:) = Q(ff:bf-1,:) + 1;
q = Q(:);


%% Test
chi3 = 5e-12;
dP = spalloc(Nc, Nc, sum(q));
ddP = spalloc(Nc, Nc, sum(q));
PNL = @(x) 3/4*chi3*(q.*abs(x)).^2 .*x;
dPNL = @(x) spdiags(3/2*chi3*(q.*abs(x)).^2 + 3/4*chi3*(q.*abs(x)).*x, 0, dP);
ddPNL = @(x) spdiags(9/2*chi3*(q.*abs(x)).* x, 0, dP);
param2D.NLtol = 1e-12; % Residual error tolerance

pathhead = sprintf("nltest/per_kerr_");
param2D.dir = pathhead;
[x, y, EzL, Ez] = fdfd2D_per_NL(omeg, thetai, eps_zz, PNL, dPNL, ddPNL, param2D);


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

% Add outline of Total-Field Region
plot([x(param2D.Lx+param2D.bx),x(param2D.Lx+param2D.bx)], ...
    [y(1),y(end)], '-k', "linewidth", 1.5)
% plot([x(end-param2D.Lx-param2D.bx),x(end-param2D.Lx-param2D.bx)], ...
%     [y(1),y(end)], '-k', "linewidth", 1.5)

% Add outline of 'box'
plot([x(ff),x(ff)], [y(1),y(end)], '--k', "linewidth", 1.5)
plot([x(bf),x(bf)], [y(1),y(end)], '--k', "linewidth", 1.5)

%quiver(x(sfr+20), y(20), lamby*1e6, lambx*1e6, '-r', "linewidth", 2.5)

freq_text = sprintf(" %.0f THz", omeg*1e-12);
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

%% Plot the linear solution too
f = figure(3);
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

% Add outline of Total-Field Region
plot([x(param2D.Lx+param2D.bx),x(param2D.Lx+param2D.bx)], ...
    [y(1),y(end)], '-k', "linewidth", 1.5)
%plot([x(end-param2D.Lx-param2D.bx),x(end-param2D.Lx-param2D.bx)], ...
%    [y(1),y(end)], '-k', "linewidth", 1.5)

% Add outline of 'box'
plot([x(ff),x(ff)], [y(1),y(end)], '--k', "linewidth", 1.5)
plot([x(bf),x(bf)], [y(1),y(end)], '--k', "linewidth", 1.5)

%quiver(x(sfr+20), y(20), lamby*1e6, lambx*1e6, '-r', "linewidth", 2.5)

freq_text = sprintf(" %.0f THz", omeg*1e-12);
title(strcat('Linear solution $E_z(x,y)$ for $\omega$ = ', freq_text));
%text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
%text(z(sfr+10), -0.9, '\epsilon = \epsilon_0');
%text(z(p+50), -0.9, '\epsilon = 5\epsilon_0')
%text(z(q+10), -0.9, '\epsilon = \epsilon_0')
%%
figfile = strcat(pathhead, "Ez_linear.png");
saveas(f, figfile)