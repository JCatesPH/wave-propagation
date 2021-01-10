%% Finite-Difference Frequency-Domain (FDFD) method
%
%   Principally from: http://www.hade.ch/docs/report_FDFD.pdf
%
%% Author: Jalen Cates
%% Date: 1/1/2020
%------------------------------------------------------------
% Physical constants
c0 = 3e8; % Speed of light
eps0 = 8.85418781762039e-12;
mu0 = 1.25663706212e-6;

% Domain parameters
nz = 1000;
dz = 1e-9;
z = [1:nz] .* dz;

% Lorentz media parameters
epsinf = eps0;
epss = 2.25 * eps0;
delta = 0.28e16;
omeg0 = 4e16;
eps = @(omeg) epsinf .+ (epss .- epsinf) .* omeg0 .^2 ./ (omeg0 .^2 .- 1i .*omeg .*delta .- omeg.^2);

% Initialize simulation values
omeg = 3.5e16;
mu = diag(mu0 .* ones(nz,1));
% Specify where boundary is
p = floor(nz/2);
q = nz - p;
perm = diag([eps0 .* ones(p,1); eps(omeg) .* ones(q,1)]);

% Source function
fsrc = @(z) exp(-1i .* omeg .* sqrt(eps0 .* mu0) .* z);
% Specify Scattered-Field and Total-Field regions
p = floor(nz/5);
q = nz - p;
Q = diag([ones(p,1); zeros(q,1); ones(p,1); zeros(q,1)]);
%perm = eps(omeg)

Dze = (diag(-1 .* ones(nz,1)) + diag(ones(nz-1,1),1));
Dzh = (diag(-1 .* ones(nz-1,1),-1) + diag(ones(nz,1)));

Dze = Dze ./ dz;
Dzh = Dzh ./ dz;

% Periodic boundary condition
Dze(end,1) = 1;
Dzh(1,end) = -1;

A = sparse([Dze, -1i .* omeg .* perm; -1i .* omeg .* mu, Dzh]);

b = (Q*A - A*Q) * [fsrc(z)'; fsrc(z .+ 0.5*dz)'];

x = A \ b;

clf;
graphics_toolkit gnuplot;
f = figure('visible','off');

plot(z,real(x(1:floor(end/2))), '-b;Re[E(\omega,z)];');
xlabel('z [m]');
hold on;
plot(z,Q(1:floor(end/2),1:floor(end/2)) * real(fsrc(z))', "-r;Re[f(z)];");
legend("show");
print("figs/test.png");

system ("wslview figs/test.png")