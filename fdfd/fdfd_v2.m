%% Finite-Difference Frequency-Domain (FDFD) method
%
%   Principally from: http://www.hade.ch/docs/report_FDFD.pdf
%
%% Author: Jalen Cates
%% Date: 12/22/2020
%------------------------------------------------------------
% Physical constants
c0 = 3e8; % Speed of light
eps0 = 8.85418781762039e-12;
mu0 = 1.25663706212e-6;

% Domain parameters
nx = 1000;
dx = 1e-9;

% Lorentz media parameters
epsinf = eps0;
epss = 2.25 * eps0;
delta = 0.28e16;
omeg0 = 4e16;
eps = @(omeg) epsinf .+ (epss .- epsinf) .* omeg0 .^2 ./ (omeg0 .^2 .- 1i .*omeg .*delta .- omeg.^2);

% Initialize simulation values
omeg = 3e16
%perm = eps(omeg)
perm = eps(omeg)
mu = mu0

x = zeros(2*nx,1); % Field vector with [E1,H1,E2,H2,...]
b = zeros(2*nx,1); % RHS vector with currents/sources
A = zeros(2*nx, 2*nx); % Matrix for differentiation
for i = 2:2:2*nx-1
  A(i,i) = 1i * mu * omeg;
  A(i,i-1) = 1 / dx;
  A(i,i+1) = -1 / dx;
  if (i<nx)
    A(i+1,i+1) = 1i * eps0 * omeg;
  else
    A(i+1,i+1) = 1i * perm * omeg;
  endif
  A(i+1,i) = 1 / dx;
  A(i+1,i+2) = -1 / dx;
  b(i+1) = sin(omeg * sqrt(eps0*mu) * i * dx) * exp(-0.5 * (i - nx / 4)^2 / 500);
endfor
A(1,1) = 1.;
A(end,end) = 1.;

%b = A * b;
x = A \ b;

clf;
graphics_toolkit gnuplot;
f = figure('visible','off');

plot([1:nx].*dx,x(1:2:end), "-b;E(\omega);")
hold on;
plot([1:nx].*dx,b(1:2:end), "-r;source;")
legend("show")
print("figs/test.png")

system ("wslview figs/test.png")