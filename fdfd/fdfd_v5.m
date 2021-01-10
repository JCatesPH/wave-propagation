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
nz = 2501;
dz = 5e-10;
z = [1:nz] .* dz;

% Lorentz media parameters
epsinf = eps0;
epss = 2.25 * eps0;
delta = 0.28e16;
omeg0 = 4e16;
eps = @(omeg) epsinf .+ (epss .- epsinf) .* omeg0 .^2 ./ (omeg0 .^2 .+ 2i .* omeg .* delta .- omeg.^2);

%---- PML absorbing boundary condition----
sig = zeros(nz,1); % initialize conductivity array
d = 200; % width of PML layer
m = 3; % polynomial order for grading sigma array (pp 292, Taflove)
neta = sqrt(mu0/eps0); 
R = 1e-12; % required reflectivity
sigma_max = -(m + 1) * log(R) / (2 * neta * d * dz);
Pright = ((1:d+1) ./ d) .^ m * sigma_max;
sig(nz-d:nz) = Pright; % lossy conductivity profile
sig(1:d+1) = fliplr(Pright);
%sigdag = sig .* mu0 ./ eps0; % Eq 7.8 Taflove, pp 275

% Initialize simulation values
mu = mu0;

% Specify Scattered-Field and Total-Field regions
sfr = floor(nz/5);
tfr = nz - sfr;
Q = diag([ones(sfr,1); zeros(tfr,1); ones(sfr,1); zeros(tfr,1)]);
%perm = eps(omeg)

Dze = (diag(-1 .* ones(nz,1)) + diag(ones(nz-1,1),1));
Dzh = (diag(-1 .* ones(nz-1,1),-1) + diag(ones(nz,1)));

Dze = Dze ./ dz;
Dzh = Dzh ./ dz;

% Periodic boundary condition
%Dze(end,1) = 1;
%Dzh(1,end) = -1;

% Freq range initialization
Nf = 10;
omegs = linspace(0, pi*6e16, Nf);
Gamma = zeros(Nf, 1);

% Specify where boundary is
p = floor(0.5 * nz);
q = floor(0.75 * nz);
r = nz - q;

graphics_toolkit gnuplot;
f = figure('visible','off');

for i = 1:Nf
  omeg = omegs(i);
  %eps1 = eps(omeg);
  eps1 = eps0;
  eps2 = eps0;
  Eco = diag([-(1i .* omeg .+ sig(1:p) ./ eps0) .* eps0 .* ones(p,1); 
    -(1i .* omeg .+ sig(p+1:q) ./ eps1) .* eps1 .* ones(q-p,1); 
    -(1i .* omeg .+ sig(q+1:end) ./ eps2) .* eps2 .* ones(r,1)]);
  Hco = diag([-(1i .* omeg .+ sig(1:p) ./ eps0) .* mu0 .* ones(p,1); 
    -(1i .* omeg .+ sig(p+1:q) ./ eps1) .* mu .* ones(q-p,1);
    -(1i .* omeg .+ sig(q+1:end) ./ eps2) .* mu0 .* ones(r,1)]);
  %kz = omeg .* sqrt(mu0 .* [(sig(1:end/2) .+ 1) .* eps0; (sig(end/2+1:end) .+ 1) .* eps1]);
  k0 = omeg .* sqrt(mu0 .* eps0);
  % Source function
  fsrc = @(z) exp(1i .* k0' .* z);
  % Create matrix to be inverted
  A = sparse([Dze, Eco; Hco, Dzh]);
  % Create RHS from source function
  b = (Q*A - A*Q) * [fsrc(z)'; -fsrc(z .+ 0.5*dz)'];
  % Solve system
  x = A \ b;
  Ez = x(1:floor(end/2));
  Vr = abs(max(Ez(d+50:sfr-10)));
  %VSWR = (1 + Vr) / (1 - Vr);
  %Gamma(i) = (VSWR - 1) / (VSWR + 1);
  Gamma(i) = Vr;

  clf;
  plot(z, real(Ez), '-b;Re[E(\omega,z)];');
  xlabel('z [m]');
  xlim([z(1) z(end)]);
  hold on;
  %lims = max(Ez(i,:));
  lims = 1;
  plot([dz*sfr,dz*sfr], [-lims,lims], '--k', "linewidth", 2.5)
  plot([z(p),z(p)], [-lims,lims], '-k', "linewidth", 2.5)
  plot([z(q),z(q)], [-lims,lims], '-k', "linewidth", 2.5)
  %plot(z,Q(1:floor(end/2),1:floor(end/2)) * real(fsrc(z))', "-r;Re[f(z)];");
  legend("hide");
  
  freqtext = sprintf("%.2e", omeg);
  text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
  text(dz*(sfr+10), -0.9, '\epsilon = \epsilon_0');
  text(z(p+10), -0.9, '\epsilon = \epsilon(\omega)')
  text(z(q+10), -0.9, '\epsilon = \epsilon_0')
  
  file_text=sprintf("figs/output%d.png",i);
  saveas (gca, file_text);
  %system ("wslview figs/test.png")
endfor

clf;
plot(omegs, Gamma);
xlabel('\omega');
ylabel('|\Gamma|');
print("figs/reflection.png");

clf;
plot(z, sig);
xlim([z(1) z(end)]);
xlabel('z [m]');
ylabel('\sigma');
print("figs/sig.png");