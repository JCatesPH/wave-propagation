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
nz = 25000;
dz = 5e-12;
z = [1:nz] .* dz;

% Lorentz media parameters
epsinf = eps0;
epss = 2.25 * eps0;
delta = 0.28e16;
omeg0 = 4e16;
eps = @(omeg) epsinf .+ (epss .- epsinf) .* omeg0 .^2 ./ (omeg0 .^2 .+ 2i .* omeg .* delta .- omeg.^2);

%---- PML absorbing boundary condition----
%       Now from Rumpf (2011)
sz = ones(nz,1); % initialize conductivity array
d = 1500; % width of PML layer
m = 2; % polynomial order for grading sigma array (pp 292, Taflove)
neta = sqrt(mu0/eps0); % Impedance
sigmax = 1;
smax = 5;
sz0 = ((1:d) ./ d) .^ m .* smax .+ 1;
sigp = sigmax .* (sin(pi .* (1:d) ./ (2 .* d))) .^ 2;
sz(nz-d+1:nz) = sz0 .* (1 .- 1i .* neta .* sigp);
sz(1:d) = fliplr(sz0 .* (1 .- 1i .* neta .* sigp));
%sigdag = sig .* mu0 ./ eps0; % Eq 7.8 Taflove, pp 275


% Specify Scattered-Field and Total-Field regions
sfr = floor(nz/5);
tfr = nz - sfr;
Qx = sparse(diag([ones(sfr,1); zeros(tfr,1)]));
%perm = eps(omeg)

dze = sparse(diag(-1 .* ones(nz,1)) + diag(ones(nz-1,1),1));
dzh = sparse(diag(-1 .* ones(nz-1,1),-1) + diag(ones(nz,1)));

% Periodic boundary condition
%Dze(end,1) = 1;
%Dzh(1,end) = -1;

% Freq range initialization
Nf = 40;
omegs = linspace(1e16, pi*6e16, Nf);
%omegs = linspace(1e13, 1e15, Nf);
Gamma = zeros(Nf, 1);

% Specify where boundary is
p = floor(0.5 * nz);
q = floor(0.75 * nz);
r = nz - q;


for i = 1:Nf
  omeg = omegs(i);
  %eps1 = eps(omeg);
  epsr0 = 1;
  epsr1 = 1;
  epsr2 = 1;
  eps_xx = sparse(diag([sz(1:p) .* epsr0 .* ones(p,1); 
            sz(p+1:q) .* epsr1 .* ones(q-p,1); 
            sz(q+1:end) .* epsr2 .* ones(r,1)]));
  mu_yy = diag(sz);
  %kz = omeg .* sqrt(mu0 .* [(sig(1:end/2) .+ 1) .* eps0; (sig(end/2+1:end) .+ 1) .* eps1]);
  k0 = omeg / c0;
  Dze = dze ./ (k0 .* dz);
  Dzh = dzh ./ (k0 .* dz);
  % Source function
  fsrc = @(z) exp(1i .* k0' .* z);
  % Create matrix to be inverted
  Ax = sparse(Dzh * inv(mu_yy) * Dze + eps_xx);
  % Create RHS from source function
  bx = sparse(Qx*Ax - Ax*Qx) * fsrc(z)';
  % Solve system
  ex = Ax \ bx;
  Vr = abs(max(ex(d+20:sfr-10)));
  %VSWR = (1 + Vr) / (1 - Vr);
  %Gamma(i) = (VSWR - 1) / (VSWR + 1);
  Gamma(i) = Vr;

  outfile = sprintf("data/ex_%02d.csv", i);
  csvwrite(outfile, [z', real(ex)]);

endfor

csvwrite("data/Gamm.csv", [omegs', Gamma]);
csvwrite("data/sz.csv", [z', sz]);
