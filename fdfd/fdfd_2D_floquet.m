%% Finite-Difference Frequency-Domain (FDFD) method
%
%   Principally from: http://www.hade.ch/docs/report_FDFD.pdf
%
%% Author: Jalen Cates
% Date: 1/1/2020
%------------------------------------------------------------
% Physical constants
c0 = 3e8; % Speed of light
eps0 = 8.85418781762039e-12;
mu0 = 1.25663706212e-6;
neta = sqrt(mu0/eps0); % vacuum impedance

% Angle of incidence 
thetai = pi/24;

% ------------------------
% % Domain specification
% --
% Domain size parameters
lx = 50; % Number of cells to split vacuum wavelengths into
ly = 50;
Lx = 40; % width of PML layer
bufsize = 100;
sfr = Lx + bufsize; % index where total-field region begins
Nx = 500 + sfr;
Ny = 501;
Nc = Nx*Ny;
% Indices for reflectance and transmittance reference planes
refp_ind = Lx + 5;
trp_ind = Nx - Lx - 5;
% Specify where boundary is
p = floor(0.5 * (Nx-sfr) + sfr);
q = floor(0.8 * (Nx-sfr) + sfr);
r = Nx - q;

% ------------------------
% % Dispersion relations
% --
% Lorentz media parameters
epsinf = eps0;
epss = 2.25 * eps0;
delta = 0.28e16;
omeg0 = 4e16;
eps = @(omeg) epsinf + (epss - epsinf) .* omeg0 .^2 ./ (omeg0^2 + 2i .* omeg .* delta - omeg.^2);

% ------------------------
% %  PML absorbing boundary condition
%       Now from Rumpf (2012)
sx = ones(Nx,1);
m = 3; % polynomial order for grading sigma array (pp 292, Taflove)
R = 1e-8; % Desired reflection coefficient
sigmax = -(m+1) * log(R) / (2 * neta * Lx);
smax = 5;
sx0 = smax .* ((1:Lx) ./ Lx) .^ m + 1;
sigxp = sigmax .* (sin(pi .* (1:Lx) ./ (2 .* Lx))) .^ 2;
sx(Nx-Lx+1:Nx) = sx0 .* (1 - 1i .* neta .* sigxp);
sx(1:Lx) = fliplr(sx0 .* (1 - 1i .* neta .* sigxp));

% Trick to copy vector Ny times
Sx = sx(:,ones(Ny,1));
Sx = Sx(:);

% Specify Scattered-Field and Total-Field regions
tfr = Nx - sfr;
qxy = [ones(sfr,1); zeros(tfr,1)];
Qxy = qxy(:, ones(Ny,1));
Qxy = Qxy(:);
Qxy = spdiags(Qxy, 0, Nc, Nc);

% Difference matrices
%dxe = sparse(diag(sparse(-1 .* ones(Nx*Ny,1))) + diag(sparse(ones(Nx*Ny-1,1)),1));
%dye = sparse(diag(sparse(-1 .* ones(Nx*Ny,1))) + diag(sparse(ones(Nx*(Ny-1),1)),Nx));
o = ones(Nc,1);
dxe = spdiags([-o, o], [0,1], Nc, Nc);
dye = spdiags([-o, o], [0,Nx], Nc, Nc);

% Freq range initialization
Nf = 3;
omegs = linspace(430e12, 750e12, Nf);
%omegs = linspace(1e16, 2e16, Nf);
lambs = 2*pi*c0 ./ omegs;
k0 = omegs ./ c0;
ki_x = k0 .* cos(thetai);
ki_y = k0 .* sin(thetai);
R = zeros(Nf, 1);
T = zeros(Nf, 1);


% Specify relative permittivity
%   NOTE: move to loop if dispersion
%eps1 = eps(omeg);
epsr0 = 1;
%epsr1 = eps(omeg) / eps0;
epsr1 = 2.25;
epsr2 = epsr1;
eps_vec = sx .* [epsr0 .* ones(p,1); 
                epsr1 .* ones(q-p,1); 
                epsr2 .* ones(r,1)];
eps_zz = eps_vec(:,ones(Ny,1));
eps_zz = spdiags(eps_zz(:), 0, Nc, Nc);

% Find permeability matrices (mu_r=1)
mu_xx = spdiags(Sx.^(-1), 0, Nc, Nc);
mu_yy = spdiags(Sx, 0, Nc, Nc);

for n = 1:Nf
  omeg = omegs(n)
  % Choose grid based on wavelengths
  dx = lambs(n) / (lx * max(sqrt([epsr0, epsr1, epsr2])));
  dy = lambs(n) / (ly * max(sqrt([epsr0, epsr1, epsr2])));
  x = [1:Nx]' .* dx;
  y = [1:Ny]' .* dy;
  Wy = Ny * dy; % Width of y for periodic BC and reflections

  % Floquet boundary condition
      %Dye = sparse(dye + diag(sparse(exp(1i * k0 * Wy) * ones(Nx,1)), -Nx*(Ny-1)));
      %Dye = sparse(dye + diag(sparse(exp(1i * ki_y(n) * Wy) * ones(Nx,1)), -Nx*(Ny-1)));
      Dye = spdiags(exp(1i * ki_y(n) * Wy)*o, -Nx*(Ny-1), dye);
  % Basic periodic condition
    %Dye = spdiags(o, -Nx*(Ny-1), dye);
  
  % Adjust operators to grid
  Dxe = dxe ./ (k0(n) .* dx);
  Dye = Dye ./ (k0(n) .* dy);
  Dxh = -(Dxe');
  Dyh = -(Dye');
  % Source function
  fsrc = planewave(x, y, ki_x(n), ki_y(n));
  % Create matrix to be inverted
  Ae = sparse(Dxh / mu_yy * Dxe + Dyh / mu_xx * Dye + eps_zz);
  % Create RHS from source function
  b = (Qxy*Ae - Ae*Qxy) * fsrc;
  % Solve system
  ez = Ae \ b;
  Ez = reshape(ez, Nx, Ny);
  [R(n), T(n)] = reflectivity(Ez(refp_ind,:), Ez(trp_ind,:), k0(n), ki_x(n), ki_y(n), y, Ny, Wy, 1, sqrt(epsr1));

  efile = sprintf("data/Ez_%02d.csv", n);
  writematrix(real(ez), efile);
  permfile = sprintf("data/eps_%02d.csv", n);
  writematrix(eps_vec, permfile);
  xfile = sprintf("data/X_%02d.csv", n);
  writematrix(x, xfile);
  yfile = sprintf("data/Y_%02d.csv", n);
  writematrix(y, yfile);
  srcfile = sprintf("data/fsrc_%02d.csv", n);
  writematrix(fsrc, srcfile);
end


writematrix([omegs', R], "data/reflectivity.csv");

% Angle of refraction from Snell's law
n1 = sqrt(epsr0)
n2 = sqrt(epsr1)
thetat = asin(n1 / n2 * sin(thetai))
R_Fresnel = ((n1*cos(thetai)-n2*cos(thetat))/(n1*cos(thetai)+n2*cos(thetat)))^2
R_mean = mean(R)
R_relerr = abs(R_mean - R_Fresnel) / R_Fresnel * 100
R_std = std(R);
T_mean = mean(T)
T_std = std(T);
Econs = R_mean + T_mean
