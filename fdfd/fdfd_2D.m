%% Finite-Difference Frequency-Domain (FDFD) method
%% Author: Jalen Cates
%% Date: 1/1/2020
%------------------------------------------------------------
% Physical constants
c0 = 3e8; % Speed of light
eps0 = 8.85418781762039e-12;
mu0 = 1.25663706212e-6;
neta = sqrt(mu0/eps0); % vacuum impedance

% Angle of incidence 
thetai = pi/6;

% ------------------------
% % Domain specification
% --
% Domain size parameters
lx = 50; % Number of cells to split vacuum wavelengths into
ly = 50;
Lx = 40; % width of PML layer
Ly = 40;
bufsize = 100;
SFx = Lx + bufsize; % index where total-field region begins
SFy = Ly + bufsize;
Nx = 500 + SFx;
Ny = 501;
% Indices for reflectance and transmittance reference planes
refp_ind = Lx + 5;
trp_ind = Nx - Lx - 5;
% Specify where boundary is
p = floor(0.5 * (Nx-SFx) + SFx);
q = floor(0.8 * (Nx-SFx) + SFx);
r = Nx - q;

% ------------------------
% % Dispersion relations
% --
% Lorentz media parameters
epsinf = eps0;
epss = 2.25 * eps0;
delta = 0.28e16;
omeg0 = 4e16;
eps = @(omeg) epsinf .+ (epss .- epsinf) .* omeg0 .^2 ./ (omeg0 .^2 .+ 2i .* omeg .* delta .- omeg.^2);

% ------------------------
% %  PML absorbing boundary condition
%       Now from Rumpf (2012)
sx = ones(Nx,1);
m = 4; % polynomial order for grading sigma array (pp 292, Taflove)
R = 1e-8; % Desired reflection coefficient
sigmax = -(m+1) * log(R) / (2 * neta * Lx);
smax = 4;
sx0 = smax .* ((1:Lx) ./ Lx) .^ m .+ 1;
sigxp = sigmax .* (sin(pi .* (1:Lx) ./ (2 .* Lx))) .^ 2;
sx(Nx-Lx+1:Nx) = sx0 .* (1 .- 1i .* neta .* sigxp);
sx(1:Lx) = fliplr(sx0 .* (1 .- 1i .* neta .* sigxp));

% Trick to copy vector Ny times
Sx = sx(:,ones(Ny,1));
Sx = Sx(:);

% Specify Scattered-Field and Total-Field regions
tfr = Nx - SFx;
qxy = [ones(SFx,1); zeros(tfr,1)];
Qxy = qxy(:, ones(Ny,1));
Qxy = Qxy(:);
Qxy = diag(sparse(Qxy));

% Difference matrices
dxe = sparse(diag(sparse(-1 .* ones(Nx*Ny,1))) + diag(sparse(ones(Nx*Ny-1,1)),1));
dye = sparse(diag(sparse(-1 .* ones(Nx*Ny,1))) + diag(sparse(ones(Nx*(Ny-1),1)),Nx));

% Freq range initialization
Nf = 2;
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
eps_zz = diag(sparse(eps_zz(:)));

% Find permeability matrices (mu_r=1)
mu_xx = diag(sparse(Sx .^ (-1)));
mu_yy = diag(sparse(Sx));

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
  Dye = sparse(dye + diag(sparse(exp(1i * ki_y(n) * Wy) * ones(Nx,1)), -Nx*(Ny-1)));
  % Adjust operators to grid
  Dxe = sparse(dxe ./ (k0(n) .* dx));
  Dye = sparse(Dye ./ (k0(n) .* dy));
  Dxh = -(Dxe');
  Dyh = -(Dye');
  % Source function
  fsrc = planewave(x, y, ki_x(n), ki_y(n));
  % Create matrix to be inverted
  Ae = sparse(Dxh / mu_yy * Dxe + Dyh / mu_xx * Dye + eps_zz);
  % Create RHS from source function
  b = sparse(Qxy*Ae - Ae*Qxy) * fsrc;
  % Solve system
  ez = Ae \ b;
  Ez = reshape(ez, Nx, Ny);
  [R(n), T(n)] = reflectivity(Ez(refp_ind,:), Ez(trp_ind,:), k0(n), ki_x(n), ki_y(n), x(refp_ind), y, Ny, Wy, 1, sqrt(epsr1));

  efile = sprintf("data/Ez_%02d.csv", n);
  csvwrite(efile, real(ez));
  permfile = sprintf("data/eps_%02d.csv", n);
  csvwrite(permfile, eps_vec);
  xfile = sprintf("data/X_%02d.csv", n);
  csvwrite(xfile, x);
  yfile = sprintf("data/Y_%02d.csv", n);
  csvwrite(yfile, y);
endfor


csvwrite("data/reflectivity.csv", [omegs', R]);

% Angle of refraction from Snell's law
n1 = sqrt(epsr0)
n2 = sqrt(epsr1)
thetat = asin(n1 / n2 * sin(thetai))
R_Fresnel = ((n1*cos(thetai)-n2*cos(thetat))/(n1*cos(thetai)+n2*cos(thetat)))^2
R_mean = mean(R)
R_relerr = abs(R_mean - R_Fresnel) / R_Fresnel
R_std = std(R);
T_mean = mean(T)
T_std = std(T);
Econs = R_mean + T_mean
