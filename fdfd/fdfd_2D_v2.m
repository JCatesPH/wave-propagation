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
neta = sqrt(mu0/eps0); % vacuum impedance

% Angle of incidence 
thetai = 0.00;

% Freq initialization
omeg = 600e12
lamb0 = 2 .* pi .* c0 ./ omeg;
lambx = lamb0 .* sin(thetai);
lamby = lamb0 .* cos(thetai);
k0 = omeg ./ c0;
ki_x = k0 .* sin(thetai);
ki_y = k0 .* cos(thetai);

% ------------------------
% % Domain specification
% --
% Domain size parameters
ldx = 20; % Number of cells to split vacuum wavelengths into
ldy = 20;
Ly = 60; % width of PML layer
bufsize = 100;
SFy = Ly + bufsize; % index where total-field region begins
Nx = 51;
Ny = 1000 + SFy;
Ng = Nx * Ny; % Total number of grid cells
% Indices for reflectance and transmittance reference planes
Tplane = Ly + 5;
Rplane = Ny - Ly - 5;

% Choose grid based on wavelengths
%dx = lambx / (ldx * max(sqrt([epsr0, epsr1, epsr2])));
dx = 1e-8;
dy = lamby / (ldy * max(sqrt([epsr0, epsr1, epsr2])));
x = [1:Nx]' .* dx;
y = [1:Ny]' .* dy;
Wx = Nx * dx; % Width of x for periodic BC and reflections

% Specify where boundaries are
p = floor(0.5 * (Ny-SFy) + SFy);
q = floor(0.8 * (Ny-SFy) + SFy);
r = Ny - q;

% ------------------------
% % Dispersion relations
% --
% Specify the relative eps
epsr0 = 1;
epsr1 = 1;
epsr2 = epsr1;

% ------------------------
% %  PML absorbing boundary condition
%       Now from Rumpf (2012)
sy = ones(Ny,1);
m = 4; % polynomial order for grading sigma array (pp 292, Taflove)
R = 1e-8; % Desired reflection coefficient
sigymax = -(m+1) * log(R) / (Ly * dy);
sigyp = sigymax .* ((1:Ly) ./ Ly) .^ m;
% -- NEW UPML --
syr = 1 - 1i * sigyp / k0;
sy(Ny-Ly+1:end) = syr;
sy(1:Ly) = fliplr(syr);

% Use PML parameters to calculate rel mu and eps
Sy = sy(:,ones(Nx,1))';
Sy = Sy(:);

% ------------------------
% %  Masking matrix
tfr = Ny - SFy;
qxy = [ones(SFy,1); zeros(tfr,1)];
Qxy = qxy(:, ones(Nx,1))';
Qxy = Qxy(:);
Qxy = spdiags(Qxy, 0, Ng, Ng);


% ------------------------
% %  Relative permittivity and permeability
%   
eps_vec = [epsr0 .* ones(p,1); 
             epsr1 .* ones(q-p,1); 
             epsr2 .* ones(r,1)];
eps_zz = eps_vec(:,ones(Nx,1))';
eps_zz = sparse(eps_zz(:));

eps_zz = spdiags(eps_zz .* Sy, 0, Ng, Ng);
% Find permeability matrices (mu_r=1)
mu_xx = spdiags(Sy, 0, Ng, Ng);
mu_yy = spdiags(Sy.^(-1), 0, Ng, Ng);

% ------------------------
% %  Difference matrices
dp = ones(Ng,1);
dm = -1 .* ones(Ng,1);
pbc = zeros(Nx,1);
pbc(1) = 1;
pbc = pbc(:, ones(Ny,1));
pbc = pbc(:);
Dxe = spdiags([dm, dp, pbc]./ (k0 .* dx), [0, 1, -Nx+1], Ng, Ng); % Third entry: 1 -> periodic BC
Dye = spdiags([dm, dp]./ (k0 .* dy), [0, Nx], Ng, Ng);

Dxh = -(Dxe');
Dyh = -(Dye');

% ------------------------
% %  Starting calculation
% Source function
fsrc = planewave(x, y, ki_x, ki_y);
% Create matrix to be inverted
Ae = sparse(Dxh / mu_yy * Dxe + Dyh / mu_xx * Dye + eps_zz);
% Create RHS from source function
b = sparse(Qxy*Ae - Ae*Qxy) * fsrc;
% Solve system
ez = Ae \ b;
Ez = reshape(ez, Nx, Ny);
[R, T] = reflectivity1(Ez(:,Rplane), Ez(:,Tplane), k0, ki_x, ki_y, x, Nx, Wx, 1, sqrt(epsr1));

efile = sprintf("data/Ez.csv");
csvwrite(efile, real(ez));
permfile = sprintf("data/eps.csv");
csvwrite(permfile, eps_vec);
xfile = sprintf("data/X.csv");
csvwrite(xfile, x);
yfile = sprintf("data/Y.csv");
csvwrite(yfile, y);
csvwrite("data/reflectivity.csv", [omeg, R, T]);

% Angle of refraction from Snell's law
n1 = sqrt(epsr0)
n2 = sqrt(epsr1)
thetai
thetat = asin(n1 / n2 * sin(thetai))
R
R_Fresnel = ((n1*cos(thetai)-n2*cos(thetat))/(n1*cos(thetai)+n2*cos(thetat)))^2
R_relerr = abs(R - R_Fresnel) / R_Fresnel
T
Econs = R + T
plotdata_2D_v2