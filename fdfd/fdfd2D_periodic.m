function [x, y, Ez] = fdfd2D(omeg,thetai,eps_zz,params)
%fdfd2D - Compute FDFD solution for given frequency and parameters.
%   Currently, the following assumptions are made:
%       - isotropic so permittivity and permeability are scalars
%       - periodic and y and PML in x
%       - space is split into two half regions in x
%
% ------------------------------------------------------------
% Author: Jalen Cates
% Date: 2/4/2021
% ------------------------------------------------------------
% Physical constants
c0 = 3e8; % Speed of light
eps0 = 8.85418781762039e-12;
mu0 = 1.25663706212e-6;
neta = sqrt(mu0/eps0); % vacuum impedance

% ------------------------
% % Free space field parameters
% --
lamb0 = 2*pi*c0 ./ omeg;
k0 = omeg ./ c0;
ki_x = k0 * cos(thetai);
ki_y = k0 * sin(thetai);
%lambx = lamb0*cos(thetai);

% ------------------------
% % Choose grid based on passed params
if isfield(params, 'lx') && isfield(params, 'ly')
    dx = lamb0 / params.lx;
    dy = lamb0 / params.ly;
elseif isfield(params, 'dx') && isfield(params, 'dy')
    dx = params.dx;
    dy = params.dy;
else
    error("Must specify lx,ly or dx,dy for domain.");
end

x = (1:params.Nx)' .* dx;
y = (1:params.Ny)' .* dy;
Wy = params.Ny * dy; % Width of y for periodic BC and reflections


% ------------------------
% % Domain calculations
% --
sfr = params.Lx + params.bufsize; % index where total-field region begins
Nc = params.Nx*params.Ny;
% Specify where boundary is
p = floor(0.5 * (params.Nx-sfr) + sfr);
q = params.Nx - p;


% ------------------------
% %  PML absorbing boundary condition
%       Now from Rumpf (2012)
sx = ones(params.Nx,1);
sx0 = 1 + params.sx_smax .* ((1:params.Lx) ./ params.Lx) .^ params.sx_m;
sigmax = -(params.sx_m+1) * log(params.sx_R) / (2 * neta * params.Lx);
%sigmax = params.sx_sigmax;
sigxp = sigmax .* sin(pi .* (1:params.Lx) ./ (2 .* params.Lx)) .^ 2;
sx(params.Nx-params.Lx+1:params.Nx) = sx0 .* (1 - 1i .* neta .* sigxp);
sx(1:params.Lx) = fliplr(sx0 .* (1 - 1i .* neta .* sigxp));



% Trick to copy vector params.Ny times
Sx = sx(:,ones(params.Ny,1));
Sx = Sx(:);

% ------------------------
% % Specify Scattered-Field and Total-Field regions
%
tfr = params.Nx - sfr;
qxy = [ones(sfr,1); zeros(tfr,1)];
Qxy = qxy(:, ones(params.Ny,1));
Qxy = Qxy(:);
Qxy = spdiags(Qxy, 0, Nc, Nc);


% ------------------------
% % Compute tensors for perms
%
eps_zz = Sx .* eps_zz(:);
eps_zz = spdiags(eps_zz, 0, Nc, Nc);

mu_xx = spdiags(Sx.^(-1), 0, Nc, Nc);
mu_yy = spdiags(Sx, 0, Nc, Nc);


% ------------------------
% % Difference matrices
%
o = ones(Nc,1);
dxe = spdiags([-o, o], [0,1], Nc, Nc);
dye = spdiags([-o, o], [0,params.Nx], Nc, Nc);

% Choose Floquet or basic periodic BC based on boolean param
if params.Floquet == 1
    Dye = spdiags(exp(1i * ki_y * Wy)*o, -params.Nx*(params.Ny-1), dye);
else
    Dye = spdiags(o, -params.Nx*(params.Ny-1), dye);
end

% Scale operators to grid
Dxe = dxe ./ (k0 .* dx);
Dye = Dye ./ (k0 .* dy);
% Find FD operators for magnetic fields
Dxh = -(Dxe');
Dyh = -(Dye');

% Source function
fsrc = planewave(x, y, ki_x, ki_y);
% Create matrix to be inverted
Ae = sparse(Dxh / mu_yy * Dxe + Dyh / mu_xx * Dye + eps_zz);
% Create RHS from source function
b = (Qxy*Ae - Ae*Qxy) * fsrc;
% Solve system
ez = Ae \ b;

% Save the information if directory is passed.
if isfield(params, 'dir')
    % Check if folder exists already. Must parse string.
    spath = split(params.dir);
    ppath = join(spath(1:end-1), "/");
    if isfolder(ppath) == 0
       mkdir ppath; 
    end
    sxfile = strcat(params.dir, "sx.csv");
    dlmwrite(sxfile, sx);
    efile = strcat(params.dir, "Ez.csv");
    dlmwrite(efile, ez);
    xfile = strcat(params.dir, "X.csv");
    dlmwrite(xfile, x);
    yfile = strcat(params.dir, "Y.csv");
    dlmwrite(yfile, y);
    paramfile = strcat(params.dir, "params.txt");
    paramtab = struct2table(params);
    paramtab = addvars(paramtab, omeg, thetai);
    writetable(paramtab, paramfile);
end

Ez = reshape(ez, params.Nx, params.Ny);

end

