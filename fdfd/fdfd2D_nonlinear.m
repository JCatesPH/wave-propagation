function [x, y, Ez] = fdfd2D_nonlinear(omeg, thetai, eps_zz, PNL, dPNL, params)
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

% ------------------------
% % Domain calculations
% --
SFx = params.Lx + params.bx; % index where total-field region begins
SFy = params.Ly + params.by;
Nc = params.Nx*params.Ny;

% ------------------------
% %  PML absorbing boundary condition
%       Now from Rumpf (2011) and Shin (2012)
sx = ones(params.Nx,1);
sx0 = 1 + params.sx_smax .* ((1:params.Lx) ./ params.Lx) .^ params.sx_m;
sigmax_x = -(params.sx_m+1) * log(params.sx_R) / (2 * neta * params.Lx);
%sigmax = params.sx_sigmax;
sigxp = sigmax_x .* sin(pi .* (1:params.Lx) ./ (2 .* params.Lx)) .^ 2;
sx(params.Nx-params.Lx+1:params.Nx) = sx0 .* (1 - 1i .* neta .* sigxp);
sx(1:params.Lx) = fliplr(sx0 .* (1 - 1i .* neta .* sigxp));

sy = ones(params.Nx,1);
sy0 = 1 + params.sy_smax .* ((1:params.Ly) ./ params.Ly) .^ params.sx_m;
sigmax_y = -(params.sy_m+1) * log(params.sy_R) / (2 * neta * params.Ly);
%sigmax = params.sx_sigmax;
sigyp = sigmax_y .* sin(pi .* (1:params.Ly) ./ (2 .* params.Ly)) .^ 2;
sy(params.Ny-params.Ly+1:params.Ny) = sy0 .* (1 - 1i .* neta .* sigyp);
sy(1:params.Ly) = fliplr(sy0 .* (1 - 1i .* neta .* sigyp));


% Trick to copy vector params.Ny times
Sx = sx(:,ones(params.Ny,1));
Sx = Sx(:);
Sy = sy(:,ones(params.Ny,1));
Sy = Sy(:);

% ------------------------
% % Specify Scattered-Field and Total-Field regions
%
qmat = ones(params.Nx, params.Ny);
qmat(SFx+1:end-SFx, SFy+1:end-SFy) = 0*qmat(SFx+1:end-SFx, SFy+1:end-SFy);
qxy = qmat(:);
Qxy = spdiags(qxy, 0, Nc, Nc);


% ------------------------
% % Compute tensors for perms
%
eps_zz = Sx .* Sy .* eps_zz(:);
eps_zz = spdiags(eps_zz, 0, Nc, Nc);

mu_xx = spdiags(Sy ./ Sx, 0, Nc, Nc);
mu_yy = spdiags(Sx ./ Sy, 0, Nc, Nc);


% ------------------------
% % Difference matrices
%
o = ones(Nc,1);
dxe = spdiags([-o, o], [0,1], Nc, Nc);
dye = spdiags([-o, o], [0,params.Nx], Nc, Nc);

% Scale operators to grid
Dxe = dxe ./ (k0 .* dx);
Dye = dye ./ (k0 .* dy);

% Find FD operators for magnetic fields
Dxh = -(Dxe');
Dyh = -(Dye');

% Source function
fsrc = planewave(x, y, ki_x, ki_y);
% Create matrix to be inverted
Ae = sparse(Dxh / mu_yy * Dxe + Dyh / mu_xx * Dye + eps_zz);
% Create RHS from source function
b = (Qxy*Ae - Ae*Qxy) * fsrc;
%spy(Ae)

% Solve system
% solveopt = optimoptions('fsolve','Display','iter', ...
%     'Algorithm', 'trust-region', 'SubproblemAlgorithm', 'cg', ...
%     'Jacobian', 'on');
% solveopt = optimoptions('fsolve','Display','iter', 'Algorithm', 'levenberg-marquardt');
ez_lin = Ae \ b;

fprintf("Linear solution found. Starting nonlinear solver.\n");
Xs = nlSolve(ez_lin, Ae, b, PNL, dPNL);
ez = Xs.sol;
fprintf("Nonlinear solver converged. \n Error: %e \n Function calls: %i\n", Xs.err, Xs.ncalls);

% Save the information if directory is passed.
if isfield(params, 'dir')
    % Check if folder exists already. Must parse string.
    spath = split(params.dir, "/");
    ppath = join(spath(1:end-1), "/");
    if isfolder(ppath) == 0
       mkdir(ppath); 
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

