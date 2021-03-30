%% Usage: [Er,Hr] = fdtd_1d (source, epsr, con, dt, dz, params.Nz, nsteps, Nr)
%  Input:
%    source [func handle] : initial E(t=0;x) 
%    epsr    [1] or [1,n] : relative permittivity
%   con     [1] or [1,n] : conductivity
%    dt               [1] : step size in t
%    dz               [1] : step size in z
%    params.Nz               [1] : number of steps in z
%    nsteps           [1] : number of time steps to make
%    Nr               [1] : interval for saved time steps
%% Author: Jalen Cates 
% Created: 2020-12-20

function Er = fdtd_NL(source, params)
   % Constants
   c0 = 3e8; % Speed of light
   eps0 = 8.85418781762039e-12;
   mu0 = 1.25663706212e-6;
   
   % Initialize variables
   t = 0;
   gaz = 1 ./ (params.epsr + (params.con .* params.dt ./ eps0) + (params.chi .*params.dt./params.tau));
   gbz = params.con .* params.dt ./ eps0;
   gcz = params.chi .* params.dt ./ params.tau;
   edt = exp(-params.dt/params.tau);
   Ez = zeros(params.Nx, params.Ny);
   Hx = zeros(params.Nx, params.Ny);
   Hy = zeros(params.Nx, params.Ny);
   Dz = zeros(params.Nx, params.Ny);
   Iz = zeros(params.Nx, params.Ny);
   Sz = zeros(params.Nx, params.Ny);
   % Initialize matrix of saved steps
   Er = zeros(params.Nx, params.Ny, ceil(params.Nsteps/params.Nsaved));
   Er(:,:,1) = Ez;
   k = 1; % Keep track of steps saved
   
   for j = 1:params.Nsteps
     % Update flux density
     Dz(2:end,2:end) = Dz(2:end,2:end) + 0.5 .* ...
         (Hy(2:end,2:end) - Hy(1:end-1,2:end) ...
         - Hx(2:end,2:end) + Hx(2:end,1:end-1));
     % Update electric field
     Ez = gaz .* (Dz - Iz - edt .* Sz);
     % Update current term
     Iz = Iz + gbz .* Ez;
     % Could add E-field source here.
     Dz = Dz + source(t); 
     % Update Debye term
     Sz = edt .* Sz + gcz .* Ez;
     % PML
     
     % Update magnetic fields
     Hx(:,1:end-1) = Hx(:,1:end-1) + 0.5 .* (Ez(:,1:end-1) - Ez(:,2:end));
     Hy(1:end-1,:) = Hy(1:end-1,:) + 0.5 .* (Ez(2:end,:) - Ez(1:end-1,:));
     
     % Decide whether to save step
     if (mod(j,params.Nsaved) == 0)
       Er(:,:,k) = Ez;
       k = k + 1;
     end
     t = t + params.dt;
   end
   Er(:,:,end) = Ez;
end
