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
   Kx = zeros(params.Nx, params.Ny);
   Ky = zeros(params.Nx, params.Ny);
   % Initialize PML parameters
   xn = 0.333*((1:params.Lx)./params.Lx).^3;
   yn = 0.333*((1:params.Ly)./params.Ly).^3;

   gi2 = ones(params.Nx, params.Ny);
   gi2(1:params.Lx,:) = repmat(fliplr(1./(xn+1)), params.Ny, 1)';
   gi2(end-params.Lx+1:end,:) = repmat(1./(xn+1), params.Ny, 1)';
   
   gj2 = ones(params.Nx, params.Ny);
   gj2(:,1:params.Ly) = repmat(fliplr(1./(yn+1)), params.Nx, 1);
   gj2(:,end-params.Ly+1:end) = repmat(1./(yn+1), params.Nx, 1);
   
   gi3 = ones(params.Nx, params.Ny);
   gi3(1:params.Lx,:) = repmat(fliplr((1-xn)./(xn+1)), params.Ny, 1)';
   gi3(end-params.Lx+1:end,:) = repmat((1-xn)./(xn+1), params.Ny, 1)';
   
   gj3 = ones(params.Nx, params.Ny);
   gj3(:,1:params.Ly) = repmat(fliplr((1-yn)./(yn+1)), params.Nx, 1);
   gj3(:,end-params.Ly+1:end) = repmat((1-yn)./(yn+1), params.Nx, 1);

   xn = 0.25*(((1:params.Lx)-0.5)./params.Lx).^3;
   yn = 0.25*(((1:params.Ly)-0.5)./params.Ly).^3;
   
   fi1 = zeros(params.Nx, params.Ny);
   fi1(1:params.Lx,:) = repmat(fliplr(xn), params.Ny, 1)';
   fi1(end-params.Lx+1:end,:) = repmat(xn, params.Ny, 1)';
   
   fj1 = zeros(params.Nx, params.Ny);
   fj1(:,1:params.Ly) = repmat(fliplr(yn), params.Nx, 1);
   fj1(:,end-params.Ly+1:end) = repmat(yn, params.Nx, 1);
   
   fi2 = ones(params.Nx, params.Ny);
   fi2(1:params.Lx,:) = repmat(fliplr(1./(xn+1)), params.Ny, 1)';
   fi2(end-params.Lx+1:end,:) = repmat(1./(xn+1), params.Ny, 1)';
   
   fj2 = ones(params.Nx, params.Ny);
   fj2(:,1:params.Ly) = repmat(fliplr(1./(yn+1)), params.Nx, 1);
   fj2(:,end-params.Ly+1:end) = repmat(1./(yn+1), params.Nx, 1);
   
   fi3 = ones(params.Nx, params.Ny);
   fi3(1:params.Lx,:) = repmat(fliplr((1-xn)./(xn+1)), params.Ny, 1)';
   fi3(end-params.Lx+1:end,:) = repmat((1-xn)./(xn+1), params.Ny, 1)';
   
   fj3 = ones(params.Nx, params.Ny);
   fj3(:,1:params.Ly) = repmat(fliplr((1-yn)./(yn+1)), params.Nx, 1);
   fj3(:,end-params.Ly+1:end) = repmat((1-yn)./(yn+1), params.Nx, 1);
   
   % Initialize matrix of saved steps
   Er = zeros(params.Nx, params.Ny, ceil(params.Nsteps/params.Nsaved));
   Er(:,:,1) = Ez;
   k = 1; % Keep track of steps saved
   
   for j = 1:params.Nsteps
     % Update flux density
     Dz(2:end,2:end) = gi3(2:end,2:end) .* gj3(2:end,2:end) .* Dz(2:end,2:end) ...
         + 0.5 .* gi2(2:end,2:end) .* gj2(2:end,2:end) ...
         .* (Hy(2:end,2:end) - Hy(1:end-1,2:end) ...
         - Hx(2:end,2:end) + Hx(2:end,1:end-1));
     % Could add E-field source here.
     Dz = Dz + source(t); 
     % Update electric field
     Ez = gaz .* (Dz - Iz - edt .* Sz);
     % Update current term
     Iz = Iz + gbz .* Ez;
     % Update Debye term
     Sz = edt .* Sz + gcz .* Ez;
     % PML
     curlEy = Ez(:,1:end-1) - Ez(:,2:end);
     Kx(:,1:end-1) = Kx(:,1:end-1) + curlEy;
     curlEx = Ez(2:end,:) - Ez(1:end-1,:);
     Ky(1:end-1,:) = Ky(1:end-1,:) + curlEx;
     
     % Update magnetic fields
     Hx(:,1:end-1) = fj3(:,1:end-1) .* Hx(:,1:end-1) + ...
         fj2(:,1:end-1) .* 0.5 .* curlEy + fi1(:,1:end-1) .* Kx(:,1:end-1);
     Hy(1:end-1,:) = fi3(1:end-1,:) .* Hy(1:end-1,:) - ...
         0.5 .* fi2(1:end-1,:) .* curlEx + fj1(1:end-1,:) .* Ky(1:end-1,:);
     % Decide whether to save step
     if (mod(j,params.Nsaved) == 0)
       Er(:,:,k) = Ez;
       k = k + 1;
     end
     t = t + params.dt;
   end
   Er(:,:,end) = Ez;
end
