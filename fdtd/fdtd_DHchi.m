%% Usage: [Er,Hr] = fdtd_1d (source, epsr, con, dt, dz, nz, nsteps, Nr)
%  Input:
%    source [func handle] : initial E(t=0;x) 
%    epsr    [1] or [1,n] : relative permittivity
%   con     [1] or [1,n] : conductivity
%    dt               [1] : step size in t
%    dz               [1] : step size in z
%    nz               [1] : number of steps in z
%    nsteps           [1] : number of time steps to make
%    Nr               [1] : interval for saved time steps
%% Author: Jalen Cates 
% Created: 2020-12-20

function [Er,Hr] = fdtd_DHchi(source, epsr, con, chi1, tau, dt, dz, nz, nsteps, Nr)
   % Constants
   c0 = 3e8; % Speed of light
   eps0 = 8.85418781762039e-12;
   mu0 = 1.25663706212e-6;
   % Initialize variables
   t = 0;
   gax = 1 ./ (epsr + (con .* dt ./ eps0) + (chi1.*dt./tau));
   gbx = con .* dt ./ eps0;
   gcx = chi1 .* dt ./ tau;
   edt = exp(-dt/tau);
   Ex = zeros(1,nz);
   Hy = zeros(1,nz);
   Dx = zeros(1,nz);
   Ix = zeros(1,nz);
   Sx = zeros(1,nz);
   % Initialize matrix of saved steps
   Er(1,:) = Ex;
   Hr(1,:) = Hy;
   k = 1; % Keep track of steps saved
   
   % Absorbing boundary condition E^(n) (0) = E^(n-2) (1)
   Ea2 = Ex(2);
   Ea1 = Ea2;
   Eb2 = Ex(end-1);
   Eb1 = Eb2;
   for j = 1:nsteps
     Dx(2:end) = Dx(2:end) + 0.5 .* (Hy(1:end-1) - Hy(2:end));
     
     Ex(1:end-1) = gax(1:end-1) .* (Dx(1:end-1) - Ix(1:end-1) - edt .* Sx(1:end-1));
     
     Ix(1:end-1) = Ix(1:end-1) + gbx(1:end-1) .* Ex(1:end-1);
     
     % Could add E-field source here.
     Dx(end/4) = Dx(end/4) + source(t); 
     
     Sx(1:end-1) = edt .* Sx(1:end-1) + gcx(1:end-1) .* Ex(1:end-1);
     % Force absorbing boundary condition
     Ex(1) = Ea2;
     Ea2 = Ea1;
     Ea1 = Ex(2);
     Ex(end) = Eb2;
     Eb2 = Eb1;
     Eb1 = Ex(end-1);
     
     Hy(1:end-2) = Hy(1:end-2) + 0.5 .* (Ex(1:end-2) - Ex(2:end-1));
     
     % Decide whether to save step
     if (mod(j,Nr) == 1)
       Er(k,:) = Ex;
       Hr(k,:) = Hy;
       k = k + 1;
     end
     t = t + dt;
   end
   Er(end+1,:) = Ex;
   Hr(end+1,:) = Hy;
end
