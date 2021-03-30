%% Usage: [Er,Hr] = fdtd_1d (Ex, Hy, epsr, con, nsteps, Nr)
%  Input:
%    Ex          [1,n] : initial E(t=0;x) 
%    Hy          [1,n] : initial H(t=0;x)
%    epsr [1] or [1,n] : relative permittivity
%    con  [1] or [1,n] : conductivity
%    dt            [1] : step size in t
%    dz            [1] : step size in z
%    nsteps        [1] : number of time steps to make
%    Nr            [1] : interval for saved time steps

%% Author: Jalen Cates <jmcates@jmcates-Surface-Pro-6>
% Created: 2020-12-20

function [Er,Hr] = fdtd_1d(Ex, Hy, epsr, con, dt, dz, nsteps, Nr)
   c0 = 3e8; % Speed of light
   eps0 = 8.85418781762039e-12;
   mu0 = 1.25663706212e-6;
   % Calculate coefficients from parameters
   ca = (1 - con .* dt ./ (2 * eps0 .* epsr)) ./ (1 + con .* dt ./ (2 * eps0 .* epsr));
   cb = c0 * dt ./ (epsr * dz) ./ (1 + con .* dt ./ (2 * eps0 .* epsr));
   cc = c0 * dt / dz;
  
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
     Ex(2:end) = ca .* Ex(2:end) + cb .* (Hy(1:end-1) - Hy(2:end));
     
     % Could add E-field source here.
      % source
     % Force absorbing boundary condition
     Ex(1) = Ea2;
     Ea2 = Ea1;
     Ea1 = Ex(2);
     Ex(end) = Eb2;
     Eb2 = Eb1;
     Eb1 = Ex(end-1);
     
     Hy(1:end-1) = Hy(1:end-1) + cc .* (Ex(1:end-1) - Ex(2:end));
     
     % Decide whether to save step
     if (mod(j,Nr) == 0)
       k = k + 1;
       Er(k,:) = Ex;
       Hr(k,:) = Hy;
     end
   end
end
