## Copyright (C) 2020 Jalen Cates
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## Usage: [Er,Hr] = fdtd_1d (source, epsr, con, dt, dz, nz, nsteps, Nr)
##  Input:
##    source [func handle] : initial E(t=0;x) 
##    epsr    [1] or [1,n] : relative permittivity
##    con     [1] or [1,n] : conductivity
##    dt               [1] : step size in t
##    dz               [1] : step size in z
##    nz               [1] : number of steps in z
##    nsteps           [1] : number of time steps to make
##    Nr               [1] : interval for saved time steps

## Author: Jalen Cates <jmcates@jmcates-Surface-Pro-6>
## Created: 2020-12-20

function Er = fdtd_plasma (source, vc, omp, dt, dz, nz, nsteps, Nr)
   % Constants
   c0 = 3e8; % Speed of light
   eps0 = 8.85418781762039e-12;
   mu0 = 1.25663706212e-6;
   % Initialize variables
   t = 0;
   Ex = zeros(1,nz);
   Hy = zeros(1,nz);
   Dx = zeros(1,nz);
   Sx = zeros(1,nz);
   sxm2 = zeros(1,nz);
   sxm1 = zeros(1,nz);
   evc = exp(-vc .* dt);
   gax = omp.^2 .* dt ./ vc;
   gax(isnan(gax)) = 0;
   check = sum(isnan(gax))
   check = sum(isnan(evc))
   % Initialize matrix of saved steps
   Er(1,:) = Ex;
   k = 1; % Keep track of steps saved
   
   % Absorbing boundary condition E^(n) (0) = E^(n-2) (1)
   Ea2 = Ex(2);
   Ea1 = Ea2;
   Eb2 = Ex(end-1);
   Eb1 = Eb2;
   for j = 1:nsteps
     Dx(2:end) = Dx(2:end) .+ 0.5 .* (Hy(1:end-1) .- Hy(2:end));
     Ex(1:end-1) = Dx(1:end-1) .- Sx(1:end-1);
     Sx = (1 .+ evc) .* sxm1 .- evc .* sxm2 .+ gax .* (1 .- evc) .* Ex;
     sxm2 = sxm1;
     sxm1 = Sx;
     
     % Could add E-field source here.
     Dx(end/10) = Dx(end/10) + source(t); 
     
     % Force absorbing boundary condition
     Ex(1) = Ea2;
     Ea2 = Ea1;
     Ea1 = Ex(2);
     Ex(end) = Eb2;
     Eb2 = Eb1;
     Eb1 = Ex(end-1);
     
     Hy(1:end-2) = Hy(1:end-2) .+ 0.5 .* (Ex(1:end-2) .- Ex(2:end-1));
     
     % Decide whether to save step
     if (mod(j,Nr) == 1)
       Er(k,:) = Ex;
       k = k + 1;
     endif
     t = t + dt;
   endfor
   Er(end+1,:) = Ex;
endfunction
