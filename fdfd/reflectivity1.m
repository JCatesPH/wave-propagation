## Copyright (C) 2021 
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

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} reflectivity (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <jmcates@DESKTOP-JQUP7RN>
## Created: 2021-01-10

function [R, T] = reflectivity1 (Er, Et, k0, kti, kli, Xp, Np, Wp, n1, n2)
  kt = kti .- 2 .* pi .* (-floor(Np/2):floor(Np/2)) ./ Wp; % Transverse components
  klr = sqrt((k0 .* n1).^2 .- abs(kt).^2); % Longitudinal components
  klt = sqrt((k0 .* n2).^2 .- abs(kt).^2); % Longitudinal components
  Sr = flipud(fftshift(fft(Er .* exp(-1i * kti * Xp)))) ./ Np;
  DEm = abs(Sr).^ 2 .* real(klr ./ kli)';
  R = sum(DEm);
  St = flipud(fftshift(fft(Et .* exp(-1i * kti * Xp)))) ./ Np;
  DEm = abs(St).^ 2 .* real(klt ./ kli)';
  T = sum(DEm);
  %R = abs(Ey(25)).^2;
endfunction
