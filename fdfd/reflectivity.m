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

function [R, T] = reflectivity (Er, Et, k0, kxi, kyi, xp, y, Ny, Wy, nref, ntr)
  kym = kyi .- 2 .* pi .* (-floor(Ny/2):floor(Ny/2)) ./ Wy; % Transverse components
  kxr = sqrt((k0 .* nref).^2 .- abs(kym).^2); % Longitudinal components
  kxt = sqrt((k0 .* ntr).^2 .- abs(kym).^2); % Longitudinal components
  Sr = flipud(fftshift(fft(Er .* exp(-1i * kyi * y)'))) ./ Ny;
  DEm = abs(Sr).^ 2 .* real(kxr ./ kxi);
  R = sum(DEm);
  St = flipud(fftshift(fft(Et .* exp(-1i * kyi * y)'))) ./ Ny;
  DEm = abs(St).^ 2 .* real(kxt ./ kxi);
  T = sum(DEm);
  %R = abs(Ey(25)).^2;
endfunction
