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
## @deftypefn {} {@var{retval} =} calcKerr (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <jmcates@DESKTOP-JQUP7RN>
## Created: 2021-01-13

function PNL = calcKerr(Ap, Am, z, kz)
  chi3 = 1e-20;
  eps0 = 8.85e-12;
  % To change to expected frequency ordering for IFFT
  %   IE: Pos freq -> (zero freq, pos freq, neg freq) 
  Aps = [0; Ap; flipud(conj(Ap(1:end-1)))]; 
  Ams = [0; Am; flipud(conj(Am(1:end-1)))];
  kzs = [0; kz; flipud(kz(1:end-1))];
  N = length(Aps); % For normalization
  Ep = ifft(Aps.*exp(-1i*kzs .* z));
  Em = ifft(Ams.*exp(1i*kzs .* z));
  PNL = fft(chi3 * eps0 * (Ep+Em).^3) / N;
  PNL = PNL(2:end/2+1);
endfunction
