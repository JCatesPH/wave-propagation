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
## @deftypefn {} {@var{retval} =} twocolorpulse (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <jmcates@DESKTOP-JQUP7RN>
## Created: 2021-01-13

function E = twocolorpulse (t, x, I0, wx, taup, omeg0, R, phi)
  % R is relative strength of fundamental and second harmonic
  c = 3e8;
  eps0 = 8.85e-12;
  E = sqrt(2*I0/(c*eps0)) * exp(-2*log(2)* (x./ wx).^2) .* ( sqrt(1-R) * exp(-2*log(2)* (t./ taup).^2) .* cos(omeg0*t) .+ sqrt(R) * exp(-8*log(2)* (t./ taup).^2) .* cos(2*omeg0 .* t .+ phi));
endfunction
