## Copyright (C) 2021 jcate
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
## @deftypefn {} {@var{retval} =} boundary (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jcate <jcate@DESKTOP-JQUP7RN>
## Created: 2021-01-14

function [Ap2, Am2] = boundary (Ap1, Am1, z, kz1, kz2, kt)
  Ap2 = zeros(length(Ap1),1);
  Am2 = zeros(length(Am1),1);
  for j = 1:length(kz1)
    %j
    if kz1(j) == 0
      M1 = [1, 1; 0, 0];
    else
      M1 = [exp(-i*kz1(j)*z), exp( i*kz1(j)*z);
        -i*(kz1(j) .+ kt.^2 ./kz1(j)).*exp(-i*kz1(j)*z), i*(kz1(j) .+ kt.^2./kz1(j)).*exp(-i*kz1(j)*z)];
    endif
    if kz2(j) == 0
      M2 = [1, 1; 0, 0];
    else
      M2 = [exp(-i*kz2(j)*z), exp( i*kz2(j)*z);
        -i*(kz2(j) .+ kt.^2 ./kz2(j)).*exp(-i*kz2(j)*z), i*(kz2(j) .+ kt.^2./kz2(j)).*exp(-i*kz2(j)*z)];
    endif
    A2 = M2 \ M1 * [Ap1(j); Am1(j)];
    Ap2(j) = A2(1);
    Am2(j) = A2(2);
  endfor
endfunction
