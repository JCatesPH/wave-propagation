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
## @deftypefn {} {@var{retval} =} customRK4 (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <jmcates@DESKTOP-JQUP7RN>
## Created: 2021-01-13

function [Apf, Amf] = RK4step(Api, Ami, omeg, kz, z, dz, P, J)
  eps0 = 8.85e-12;
  c = 3e8;
  
  kp = zeros(length(Api),4);
  km = zeros(length(Ami),4);
  fp = @(zj, Apj, Amj) -i*omeg.^2 ./ (2*kz*eps0*c^2) .* P(Apj,Amj,zj,kz) .* exp(i * kz * zj); %- omeg./(2*kz*eps0*c^2) .* J(z,Ap) .* exp(i * kz .* z);
  fm = @(zj, Apj, Amj) i*omeg.^2 ./ (2*kz*eps0*c^2) .* P(Apj,Amj,zj,kz) .* exp(-i * kz * zj);
  
  kp(:,1) = fp(z, Api, Ami);
  km(:,1) = fm(z, Api, Ami);
  
  kp(:,2) = fp(z+dz/2, Api.+dz.*kp(:,1)./2, Ami.+dz.*km(:,1)./2);
  km(:,2) = fm(z+dz/2, Api.+dz.*kp(:,1)./2, Ami.+dz.*km(:,1)./2);

  kp(:,3) = fp(z+dz/2, Api.+dz.*kp(:,2)./2, Ami.+dz.*km(:,2)./2);
  km(:,3) = fm(z+dz/2, Api.+dz.*kp(:,2)./2, Ami.+dz.*km(:,2)./2);

  kp(:,4) = fp(z+dz, Api.+dz.*kp(:,3), Ami.+dz.*km(:,3));
  km(:,4) = fm(z+dz, Api.+dz.*kp(:,3), Ami.+dz.*km(:,3));

  kp(isnan(kp)) = 0.0;
  km(isnan(km)) = 0.0;
  
  Apf = Api + dz.* kp * [1; 2; 2; 1] ./ 6;
  Amf = Ami + dz.* km * [1; 2; 2; 1] ./ 6;
endfunction
