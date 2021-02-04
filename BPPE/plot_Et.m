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
## @deftypefn {} {@var{retval} =} plot_Ex (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author:  <jmcates@DESKTOP-JQUP7RN>
## Created: 2021-01-13

function f = plot_Et (E, t, omeg, filename)
  Ew = abs(fft(E));
  %graphics_toolkit gnuplot;
  %f = figure('visible','off');
  subplot(2, 1, 1);
    plot(t*1e15, real(E));
    xlabel("t [fs]");
    ylabel("E(z_0,t)");
    %title("E(z_0,t)");
   subplot(2, 1, 2);
    plot(omeg*1e-12, Ew(2:end/2+1).^2);
    xlabel('\omega [THz]');
    ylabel('|E(z_0,\omega)|^2');
    %title('|E(z_0,\omega)|^2');
  print(["figs/" filename]);  
endfunction
