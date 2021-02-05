%% Author:  <jmcates@DESKTOP-JQUP7RN>
% Created: 2021-01-10

function Fxy = planewave(xvals, yvals, kx, ky)
  [X,Y] = meshgrid(xvals,yvals);
  Fxy = exp(1i .* (kx .* X + ky .* Y));
  Fxy = Fxy';
  Fxy = Fxy(:);
end
