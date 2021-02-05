%% Author:  <jmcates@DESKTOP-JQUP7RN>
% Created: 2021-01-10

function [R, T] = reflectivity (Er, Et, k0, kxi, kyi, y, Ny, Wy, nref, ntr)
  kym = kyi - 2 .* pi .* (-floor(Ny/2):floor(Ny/2)) ./ Wy; % Transverse components
  kxr = sqrt((k0 .* nref).^2 - abs(kym).^2); % Longitudinal components
  kxt = sqrt((k0 .* ntr).^2 - abs(kym).^2); % Longitudinal components
  Sr = flipud(fftshift(fft(Er .* exp(-1i * kyi * y)'))) ./ Ny;
  DEm = abs(Sr).^ 2 .* real(kxr ./ kxi);
  R = sum(DEm);
  St = flipud(fftshift(fft(Et .* exp(-1i * kyi * y)'))) ./ Ny;
  DEm = abs(St).^ 2 .* real(kxt ./ kxi);
  T = sum(DEm);
  %R = abs(Ey(25)).^2;
end
