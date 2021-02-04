%% Finite-Difference Frequency-Domain (FDFD) method
%
%   Principally from: http://www.hade.ch/docs/report_FDFD.pdf
%
%% Author: Jalen Cates
%% Date: 1/1/2020
%------------------------------------------------------------
c = 3e8;
eps0 = 8.85e-12;
% % --- Initial Beam Parameters --- % %
%  For 1d case, set x=0,R=0
I0 = 1e16;
wx = 1;
taup = 60e-15;
%omeg0 = 600e12;
omeg0 = 300e12;
phi = 0;
R = 0;

% % --- Domain Parameters --- % %
Nx = 2; % Number of x points
Nt = 3000;  % Number of t points (actually N+1 points)
dx = 1e-8; % Cell size
dt = 1e-15;
#x = (-floor(Nx/2):floor(Nx/2))' .* dx; % Initialize domain variables
x = 0;
kt = 0; % Transverse wave vectors
t = (-floor(Nt/2):floor(Nt/2))' .* dt;
omeg = (-floor(Nt/2):floor(Nt/2))' ./ (Nt*dt);
%omeg = (-floor(Nt/2):floor(Nt/2))' .* omeg0 ./ Nt + omeg0;

% % --- Dispersion relations --- % %
n0 = 1; % start in vacuum
n1 = 1.5; % Index of refraction in second region (could be function)
% Currently chi3 is in Kerr function

kz0 = sqrt(n0.^2 .* omeg.^2 ./c^2 - kt.^2);
kz1 = sqrt(n1.^2 .* omeg.^2 ./c^2 - kt.^2);
kz2 = sqrt(n0.^2 .* omeg.^2 ./c^2 - kt.^2);

% % --- Domain z
Nz = 1e3; % Number of steps in z
z1 = 3e-6; % Boundary 1
z2 = 5.5e-6;
Z = 10e-6; % End of domain
z = linspace(0, Z, Nz)';
dz = Z / Nz;




% Calculate initial beam in real and frequency space
E0 = twocolorpulse(t, x, I0, wx, taup, omeg0, R, phi);
Ap0 = fftshift(fft(E0));
Am0 = fftshift(fft(0.4*E0));

figure(1);
plot_Et(E0, t, omeg, "E0.png");

[Ap1, Am1] = boundary(Ap0, Am0, z1, kz0, kz1, kt);

##Ep1 = ifft(ifftshift(Ap1 .* exp(i*kz1*z1)));
##Em1 = ifft(ifftshift(Am1 .* exp(-i*kz1*z1)));
##figure(2);
##plot(t, Ep1, '-b;E_+(z_1);', t, Em1, '-r;E_-(z_1);');

[Ap2, Am2] = boundary(Ap1, Am1, z2, kz1, kz2, kt);

##Ep2 = ifft(ifftshift(Ap2 .* exp(i*kz2*z2)));
##Em2 = ifft(ifftshift(Am2 .* exp(-i*kz2*z2)));
##figure(3);
##plot(t, Ep2, '-b;E_+(z_2);', t, Em2, '-r;E_-(z_2);');

[Am1, Ap1] = boundary(Am2, Ap2, z2, kz2, kz1, kt);

Ep1 = ifft(ifftshift(Ap1 .* exp(i*kz1*z1)));
Em1 = ifft(ifftshift(Am1 .* exp(-i*kz1*z1)));
figure(2);
plot(t, Ep1, '-b;E_+(z_1);', t, Em1, '-r;E_-(z_1);');

[Am0, Ap0] = boundary(Am1, Ap1, z1, kz1, kz0, kt);

##Ep0 = ifft(ifftshift(Ap0f));
##Em0 = ifft(ifftshift(Am0f));
##figure(5);
##plot(t, Ep0, '-b;E_+(z_0);', t, Em0, '-r;E_-(z_0);');
j0 = floor(z1/dz);
P = @calcKerr;
[Ap, Am] = RK4step(Ap1, Am1, omeg, kz1, z(j0), dz, P, 0);

Ns = 200;
for j = 1:Ns
  zj = z(j0+j);
  [Ap, Am] = RK4step(Ap, Am, omeg, kz1, zj, dz, P, 0);
endfor
##for j = 2:300
##  [Ap, Am] = RK4step(Ap, Am, omeg, kz0, z(j), dz, P, 0);
##endfor

Ep = ifft(ifftshift(Ap).* exp(i*kz1*zj));
Em = ifft(ifftshift(Am).* exp(-i*kz1*zj));
figure(3);
plot(t, Ep, '-b;E_+;', t, Em, '-r;E_-;');

figure(4);
plot(t, real(Ep+Em));
