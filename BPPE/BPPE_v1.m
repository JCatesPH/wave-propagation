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
I0 = 1e17;
wx = 1;
taup = 50e-15;
%omeg0 = 600e12;
omeg0 = 800e12;
phi = 0;
R = 0;

% % --- Domain Parameters --- % %
Nx = 2; % Number of x points
Nt = 5000;  % Number of t points (actually N+1 points)
dx = 1e-8; % Cell size
dt = 1e-16;
#x = (-floor(Nx/2):floor(Nx/2))' .* dx; % Initialize domain variables
x = 0;
kt = 0; % Transverse wave vectors
t = [-(floor(Nt/2):-1:1), (1:floor(Nt/2))]' .* dt;
omeg = 2*pi*(1:floor(Nt/2))' ./ (Nt*dt);
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
Z = 20e-6; % End of domain
z = linspace(0, Z, Nz)';
dz = Z / Nz;
zj1 = 300; % Boundary 1
z1 = z(zj1);
zj2 = 500; % Boundary 2
z2 = z(zj2);


% Calculate initial beam in real and frequency space
%   NOTE: Normalization 1/Nt comes from Parseval's theorem
 %    Moved normalization to plotting only, as fft in octave is normalized in nonsymmetric way
E0 = twocolorpulse(t, x, I0, wx, taup, omeg0, R, phi);
Ap = fft(E0); 
Ap = Ap(2:end/2+1);
Am = fft(0*E0);
Am = Am(2:end/2+1);

[pmax, iw] = max(1.1*abs(Ap/Nt).^2);
omeg0 = omeg(iw); % Central freq of pulse

% Plot initial pulse and spectral power
figure(1);
plot_Et(E0, t, omeg, "E0.png");


figure(2);
l = 1;
clf;
P = @calcKerr;
for j = 1:zj1
  zj = z(j);
  if (mod(j,5)==1)
    semilogy(omeg/omeg0, abs(Ap/Nt).^2, '-b;|A_p(z,\omega)|^2;', omeg/omeg0, abs(Am/Nt).^2, '-r;|A_m(z,\omega)|^2;');
    #ylabel('|A_+|^2');
    xlabel('\omega [\omega_0]');
    grid on;
    xlim([0, 5]);
    ylim([1, pmax]);
    titletext = sprintf("Power spectrum at z = %.2f", zj*1e6);
    title([titletext ' {\mu}m']);
    filename = sprintf("figs/Ap_%u.png", l);
    print(filename);
    l = l + 1;
  endif
  [Ap, Am] = RK4step(Ap, Am, omeg, kz0, zj, dz, P, 0);
endfor
[Ap1, Am1] = boundary(Ap, Am, z1, kz0, kz1, kt);

##for j = zj1+1:zj2
##  zj = z(j);
##  if (mod(j,50)==1)
##    clf;
##    semilogy(omeg, abs(Ap).^2, '-b;|A_+|^2;', omeg, abs(Am).^2, '-r;|A_-|^2;');
##    xlabel('\omega');
##    #ylabel('|A_+|^2');
##    ylim([1, pmax]);
##    titletext = sprintf("Spectral Power at %.2e second region", zj);
##    title(titletext);
##    filename = sprintf("figs/Ap_%u.png", j);
##    print(filename);
##  endif
##  [Ap, Am] = RK4step(Ap, Am, omeg, kz0, zj, dz, P, 0);
##endfor