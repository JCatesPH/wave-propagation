%% Finite-Difference Frequency-Domain (FDFD) method
%
%   Principally from: http://www.hade.ch/docs/report_FDFD.pdf
%
%% Author: Jalen Cates
%% Date: 1/1/2020
%------------------------------------------------------------
c = 3e8;
eps0 = 8.85e-12;
% % --- Domain Parameters --- % %
Nx = 2; % Number of x points
Nt = 200;  % Number of t points (actually N+1 points)
dx = 1e-8; % Cell size
dt = 1e-9;
#x = (-floor(Nx/2):floor(Nx/2))' .* dx; % Initialize domain variables
x = 0;
kt = 0; % Transverse wave vectors
t = (-floor(Nt/2):floor(Nt/2))' .* dt;
omeg = (-floor(Nt/2):floor(Nt/2))' ./ dt;

% % --- Boundaries, z_i
n0 = 1; % start in vacuum
z1 = 10e-6;
n1 = 1.5; % Index of refraction in second region (could be function)
kz0 = sqrt(n0.^2 .* omeg.^2 ./c^2 - kt.^2)

% % --- Initial Beam Parameters --- % %
%  For 1d case, set x=0,R=0
I0 = 1;
wx = 1;
taup = 1e-7;
%omeg0 = 600e12;
omeg0 = 0;
phi = 0;
R = 0;

E0 = twocolorpulse(t, x, I0, wx, taup, omeg0, R, phi);
A0 = flipud(fftshift(fft(E0)));

plot_Et(E0, t, omeg, "E0.png");
plot_Et(A0, omeg, t, "A0.png");

Ap = @(z,Ai) -1i*omeg.^2 ./(2*kz*eps0*c^2) .* calcKerr(Ai) .* exp(1i * kz * z);
odestruct = odeset("AbsTol", "1e-6");

