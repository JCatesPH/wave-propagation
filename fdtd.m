%% Finite-Difference Time Domain algorithm (1D)
% 
% For testing of more sophisticated algorithms and conceptual understanding.
%   Principally from "Electromagnetic Simulation using the FDTD method" by Dennis M. Sullivan
%% 

% Use change of variables: ^E = sqrt(eps0/mu0) E

% Beam parameters
% ---------------
c0 = 3e8; % Speed of light
bw = 5000 % Beamwidth
t0 = 40.0; 
taup = 12;

% Integration parameters
% ---------------
nx = 1000 % Number of x-grid points
dx = 1 % Size of x-step
dt = dx / (2 * c0) % Size of t-step

E0 = @(x) exp(-x.^2 / bw) .* cos(100*x); % Initial pulse

x = [-nx/2:nx/2] .* dx;

Ei = E0(x);
Hi = zeros(1,length(Ei));

plot(x, Ei);

[Ex, Hy] = fdtd_1d(Ei, Hi, 10000, 100);


for n = 1:rows(Ex)
  plot(x, Ex(n,:));
  xlabel("x")
  ylabel("E(x)")
  ylim([-1.0 1.0])
  file_text=sprintf("figs/output%d.png",n)
  saveas (gca, file_text);
endfor
