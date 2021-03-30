%% Finite-Difference Time Domain algorithm (1D)
% 
% For testing of more sophisticated algorithms and conceptual understanding.
%   Principally from "Electromagnetic Simulation using the FDTD method" by Dennis M. Sullivan
%% 

% Use change of variables: ^E = sqrt(eps0/mu0) E

% Beam parameters
% ---------------
bw = 1e2 % Beamwidth
t0 = 40.0; 
taup = 12;

% Integration parameters
% ---------------
nz = 4000 % Number of x-grid points
dz = 0.1 % Size of x-step

epsr = ones(1,nz);
con = zeros(1,nz);
epsr(nz/2:end) = 12.4 .* epsr(nz/2:end); % Relative permittivity of GaAs in second half of domain
con(nz/2:end) = 5e-4 + con(nz/2:end); % Conductivity of GaAs in second half of domain

c0 = 3e8;
dt = dz / (2 * c0) % Size of t-step
NMAX = 10000; % Number of time steps
Nsaved = 200; % Interval between saves

E0 = @(z) exp(-z.^2 / bw); % Initial pulse
z = [-nz/4:3*nz/4] .* dz;


Ei = E0(z);
Hi = zeros(1,length(Ei));

plot(z, Ei);

[Ex, Hy] = fdtd_1d(Ei, Hi, epsr, con, dt, dz, NMAX, Nsaved);

clf;
graphics_toolkit gnuplot
f = figure('visible','off')
dimEx = size(Ex);

for n = 1:dimEx(1)
  % Plot the E-field
  subplot(2,1,1);
    plot(z, Ex(n,:));
    xlabel("z");
    ylabel('E_x(z)');
    ylim([-1.0 1.0]);
    grid on;
  subplot(2,1,2);
    plot(z, Hy(n,:));
    xlabel("z");
    ylabel('H_y(z)');
    ylim([-1.0 1.0]);
    grid on;
  file_text=sprintf("figs/output%d.png",n)
  saveas (gca, file_text);
end
