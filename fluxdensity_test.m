%% Finite-Difference Time Domain algorithm (1D)
% 
% For testing of more sophisticated algorithms and conceptual understanding.
%   Principally from "Electromagnetic Simulation using the FDTD method" by Dennis M. Sullivan
%% 

% Use change of variables: ^E = sqrt(eps0/mu0) E

% Beam parameters
% ---------------
c0 = 3e8;
bw = 1e2 % Beamwidth
t0 = 50.0; % Central time of pulse
taup = 20; % Pulse duration
E0 = @(t) exp(-0.5 .* ((t .- t0)./ taup).^2 ); % Source pulse

% Integration parameters
% ---------------
nz = 1000 % Number of x-grid points
dz = 0.005 % Size of x-step

epsr = ones(1,nz);
con = zeros(1,nz);
chi = zeros(1,nz);
epsr(nz/2:3*nz/4) = 2 .* epsr(nz/2:3*nz/4); % Relative permittivity of 2 in second half of domain
con(nz/2:3*nz/4) = 0.01 .+ con(nz/2:3*nz/4); % Conductivity of 0.01 in second half of domain
chi(nz/2:3*nz/4) = 2 .+ con(nz/2:3*nz/4); 

dt = dz / (2 * c0) % Size of t-step
NMAX = 2000; % Number of time steps
Nsaved = 100; % Interval between saves


z = [-nz/2:nz/2-1] .* dz;

[Ex, Hy] = fdtd_DH(E0, epsr, con, dt, dz, nz, NMAX, Nsaved);

clf;
graphics_toolkit gnuplot;
f = figure('visible','off');

for n = 1:rows(Ex)
  % Plot the E-field
  subplot(2,1,1);
    plot(z, Ex(n,:));
    xlabel("z");
    ylabel('E_x(z)');
    xlim([z(1) z(end)]);
    ylim([-1.0 1.0]);
    grid on;
  subplot(2,1,2);
    plot(z, Hy(n,:));
    xlabel("z");
    ylabel('H_y(z)');
    xlim([z(1) z(end)]);
    ylim([-1.0 1.0]);
    grid on;
  file_text=sprintf("figs/output%d.png",n)
  saveas (gca, file_text);
endfor
