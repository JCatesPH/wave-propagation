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
t0 = 1e-9; % Central time of pulse
tp = 1e-10; % Pulse duration
E0 = @(t) exp(-0.5 .* ((t .- t0)./ tp).^2 ); % Source pulse
tau = 1e-9; % Relaxation time of material
epsr1 = 2.; % Rel Perm in second material
con1 = 0.001; % Conductivity in second material
chi1 = 1.; % Linear susceptibility in second material

% Integration parameters
% ---------------
nz = 1000 % Number of x-grid points
dz = 5e-3 % Size of x-step
z1 = floor(nz/2);
%z2 = floor(3*nz/4);
z2 = nz;

epsr = ones(1,nz);
con = zeros(1,nz);
chi = zeros(1,nz);
epsr(z1:z2) = epsr1 .* epsr(z1:z2); % Rel permittivity in second material
con(z1:z2) = con1 .+ con(z1:z2); % Conductivity in second material
chi(z1:z2) = chi1 .+ con(z1:z2); 

dt = dz / (2 * c0) % Size of t-step
NMAX = 2000; % Number of time steps
Nsaved = 40; % Interval between saves


z = [-nz/2:nz/2-1] .* dz;

[Ex, Hy] = fdtd_DHchi(E0, epsr, con, chi, tau, dt, dz, nz, NMAX, Nsaved);

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
