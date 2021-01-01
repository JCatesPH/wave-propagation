%% Finite-Difference Time Domain algorithm (1D)
% 
% For testing of more sophisticated algorithms and conceptual understanding.
%   Principally from "Electromagnetic Simulation using the FDTD method" by Dennis M. Sullivan
%% 

% Use change of variables: ^E = sqrt(eps0/mu0) E

% Beam parameters
% ---------------
t0 = 1e-16;
tp = 1e-15;
freq = 4000e12;
omega = freq .* 2 .* pi ;
func = @(t) sin(omega .* t) .* exp(-0.5 .* ((t-t0)/tp).^2); % Initial pulse

% Integration parameters
% ---------------
nz = 6000 % Number of x-grid points
dz = 1e-9; % Size of x-step
c0 = 3e8;
dt = dz / (2 * c0) % Size of t-step
NMAX = 1000; % Number of time steps
Nsaved = 100; % Interval between saves

z = ([1:nz] .- 1) .* dz;
kz = [ -(ceil((nz-1)/2):-1:1), 0, (1:floor((nz-1)/2)) ] ./ (nz * dz);

% Material parameters
% ---------------
vc = 57e12; % Electron collision frequency
fp = 2000e12; % Plasma resonance frequency

vcv = [0 .* (1:nz/2), vc .* (nz/2+1:nz)];
omp = [0 .* (1:nz/2), 2 .* pi .* fp .* (nz/2+1:nz)]; % Plasma frequency in domains

% ---------------
% Integration
Ex = fdtd_plasma (func, vcv, omp, dt, dz, nz, NMAX, Nsaved);

% ---------------
% Plotting
clf;
graphics_toolkit gnuplot
f = figure('visible','off')

for n = 1:rows(Ex)
  Ehf = fftshift(fft(Ex(n,:)));
  % Plot the E-field
  subplot(2,1,1);
    plot(z, Ex(n,:));
    xlabel("z");
    ylabel('E_x(z)');
    ylim([-1.0 1.0]);
    grid on;
  subplot(2,1,2);
    plot(kz, abs(Ehf).^2);
    xlabel('k_z');
    ylabel('E_x(k_z)');
    grid on;
  file_text=sprintf("figs/output%d.png",n)
  saveas (gca, file_text);
endfor
