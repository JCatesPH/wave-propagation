%% Finite-Difference Time Domain algorithm (1D)
% 
% For testing of more sophisticated algorithms and conceptual understanding.
%   Principally from "Electromagnetic Simulation using the FDTD method" by Dennis M. Sullivan
%% 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
c0 = 3e8;


%% Integration parameters
% ---------------
Nx = 750; % Number of x-grid points
Ny = 750; % Number of y-grid points
dx = 0.005; % Size of x-step
dy = 0.005; % Size of y-step

x = [-Nx/2:Nx/2-1] .* dx;
y = [-Ny/2:Ny/2-1] .* dy;
[X,Y] = meshgrid(x,y);

dt = dx / (2 * c0); % Size of t-step
NMAX = 1000; % Number of time steps
Nsaved = 50; % Interval between saves

%% Domain Parameters
tau = 1e-11; % Relaxation time of material
epsr1 = 1; % Rel Perm in second material
con1 = 0.0; % Conductivity in second material
chi1 = 0.0; % Linear susceptibility in second material

% Define 'box' of second material
ff = 250;
bf = 500;
ls = 250;
rs = 500;

epsr = ones(Nx,Ny); % Rel permittivity in second material
con = zeros(Nx,Ny); % Conductivity in second material
chi = zeros(Nx,Ny); % Dispersion term
epsr(ff:bf-1,ls:rs) = epsr1 .* epsr(ff:bf-1,ls:rs);
epsr(bf,ls:rs-1) = epsr1 .* epsr(bf,ls:rs-1);

%% Beam parameters
% ---------------
% Use change of variables: ^E = sqrt(eps0/mu0) E
E0 = 5; % Peak amplitude
w0 = 5e-1; % Beamwidth
t0 = 4e-10; % Central time of pulse
tp = 2e-10; % Pulse duration

src = @(t) gaussianSource(t, X, Y, E0, w0, t0, tp);

%% Check beam
% t = (1:NMAX).*dt;
% omeg = 2*pi/(t(end)-t(1)) * ((1:NMAX) - NMAX/2)';
% 
% E0t = zeros(NMAX, Nx, Ny);
% for n = 1:NMAX
%     E0t(n,:,:) = src(t(n));
% end
% E0w = fft(E0t) ./ NMAX;
% 
% figure(3);
% clf;
% %surf(X, Y, real(E0w));
% colormap hot;
% [M,c] = contourf(x, omeg, real(fftshift(E0w(:,:,floor(end/2)), 1)), 20);
% set(c,'LineColor','none');
% xlabel('$r$');
% ylabel('$\omega$')
% colorbar();

%% Create parameter structure and begin simulation
params.epsr = epsr;
params.con = con;
params.chi = chi;
params.tau = tau;
params.dt = dt;
params.Nx = Nx;
params.Ny = Ny;
params.Nsteps = NMAX;
params.Nsaved = Nsaved

Ez = fdtd_NL(src, params);

fprintf('Simulation finished. Starting plotting and logging.\n');

%%
f = figure(4);
clf;
dimEz = size(Ez);
Ezt = Ez(:,:,1);
nlevels = 30;
colormap parula;
[C,h] = contourf(X, Y, Ezt);
set(h,'LineColor','none');
set(h,'ZDataSource','Ezt');
xlabel('$x$');
ylabel('$y$');
caxis([-1 1]);
colorbar();

for n = 1:1:dimEz(3)
    Ezt = Ez(:,:,n);
    refreshdata
    file_text=sprintf("figs/output%d.png",n);
    saveas(f, file_text);
end
%%
% V = VideoWriter('figs/newfile.avi','Motion JPEG AVI');
% V.FrameRate = 3;
% open(V);
% 
% fig = figure(5);
% clf;
% axis tight manual
% ax = gca;
% ax.XLabel.String = '$z$';
% ax.YLabel.String = '$E_x(z)$';
% ax.XLim = [z(1) z(end)];
% ax.YLim = [-1 1];
% ax.NextPlot = 'replaceChildren';
% %fig.Visible = 'off';
% 
% P = plot(z, Ex(1,:));
% 
% dimEx = size(Ex); 
% F(dimEx(1)) = struct('cdata',[],'colormap',[]);
% 
% for n = 1:dimEx(1)
%     set(P, 'YData', Ex(n,:));
%     drawnow
%     F(n) = getframe(fig);
%     writeVideo(V, F(n));
% end
% 
% close(V);
% %%
% fig = figure(6);
% clf;
% movie(fig,F,5,2)
