%%
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%%
Nx = 450;
Ny = 450;
dx = 2e-7;
dy = 2e-7;
x = (1:Nx)' .* dx;
y = (1:Ny)' .* dy;


%% Define the wavevectors and compute source.
Nf = 31;
omeg0 = 590e12;
omegs = linspace(430e12, 750e12, Nf);
sw = 100e12;
phi = pi/4;
k0 = omegs ./ 3e8;
kx = k0;
ky = 0*k0;
x0 = 400e-7;

A = exp(-4*log(2)*(omegs-omeg0).^2 / sw^2) .* exp(1i*k0.*x0);
Ez = A(1)*planewave(x, y, kx(1), ky(1));
Ez = reshape(Ez, Nx, Ny);

% f = figure
% clf
% [M,c] = contourf(x*1e6, y*1e6, real(Ez)', 64);
% set(c,'LineColor','none')
% xlabel('x [{\mu}m]');
% ylabel('y [{\mu}m]');
% colorbar("EastOutside");

for j = 2:length(omegs)
    Fxy = planewave(x, y, kx(j), ky(j));
    Fxy = reshape(Fxy, Nx, Ny);
    
%     f = figure
%     clf;
%     plot(x*1e6, real(Fxy(:,1)))
%     xlabel('x [{\mu}m]');
    
    Ez = Ez + A(j)*Fxy;
end
%%
% f = figure
% clf
% [M,c] = contourf(xsam*1e6, ysam*1e6, real(Ez)', 64);
% set(c,'LineColor','none')
% xlabel('x [{\mu}m]');
% ylabel('y [{\mu}m]');
% colorbar("EastOutside");

f = figure(1)
clf;
plot(x*1e6, real(Ez(:,1)))
xlabel('x [{\mu}m]');
title('Total field E(x)');

f = figure(2)
clf;
plot(omegs, abs(A))
title('Spectral amplitudes A(k)')