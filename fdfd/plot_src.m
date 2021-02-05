%%
R = readmatrix("data/reflectivity.csv");
omegs = R(:,1);

x1 = readmatrix("data/X_01.csv");
x2 = readmatrix("data/X_02.csv");
x3 = readmatrix("data/X_03.csv");
y1 = readmatrix("data/Y_01.csv");
y2 = readmatrix("data/Y_02.csv");
y3 = readmatrix("data/Y_03.csv");
fsrc1 = readmatrix("data/fsrc_01.csv");
fsrc2 = readmatrix("data/fsrc_02.csv");
fsrc3 = readmatrix("data/fsrc_03.csv");

%% Define the grid for summed source.
Nx = 640;
Ny = 501;
xmin = min([min(x1); min(x2); min(x3)]);
xmax = max([max(x1); max(x2); max(x3)]);
ymin = min([min(y1); min(y2); min(y3)]);
ymax = max([max(y1); max(y2); max(y3)]);

xsam = linspace(xmin, xmax, Nx);
ysam = linspace(ymin, ymax, Ny);

%% Define the wavevectors and compute source.
thetai = pi/24;
k0 = omegs ./ 3e8;
kx = k0 * cos(thetai);
ky = k0 * sin(thetai);

A = [1; 3; 1];
Ez = A(1)*planewave(xsam, ysam, kx(1), ky(1));
Ez = reshape(Ez, length(xsam), length(ysam));

f = figure
clf
[M,c] = contourf(xsam*1e6, ysam*1e6, real(Ez)', 64);
set(c,'LineColor','none')
xlabel('x [{\mu}m]');
ylabel('y [{\mu}m]');
colorbar("EastOutside");

for j = 2:length(omegs)
    Fxy = planewave(xsam, ysam, kx(j), ky(j));
    Fxy = reshape(Fxy, length(xsam), length(ysam));
    
    f = figure
    clf
    [M,c] = contourf(xsam*1e6, ysam*1e6, real(Fxy)', 64);
    set(c,'LineColor','none')
    xlabel('x [{\mu}m]');
    ylabel('y [{\mu}m]');
    colorbar("EastOutside");
    
    Ez = Ez + A(j)*Fxy;
end
%%
f = figure
clf
[M,c] = contourf(xsam*1e6, ysam*1e6, real(Ez)', 64);
set(c,'LineColor','none')
xlabel('x [{\mu}m]');
ylabel('y [{\mu}m]');
colorbar("EastOutside");