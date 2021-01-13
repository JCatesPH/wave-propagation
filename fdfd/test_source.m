% Domain parameters
nx = 5;
ny = 5;
dx = 1e-10;
dy = 1e-10;
x = [1:nx]' .* dx;
y = [1:ny]' .* dy;
[X,Y] = meshgrid(x, y);

omeg = 1e18;
k0 = omeg / c0;
ki_x = k0;
ki_y = 0;

% Source function
%fsrc = @(x,y) exp(1i .* (ki_x .* x + ki_y .* y));
fsrc = planewave(x, y, ki_x, ki_y)

graphics_toolkit gnuplot;
f = figure('visible','off');

clf;
contourf(x, y, real(fsrc));
xlabel('x');
ylabel('y');
print("figs/source_2D.png");