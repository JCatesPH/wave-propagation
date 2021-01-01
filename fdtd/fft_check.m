z = [-nz/2:nz/2-1] .* dz;

Eh0 = fftshift(fft(Ex(7,:)));
Ehf = fftshift(fft(Ex(38,:)));
kz = [ -(ceil((nz-1)/2):-1:1), 0, (1:floor((nz-1)/2)) ] ./ (nz * dz);

clf;
graphics_toolkit qt;
f = figure('visible','on');
ax(1) = subplot(2,1,1);
  plot(kz, abs(Eh0).^2);
  xlim([-20 20]);
  title(ax(1),"Initial pulse");
  xlabel('k_z');
  ylabel('E_x(k_z)');
ax(2) = subplot(2,1,2);
  plot(kz, abs(Ehf).^2);
  xlim([-20 20]);
  title(ax(2),"Reflected and transmitted pulse");
  xlabel('k_z');
  ylabel('E_x(k_z)');
  
  
print("figs/transmission_Eh_case1.svg")

