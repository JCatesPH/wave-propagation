% To probe how Octave does FFT
N = 50;
x = (-floor(N/2):floor(N/2))';
kx = (-floor(N/2):floor(N/2))';

mu = 0;
sig = N/10;
f = @(x) exp(-0.5*(x-mu).^2 / sig.^2);
F = f(x);
G = fftshift(fft(F));

clf;
subplot(2, 1, 1);
  plot(x, F);
subplot(2, 1, 2);
  plot(kx, real(G), "-r", kx, imag(G), "-b", kx, abs(G), "--k");
