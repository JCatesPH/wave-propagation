N = 1001;
dt = 1e-1;
t = [(0:N) .* dt];

func = @(t) sin(2 * pi .* t) .+ sin(4 .* pi .* t);
%func = @(t) exp(-0.5 .* (t .- t(floor(N/4))).^2);
yt = func(t);

yf = fftshift(fft(yt));
freqs = [ -(ceil((N-1)/2):-1:1), 0, (1:floor((N-1)/2)) ] ./ (N * dt);

clf;
subplot(2,1,1);
  plot(t, yt);
  xlabel("t");
  ylabel('y(t)');
subplot(2,1,2);
  plot(freqs, abs(yf(1:end-1)));
  xlabel("freqs");
  ylabel('y(f)');
