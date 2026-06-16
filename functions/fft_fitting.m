function dby = fft_fitting(a,X)

fs = 1000;
T = 10*fs;

% weights = [a(1)*exp(1i*a(4)) a(2)*exp(1i*a(5)) a(3)*exp(1i*a(6))];
% weights = [a(1)*exp(1i*a(3)) a(2)*exp(1i*a(4))];
weights = [a(1)*exp(1i*a(3)) a(2)];
% weights = [a(1) a(2)];
fy = fft(real(weights*X));
dby = 10*log10(abs(fy).^2/(fs*T));
dby = dby(10:250);