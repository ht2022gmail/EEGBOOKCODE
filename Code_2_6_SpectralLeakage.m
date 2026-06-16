clear all; close all; clc;

%% Demostrating spectral leakage:
% Define the time series:
srate = 1000;
f = 10;
t = -100:1/srate:100;
Ntime = length(t);
x = cos(2*pi*f*t);

% Cut the 10 sec time window:
Theta = (t>=-2.5 & t<=2.5);
NtimeTheta = sum(Theta);
xobs = Theta.*x;

% FFT:
f = [0:(Ntime-1)]/Ntime*srate;
fx = fft(x);
fxobs = fft(xobs);

%% Visualing the results:
figure(1); clf; hold on;
ff = fill(2.5*[-1 1 1 -1 -1], 1.1*[-1 -1 1 1 -1], 'w');
set(ff, 'EdgeColor', 'r');
plot(t, x, 'k', 'LineWidth', 2);
plot(t(t>=-2.5 & t<=2.5), x(t>=-2.5 & t<=2.5), 'r', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('EEG (\mu V)');
xlim([-5 5]); ylim([-1.2 1.2])
set(gca, 'plotboxaspectratio', [3 1 1]);
set(gca, 'XTick', [-10:10], 'ytick', [-1:1]);
set(gcf, 'color', 'w')

figure(3); clf; hold on;
plot(f, 10*log(abs(fx).^2/(srate*Ntime)), 'k', 'LineWidth', 2);
xlim([9 11]); ylim([-200 50])
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
set(gca, 'plotboxaspectratio', [1.5 1 1]);
set(gca, 'XTick', [9:11], 'ytick', [-200:50:50]);
set(gcf, 'color', 'w')
box on

figure(4); clf;
plot(f, 10*log(abs(fxobs).^2/(srate*NtimeTheta)), 'r', 'LineWidth', 2);
xlim([9 11]); ylim([-200 50])
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
set(gca, 'plotboxaspectratio', [1.5 1 1]);
set(gca, 'XTick', [9:11], 'ytick', [-200:50:50]);
set(gcf, 'color', 'w')


