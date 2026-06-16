clear all; close all; clc;

sigma = 0.3;
fs = 1000;
dt = 1/fs;
t=-5:dt:5;
N = length(t);

%%
w_rect = double(t>=-1/2 & t<=1/2);
w_hann = 1/2*(1-cos(2*pi*(t+1/2))).*w_rect;
w_hamming = (0.54-0.46*cos(2*pi*(t+1/2))).*w_rect;
w_gauss = exp(-t.^2/sigma^2).*w_rect;

figure(1); clf; hold on;
plot(t, w_rect, 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-2 2]); set(gcf, 'color', 'w');
xlabel('Time (s)');

figure(2); clf; hold on;
plot(t, w_hann, 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-2 2]); set(gcf, 'color', 'w');
xlabel('Time (s)');

figure(3); clf; hold on;
plot(t, w_hamming, 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-2 2]); set(gcf, 'color', 'w');
xlabel('Time (s)');

figure(4); clf; hold on;
plot(t, w_gauss, 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-2 2]); set(gcf, 'color', 'w');
xlabel('Time (s)');

%%

p_rect = 10*log10(abs(fft(w_rect)).^2/fs^2);
p_hann = 10*log10(abs(fft(w_hann)).^2/fs^2);
p_hamming = 10*log10(abs(fft(w_hamming)).^2/fs^2);
p_gauss = 10*log10(abs(fft(w_gauss)).^2/fs^2);

f = (0:(N-1))/N*fs - fs/2;

figure(5); clf; hold on;
plot(f, fftshift(p_rect), 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-20 20]); ylim([-100 0]); set(gcf, 'color', 'w');
xlabel('frequency (s)');
ylabel('log(PSD) (dB)');

figure(6); clf; hold on;
plot(f, fftshift(p_hann), 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-20 20]); ylim([-100 0]); set(gcf, 'color', 'w');
xlabel('frequency (s)');
ylabel('log(PSD) (dB)');

figure(7); clf; hold on;
plot(f, fftshift(p_hamming), 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-20 20]); ylim([-100 0]); set(gcf, 'color', 'w');
xlabel('frequency (s)');
ylabel('log(PSD) (dB)');

figure(8); clf; hold on;
plot(f, fftshift(p_gauss), 'k');
set(gca, 'plotboxaspectratio', [2 1 1])
xlim([-20 20]); ylim([-100 0]); set(gcf, 'color', 'w');
xlabel('frequency (s)');
ylabel('log(PSD) (dB)');
