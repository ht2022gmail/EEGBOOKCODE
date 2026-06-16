clear all; close all; clc;
addpath('./functions/');

fs = 1000; % sampling frequency
dt = 1/1000; % sampling time
f0 = 10; % signal frequency
T = 10; % duration
s = 0.4; % skewness parameter

t = 0:dt:T;
t0 = 0:dt:(1/f0-dt);
f = [0:length(t)-1]/length(t)*fs;

s0 = [0.3  0.7];
clrs = {'r-', 'k:', 'b', 'm'};

figure(1); clf; hold on;
figure(2); clf; hold on;

for nn=1:length(s0)
    s = s0(nn);
    % creating skewed waveform:
    x0 = zeros(size(t0));
    idx = dsearchn(t0', s/f0);
    f1 = f0/(2*s);
    f2 = f0/(2*(1-s));
    x0(1:idx) = sin(2*pi*f1*t0(1:idx)-pi/2);
    x0(idx+1:end) = cos(2*pi*f2*(t0(idx+1:end)-t0(idx)));
    
    repeat = ceil(length(t)/length(t0));
    x = repmat(x0, 1, repeat);
    x = x(1:length(t));
    
    % FFT:
    fx = fft(x);
    px = abs(fx).^2/(fs*length(t));
    dbx = 10*log10(px);
    [peaks, loc_peaks] = findpeaks(dbx);

    figure(1);
    plot(t, x, clrs{nn}, 'linewidth', 2);
    xlim([1 1.4]); ylim([-1.5 1.5]);

    figure(2);
    plot(f, dbx, clrs{nn}, 'linewidth', 2);
    xlim([0 35]); ylim([-80 15]);
end

figure(1); 
set(gca, 'PlotBoxAspectRatio', [2 1 1], 'xtick', [1:0.1:1.4], 'YTick', [-1:0.5:1]);
xlabel('Time (s)', 'FontSize', 14);
ylabel('EEG (\mu V)', 'FontSize', 14);
set(gcf, 'Color', 'w');

figure(2); 
set(gca, 'PlotBoxAspectRatio', [2 1 1], 'xtick', 0:5:35, 'YTick', -80:20:0);
xlabel('Frequency (Hz)', 'FontSize', 14);
ylabel('10 log_{10} (power)', 'Fontsize', 14)
set(gcf, 'Color', 'w');
%% fitting
% x1 = [exp(1i*(2*pi*f0*t)); exp(1i*(2*pi*2*f0*t)); exp(1i*(2*pi*3*f0*t))];
x1 = [exp(1i*(2*pi*f0*t)); exp(1i*(2*pi*2*f0*t))];
fx1 = fft(x1')';
px1 = abs(fx1).^2/(fs*length(t));
dbx1 = 10*log10(px1);
a = nlinfit(x1, dbx(10:250), @fft_fitting, [1 0 0]);
% weights = [a(1)*exp(1i*a(4)) a(2)*exp(1i*a(5)) a(3)*exp(1i*a(6))];
weights = [a(1)*exp(1i*a(3)) a(2)];
% weights = [a(1) a(2)];
x_fit = real(weights*x1);
fx_fit = fft(x_fit);
px_fit = abs(fx_fit).^2/(fs*length(t));
dbx_fit = 10*log10(px_fit);

figure(3); clf; hold on;
plot(t, x, 'k', 'LineWidth', 2);
% plot(t, -x_fit, 'r', 'LineWidth', 2)
xlim([1 1.4]); ylim([-2 2]);
set(gcf, 'Color', 'w');set(gca, 'PlotBoxAspectRatio', [2 1 1], 'xtick', [1:0.1:1.4], 'YTick', [-1:0.5:1]);
xlabel('Time (s)', 'FontSize', 14);
ylabel('EEG (\mu V)', 'FontSize', 14);
set(gcf, 'Color', 'w');

figure(4); clf; hold on;
plot(t, -real(weights(1)*x1(1,:)), 'r', 'LineWidth', 2);
plot(t, -real(weights(2)*x1(2,:)), 'r:', 'LineWidth', 2)
xlim([1 1.4]); ylim([-2 2]);
set(gcf, 'Color', 'w');set(gca, 'PlotBoxAspectRatio', [2 1 1], 'xtick', [1:0.1:1.4], 'YTick', [-1:0.5:1]);
xlabel('Time (s)', 'FontSize', 14);
ylabel('EEG (\mu V)', 'FontSize', 14);
set(gcf, 'Color', 'w');

figure(5); clf; hold on;
plot(f, dbx, 'k', 'LineWidth', 2);
plot(f, dbx_fit, 'r', 'LineWidth', 2)
xlim([0 25]); ylim([-100 20]);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'xtick', 0:5:35, 'YTick', -80:20:0);
xlabel('Frequency (Hz)', 'FontSize', 14);
ylabel('10 log_{10} (power)', 'Fontsize', 14)
set(gcf, 'Color', 'w');

