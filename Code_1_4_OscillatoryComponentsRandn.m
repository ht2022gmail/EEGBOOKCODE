clear all; close all; clc;
load('./data/sampledataEEG.mat');
EEG  = ALLEEG(1); % loading continuous data

%% Fourier transform
X = EEG.data;
[Nelectrodes, Ntime, Ntrials] = size(X);

f = [0:(Ntime-1)]/Ntime*EEG.srate;

% Select an electrode for example:
elec_name = 'Oz';
elec = find(strcmp({EEG.chanlocs.labels}, elec_name));
x = squeeze(X(elec,:,:)); % time * trials
mu = mean(x);
sd = std(x);
xrandn = randn(1, length(x))*sd;

fx = fft(xrandn);
Xpower = 10*log10(1/(EEG.srate*Ntime)*abs((fx).^2)*2);

%% Welch method:
win = EEG.srate * 3;    % 3-second window
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT
[px_welch, f_welch, cl_welch] = ...
    pwelch(xrandn, win, overlap, nfft, EEG.srate, 'ConfidenceLevel', 0.95);

figure(1); clf; hold on

freqs = [1 4 8 13 20 50];
clrs = {'r', 'g', 'b', 'c', 'm'};
for nn=1:5
idx1 = dsearchn(f_welch, freqs(nn));
idx2 = dsearchn(f_welch, freqs(nn+1));
xx = [f_welch(idx1:idx2); f_welch(idx2:(-1):idx1); f_welch(idx1)];
yy = [10*log10(px_welch(idx1:idx2)); -20*ones(length(idx1:idx2),1); 10*log10(px_welch(idx1))];
ff = fill(xx, yy, clrs{nn});
set(ff, 'FaceAlpha', 0.6)
end

plot(f_welch, 10*log10(px_welch), 'k', 'linewidth', 2);
xlabel('Frequency (Hz)');
ylabel('10 log_{10} (power)');
set(gca, 'PlotBoxAspectRatio', [2 1 1]);
axis([0 50 -10 20])

set(gcf, 'color', 'w')

%% Bandpass filtering:
x_delta = bandpass(xrandn, [1 4], EEG.srate);
x_theta = bandpass(xrandn, [4 8], EEG.srate);
x_alpha = bandpass(xrandn, [8 13], EEG.srate);
x_beta = bandpass(xrandn, [13 30], EEG.srate);
x_gamma = bandpass(xrandn, [30 50], EEG.srate);

t = EEG.times/1000;

figure(2); clf; 
subplot(5,1,1); plot(t, x_delta, 'r', 'LineWidth', 2); 
% xlabel('Time (s)'); 
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [8 1 1]);
xlim([10 13]);

subplot(5,1,2); plot(t, x_theta, 'g', 'LineWidth', 2);
% xlabel('Time (s)'); 
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [8 1 1]);
xlim([10 13]);

subplot(5,1,3); plot(t, x_alpha, 'b', 'LineWidth', 2);
% xlabel('Time (s)'); 
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [8 1 1]);
xlim([10 13]);

subplot(5,1,4); plot(t, x_beta, 'c', 'LineWidth', 2);
% xlabel('Time (s)'); 
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [8 1 1]);
xlim([10 13]);

subplot(5,1,5); plot(t, x_gamma, 'm', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [8 1 1]);
xlim([10 13]);

set(gcf, 'color', 'w')

