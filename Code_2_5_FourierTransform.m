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
x = x - mean(x);

fx = fft(x);
xpow = abs((fx)).^2/Ntime/EEG.srate;
xpow_decibel = 10*log10(xpow);
% 10*log10();


%% Visualizing the results:
figure(1); clf; 
plot(EEG.times/1000, x, 'k-');
xlabel('Time (s)');
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [1.3 1 1], 'xtick', 100:110, 'ytick', -60:20:60);
xlim([100 110])
set(gcf, 'color', 'w');

figure(2); clf; hold on;
plot(f, xpow, 'k');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
xlabel('Frequency (Hz)');
ylabel('EEG power (\mu V^2)');
set(gca, 'PlotBoxAspectRatio', [1.3 1 1]);
xlim([0 40]); ylim([0 200])
set(gcf, 'color', 'w');

figure(3); clf; hold on;
plot(f, xpow_decibel, 'k');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')
xlabel('Frequency (Hz)');
ylabel('10 log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.3 1 1]);
xlim([0 40]); ylim([-30 30])
set(gcf, 'color', 'w');


