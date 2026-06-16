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
Xchan = squeeze(X(elec,:,:)); % time * trials

fXchan = fft(Xchan);
Xpower = 10*log10(1/(EEG.srate*Ntime)*abs((fXchan).^2)*2);

%% Welch method:
win = EEG.srate * 3;    % 3-second window
overlap = 0.5*win;  % 50% overlap 
nfft = win;         % # samples for FFT
[px_welch, f_welch, cl_welch] = ...
    pwelch(Xchan, win, overlap, nfft, EEG.srate, 'ConfidenceLevel', 0.95);

%% Visualizing the results:
figure(1); clf; hold on;
plot(f, Xpower, 'k');
xlabel('Frequency (Hz)');
ylabel('Power (dB)')

plot(f_welch, 10*log10(px_welch), 'r', 'linewidth', 2);
plot(f_welch, 10*log10(px_welch+cl_welch(:,2)), 'r--', 'linewidth', 2);
plot(f_welch, 10*log10(px_welch-cl_welch(:,2)), 'r--', 'linewidth', 2);
xlabel('Frequency (Hz)');
ylabel('log_{10} (EEG power)');
set(gca, 'PlotBoxAspectRatio', [1.2 1 1]);
xlim([0 40]); ylim([-30 30])
set(gcf, 'color', 'w');

figure(2); clf; 
plot(EEG.times/1000, Xchan, 'k-');
xlabel('Time (s)');
ylabel('EEG (\mu V)');
set(gca, 'PlotBoxAspectRatio', [1.4 1 1], 'xtick', 100:0.2:101, 'ytick', -50:10:60);
xlim([100 101])
set(gcf, 'color', 'w');
